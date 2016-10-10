import shutil
import sys

from defaults import *
from vcf_utils import *
from sv_interval import SVInterval, get_gaps_file
from pindel_reader import PindelReader
from breakdancer_reader import BreakDancerReader
from breakseq_reader import BreakSeqReader
from cnvnator_reader import CNVnatorReader
from softclip_reader import SoftClipReader
from generate_sv_intervals import parallel_generate_sc_intervals
from run_spades import run_spades_parallel
from age import run_age_parallel
from generate_final_vcf import convert_metasv_bed_to_vcf
from fasta_utils import get_contigs
from genotype import parallel_genotype_intervals
from process_sv_calls import load_sv_callers_files,merge_sv_callers_files
from _version import __version__

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


def create_dirs(dirlist):
    for dirname in dirlist:
        if not os.path.isdir(dirname):
            logger.info("Creating directory %s" % (dirname))
            os.makedirs(dirname)


def run_metasv(args):
    logger.info("Running MetaSV %s" % __version__)
    logger.info("Command-line %s" % (" ".join(sys.argv)))
    logger.info("Arguments are " + str(args))
    
    
    # Check if there is work to do
    if not (args.pindel_vcf + args.breakdancer_vcf + args.breakseq_vcf + args.cnvnator_vcf +
            args.pindel_native + args.breakdancer_native + args.breakseq_native + args.cnvnator_native +
            args.manta_vcf + args.lumpy_vcf + args.cnvkit_vcf, args.wham_vcf):
        logger.warning("Nothing to merge since no SV file specified")

    # Simple check for arguments
    if not args.disable_assembly:
        if not args.spades:
            logger.error("Spades executable not specified")
            return os.EX_USAGE

        if not args.age:
            logger.error("AGE executable not specified")
            return os.EX_USAGE

    
    enable_flags = ENABLE_FLAGS[args.restart_step]

    
    # Create the directories for working
    bedtools_tmpdir = os.path.join(args.workdir, "bedtools")
    create_dirs([args.workdir, args.outdir, bedtools_tmpdir])

    # Reference handling
    if not os.path.isfile(args.reference + ".fai"):
        logger.error("Reference file %s is not indexed" % (args.reference))
        return 1

    fasta_handle = pysam.Fastafile(args.reference) if os.path.isfile(args.reference) else None
    contigs = get_contigs(args.reference)
    include_intervals = sorted(
        [SVInterval(contig.name, 0, contig.length, contig.name, "include", length=contig.length) for contig in contigs])

    # Generate the list of contigs to process
    contig_whitelist = set(args.chromosomes) if args.chromosomes else set([contig.name for contig in contigs])
    if args.keep_standard_contigs:
        contig_whitelist &= set(
            [str(i) for i in xrange(1, 23)] + ["chr%d" % (i) for i in xrange(1, 23)] + ["X", "Y", "MT", "chrX", "chrY",
                                                                                        "chrM"])
    logger.info("Only SVs on the following contigs will be reported: %s" % (sorted(list(contig_whitelist))))

    # Load the intervals from different files
    vcf_name_list = [("CNVnator", args.cnvnator_vcf), ("Pindel", args.pindel_vcf),
                     ("BreakDancer", args.breakdancer_vcf),
                     ("BreakSeq", args.breakseq_vcf), ("HaplotypeCaller", args.gatk_vcf),
                     ("Lumpy", args.lumpy_vcf), ("Manta", args.manta_vcf), ("CNVkit", args.cnvkit_vcf),
                     ("WHAM", args.wham_vcf)]
    native_name_list = [("CNVnator", args.cnvnator_native, CNVnatorReader),
                        ("Pindel", args.pindel_native, PindelReader),
                        ("BreakSeq", args.breakseq_native, BreakSeqReader),
                        ("BreakDancer", args.breakdancer_native, BreakDancerReader),
                        ("SoftClip", args.softclip_native, SoftClipReader)]


    gap_intervals = []
    if args.filter_gaps:
        gaps = args.gaps if args.gaps else get_gaps_file(contig_whitelist)
        gap_intervals = sorted(load_gap_intervals(gaps))

    loaded_pickle = load_sv_callers_files(native_name_list=native_name_list, vcf_name_list=vcf_name_list,
                                          workdir=args.workdir, outdir=args.outdir, sample=args.sample,
                                          contigs=contigs, fasta_handle=fasta_handle, mean_read_length=args.mean_read_length, isize_sd=args.isize_sd,
                                          gap_intervals=gap_intervals, contig_whitelist=contig_whitelist, include_intervals=include_intervals,
                                          wiggle=args.wiggle, inswiggle=args.inswiggle, svs_to_report=args.svs_to_report, overlap_ratio=args.overlap_ratio,
                                          minsvlen=args.minsvlen, maxsvlen=args.maxsvlen, 
                                          enable_per_tool_output=args.enable_per_tool_output, enabled=enable_flags[STEP_LOAD])

    merged_bed, preasm_vcf = merge_sv_callers_files(loaded_pickle=loaded_pickle, sample=args.sample, 
                                                    workdir=args.workdir, contigs=contigs, fasta_handle=fasta_handle, overlap_ratio=args.overlap_ratio, 
                                                    minsvlen=args.minsvlen, maxsvlen=args.maxsvlen, wiggle=args.wiggle, inswiggle=args.inswiggle,
                                                    enabled=enable_flags[STEP_MERGE])

    final_vcf = os.path.join(args.outdir, "variants.vcf")
    # Run assembly here
    if not args.disable_assembly:
        logger.info("Running assembly")

        spades_tmpdir = os.path.join(args.workdir, "spades")
        age_tmpdir = os.path.join(args.workdir, "age")

        create_dirs([spades_tmpdir, age_tmpdir])

        assembly_bed = merged_bed

        logger.info("Will run assembly now")

        assembled_fasta, ignored_bed = run_spades_parallel(bams=args.bams, spades=args.spades, spades_options=args.spades_options, bed=assembly_bed,
                                                           work=spades_tmpdir, pad=args.assembly_pad,
                                                           nthreads=args.num_threads,
                                                           chrs=list(contig_whitelist),
                                                           max_interval_size=args.spades_max_interval_size,
                                                           timeout=args.spades_timeout,
                                                           svs_to_assemble=args.svs_to_assemble,
                                                           stop_on_fail=args.stop_spades_on_fail,
                                                           max_read_pairs=args.extraction_max_read_pairs,
                                                           assembly_max_tools=args.assembly_max_tools, enabled=enable_flags[STEP_SPADES_ASSEMBLY])
        breakpoints_bed = run_age_parallel(intervals_bed=assembly_bed, reference=args.reference,
                                           assembly=assembled_fasta,
                                           pad=args.assembly_pad, age=args.age, timeout=args.age_timeout, chrs=list(contig_whitelist),
                                           nthreads=args.num_threads,
                                           min_contig_len=AGE_MIN_CONTIG_LENGTH, min_del_subalign_len=args.min_del_subalign_len,
                                           min_inv_subalign_len=args.min_inv_subalign_len,
                                           age_window=args.age_window,
                                           age_workdir=age_tmpdir, enabled=enable_flags[STEP_AGE_ALIGNMENT])

        final_bed = os.path.join(args.workdir, "final.bed")
        if breakpoints_bed:
            if ignored_bed:
                pybedtools.BedTool(breakpoints_bed) \
                    .cat(pybedtools.BedTool(ignored_bed), postmerge=False) \
                    .sort().saveas(final_bed)
            else:
                pybedtools.BedTool(breakpoints_bed).saveas(final_bed)
        elif ignored_bed:
            pybedtools.BedTool(ignored_bed).sort().saveas(final_bed)
        else:
            final_bed = None

        genotyped_bed = parallel_genotype_intervals(final_bed, args.bams,
                                                    workdir=os.path.join(args.workdir, "genotyping"),
                                                    nthreads=args.num_threads, chromosomes=list(contig_whitelist),
                                                    window=args.gt_window, isize_mean=args.isize_mean,
                                                    isize_sd=args.isize_sd,
                                                    normal_frac_threshold=args.gt_normal_frac, 
                                                    enabled=enable_flags[STEP_GENOTYPE])

        logger.info("Output final VCF file")

        convert_metasv_bed_to_vcf(bedfile=genotyped_bed, vcf_out=final_vcf, workdir=args.workdir, 
                                  sample=args.sample, reference=args.reference, pass_calls=False, 
                                  enabled =enable_flags[STEP_GEN_VCF])
    else:
        shutil.copy(preasm_vcf, final_vcf)
        pysam.tabix_index(final_vcf, force=True, preset="vcf")

    logger.info("Clean up pybedtools")

    pybedtools.cleanup(remove_all=True)

    logger.info("All Done!")

    return os.EX_OK
