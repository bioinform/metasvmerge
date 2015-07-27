#!/usr/bin/env python

from collections import defaultdict
import shutil

from defaults import *
from vcf_utils import *
from sv_interval import SVInterval, get_gaps_file, interval_overlaps_interval_list, merge_intervals
from pindel_reader import PindelReader
from breakdancer_reader import BreakDancerReader
from breakseq_reader import BreakSeqReader
from cnvnator_reader import CNVnatorReader
from generate_sv_intervals import parallel_generate_sc_intervals
from run_spades import run_spades_parallel
from run_age import run_age_parallel
from generate_final_vcf import convert_metasv_bed_to_vcf
from fasta_utils import get_contigs
from genotype import parallel_genotype_intervals

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


def create_dirs(dirlist):
    for dirname in dirlist:
        if not os.path.isdir(dirname):
            logger.info("Creating directory %s" % (dirname))
            os.makedirs(dirname)


def run_metasv(args):
    sample = args.sample
    reference = args.reference
    pindel_vcf=args.pindel_vcf
    pindel_native=args.pindel_native
    breakdancer_vcf=args.breakdancer_vcf
    breakdancer_native=args.breakdancer_native
    breakseq_vcf=args.breakseq_vcf
    breakseq_native=args.breakseq_native
    cnvnator_vcf=args.cnvnator_vcf
    cnvnator_native=args.cnvnator_native
    gatk_vcf=args.gatk_vcf
    gaps=args.gaps
    filter_gaps=args.filter_gaps
    keep_standard_contigs=args.keep_standard_contigs
    wiggle=args.wiggle
    overlap_ratio=args.overlap_ratio
    workdir=args.workdir
    outdir=args.outdir
    boost_ins=args.boost_ins
    bam=args.bam
    chromosomes=args.chromosomes
    num_threads=args.num_threads
    spades=args.spades
    age=args.age
    disable_assembly=args.disable_assembly
    svs_to_assemble=args.svs_to_assemble
    asm_max_size=args.spades_max_interval_size
    minsvlen=args.minsvlen
    maxsvlen=args.maxsvlen
    inswiggle=args.inswiggle
    enable_per_tool_output=args.enable_per_tool_output
    min_support=args.min_ins_support
    min_support_frac=args.min_ins_support_frac
    max_intervals=args.max_ins_intervals
    stop_spades_on_fail=args.stop_spades_on_fail
    gt_window=args.gt_window
    gt_normal_frac=args.gt_normal_frac
    isize_mean=args.isize_mean
    isize_sd=args.isize_sd
    extraction_max_read_pairs=args.extraction_max_read_pairs
    svs_to_report=args.svs_to_report
    min_mapq=args.min_mapq
    min_avg_base_qual=args.min_avg_base_qual
    min_soft_clip=args.min_soft_clip
    max_soft_clip=args.max_soft_clip
    max_nm=args.max_nm
    min_matches=args.min_matches

    # Check if there is work to do
    if not (
                                        pindel_vcf + breakdancer_vcf + breakseq_vcf + cnvnator_vcf + pindel_native + breakdancer_native + breakseq_native + cnvnator_native):
        logger.error("Nothing to merge since no SV file specified")

    # Create the directories for working
    bedtools_tmpdir = os.path.join(workdir, "bedtools")
    create_dirs([workdir, outdir, bedtools_tmpdir])

    # Reference handling
    if not os.path.isfile(reference + ".fai"):
        logger.error("Reference file %s is not indexed" % (reference))
        return 1

    fasta_handle = pysam.Fastafile(reference) if os.path.isfile(reference) else None
    contigs = get_contigs(reference)
    include_intervals = sorted(
        [SVInterval(contig.name, 0, contig.length, contig.name, "include", length=contig.length) for contig in contigs])

    # Generate the list of contigs to process
    contig_whitelist = set(chromosomes) if chromosomes else set([contig.name for contig in contigs])
    if keep_standard_contigs:
        contig_whitelist &= set(
            [str(i) for i in xrange(1, 23)] + ["chr%d" % (i) for i in xrange(1, 23)] + ["X", "Y", "MT", "chrX", "chrY",
                                                                                        "chrM"])
    logger.info("Only SVs on the following contigs will be reported: %s" % (sorted(list(contig_whitelist))))

    # Load the intervals from different files
    vcf_name_list = [("CNVnator", cnvnator_vcf), ("Pindel", pindel_vcf), ("BreakDancer", breakdancer_vcf),
                     ("BreakSeq", breakseq_vcf), ("HaplotypeCaller", gatk_vcf)]
    native_name_list = [("CNVnator", cnvnator_native, CNVnatorReader),
                        ("Pindel", pindel_native, PindelReader),
                        ("BreakSeq", breakseq_native, BreakSeqReader),
                        ("BreakDancer", breakdancer_native, BreakDancerReader)]

    tools = []
    intervals = {}
    sv_types = set()

    gap_intervals = []
    if filter_gaps:
        if not gaps: gaps = get_gaps_file(contig_whitelist)
        gap_intervals = sorted(load_gap_intervals(gaps))

    # Handles native input
    logger.info("Load native files")
    for toolname, nativename, svReader in native_name_list:
        # If no native file is given, ignore the tool
        if not nativename: continue

        tools.append(toolname)
        intervals[toolname] = defaultdict(list)

        for native_file in nativename:
            for record in svReader(native_file, svs_to_report=svs_to_report):
                interval = record.to_sv_interval()

                if not interval:
                    # This is the case for SVs we want to skip
                    continue
                if not interval_overlaps_interval_list(interval, gap_intervals) and interval.chrom in contig_whitelist:

                    # Check length
                    if interval.length < minsvlen:
                        continue

                    # Set wiggle
                    if interval.sv_type == "INS":
                        interval.wiggle = max(inswiggle, wiggle)
                    else:
                        interval.wiggle = wiggle

                    intervals[toolname][interval.sv_type].append(interval)

        sv_types |= set(intervals[toolname].keys())

    # Handles the VCF input cases, we will just deal with these cases
    logger.info("Load VCF files")
    for toolname, vcfname in vcf_name_list:
        # If no VCF is given, ignore the tool
        if not vcfname:
            continue

        tools.append(toolname)
        intervals[toolname] = {}

        vcf_list = []
        for vcffile in vcfname:
            if os.path.isdir(vcffile):
                logger.info("Will load from per-chromosome VCFs from directory %s for tool %s" % (vcffile, toolname))
                vcf_list += [os.path.join(vcffile, "%s.vcf.gz" % contig.name) for contig in contigs if
                             (not contig_whitelist or contig.name in contig_whitelist)]
            else:
                vcf_list.append(vcffile)

        for vcffile in vcf_list:
            load_intervals(vcffile, intervals[toolname], gap_intervals, include_intervals, toolname, contig_whitelist,
                           minsvlen=minsvlen, wiggle=wiggle, inswiggle=inswiggle, svs_to_report=svs_to_report)
        sv_types |= set(intervals[toolname].keys())

    logger.info("SV types are %s" % (str(sv_types)))
    tool_merged_intervals = {}
    final_intervals = []

    # This will just output per-tool VCFs, no intra-tool merging is done yet
    if enable_per_tool_output:
        logger.info("Output per-tool VCFs")
        for toolname in intervals:
            tool_out = os.path.join(outdir, "%s.vcf" % (toolname.lower()))

            logger.info("Outputting single tool VCF for %s" % (str(toolname)))
            vcf_template_reader = vcf.Reader(open(os.path.join(mydir, "resources/template.vcf"), "r"))
            vcf_template_reader.samples = [sample]

            intervals_tool = []
            tool_out_fd = open(tool_out, "w")
            vcf_writer = vcf.Writer(tool_out_fd, vcf_template_reader)
            chr_intervals_tool = {contig.name: [] for contig in contigs}
            for sv_type in sv_types:
                if sv_type in intervals[toolname]:
                    intervals_tool.extend([copy.deepcopy(interval) for interval in intervals[toolname][sv_type]])
            for interval in intervals_tool:
                # Marghoob says that this is just to fill-in some metadata
                interval.do_validation(overlap_ratio)

                interval.fix_pos()
                chr_intervals_tool[interval.chrom].append(interval)

            for contig in contigs:
                chr_intervals_tool[contig.name].sort()
                for interval in chr_intervals_tool[contig.name]:
                    vcf_record = interval.to_vcf_record(fasta_handle, sample)
                    if vcf_record is not None:
                        vcf_writer.write_record(vcf_record)
            tool_out_fd.close()
            vcf_writer.close()
            logger.info("Indexing single tool VCF for %s" % (str(toolname)))
            pysam.tabix_index(tool_out, force=True, preset="vcf")


    # Do merging here
    logger.info("Do merging")
    for sv_type in sv_types:
        logger.info("Processing SVs of type %s" % sv_type)
        tool_merged_intervals[sv_type] = []

        # Do the intra-tool merging
        logger.info("Intra-tool Merging SVs of type %s" % sv_type)
        for tool in tools:
            logger.debug("Is %s in tool keys? %s" % (sv_type, str(intervals[tool].keys())))
            if sv_type not in intervals[tool]:
                logger.debug("%s not in tool %s" % (sv_type, tool))
                continue
            logger.info("First level merging for %s for tool %s" % (sv_type, tool))
            tool_merged_intervals[sv_type] += merge_intervals(intervals[tool][sv_type])

        # Do the inter-tool merging
        logger.info("Inter-tool Merging SVs of type %s" % sv_type)
        merged_intervals = merge_intervals(tool_merged_intervals[sv_type])

        # Intervals which overlap well with merged_intervals
        intervals1 = []
        # Intervals which do not overlap well with merged_intervals.
        # Used to filter out small intervals which got merged with large intervals
        intervals2 = []

        logger.info("Checking overlaps SVs of type %s" % sv_type)
        for interval in tool_merged_intervals[sv_type]:
            if interval_overlaps_interval_list(interval, merged_intervals, overlap_ratio, overlap_ratio):
                intervals2.append(interval)
            else:
                intervals1.append(interval)
        final_intervals.extend(merge_intervals(intervals1) + merge_intervals(intervals2))

    final_chr_intervals = {contig.name: [] for contig in contigs}
    for interval in final_intervals:
        if minsvlen <= interval.length <= maxsvlen:
            interval.do_validation(overlap_ratio)
            interval.fix_pos()
            final_chr_intervals[interval.chrom].append(interval)

    # This is the merged VCF without assembly, ok for deletions at this point
    logger.info("Output merged VCF without assembly ")
    vcf_template_reader = vcf.Reader(open(os.path.join(mydir, "resources/template.vcf"), "r"))
    vcf_template_reader.samples = [sample]
    preasm_vcf = os.path.join(workdir, "pre_asm.vcf")
    vcf_fd = open(preasm_vcf, "w")
    vcf_writer = vcf.Writer(vcf_fd, vcf_template_reader)

    final_stats = {}

    bed_intervals = []
    for contig in contigs:
        final_chr_intervals[contig.name].sort()
        for interval in final_chr_intervals[contig.name]:
            vcf_record = interval.to_vcf_record(fasta_handle)
            if vcf_record is not None:
                key = (interval.sv_type, "PASS" if interval.is_validated else "LowQual",
                       "PRECISE" if interval.is_precise else "IMPRECISE", tuple(sorted(list(interval.sources))))
                if key not in final_stats:
                    final_stats[key] = 0
                final_stats[key] += 1
                vcf_writer.write_record(vcf_record)
            bed_interval = interval.to_bed_interval(sample)
            if bed_interval is not None:
                bed_intervals.append(bed_interval)
    vcf_fd.close()
    vcf_writer.close()

    # Also save a BED file representation of the merged variants without assembly
    merged_bed = None
    if bed_intervals:
        merged_bed = os.path.join(workdir, "metasv.bed")
        pybedtools.BedTool(bed_intervals).saveas(merged_bed)

    for key in sorted(final_stats.keys()):
        logger.info(str(key) + ":" + str(final_stats[key]))

    final_vcf = os.path.join(outdir, "variants.vcf")

    # Run assembly here
    if not disable_assembly:
        logger.info("Running assembly")
        if spades is None:
            logger.error("Spades executable not specified")
            return 1

        if age is None:
            logger.error("AGE executable not specified")
            return 1

        spades_tmpdir = os.path.join(workdir, "spades")
        age_tmpdir = os.path.join(workdir, "age")

        create_dirs([spades_tmpdir, age_tmpdir])

        assembly_bed = merged_bed

        # this does the improved assembly location finder with softclipped reads
        if boost_ins and "INS" in svs_to_assemble:
            logger.info("Generating intervals for insertions")
            assembly_bed = parallel_generate_sc_intervals([bam.name], list(contig_whitelist), merged_bed, workdir,
                                                          num_threads=num_threads, min_support=min_support,
                                                          min_support_frac=min_support_frac,
                                                          max_intervals=max_intervals, min_mapq=min_mapq,
                                                          min_avg_base_qual=min_avg_base_qual,
                                                          min_soft_clip=min_soft_clip, max_soft_clip=max_soft_clip,
                                                          max_nm=max_nm, min_matches=min_matches)
            logger.info("Generated intervals for assembly in %s" % assembly_bed)

        logger.info("Will run assembly now")

        assembled_fasta, ignored_bed = run_spades_parallel(bam=bam.name, spades=spades, bed=assembly_bed,
                                                           work=spades_tmpdir, pad=SPADES_PAD, nthreads=num_threads,
                                                           chrs=list(contig_whitelist),
                                                           max_interval_size=asm_max_size,
                                                           svs_to_assemble=svs_to_assemble,
                                                           stop_on_fail=stop_spades_on_fail,
                                                           max_read_pairs=extraction_max_read_pairs)
        breakpoints_bed = run_age_parallel(intervals_bed=assembly_bed, reference=reference, assembly=assembled_fasta,
                                           pad=AGE_PAD, age=age, chrs=list(contig_whitelist), nthreads=num_threads,
                                           min_contig_len=AGE_MIN_CONTIG_LENGTH, age_workdir=age_tmpdir)

        final_bed = os.path.join(workdir, "final.bed")
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

        genotyped_bed = parallel_genotype_intervals(final_bed, bam.name, workdir=os.path.join(workdir, "genotyping"),
                                                    nthreads=num_threads, chromosomes=list(contig_whitelist),
                                                    window=gt_window, isize_mean=isize_mean, isize_sd=isize_sd,
                                                    normal_frac_threshold=gt_normal_frac)

        logger.info("Output final VCF file")

        convert_metasv_bed_to_vcf(bedfile=genotyped_bed, vcf_out=final_vcf, sample=sample, pass_calls=False)
    else:
        shutil.copy(preasm_vcf, final_vcf)
        pysam.tabix_index(final_vcf, force=True, preset="vcf")

    logger.info("Clean up pybedtools")

    pybedtools.cleanup(remove_all=True)

    logger.info("All Done!")
