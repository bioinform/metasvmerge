import os

import pybedtools
import pickle
import logging
from collections import defaultdict
from sv_interval import interval_overlaps_interval_list, merge_intervals, merge_intervals_recursively
from vcf_utils import *

from defaults import *
from cluster_interval import cluster_intervals_parallel



FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


def load_sv_callers_files(native_name_list=[], vcf_name_list=[], workdir=None, outdir=None, sample=None,
                      contigs=[], fasta_handle=None,
                      mean_read_length=MEAN_READ_LENGTH, isize_sd=ISIZE_SD,
                      gap_intervals=[], contig_whitelist=[], include_intervals=[], 
                      wiggle=WIGGLE, inswiggle=INS_WIGGLE, svs_to_report=SVS_SUPPORTED, overlap_ratio=OVERLAP_RATIO,
                      minsvlen=MIN_SV_LENGTH, maxsvlen=MAX_SV_LENGTH, enable_per_tool_output=False, enabled=True):

    logger = logging.getLogger(load_sv_callers_files.__name__)

    if workdir and not os.path.isdir(workdir):
        os.makedirs(workdir)


    if not enabled:
        loaded_pickle=os.path.join(workdir, "loaded_calls.pickle")
        logger.info("Skipped loading STEP. Will use %s"%(loaded_pickle))
        return loaded_pickle


    # Handles native input
    logger.info("Load native files")
    tools = []
    intervals = {}
    sv_types = set()
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
                BD_min_inv_len = mean_read_length+4*isize_sd
                if toolname=="BreakDancer" and interval.sv_type == "INV" and  abs(interval.length)< BD_min_inv_len:
                    #Filter BreakDancer artifact INVs with size < readlength+4*isize_sd
                    continue
                if not interval_overlaps_interval_list(interval, gap_intervals) and interval.chrom in contig_whitelist:
                
                    # Check length
                    if abs(interval.length) < minsvlen and interval.sv_type not in  ["INS", "ITX", "CTX"]:
                        continue

                    if 0 < abs(interval.length) < minsvlen and interval.sv_type == "INS":
                        continue


                    # Set wiggle
                    if interval.sv_type not in ["ITX","CTX"]:
                        interval.wiggle = max(inswiggle if interval.sv_type == "INS" else 0, wiggle)
                    else:
                        interval.wiggle = TX_WIGGLE
                
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
                           minsvlen=minsvlen, wiggle=wiggle, inswiggle=inswiggle,
                           svs_to_report=svs_to_report, maxsvlen=maxsvlen)
        sv_types |= set(intervals[toolname].keys())

    logger.info("SV types are %s" % (str(sv_types)))

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


        
    loaded_pickle=os.path.join(workdir, "loaded_calls.pickle")
    pickle.dump([intervals, sv_types, tools],open(loaded_pickle,"w"))
    return loaded_pickle
    
    
    
def merge_sv_callers_files_old(loaded_pickle=None, sample=None, workdir=None, contigs=[], fasta_handle=None,
                           overlap_ratio=OVERLAP_RATIO, minsvlen=MIN_SV_LENGTH, 
                           maxsvlen=MAX_SV_LENGTH, wiggle=WIGGLE, inswiggle=INS_WIGGLE, enabled=True):
                           
    logger = logging.getLogger(merge_sv_callers_files.__name__)

    if workdir and not os.path.isdir(workdir):
        os.makedirs(workdir)

    if not enabled:
        merged_bed = os.path.join(workdir, "metasv.bed")
        preasm_vcf = os.path.join(workdir, "pre_asm.vcf")
        logger.info("Skipped merging STEP. Will use %s and %s"%(merged_bed,preasm_vcf))
        return merged_bed, preasm_vcf


    if not loaded_pickle or not os.path.isfile(loaded_pickle):
        logger.warning("No loaded_pickle file specified.")
        return None, None

    intervals, sv_types, tools = pickle.load(open(loaded_pickle,"r"))

    # Do merging here
    logger.info("Do merging")

    tool_merged_intervals = {}
    final_intervals = []
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
        final_intervals.extend(merge_intervals_recursively(tool_merged_intervals[sv_type],overlap_ratio))

    final_chr_intervals = {contig.name: [] for contig in contigs}
    for interval in final_intervals:
        interval.do_validation(overlap_ratio)
        interval.fix_precise_coords()
        interval.fix_pos()
        if minsvlen <= abs(interval.length) <= maxsvlen or interval.sv_type in ["ITX", "CTX"] or (interval.length==0 and interval.sv_type == "INS"):
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

    return merged_bed, preasm_vcf
    
    
    
def merge_sv_callers_files(loaded_pickle=None, sample=None, workdir=None, contigs=[], fasta_handle=None,
                           overlap_ratio=OVERLAP_RATIO, minsvlen=MIN_SV_LENGTH, 
                           maxsvlen=MAX_SV_LENGTH, wiggle=WIGGLE, inswiggle=INS_WIGGLE, enabled=True):
                           
    logger = logging.getLogger(merge_sv_callers_files.__name__)

    if workdir and not os.path.isdir(workdir):
        os.makedirs(workdir)

    if not enabled:
        merged_bed = os.path.join(workdir, "metasv.bed")
        preasm_vcf = os.path.join(workdir, "pre_asm.vcf")
        logger.info("Skipped merging STEP. Will use %s and %s"%(merged_bed,preasm_vcf))
        return merged_bed, preasm_vcf

    merge_workdir = os.path.join(workdir, "merge")
    if  merge_workdir and not os.path.isdir(merge_workdir):
        os.makedirs(merge_workdir)

    if not loaded_pickle or not os.path.isfile(loaded_pickle):
        logger.warning("No loaded_pickle file specified.")
        return None, None

    intervals, sv_types, tools = pickle.load(open(loaded_pickle,"r"))

    # Do merging here
    logger.info("Do merging")

    tool_merged_intervals = {}
    final_intervals = []
    for sv_type in sv_types:
        logger.info("Processing SVs of type %s" % sv_type)
        tool_merged_intervals[sv_type] = []

        cluster_wiggle = 0
        cluster_padding = 0
        if sv_type == "INS":
            cluster_padding = inswiggle
        elif sv_type == "ITX":
            cluster_padding = TX_WIGGLE
        else:
            cluster_padding = wiggle

        # Do the intra-tool merging
        logger.info("Intra-tool Merging SVs of type %s" % sv_type)
        for tool in tools:
            logger.debug("Is %s in tool keys? %s" % (sv_type, str(intervals[tool].keys())))
            if sv_type not in intervals[tool]:
                logger.debug("%s not in tool %s" % (sv_type, tool))
                continue

            logger.info("Making bedfile for %s for tool %s" % (sv_type, tool))
            sv_bedfile = os.path.join(merge_workdir, "%s_%s.bed" % (tool,sv_type))
            sv_bed=pybedtools.BedTool(map(lambda x:x[1].to_cluster_bed_interval(sample,x[0],cluster_padding),
                                          enumerate(intervals[tool][sv_type]))).sort().saveas(sv_bedfile)
            intra_clustered_sv_bedfile = os.path.join(merge_workdir, "%s_%s_clustered.bed" % (tool,sv_type))
            intra_clustered_sv_bed = pybedtools.BedTool(cluster_intervals_parallel(sv_bed, merge_workdir, wiggle=cluster_wiggle, overlap=overlap_ratio, 
                                        output_cluster_size=False, nthreads=1)).saveas(intra_clustered_sv_bedfile)
            logger.info("First level merging for %s for tool %s using wiggle %s and padding %s" % (sv_type, tool, cluster_wiggle, cluster_padding))
            tool_merged_intervals[sv_type] += merge_clustered_intervals(intra_clustered_sv_bed,intervals[tool][sv_type])

        # Do the inter-tool merging
        sv_bedfile = os.path.join(merge_workdir, "all_%s.bed" % (sv_type))
        sv_bed=pybedtools.BedTool(map(lambda x:x[1].to_cluster_bed_interval(sample,x[0],cluster_padding),
                                      enumerate(tool_merged_intervals[sv_type]))).sort().saveas(sv_bedfile)
        clustered_sv_bedfile = os.path.join(merge_workdir, "all_%s_clustered.bed" % (sv_type))
        clustered_sv_bed = pybedtools.BedTool(cluster_intervals_parallel(sv_bed, merge_workdir, wiggle=cluster_wiggle, overlap=overlap_ratio, 
                                    output_cluster_size=False, nthreads=1)).saveas(clustered_sv_bedfile)
        logger.info("Inter-tool Merging SVs of type %s using wiggle %s and padding %s" % (sv_type,cluster_wiggle, cluster_padding))
        final_intervals.extend(merge_clustered_intervals(clustered_sv_bed,tool_merged_intervals[sv_type]))

    final_chr_intervals = {contig.name: [] for contig in contigs}
    for interval in final_intervals:
        interval.do_validation(overlap_ratio)
        interval.fix_precise_coords()
        interval.fix_pos()
        if minsvlen <= abs(interval.length) <= maxsvlen or interval.sv_type in ["ITX", "CTX"] or (interval.length==0 and interval.sv_type == "INS"):
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

    return merged_bed, preasm_vcf    