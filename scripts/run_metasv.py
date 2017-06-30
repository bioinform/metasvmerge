#!/usr/bin/env python

import sys
import argparse
from metasv.main import run_metasv
from metasv.defaults import *
from metasv._version import __version__

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge SVs from multiple tools for accurate SV calling",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    input_parser = parser.add_argument_group("Input data options")
    input_parser.add_argument("--sample", metavar="Sample", help="Sample name", required=True)
    input_parser.add_argument("--pindel_vcf", nargs="+", metavar="pindel_vcf", help="VCF file or dir for Pindel VCFs",
                              default=[])
    input_parser.add_argument("--pindel_native", nargs="+", metavar="File list", help="Pindel native files",
                              default=[])
    input_parser.add_argument("--breakdancer_vcf", nargs="+", metavar="breakdancer_vcf",
                              help="VCF file or dir for BreakDancer VCFs", default=[])
    input_parser.add_argument("--breakdancer_native", nargs="+", metavar="File list", help="BreakDancer native files",
                              default=[])
    input_parser.add_argument("--breakseq_vcf", nargs="+", metavar="breakseq_vcf",
                              help="VCF file or dir for BreakSeq VCFs", default=[])
    input_parser.add_argument("--breakseq_native", nargs="+", metavar="breakseq_native",
                              help="BreakSeq native GFF files", default=[])
    input_parser.add_argument("--cnvnator_vcf", nargs="+", metavar="cnvnator_vcf",
                              help="VCF file or dir for CNVnator VCFs", default=[])
    input_parser.add_argument("--cnvnator_native", nargs="+", metavar="File list", help="CNVnator native files",
                              default=[])
    input_parser.add_argument("--softclip_native", nargs="+", metavar="File list", help="SoftClip native files",
                              default=[])
    input_parser.add_argument("--gatk_vcf", nargs="+", metavar="file", help="VCF file or dir for gatk VCFs",
                              default=[])
    input_parser.add_argument("--manta_vcf", nargs="+", help="VCF file or dir for Manta VCFs",
                              default=[])
    input_parser.add_argument("--lumpy_vcf", nargs="+", help="VCF file or dir for Lumpy VCFs",
                              default=[])
    input_parser.add_argument("--cnvkit_vcf", nargs="+", help="VCF file or dir for CNVkit VCFs",
                              default=[])
    input_parser.add_argument("--wham_vcf", nargs="+", help="VCF file or dir for WHAM VCFs",
                              default=[])
                              
    input_parser.add_argument("--mean_read_length", type=float, default=MEAN_READ_LENGTH, help="Mean read length")

    reference_parser = parser.add_argument_group("Reference options")
    reference_parser.add_argument("--reference", metavar="reference", help="Reference file", required=True)
    reference_parser.add_argument("--chromosomes",
                                  help="Chromosome list to process. If unspecified, then all chromosomes will be considered.",
                                  nargs="+", default=[])
    reference_parser.add_argument("--gaps", metavar="gaps", help="Gap bed file", default=None)
    reference_parser.add_argument("--filter_gaps", help="Filter out gaps", action="store_true")
    reference_parser.add_argument("--keep_standard_contigs", action="store_true",
                                  help="Keep only the major contigs + MT")

    bam_parser = parser.add_argument_group("Input BAM options")
    bam_parser.add_argument("--bams", nargs="+", help="BAMs", default=[])
    bam_parser.add_argument("--isize_mean", type=float, default=ISIZE_MEAN, help="Insert size mean")
    bam_parser.add_argument("--isize_sd", type=float, default=ISIZE_SD, help="Insert size standard deviation")

    merging_parser = parser.add_argument_group("Tool output merging options")
    merging_parser.add_argument("--wiggle", help="Wiggle for interval overlap", default=WIGGLE, type=int)
    merging_parser.add_argument("--inswiggle", help="Wiggle for insertions, overides wiggle", default=INS_WIGGLE,
                                type=int)
    merging_parser.add_argument("--minsvlen", help="Minimum length acceptable to be an SV", default=MIN_SV_LENGTH,
                                type=int)
    merging_parser.add_argument("--maxsvlen", help="Maximum length SV to report", default=MAX_SV_LENGTH,
                                type=int)
    merging_parser.add_argument("--overlap_ratio", help="Reciprocal overlap ratio", default=OVERLAP_RATIO, type=float)

    as_parser = parser.add_argument_group("Assembly options")
    as_parser.add_argument("--spades", help="Path to SPAdes executable")
    as_parser.add_argument("--spades_options", help="Options for SPAdes", default="")
    as_parser.add_argument("--spades_timeout", help="Maximum time (in seconds) for running SPAdes on an interval", default=SPADES_TIMEOUT, type=int)
    as_parser.add_argument("--disable_assembly", action="store_true", help="Disable assembly")
    as_parser.add_argument("--svs_to_assemble", nargs="+", help="SVs to assemble", default=["INS", "INV", "DUP"],
                           choices=SVS_ASSEMBLY_SUPPORTED)
    as_parser.add_argument("--extraction_max_read_pairs", type=int, default=EXTRACTION_MAX_READ_PAIRS,
                           help="Maximum number of pairs to extract for assembly")
    as_parser.add_argument("--spades_max_interval_size", type=int, default=SPADES_MAX_INTERVAL_SIZE,
                           help="Maximum SV length for assembly")
    as_parser.add_argument("--assembly_max_tools", type=int, default=ASSEMBLY_MAX_TOOLS,
                           help="Skip assembly if more than this many tools support a call (default 1)")
    as_parser.add_argument("--assembly_pad", type=int, default=SPADES_PAD,
                           help="Padding base pairs to use for assembling breakpoint with Spades and AGE")
    as_parser.add_argument("--stop_spades_on_fail", action="store_true", help="Abort on SPAdes failure")
    as_parser.add_argument("--age", help="Path to AGE executable")
    as_parser.add_argument("--age_timeout", help="Maximum time (in seconds) for running AGE on an assembled contig", default=AGE_TIMEOUT, type=int)
    as_parser.add_argument("--min_inv_subalign_len", help="Minimum length of inversion sub-alginment", type=int,
                        default=MIN_INV_SUBALIGN_LENGTH)
    as_parser.add_argument("--min_del_subalign_len", help="Minimum length of deletion sub-alginment", type=int,
                        default=MIN_DEL_SUBALIGN_LENGTH)
    as_parser.add_argument("--min_dup_subalign_len", help="Minimum length of duplication sub-alginment", type=int,
                        default=MIN_DUP_SUBALIGN_LENGTH)
    as_parser.add_argument("--age_window", help="Window size for AGE to merge nearby breakpoints", type=int,
                        default=AGE_WINDOW_SIZE)

    gt_parser = parser.add_argument_group("Genotyping options")
    gt_parser.add_argument("--gt_window", type=int, default=GT_WINDOW, help="Window for genotyping")
    gt_parser.add_argument("--gt_normal_frac", type=float, default=GT_NORMAL_FRAC,
                           help="Min. fraction of reads supporting reference for genotyping")

    out_parser = parser.add_argument_group("Output options")
    out_parser.add_argument("--svs_to_report", nargs="+", help="SV types to report", default=SVS_SUPPORTED,
                            choices=SVS_SUPPORTED)
    out_parser.add_argument("--enable_per_tool_output", action="store_true",
                            help="Enable output of merged SVs for individual tools")

    work_parser = parser.add_argument_group("Running environment options")
    work_parser.add_argument("--workdir", help="Scratch directory for working", default="work")
    work_parser.add_argument("--num_threads", help="Number of threads to use", type=int, default=1)
    work_parser.add_argument("--outdir", help="Output directory", required=True)
    work_parser.add_argument("--restart_step", help="Step from where to restart", choices=STEPS, default=STEP_LOAD)

    other_parser = parser.add_argument_group("Other options")
    other_parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    args.svs_to_assemble = set(args.svs_to_assemble) & set(args.svs_to_report)
    sys.exit(run_metasv(args))
