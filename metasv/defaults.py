MIN_SV_LENGTH = 50
MAX_SV_LENGTH = 1000000
OVERLAP_RATIO = 0.5
WIGGLE = 100
INS_WIGGLE = 100
TX_WIGGLE = 5
SVS_SUPPORTED = set(["DEL", "DUP", "INS", "INV", "ITX", "CTX"])
SVS_ASSEMBLY_SUPPORTED = set(["DEL", "INS" , "INV", "DUP"])
SVS_SOFTCLIP_SUPPORTED = set(["DEL", "INS" , "INV", "DUP"])
MEAN_READ_LENGTH=100

# For generating candidate intervals for insertion assembly
MIN_SUPPORT_SC_ONLY = 2
MIN_SUPPORT_INS = 15
MIN_SUPPORT_DEL = 10 
MIN_SUPPORT_INV = 10 
MIN_SUPPORT_DUP = 10 
MIN_SUPPORT_FRAC_INS = 0.05
MIN_SUPPORT_FRAC_DEL = 0.04
MIN_SUPPORT_FRAC_INV = 0.015
MIN_SUPPORT_FRAC_DUP = 0.015
MEAN_READ_COVERAGE=50
MIN_INS_COVERAGE_FRAC=0.5
MAX_INS_COVERAGE_FRAC=1.5
MAX_INTERVALS = 10000
SC_PAD = 500
SC_MIN_SOFT_CLIP = 20
SC_MIN_AVG_BASE_QUAL = 20
SC_MIN_MAPQ = 5
SC_MAX_NM = 10
SC_MIN_MATCHES = 50


ISIZE_MIN = 250
ISIZE_MAX = 450
ISIZE_MEAN = 350.0
ISIZE_SD = 50.0

ASM_FULL = "run"
ASM_DEFER = "defer"
ASM_DISABLE = "disable"
ASM_SLICED = "parallel"
ASM_MERGE = "merge"
ASM_RUN_MODES = set([ASM_FULL, ASM_DEFER, ASM_DISABLE, ASM_SLICED, ASM_MERGE])

# For assembly read-extraction
EXTRACTION_MAX_READ_PAIRS = 10000
EXTRACTION_MAX_NM = 5
EXTRACTION_MAX_INTERVAL_TRUNCATION = 10000
EXTRACTION_TRUNCATION_PAD = 4000

# For running SPAdes
ASSEMBLY_MAX_TOOLS = 1
SPADES_TIMEOUT = 300  # in seconds
SPADES_PAD = 500
SPADES_MAX_INTERVAL_SIZE = 50000
SPADES_MAX_INTERVAL_SIZE_2BP = 1000000

# For running AGE
AGE_TIMEOUT = 300  # in seconds
AGE_MIN_CONTIG_LENGTH = 200
AGE_PAD = 500
AGE_MAX_REGION_LENGTH = 1000000
AGE_MAX_INTERVAL_TRUNCATION = 10000
AGE_TRUNCATION_PAD = 2000
AGE_DIST_TO_BP = 400
MIN_INV_SUBALIGN_LENGTH = 50
MIN_DEL_SUBALIGN_LENGTH = 50
AGE_WINDOW_SIZE = 20
AGE_DIST_TO_EXP_BP_INS = 25

# For genotyping
GT_WINDOW = 100
GT_NORMAL_FRAC = 0.05