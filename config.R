# R/config.R

# ── Paths ─────────────────────────────────────────────────────────────────────
PATH_DATA_RAW       <- file.path("data", "raw")
PATH_DATA_DERIVED   <- file.path("data", "derived")
PATH_DATA_PROCESSED <- PATH_DATA_DERIVED        # RNA-seq alias

PATH_OUTPUT         <- "output"
PATH_OUTPUT_TABLES  <- PATH_OUTPUT
PATH_OUTPUT_FIGS    <- PATH_OUTPUT
PATH_OUTPUT_FIGURES <- PATH_OUTPUT
PATH_OUTPUT_REPORT  <- PATH_OUTPUT

# ── Analysis constants ────────────────────────────────────────────────────────
ZERO_CUTOFF       <- 0.30
PSEUDOCOUNT       <- 1e-6
ALPHA_Q           <- 0.05
TREND_ALPHA       <- 0.20
PBS_LEVEL         <- "PBS_Control"  # Changed from "PBS Control"
EXCLUDE_EXPOSURES <- c("Untreated_Control", "Peat_5")

ALL_HIGH_DOSE_EXPOSURES <- c("Pine_25", "Peat_25", "Eucalyptus_25", "RedOak_25")

EXPOSURE_COLORS_DEEP <- c(
  "Pine_25"       = "#1B4F8A",
  "Peat_25"       = "gold4",
  "Eucalyptus_25" = "darkorchid3",
  "RedOak_25"     = "deeppink2"
)
EXPOSURE_COLORS_LIGHT <- c(
  "Pine_25"       = "#D6E8F7",
  "Peat_25"       = "#FFF3B0",
  "Eucalyptus_25" = "#F0E6FA",
  "RedOak_25"     = "#FFE6F0"
)

# ── RNA-seq constants ─────────────────────────────────────────────────────────
ADJ_P_CUTOFF            <- 0.05
LOG2FC_CUTOFF           <- 1.0
FDR_THRESHOLD           <- 0.25
DREAM_N_CORES           <- 4
GENE_BACKGROUND_THRESHOLD <- 10

color_exposure <- c(
  "PBS_Control"      = "grey50",
  "Peat_5"           = "gold1",
  "Peat_25"          = "gold4",
  "Eucalyptus_5"     = "mediumpurple1",
  "Eucalyptus_25"    = "darkorchid3",
  "Pine_5"           = "steelblue1",
  "Pine_25"          = "#1B4F8A",
  "RedOak_5"         = "lightpink",
  "RedOak_25"        = "deeppink2",
  "Untreated_Control"= "forestgreen"
)

highlighted_genes <- c("SPDEF", "IL1A", "IL1B", "SRXN1",
                       "CYP1A1", "TXNRD1", "TNFA", "CXCL8", "SOX9", 
                       "HIC1", "PDIA2", "PTGER2", "GPR68", "SFRP2", "F2RL3")

priority_pathways <- c(
  "REACTIVE_OXYGEN_SPECIES_PATHWAY", "ANDROGEN_RESPONSE",
  "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE",
  "IL6_JAK_STAT3_SIGNALING", "TNFA_SIGNALING_VIA_NFKB",
  "TGF_BETA_SIGNALING", "OXIDATIVE_PHOSPHORYLATION",
  "G2M_CHECKPOINT", "GLYCOLYSIS", "XENOBIOTIC_METABOLISM",
  "WNT_BETA_CATENIN_SIGNALING", "PI3K_AKT_MTOR_SIGNALING",
  "CHOLESTEROL_HOMEOSTASIS", "NOTCH_SIGNALING"
)

# ── Factor levels ─────────────────────────────────────────────────────────────
SEX_LEVELS      <- c("M", "F")
CELLTYPE_LEVELS <- c("SAE", "LAE")  # These match your metadata
HORMONE_LEVELS  <- c("NONE", "Estradiol")

# ── Cytokine LLOD reference ───────────────────────────────────────────────────
cytokine_llod <- data.frame(
  Analyte = c(
    "EOTAXIN","EOTAXIN3","GMCSF","IFNY","IL1A","IL1B","IL2","IL4","IL5","IL6",
    "IL7","IL10","IL12P40","IL12P70","IL13","IL15","IL16","IL17","IP10","MCP1",
    "MCP4","MDC","MIP1A","MIP1B","TARC","TNFA","TNFB","VEGF","IL8"
  ),
  LLOD = c(
    0.20, 1.44, 0.16, 0.37, 0.09, 0.05, 0.09, 0.02, 0.14, 0.06,
    0.12, 0.07, 0.04, 0.33, 0.11, 0.24, 0.15, 2.83, 0.31, 0.12,
    0.16, 0.05, 1.00, 0.05, 0.09, 0.11, 0.04, 0.08, 1.12
  ),
  stringsAsFactors = FALSE
)

# Legacy LLOQ table (used in older ART chunks)
cytokine_lloq <- data.frame(
  Analyte = c(
    "EOTAXIN","EOTAXIN3","GMCSF","IFNY","IL1A","IL1B","IL2","IL4","IL5","IL6","IL7",
    "IL8","IL10","IL12P40","IL12P70","IL13","IL15","IL16","IL17","IP10","MCP1",
    "MCP4","MDC","MIP1A","MIP1B","TARC","TNFA","TNFB","VEGF"
  ),
  LLOQ = c(
    2.14, 11.1, 0.842, 1.76, 2.85, 0.646, 0.89, 0.218, 4.41, 0.633, 0.851,
    0.591, 0.298, 1.32,  1.22, 4.21, 0.774, 19.1, 3.19,  1.29, 1.45,  0.518,
    6.43, 0.357, 1.20,   2.14, 0.69, 0.465, 7.7
  ),
  stringsAsFactors = FALSE
)

cat("✓ Configuration loaded\n")
