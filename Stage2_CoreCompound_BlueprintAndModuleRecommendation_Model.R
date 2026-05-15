############################################################
## ACMPI-Nano v1.2
## Stage 2: Ideal Nanoformulation Blueprint Model
## PERFECT VERSION | TRUE MODEL
##
## Purpose:
##   Translate standardized compound-level physicochemical,
##   ADMET, druggability and toxicity features into an
##   interpretable ideal nanoformulation blueprint.
##
## Key upgrades in v1.2:
##   - Reads ONLY the original master input table.
##   - Does NOT read Stage 1 outputs.
##   - Separates carrier backbones from functional modules.
##   - Adds multivalent assembly drive.
##   - Adds cyclodextrin cavity-match features.
##   - Adds formulation-priority gate.
##   - Prevents low-need compounds from being over-interpreted.
##   - Adds class definition table, contribution table, QC table.
##   - Saves all tables and high-resolution PNG/PDF figures.
############################################################

rm(list = ls())

############################################################
## 0. User settings
############################################################

target_dir <- "/media/desk16/iy15915/中药之开创/Baicalein"
input_file <- file.path(target_dir, "Baicalein_Master_Model_Input.csv")

output_prefix <- "ACMPI_Nano_Stage2_v1p2_IdealNanoformulationBlueprint"

FIG_DPI <- 600
N_SENSITIVITY <- 1000
WEIGHT_JITTER <- 0.20
SET_SEED <- 20260427

set.seed(SET_SEED)

############################################################
## 1. Required packages
############################################################

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "openxlsx",
  "scales",
  "stringr",
  "grid"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(scales)
library(stringr)
library(grid)

############################################################
## 2. Output folders
############################################################

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_dir <- file.path(
  target_dir,
  paste0("ACMPI_Nano_Stage2_v1p2_Output_", timestamp)
)

table_dir <- file.path(out_dir, "01_Tables")
fig_png_dir <- file.path(out_dir, "02_Figures_PNG")
fig_pdf_dir <- file.path(out_dir, "03_Figures_PDF")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_pdf_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 3. Read input table
############################################################

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dat <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

if (nrow(dat) < 1) {
  stop("Input table is empty.")
}

############################################################
## 4. Required columns
############################################################

required_cols <- c(
  "compound", "PubChem_CID", "SMILES", "InChIKey",
  
  "Molecular_Weight", "MolLogP", "TPSA", "HBD", "HBA",
  "Rotatable_Bonds", "Ring_Count", "Aromatic_Rings",
  "Fraction_Csp3", "Formal_Charge", "Molar_Refractivity",
  
  "Consensus_LogP", "ESOL_LogS", "ESOL_Class",
  "ADMETlab_logS", "ADMETlab_logD",
  
  "GI_absorption", "Caco2", "PAMPA", "F30",
  
  "BBB_permeant", "ADMETlab_BBB", "Pgp_substrate",
  "pgp_sub", "pgp_inh",
  
  "PPB", "Fu", "logVDss",
  
  "CYP_inhibitor_count", "t_half", "cl_plasma",
  
  "Lipinski_violations", "Veber_violations", "Bioavailability_Score",
  "PAINS_alerts", "Brenk_alerts", "Synthetic_Accessibility", "QED",
  
  "LD50", "Toxicity_Class", "DILI", "Ames",
  "Carcinogenicity", "hERG",
  "ProTox_Organ_Active_Count",
  "ProTox_Endpoint_Active_Count",
  "ProTox_Metabolism_Active_Count"
)

missing_cols <- setdiff(required_cols, colnames(dat))

if (length(missing_cols) > 0) {
  stop(
    "The input table is missing required columns:\n",
    paste(missing_cols, collapse = "\n")
  )
}

dat <- dat[, required_cols]

############################################################
## 5. Helper functions
############################################################

to_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

clip01 <- function(x) {
  pmax(0, pmin(1, x))
}

risk_high <- function(x, good, bad) {
  clip01((to_num(x) - good) / (bad - good))
}

risk_low <- function(x, good, bad) {
  clip01((good - to_num(x)) / (good - bad))
}

benefit_high <- function(x, low, high) {
  clip01((to_num(x) - low) / (high - low))
}

benefit_low <- function(x, high, low) {
  clip01((high - to_num(x)) / (high - low))
}

benefit_mid <- function(x, lower, upper) {
  x <- to_num(x)
  span <- upper - lower
  
  out <- ifelse(
    x >= lower & x <= upper,
    1,
    ifelse(
      x < lower,
      clip01((x - (lower - span)) / span),
      clip01(((upper + span) - x) / span)
    )
  )
  
  clip01(out)
}

row_weighted_mean <- function(mat, weights) {
  mat <- as.matrix(mat)
  weights <- as.numeric(weights)
  
  apply(mat, 1, function(x) {
    keep <- !is.na(x) & !is.na(weights)
    if (!any(keep)) return(NA_real_)
    sum(x[keep] * weights[keep]) / sum(weights[keep])
  })
}

row_mean_safe <- function(mat) {
  mat <- as.matrix(mat)
  
  apply(mat, 1, function(x) {
    x <- as.numeric(x)
    if (all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
  })
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

class_score <- function(x, levels_map, default = 0.5) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(default, length(x))
  
  for (nm in names(levels_map)) {
    out[x == tolower(nm)] <- levels_map[[nm]]
  }
  
  out
}

clean_label <- function(x) {
  x <- gsub("_Score$", "", x)
  x <- gsub("_Requirement$", "", x)
  x <- gsub("_Fit$", "", x)
  x <- gsub("_", " ", x)
  x
}

wrap_label <- function(x, width = 24) {
  stringr::str_wrap(x, width = width)
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

normalize_weights <- function(w) {
  w <- as.numeric(w)
  w / sum(w)
}

jitter_weights <- function(w, jitter = 0.20) {
  w2 <- w * runif(length(w), min = 1 - jitter, max = 1 + jitter)
  normalize_weights(w2)
}

getv <- function(one_row, nm) {
  if (!nm %in% colnames(one_row)) {
    stop("Missing required internal feature: ", nm)
  }
  
  val <- suppressWarnings(as.numeric(one_row[[nm]][1]))
  
  if (length(val) == 0) {
    stop("Empty internal feature: ", nm)
  }
  
  val
}

weighted_sum_strict <- function(values, weights, carrier_name = "unknown") {
  values <- as.numeric(values)
  weights <- as.numeric(weights)
  
  if (length(values) != length(weights)) {
    stop("Length mismatch in scoring: ", carrier_name)
  }
  
  keep <- !is.na(values) & !is.na(weights)
  
  if (!any(keep)) {
    stop(
      "All scoring components are NA for: ",
      carrier_name,
      ". Please check QC output."
    )
  }
  
  sum(values[keep] * weights[keep]) / sum(weights[keep])
}

############################################################
## 6. Internal Stage 1 diagnostic scores
############################################################

score_solubility_defect <- function(df) {
  s_esol_logS <- risk_low(df$ESOL_LogS, good = -3, bad = -7)
  s_admet_logS <- risk_low(df$ADMETlab_logS, good = -3, bad = -7)
  s_cons_logP <- risk_high(df$Consensus_LogP, good = 2.5, bad = 6)
  s_logD <- risk_high(df$ADMETlab_logD, good = 2.5, bad = 6)
  
  s_class <- class_score(
    df$ESOL_Class,
    levels_map = c(
      "Highly soluble" = 0.00,
      "Very soluble" = 0.00,
      "Soluble" = 0.15,
      "Moderately soluble" = 0.40,
      "Poorly soluble" = 0.75,
      "Insoluble" = 1.00
    ),
    default = 0.50
  )
  
  row_weighted_mean(
    cbind(s_esol_logS, s_admet_logS, s_cons_logP, s_logD, s_class),
    weights = c(0.30, 0.25, 0.20, 0.10, 0.15)
  )
}

score_permeability_defect <- function(df) {
  s_tpsa <- risk_high(df$TPSA, good = 75, bad = 140)
  
  s_gi <- ifelse(
    tolower(trimws(df$GI_absorption)) == "high",
    0.10,
    0.85
  )
  
  s_caco2 <- risk_low(df$Caco2, good = -4.5, bad = -6.0)
  s_pampa <- risk_low(df$PAMPA, good = 0.70, bad = 0.20)
  s_f30 <- risk_low(df$F30, good = 0.80, bad = 0.20)
  
  row_weighted_mean(
    cbind(s_tpsa, s_gi, s_caco2, s_pampa, s_f30),
    weights = c(0.20, 0.20, 0.20, 0.20, 0.20)
  )
}

score_barrier_transporter <- function(df) {
  s_bbb_class <- ifelse(
    tolower(trimws(df$BBB_permeant)) == "yes",
    0.10,
    0.90
  )
  
  s_bbb_prob <- risk_low(df$ADMETlab_BBB, good = 0.70, bad = 0.10)
  
  s_pgp_class <- ifelse(
    tolower(trimws(df$Pgp_substrate)) == "yes",
    0.90,
    0.10
  )
  
  s_pgp_sub <- risk_high(df$pgp_sub, good = 0.20, bad = 0.80)
  s_pgp_inh <- risk_high(df$pgp_inh, good = 0.20, bad = 0.80)
  
  row_weighted_mean(
    cbind(s_bbb_class, s_bbb_prob, s_pgp_class, s_pgp_sub, s_pgp_inh),
    weights = c(0.30, 0.25, 0.20, 0.15, 0.10)
  )
}

score_distribution_exposure <- function(df) {
  s_ppb <- risk_high(df$PPB, good = 80, bad = 99)
  s_fu <- risk_low(df$Fu, good = 0.30, bad = 0.02)
  s_vd <- risk_high(abs(to_num(df$logVDss)), good = 0.8, bad = 2.0)
  
  row_weighted_mean(
    cbind(s_ppb, s_fu, s_vd),
    weights = c(0.45, 0.35, 0.20)
  )
}

score_metabolic_liability <- function(df) {
  s_cyp <- clip01(to_num(df$CYP_inhibitor_count) / 5)
  s_half <- risk_low(df$t_half, good = 4, bad = 0.5)
  s_clearance <- risk_high(df$cl_plasma, good = 5, bad = 15)
  
  row_weighted_mean(
    cbind(s_cyp, s_half, s_clearance),
    weights = c(0.45, 0.25, 0.30)
  )
}

score_nano_assembly_suitability <- function(df) {
  s_mw <- benefit_mid(df$Molecular_Weight, lower = 150, upper = 700)
  s_logp <- benefit_mid(df$MolLogP, lower = 1.5, upper = 5.0)
  s_hbd <- benefit_mid(df$HBD, lower = 1, upper = 5)
  s_hba <- benefit_mid(df$HBA, lower = 2, upper = 10)
  s_rot <- benefit_low(df$Rotatable_Bonds, high = 10, low = 0)
  s_ring <- benefit_high(df$Ring_Count, low = 1, high = 4)
  s_arom <- benefit_high(df$Aromatic_Rings, low = 1, high = 4)
  s_fsp3 <- benefit_low(df$Fraction_Csp3, high = 0.8, low = 0.0)
  s_charge <- ifelse(to_num(df$Formal_Charge) == 0, 1.00, 0.60)
  s_mr <- benefit_mid(df$Molar_Refractivity, lower = 40, upper = 130)
  
  row_weighted_mean(
    cbind(s_mw, s_logp, s_hbd, s_hba, s_rot, s_ring, s_arom, s_fsp3, s_charge, s_mr),
    weights = c(0.10, 0.15, 0.10, 0.10, 0.10, 0.10, 0.15, 0.10, 0.05, 0.05)
  )
}

score_druggability_alert_requirement <- function(df) {
  s_lip <- clip01(to_num(df$Lipinski_violations) / 2)
  s_veb <- clip01(to_num(df$Veber_violations) / 1)
  s_bio <- 1 - clip01(to_num(df$Bioavailability_Score))
  s_pains <- clip01(to_num(df$PAINS_alerts) / 2)
  s_brenk <- clip01(to_num(df$Brenk_alerts) / 3)
  s_sa <- risk_high(df$Synthetic_Accessibility, good = 3, bad = 8)
  s_qed <- risk_low(df$QED, good = 0.70, bad = 0.20)
  
  row_weighted_mean(
    cbind(s_lip, s_veb, s_bio, s_pains, s_brenk, s_sa, s_qed),
    weights = c(0.15, 0.10, 0.15, 0.20, 0.15, 0.10, 0.15)
  )
}

score_safety_control_requirement <- function(df) {
  s_ld50 <- risk_low(df$LD50, good = 5000, bad = 50)
  s_class <- clip01((6 - to_num(df$Toxicity_Class)) / 5)
  s_dili <- clip01(to_num(df$DILI))
  s_ames <- clip01(to_num(df$Ames))
  s_carc <- clip01(to_num(df$Carcinogenicity))
  s_herg <- clip01(to_num(df$hERG))
  s_org <- clip01(to_num(df$ProTox_Organ_Active_Count) / 5)
  s_end <- clip01(to_num(df$ProTox_Endpoint_Active_Count) / 8)
  s_met <- clip01(to_num(df$ProTox_Metabolism_Active_Count) / 6)
  
  row_weighted_mean(
    cbind(s_ld50, s_class, s_dili, s_ames, s_carc, s_herg, s_org, s_end, s_met),
    weights = c(0.10, 0.10, 0.15, 0.12, 0.12, 0.12, 0.12, 0.10, 0.07)
  )
}

add_stage1_internal_scores <- function(df) {
  out <- df
  
  out$Solubility_Defect_Score <- score_solubility_defect(out)
  out$Permeability_Defect_Score <- score_permeability_defect(out)
  out$Barrier_Transporter_Score <- score_barrier_transporter(out)
  out$Distribution_Exposure_Score <- score_distribution_exposure(out)
  out$Metabolic_Liability_Score <- score_metabolic_liability(out)
  out$Nano_Assembly_Suitability_Score <- score_nano_assembly_suitability(out)
  out$Druggability_Alert_Requirement <- score_druggability_alert_requirement(out)
  out$Safety_Control_Requirement <- score_safety_control_requirement(out)
  
  out$Nano_Optimization_Need <- 
    0.22 * out$Solubility_Defect_Score +
    0.18 * out$Permeability_Defect_Score +
    0.15 * out$Barrier_Transporter_Score +
    0.10 * out$Distribution_Exposure_Score +
    0.15 * out$Metabolic_Liability_Score +
    0.20 * out$Safety_Control_Requirement
  
  out$Delivery_Defect_Burden <- 
    0.25 * out$Solubility_Defect_Score +
    0.25 * out$Permeability_Defect_Score +
    0.20 * out$Barrier_Transporter_Score +
    0.15 * out$Distribution_Exposure_Score +
    0.15 * out$Metabolic_Liability_Score
  
  out$Nanoformulation_Value_Index <- clip01(
    0.40 * out$Nano_Optimization_Need +
      0.45 * out$Nano_Assembly_Suitability_Score +
      0.15 * (1 - out$Druggability_Alert_Requirement)
  )
  
  out$Nanoformulation_Development_Value <- dplyr::case_when(
    out$Nanoformulation_Value_Index >= 0.70 ~ "High NDV",
    out$Nanoformulation_Value_Index >= 0.45 ~ "Moderate NDV",
    TRUE ~ "Low NDV"
  )
  
  out
}

############################################################
## 7. Stage 2 molecular design features
############################################################

add_stage2_design_features <- function(df) {
  out <- df
  
  out$Molecular_Size_Fit <- benefit_mid(out$Molecular_Weight, lower = 100, upper = 700)
  out$Balanced_LogP <- benefit_mid(out$MolLogP, lower = 1.5, upper = 5.0)
  out$Balanced_Consensus_LogP <- benefit_mid(out$Consensus_LogP, lower = 1.5, upper = 5.0)
  out$Lipophilicity_Fit <- benefit_mid(out$Consensus_LogP, lower = 2.0, upper = 6.0)
  
  out$Hydrophobic_Core_Fit <- row_mean_safe(
    cbind(
      benefit_mid(out$MolLogP, lower = 2.0, upper = 6.0),
      benefit_mid(out$Consensus_LogP, lower = 2.0, upper = 6.0),
      risk_low(out$ESOL_LogS, good = -3, bad = -7)
    )
  )
  
  out$Aromaticity <- benefit_high(out$Aromatic_Rings, low = 1, high = 4)
  out$Ring_Structure <- benefit_high(out$Ring_Count, low = 1, high = 4)
  out$Rigidity <- benefit_low(out$Rotatable_Bonds, high = 10, low = 0)
  
  out$Hydrogen_Bonding <- row_mean_safe(
    cbind(
      benefit_mid(out$HBD, lower = 1, upper = 5),
      benefit_mid(out$HBA, lower = 2, upper = 10)
    )
  )
  
  out$Molar_Refractivity_Fit <- benefit_mid(out$Molar_Refractivity, lower = 40, upper = 130)
  out$Charge_Compatibility <- ifelse(to_num(out$Formal_Charge) == 0, 1.00, 0.60)
  out$Druggability_Manageability <- 1 - out$Druggability_Alert_Requirement
  
  out$Controlled_Release_Need <- row_mean_safe(
    cbind(
      out$Metabolic_Liability_Score,
      out$Safety_Control_Requirement
    )
  )
  
  ##########################################################
  ## v1.2: refined formulation-discriminating features
  ##########################################################
  
  out$Multivalent_Assembly_Drive <- row_mean_safe(
    cbind(
      benefit_high(out$Aromatic_Rings, low = 1.5, high = 3.0),
      benefit_mid(out$HBD, lower = 2, upper = 5),
      benefit_mid(out$HBA, lower = 4, upper = 10),
      benefit_low(out$Rotatable_Bonds, high = 6, low = 0),
      benefit_high(out$Ring_Count, low = 2, high = 4),
      benefit_mid(out$Molar_Refractivity, lower = 65, upper = 130)
    )
  )
  
  out$Cyclodextrin_Cavity_Match <- row_mean_safe(
    cbind(
      benefit_mid(out$Molecular_Weight, lower = 120, upper = 420),
      benefit_mid(out$Aromatic_Rings, lower = 1, upper = 2),
      benefit_mid(out$MolLogP, lower = 1.0, upper = 4.0),
      benefit_low(out$Rotatable_Bonds, high = 6, low = 0),
      benefit_mid(out$Molar_Refractivity, lower = 45, upper = 95)
    )
  )
  
  out$Small_Aromatic_Inclusion_Preference <- row_mean_safe(
    cbind(
      benefit_mid(out$Molecular_Weight, lower = 120, upper = 300),
      benefit_mid(out$Aromatic_Rings, lower = 1, upper = 2),
      benefit_mid(out$Ring_Count, lower = 1, upper = 3),
      benefit_mid(out$MolLogP, lower = 1.0, upper = 3.5),
      benefit_low(out$Rotatable_Bonds, high = 5, low = 0)
    )
  )
  
  ##########################################################
  ## v1.2: formulation-priority gate
  ##########################################################
  
  out$Formulation_Priority_Index <- clip01(
    0.45 * out$Nano_Optimization_Need +
      0.25 * out$Safety_Control_Requirement +
      0.20 * out$Metabolic_Liability_Score +
      0.10 * out$Barrier_Transporter_Score
  )
  
  out$Formulation_Priority_Class <- dplyr::case_when(
    out$Formulation_Priority_Index >= 0.66 ~ "High formulation priority",
    out$Formulation_Priority_Index >= 0.40 ~ "Moderate formulation priority",
    TRUE ~ "Low formulation priority"
  )
  
  out
}

############################################################
## 8. Predefined class definitions
############################################################

make_class_definition_table <- function() {
  data.frame(
    Class_Type = c(
      rep("Carrier backbone", 5),
      rep("Functional module", 4)
    ),
    Model_Defined_Class = c(
      "Self-assembled nanoformulation",
      "Cyclodextrin inclusion system",
      "Polymeric micelle",
      "Lipid nanoparticle",
      "Polymeric controlled-release nanoparticle",
      "Barrier-aware ligand module",
      "Biomimetic surface module",
      "Exposure-reducing stabilization module",
      "Stimuli-responsive release module"
    ),
    Design_Role = c(
      "Primary carrier backbone based on multivalent non-covalent molecular assembly",
      "Inclusion/stabilization backbone for small-to-medium aromatic guest compounds",
      "Amphiphilic carrier backbone for hydrophobic-core loading",
      "Lipid-compatible carrier backbone for lipophilic partitioning and exposure modulation",
      "Polymeric matrix backbone for controlled release and exposure smoothing",
      "Conditional targeting module for barrier or transporter-related delivery limitations",
      "Conditional surface module for biomimetic delivery, immune shielding, or tissue selectivity",
      "Conditional stabilization module for reducing non-target exposure",
      "Conditional release module for metabolic or microenvironment-responsive control"
    ),
    Main_Input_Basis = c(
      "Multivalent assembly drive; hydrogen bonding; rigidity; aromaticity; balanced LogP; ring structure; molecular size; molar refractivity; charge compatibility",
      "Cyclodextrin cavity match; small aromatic inclusion preference; molecular size; balanced LogP; aromaticity; solubility defect; hydrogen bonding",
      "Hydrophobic-core fit; solubility defect; aromaticity; molecular size; hydrogen bonding; safety-control requirement",
      "Lipophilicity fit; hydrophobic-core fit; solubility defect; molecular size; distribution/exposure score; safety-control requirement",
      "Metabolic liability; safety-control requirement; controlled-release need; distribution/exposure score; molecular size; druggability manageability",
      "Barrier/transporter score; safety-control requirement; nano-optimization need; distribution/exposure score",
      "Barrier/transporter score; safety-control requirement; distribution/exposure score; metabolic liability",
      "Safety-control requirement; distribution/exposure score; metabolic liability; druggability manageability",
      "Metabolic liability; safety-control requirement; nano-optimization need; solubility defect"
    ),
    Interpretation_Note = c(
      "A high score indicates that the compound has multivalent structural features compatible with self-assembly or co-assembly, not merely a small aromatic scaffold.",
      "A high score indicates that the compound is structurally compatible with cyclodextrin-style cavity inclusion or stabilization.",
      "A high score indicates compatibility with micellar hydrophobic-core loading and formulation stabilization.",
      "A high score indicates compatibility with lipid-phase partitioning and lipid-based exposure modulation.",
      "A high score indicates that controlled-release polymeric delivery is a suitable design route.",
      "This is not a competing carrier backbone; it is a conditional surface or targeting module.",
      "This is not a competing carrier backbone; it is a conditional surface functionalization strategy.",
      "This is not a competing carrier backbone; it is a conditional safety/exposure control module.",
      "This is not a competing carrier backbone; it is a conditional release-control module."
    ),
    stringsAsFactors = FALSE
  )
}

class_definition_table <- make_class_definition_table()

############################################################
## 9. Carrier backbone scoring
############################################################

carrier_group_scores_one <- function(one_row) {
  
  Molecular_Size_Fit <- getv(one_row, "Molecular_Size_Fit")
  Balanced_LogP <- getv(one_row, "Balanced_LogP")
  Lipophilicity_Fit <- getv(one_row, "Lipophilicity_Fit")
  Hydrophobic_Core_Fit <- getv(one_row, "Hydrophobic_Core_Fit")
  Aromaticity <- getv(one_row, "Aromaticity")
  Ring_Structure <- getv(one_row, "Ring_Structure")
  Rigidity <- getv(one_row, "Rigidity")
  Hydrogen_Bonding <- getv(one_row, "Hydrogen_Bonding")
  Molar_Refractivity_Fit <- getv(one_row, "Molar_Refractivity_Fit")
  Charge_Compatibility <- getv(one_row, "Charge_Compatibility")
  Druggability_Manageability <- getv(one_row, "Druggability_Manageability")
  Controlled_Release_Need <- getv(one_row, "Controlled_Release_Need")
  
  Solubility_Defect_Score <- getv(one_row, "Solubility_Defect_Score")
  Distribution_Exposure_Score <- getv(one_row, "Distribution_Exposure_Score")
  Metabolic_Liability_Score <- getv(one_row, "Metabolic_Liability_Score")
  Safety_Control_Requirement <- getv(one_row, "Safety_Control_Requirement")
  
  Multivalent_Assembly_Drive <- getv(one_row, "Multivalent_Assembly_Drive")
  Cyclodextrin_Cavity_Match <- getv(one_row, "Cyclodextrin_Cavity_Match")
  Small_Aromatic_Inclusion_Preference <- getv(one_row, "Small_Aromatic_Inclusion_Preference")
  
  list(
    Self_Assembled_Nanoformulation = c(
      Primary = safe_mean(c(
        Multivalent_Assembly_Drive,
        Hydrogen_Bonding,
        Rigidity
      )),
      Secondary = safe_mean(c(
        Aromaticity,
        Balanced_LogP,
        Ring_Structure
      )),
      Modifier = safe_mean(c(
        Molecular_Size_Fit,
        Molar_Refractivity_Fit,
        Charge_Compatibility
      ))
    ),
    
    Cyclodextrin_Inclusion_System = c(
      Primary = safe_mean(c(
        Cyclodextrin_Cavity_Match,
        Small_Aromatic_Inclusion_Preference
      )),
      Secondary = safe_mean(c(
        Molecular_Size_Fit,
        Balanced_LogP,
        Aromaticity
      )),
      Modifier = safe_mean(c(
        Solubility_Defect_Score,
        Hydrogen_Bonding
      ))
    ),
    
    Polymeric_Micelle = c(
      Primary = safe_mean(c(
        Hydrophobic_Core_Fit,
        Solubility_Defect_Score
      )),
      Secondary = safe_mean(c(
        Aromaticity,
        Molecular_Size_Fit
      )),
      Modifier = safe_mean(c(
        Hydrogen_Bonding,
        Safety_Control_Requirement
      ))
    ),
    
    Lipid_Nanoparticle = c(
      Primary = safe_mean(c(
        Lipophilicity_Fit,
        Hydrophobic_Core_Fit
      )),
      Secondary = safe_mean(c(
        Solubility_Defect_Score,
        Molecular_Size_Fit
      )),
      Modifier = safe_mean(c(
        Distribution_Exposure_Score,
        Safety_Control_Requirement
      ))
    ),
    
    Polymeric_Controlled_Release_Nanoparticle = c(
      Primary = safe_mean(c(
        Metabolic_Liability_Score,
        Safety_Control_Requirement
      )),
      Secondary = safe_mean(c(
        Controlled_Release_Need,
        Distribution_Exposure_Score
      )),
      Modifier = safe_mean(c(
        Molecular_Size_Fit,
        Druggability_Manageability
      ))
    )
  )
}

carrier_base_weights <- list(
  Self_Assembled_Nanoformulation = c(Primary = 0.50, Secondary = 0.30, Modifier = 0.20),
  Cyclodextrin_Inclusion_System = c(Primary = 0.50, Secondary = 0.30, Modifier = 0.20),
  Polymeric_Micelle = c(Primary = 0.50, Secondary = 0.30, Modifier = 0.20),
  Lipid_Nanoparticle = c(Primary = 0.50, Secondary = 0.25, Modifier = 0.25),
  Polymeric_Controlled_Release_Nanoparticle = c(Primary = 0.45, Secondary = 0.30, Modifier = 0.25)
)

compute_carrier_scores_one <- function(one_row, jitter = NULL) {
  
  group_scores <- carrier_group_scores_one(one_row)
  
  out <- data.frame(
    carrier = character(0),
    score = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (carrier_name in names(group_scores)) {
    
    g <- group_scores[[carrier_name]]
    w <- carrier_base_weights[[carrier_name]]
    
    if (!all(names(w) %in% names(g))) {
      stop("Weight names do not match group score names for carrier: ", carrier_name)
    }
    
    g <- g[names(w)]
    
    if (!is.null(jitter)) {
      w <- jitter_weights(w, jitter = jitter)
    } else {
      w <- normalize_weights(w)
    }
    
    score <- weighted_sum_strict(
      values = g,
      weights = w,
      carrier_name = carrier_name
    )
    
    out <- rbind(
      out,
      data.frame(
        carrier = carrier_name,
        score = clip01(score),
        stringsAsFactors = FALSE
      )
    )
  }
  
  out
}

compute_carrier_scores_all <- function(df) {
  
  res <- list()
  
  for (i in seq_len(nrow(df))) {
    one_row <- df[i, , drop = FALSE]
    tmp <- compute_carrier_scores_one(one_row)
    tmp$compound <- one_row$compound[1]
    tmp$PubChem_CID <- one_row$PubChem_CID[1]
    tmp$InChIKey <- one_row$InChIKey[1]
    res[[i]] <- tmp
  }
  
  bind_rows(res) %>%
    dplyr::select(compound, PubChem_CID, InChIKey, carrier, score) %>%
    dplyr::group_by(compound) %>%
    dplyr::mutate(
      rank = rank(-score, ties.method = "first"),
      carrier_label = clean_label(carrier)
    ) %>%
    dplyr::ungroup()
}

make_carrier_feature_contribution <- function(df) {
  
  res <- list()
  
  for (i in seq_len(nrow(df))) {
    
    one_row <- df[i, , drop = FALSE]
    gs <- carrier_group_scores_one(one_row)
    
    tmp <- bind_rows(lapply(names(gs), function(nm) {
      data.frame(
        compound = one_row$compound[1],
        carrier = nm,
        carrier_label = clean_label(nm),
        group = names(gs[[nm]]),
        group_score = as.numeric(gs[[nm]]),
        group_weight = as.numeric(carrier_base_weights[[nm]][names(gs[[nm]])]),
        stringsAsFactors = FALSE
      )
    }))
    
    tmp$weighted_contribution <- tmp$group_score * tmp$group_weight
    tmp$final_score_from_groups <- ave(
      tmp$weighted_contribution,
      tmp$compound,
      tmp$carrier,
      FUN = sum
    )
    
    res[[i]] <- tmp
  }
  
  bind_rows(res)
}

############################################################
## 10. Functional module scoring
############################################################

functional_module_scores_one <- function(one_row) {
  
  Barrier_Transporter_Score <- getv(one_row, "Barrier_Transporter_Score")
  Safety_Control_Requirement <- getv(one_row, "Safety_Control_Requirement")
  Nano_Optimization_Need <- getv(one_row, "Nano_Optimization_Need")
  Distribution_Exposure_Score <- getv(one_row, "Distribution_Exposure_Score")
  Metabolic_Liability_Score <- getv(one_row, "Metabolic_Liability_Score")
  Solubility_Defect_Score <- getv(one_row, "Solubility_Defect_Score")
  Druggability_Manageability <- getv(one_row, "Druggability_Manageability")
  
  modules <- list(
    Barrier_Aware_Ligand_Module = c(
      Barrier = Barrier_Transporter_Score,
      Safety = Safety_Control_Requirement,
      Optimization = Nano_Optimization_Need,
      Distribution = Distribution_Exposure_Score
    ),
    
    Biomimetic_Surface_Module = c(
      Barrier = Barrier_Transporter_Score,
      Safety = Safety_Control_Requirement,
      Distribution = Distribution_Exposure_Score,
      Metabolic = Metabolic_Liability_Score
    ),
    
    Exposure_Reducing_Stabilization_Module = c(
      Safety = Safety_Control_Requirement,
      Distribution = Distribution_Exposure_Score,
      Metabolic = Metabolic_Liability_Score,
      Manageability = Druggability_Manageability
    ),
    
    Stimuli_Responsive_Release_Module = c(
      Metabolic = Metabolic_Liability_Score,
      Safety = Safety_Control_Requirement,
      Optimization = Nano_Optimization_Need,
      Solubility = Solubility_Defect_Score
    )
  )
  
  weights <- list(
    Barrier_Aware_Ligand_Module = c(Barrier = 0.50, Safety = 0.25, Optimization = 0.15, Distribution = 0.10),
    Biomimetic_Surface_Module = c(Barrier = 0.35, Safety = 0.25, Distribution = 0.20, Metabolic = 0.20),
    Exposure_Reducing_Stabilization_Module = c(Safety = 0.45, Distribution = 0.25, Metabolic = 0.20, Manageability = 0.10),
    Stimuli_Responsive_Release_Module = c(Metabolic = 0.40, Safety = 0.25, Optimization = 0.20, Solubility = 0.15)
  )
  
  out <- data.frame(
    module = character(0),
    score = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (module_name in names(modules)) {
    g <- modules[[module_name]]
    w <- weights[[module_name]]
    g <- g[names(w)]
    
    score <- weighted_sum_strict(
      values = g,
      weights = normalize_weights(w),
      carrier_name = module_name
    )
    
    out <- rbind(
      out,
      data.frame(
        module = module_name,
        score = clip01(score),
        stringsAsFactors = FALSE
      )
    )
  }
  
  out
}

compute_functional_module_scores_all <- function(df) {
  
  res <- list()
  
  for (i in seq_len(nrow(df))) {
    one_row <- df[i, , drop = FALSE]
    tmp <- functional_module_scores_one(one_row)
    tmp$compound <- one_row$compound[1]
    tmp$PubChem_CID <- one_row$PubChem_CID[1]
    tmp$InChIKey <- one_row$InChIKey[1]
    res[[i]] <- tmp
  }
  
  bind_rows(res) %>%
    dplyr::select(compound, PubChem_CID, InChIKey, module, score) %>%
    dplyr::group_by(compound) %>%
    dplyr::mutate(
      rank = rank(-score, ties.method = "first"),
      module_label = clean_label(module),
      recommendation = dplyr::case_when(
        score >= 0.66 ~ "Strongly recommended",
        score >= 0.45 ~ "Recommended",
        score >= 0.30 ~ "Optional",
        TRUE ~ "Low priority"
      )
    ) %>%
    dplyr::ungroup()
}

############################################################
## 11. QC table
############################################################

make_carrier_feature_qc <- function(df) {
  
  qc_numeric_cols <- c(
    "Molecular_Size_Fit",
    "Balanced_LogP",
    "Balanced_Consensus_LogP",
    "Lipophilicity_Fit",
    "Hydrophobic_Core_Fit",
    "Aromaticity",
    "Ring_Structure",
    "Rigidity",
    "Hydrogen_Bonding",
    "Molar_Refractivity_Fit",
    "Charge_Compatibility",
    "Druggability_Manageability",
    "Controlled_Release_Need",
    "Multivalent_Assembly_Drive",
    "Cyclodextrin_Cavity_Match",
    "Small_Aromatic_Inclusion_Preference",
    "Solubility_Defect_Score",
    "Permeability_Defect_Score",
    "Barrier_Transporter_Score",
    "Distribution_Exposure_Score",
    "Metabolic_Liability_Score",
    "Safety_Control_Requirement",
    "Nano_Optimization_Need",
    "Formulation_Priority_Index"
  )
  
  qc_class_cols <- c(
    "Formulation_Priority_Class"
  )
  
  qc_cols <- c(qc_numeric_cols, qc_class_cols)
  missing_internal <- setdiff(qc_cols, colnames(df))
  
  if (length(missing_internal) > 0) {
    stop(
      "Missing internal Stage 2 features:\n",
      paste(missing_internal, collapse = "\n")
    )
  }
  
  numeric_qc <- df %>%
    dplyr::select(compound, dplyr::all_of(qc_numeric_cols)) %>%
    tidyr::pivot_longer(
      cols = -compound,
      names_to = "feature",
      values_to = "value_numeric"
    ) %>%
    dplyr::mutate(
      value_raw = as.character(value_numeric),
      value_type = "numeric",
      is_missing = is.na(value_numeric),
      feature_label = clean_label(feature)
    )
  
  class_qc <- df %>%
    dplyr::select(compound, dplyr::all_of(qc_class_cols)) %>%
    tidyr::pivot_longer(
      cols = -compound,
      names_to = "feature",
      values_to = "value_raw"
    ) %>%
    dplyr::mutate(
      value_numeric = NA_real_,
      value_type = "character",
      is_missing = is.na(value_raw) | value_raw == "",
      feature_label = clean_label(feature)
    )
  
  dplyr::bind_rows(numeric_qc, class_qc) %>%
    dplyr::select(
      compound,
      feature,
      feature_label,
      value_type,
      value_raw,
      value_numeric,
      is_missing
    )
}

############################################################
## 12. Weight sensitivity analysis
############################################################

run_weight_sensitivity_one <- function(one_row, n_iter = 1000, jitter = 0.20) {
  
  baseline <- compute_carrier_scores_one(one_row, jitter = NULL)
  
  all_iter <- vector("list", n_iter)
  
  for (b in seq_len(n_iter)) {
    tmp <- compute_carrier_scores_one(one_row, jitter = jitter)
    tmp$iteration <- b
    all_iter[[b]] <- tmp
  }
  
  all_iter <- bind_rows(all_iter)
  
  top_tab <- all_iter %>%
    dplyr::group_by(iteration) %>%
    dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::count(carrier, name = "top_n") %>%
    dplyr::mutate(
      top_rank_probability = top_n / n_iter
    )
  
  score_tab <- all_iter %>%
    dplyr::group_by(carrier) %>%
    dplyr::summarise(
      sensitivity_mean = mean(score, na.rm = TRUE),
      sensitivity_sd = sd(score, na.rm = TRUE),
      sensitivity_p05 = quantile(score, 0.05, na.rm = TRUE),
      sensitivity_p95 = quantile(score, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
  
  baseline %>%
    dplyr::rename(baseline_score = score) %>%
    dplyr::left_join(score_tab, by = "carrier") %>%
    dplyr::left_join(
      top_tab[, c("carrier", "top_rank_probability")],
      by = "carrier"
    ) %>%
    dplyr::mutate(
      top_rank_probability = ifelse(is.na(top_rank_probability), 0, top_rank_probability),
      carrier_label = clean_label(carrier)
    )
}

run_weight_sensitivity_all <- function(df, n_iter = 1000, jitter = 0.20) {
  
  res <- list()
  
  for (i in seq_len(nrow(df))) {
    one_row <- df[i, , drop = FALSE]
    
    tmp <- run_weight_sensitivity_one(
      one_row,
      n_iter = n_iter,
      jitter = jitter
    )
    
    tmp$compound <- one_row$compound[1]
    tmp$PubChem_CID <- one_row$PubChem_CID[1]
    tmp$InChIKey <- one_row$InChIKey[1]
    
    res[[i]] <- tmp
  }
  
  bind_rows(res) %>%
    dplyr::select(
      compound, PubChem_CID, InChIKey,
      carrier, carrier_label,
      baseline_score,
      sensitivity_mean,
      sensitivity_sd,
      sensitivity_p05,
      sensitivity_p95,
      top_rank_probability
    )
}

############################################################
## 13. Blueprint logic
############################################################

make_blueprint_one <- function(one_row, carrier_scores_one, module_scores_one) {
  
  ranked_carrier <- carrier_scores_one %>%
    dplyr::arrange(dplyr::desc(score))
  
  ranked_module <- module_scores_one %>%
    dplyr::arrange(dplyr::desc(score))
  
  primary <- ranked_carrier$carrier[1]
  secondary <- ranked_carrier$carrier[2]
  top_module <- ranked_module$module[1]
  top_carrier_score <- ranked_carrier$score[1]
  
  barrier_high <- one_row$Barrier_Transporter_Score >= 0.50
  metabolic_high <- one_row$Metabolic_Liability_Score >= 0.50
  safety_high <- one_row$Safety_Control_Requirement >= 0.50
  assembly_high <- one_row$Nano_Assembly_Suitability_Score >= 0.66
  sol_high <- one_row$Solubility_Defect_Score >= 0.50
  perm_high <- one_row$Permeability_Defect_Score >= 0.50
  
  formulation_priority_index <- as.numeric(one_row$Formulation_Priority_Index[1])
  formulation_priority_class <- as.character(one_row$Formulation_Priority_Class[1])
  
  ideal_size <- dplyr::case_when(
    barrier_high ~ "50-120 nm",
    primary %in% c("Self_Assembled_Nanoformulation", "Polymeric_Micelle", "Cyclodextrin_Inclusion_System") ~ "60-150 nm",
    primary %in% c("Lipid_Nanoparticle") ~ "70-160 nm",
    primary %in% c("Polymeric_Controlled_Release_Nanoparticle") ~ "80-200 nm",
    TRUE ~ "60-180 nm"
  )
  
  ideal_pdi <- dplyr::case_when(
    assembly_high ~ "<=0.20 preferred; <=0.25 acceptable",
    TRUE ~ "<=0.25 acceptable"
  )
  
  zeta <- dplyr::case_when(
    barrier_high ~ "near-neutral to mildly negative (-5 to -20 mV)",
    safety_high ~ "near-neutral or mildly negative (-5 to -25 mV)",
    TRUE ~ "mildly negative (-10 to -30 mV)"
  )
  
  release_strategy <- dplyr::case_when(
    metabolic_high & safety_high ~ "sustained release with optional pH/ROS-responsive release",
    metabolic_high ~ "sustained release to reduce peak exposure and metabolic burden",
    safety_high ~ "controlled release to reduce non-target systemic exposure",
    barrier_high ~ "stable loading with target-site responsive release",
    TRUE ~ "stable encapsulation with moderate sustained release"
  )
  
  encapsulation_driver <- dplyr::case_when(
    primary == "Self_Assembled_Nanoformulation" ~
      "multivalent hydrogen-bonding and pi-pi stacking assisted self-assembly",
    primary == "Cyclodextrin_Inclusion_System" ~
      "cyclodextrin cavity inclusion and solubility/stability enhancement",
    primary == "Polymeric_Micelle" ~
      "hydrophobic-core loading with hydrogen-bond-assisted stabilization",
    primary == "Lipid_Nanoparticle" ~
      "lipophilic partitioning into lipid phase",
    primary == "Polymeric_Controlled_Release_Nanoparticle" ~
      "polymeric matrix entrapment and controlled release",
    TRUE ~
      "mixed non-covalent loading"
  )
  
  co_loading_need <- dplyr::case_when(
    one_row$Nano_Optimization_Need >= 0.66 ~ "High",
    barrier_high | metabolic_high | safety_high ~ "Moderate",
    TRUE ~ "Low"
  )
  
  final_interpretation <- dplyr::case_when(
    formulation_priority_index < 0.40 & top_carrier_score >= 0.60 ~
      "Structurally compatible but low-priority nanoformulation need",
    
    formulation_priority_index >= 0.66 & top_carrier_score >= 0.60 ~
      "High-priority nanoformulation development candidate",
    
    formulation_priority_index >= 0.40 & top_carrier_score >= 0.60 ~
      "Recommended nanoformulation development candidate",
    
    formulation_priority_index < 0.40 & top_carrier_score < 0.60 ~
      "Limited nanoformulation development priority",
    
    TRUE ~
      "Conditional nanoformulation development candidate"
  )
  
  efficacy_logic <- paste(
    c(
      if (assembly_high) "use intrinsic assembly features to improve loading and structural stability",
      if (sol_high) "improve apparent solubility",
      if (perm_high) "support absorption or membrane transport",
      if (barrier_high) "improve tissue or barrier-aware delivery",
      "increase local effective exposure through nanoscale delivery"
    ),
    collapse = "; "
  )
  
  toxicity_logic <- paste(
    c(
      if (metabolic_high) "use sustained release to reduce peak concentration and metabolic stress",
      if (safety_high) "reduce non-target systemic exposure",
      if (barrier_high) "use targeting or biomimetic surface design to improve delivery selectivity",
      "avoid unnecessary exposure escalation"
    ),
    collapse = "; "
  )
  
  data.frame(
    compound = one_row$compound[1],
    PubChem_CID = one_row$PubChem_CID[1],
    InChIKey = one_row$InChIKey[1],
    Primary_Nanoformulation_Backbone = clean_label(primary),
    Primary_Carrier_Fit_Score = round(top_carrier_score, 4),
    Secondary_Nanoformulation_Module = clean_label(secondary),
    Conditional_Functional_Module = clean_label(top_module),
    Ideal_Particle_Size = ideal_size,
    Ideal_PDI = ideal_pdi,
    Ideal_Zeta_Potential = zeta,
    Release_Strategy = release_strategy,
    Encapsulation_Driver = encapsulation_driver,
    Co_loading_Need = co_loading_need,
    Formulation_Priority_Index = round(formulation_priority_index, 4),
    Formulation_Priority_Class = formulation_priority_class,
    Final_Development_Interpretation = final_interpretation,
    Efficacy_Enhancement_Logic = efficacy_logic,
    Toxicity_Reduction_Logic = toxicity_logic,
    stringsAsFactors = FALSE
  )
}

make_blueprint_all <- function(df, carrier_scores, functional_module_scores) {
  
  res <- list()
  
  for (i in seq_len(nrow(df))) {
    one_row <- df[i, , drop = FALSE]
    
    one_carriers <- carrier_scores %>%
      dplyr::filter(compound == one_row$compound[1])
    
    one_modules <- functional_module_scores %>%
      dplyr::filter(compound == one_row$compound[1])
    
    res[[i]] <- make_blueprint_one(one_row, one_carriers, one_modules)
  }
  
  bind_rows(res)
}

make_design_logic_table <- function(df) {
  
  res <- list()
  
  for (i in seq_len(nrow(df))) {
    
    one <- df[i, , drop = FALSE]
    
    tmp <- data.frame(
      compound = one$compound[1],
      Design_Dimension = c(
        "Assembly utilization",
        "Cyclodextrin inclusion discrimination",
        "Metabolic control",
        "Barrier-aware delivery",
        "Exposure modulation",
        "Safety-oriented design",
        "Formulation-priority gate",
        "Co-loading rationale"
      ),
      Triggering_Evidence = c(
        paste0("Nano-assembly suitability score = ", round(one$Nano_Assembly_Suitability_Score, 3)),
        paste0(
          "Cyclodextrin cavity match = ", round(one$Cyclodextrin_Cavity_Match, 3),
          "; multivalent assembly drive = ", round(one$Multivalent_Assembly_Drive, 3)
        ),
        paste0("Metabolic liability score = ", round(one$Metabolic_Liability_Score, 3)),
        paste0("Barrier/transporter score = ", round(one$Barrier_Transporter_Score, 3)),
        paste0("Distribution/exposure score = ", round(one$Distribution_Exposure_Score, 3)),
        paste0("Safety-control requirement = ", round(one$Safety_Control_Requirement, 3)),
        paste0(
          "Formulation priority index = ", round(one$Formulation_Priority_Index, 3),
          " (", one$Formulation_Priority_Class, ")"
        ),
        paste0("Nano-optimization need = ", round(one$Nano_Optimization_Need, 3))
      ),
      Design_Action = c(
        "Prioritize self-assembly only when multivalent assembly drive supports stable non-covalent organization.",
        "Use cyclodextrin inclusion when small-aromatic cavity matching exceeds general self-assembly drive.",
        "Use sustained release to reduce Cmax and metabolic interaction burden.",
        "Consider ligand-modified or biomimetic surfaces when barrier targeting is relevant.",
        "Tune particle size, surface charge and release to improve effective exposure.",
        "Avoid exposure escalation; emphasize targeted or controlled delivery.",
        "Do not over-interpret carrier compatibility when formulation priority is low.",
        "Consider co-loading when barrier, metabolic or safety-control needs are moderate to high."
      ),
      stringsAsFactors = FALSE
    )
    
    res[[i]] <- tmp
  }
  
  bind_rows(res)
}

############################################################
## 14. Run Stage 2 model
############################################################

stage2 <- dat %>%
  add_stage1_internal_scores() %>%
  add_stage2_design_features()

carrier_feature_qc <- make_carrier_feature_qc(stage2)

if (any(carrier_feature_qc$is_missing)) {
  warning("Some internal carrier-fit features are missing. Check CarrierFeatureQC output.")
}

carrier_scores <- compute_carrier_scores_all(stage2)

if (all(carrier_scores$score == 0, na.rm = TRUE)) {
  stop(
    "All carrier-fit scores are zero. This indicates a scoring failure. ",
    "Please inspect CarrierFeatureQC and internal feature values."
  )
}

carrier_feature_contribution <- make_carrier_feature_contribution(stage2)
functional_module_scores <- compute_functional_module_scores_all(stage2)

sensitivity_table <- run_weight_sensitivity_all(
  stage2,
  n_iter = N_SENSITIVITY,
  jitter = WEIGHT_JITTER
)

blueprint_table <- make_blueprint_all(
  stage2,
  carrier_scores,
  functional_module_scores
)

design_logic_table <- make_design_logic_table(stage2)

############################################################
## 15. Summary tables
############################################################

stage2_summary <- carrier_scores %>%
  dplyr::filter(rank <= 3) %>%
  dplyr::arrange(compound, rank) %>%
  dplyr::mutate(
    recommendation_tier = dplyr::case_when(
      rank == 1 ~ "Primary backbone",
      rank == 2 ~ "Secondary backbone",
      rank == 3 ~ "Tertiary backbone",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::select(
    compound, PubChem_CID, InChIKey,
    recommendation_tier,
    carrier, carrier_label,
    carrier_fit_score = score,
    rank
  ) %>%
  dplyr::mutate(
    carrier_fit_score = round(carrier_fit_score, 4)
  )

carrier_scores$score <- round(carrier_scores$score, 4)
functional_module_scores$score <- round(functional_module_scores$score, 4)
carrier_feature_contribution$group_score <- round(carrier_feature_contribution$group_score, 4)
carrier_feature_contribution$weighted_contribution <- round(carrier_feature_contribution$weighted_contribution, 4)
carrier_feature_contribution$final_score_from_groups <- round(carrier_feature_contribution$final_score_from_groups, 4)

############################################################
## 16. Plot settings
############################################################

theme_acmpi <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "#2C2C2C"),
      plot.title = element_text(face = "bold", size = base_size + 4, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "#555555", hjust = 0, lineheight = 1.15),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 1, color = "#333333"),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      plot.margin = margin(24, 64, 24, 64),
      legend.position = "right"
    )
}

carrier_palette <- c(
  "Self Assembled Nanoformulation" = "#8DD3C7",
  "Cyclodextrin Inclusion System" = "#BEBADA",
  "Polymeric Micelle" = "#FDB462",
  "Lipid Nanoparticle" = "#80B1D3",
  "Polymeric Controlled Release Nanoparticle" = "#FB8072"
)

module_palette <- c(
  "Barrier Aware Ligand Module" = "#B3DE69",
  "Biomimetic Surface Module" = "#CCEBC5",
  "Exposure Reducing Stabilization Module" = "#FCCDE5",
  "Stimuli Responsive Release Module" = "#BC80BD"
)

section_palette <- c(
  "Carrier architecture" = "#8DD3C7",
  "Material parameters" = "#FDB462",
  "Design strategy" = "#BEBADA",
  "Development priority" = "#CCEBC5"
)

save_plot_dual <- function(plot, png_path, pdf_path, width, height, dpi = 600) {
  
  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  
  tryCatch(
    {
      ggsave(
        filename = pdf_path,
        plot = plot,
        width = width,
        height = height,
        device = cairo_pdf,
        bg = "white",
        limitsize = FALSE
      )
    },
    error = function(e) {
      ggsave(
        filename = pdf_path,
        plot = plot,
        width = width,
        height = height,
        device = "pdf",
        bg = "white",
        limitsize = FALSE
      )
    }
  )
}

############################################################
## 17. Figure 4: carrier-fit heatmap
############################################################

carrier_heatmap_df <- carrier_scores %>%
  dplyr::mutate(
    carrier_label_wrapped = wrap_label(carrier_label, width = 18),
    compound_wrapped = wrap_label(compound, width = 18)
  )

p_carrier_heatmap <- ggplot(
  carrier_heatmap_df,
  aes(x = carrier_label_wrapped, y = compound_wrapped, fill = score)
) +
  geom_tile(
    color = "white",
    linewidth = 0.8,
    width = 0.92,
    height = 0.82
  ) +
  geom_text(
    aes(label = sprintf("%.2f", score)),
    size = 4.2,
    fontface = "bold",
    color = "#2C2C2C"
  ) +
  scale_fill_gradientn(
    colors = c("#F7FBFF", "#DDEFD8", "#8DD3C7", "#4DAF9A"),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    name = "Carrier-fit\nscore"
  ) +
  labs(
    title = "Carrier-backbone fit landscape",
    subtitle = "Predefined carrier backbones scored by refined hierarchical compatibility rules",
    x = NULL,
    y = NULL
  ) +
  theme_acmpi(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 11),
    plot.margin = margin(24, 86, 48, 44)
  )

############################################################
## 18. Figure 4b: carrier bar plot
############################################################

make_carrier_bar_plot <- function(compound_name) {
  
  tmp <- carrier_scores %>%
    dplyr::filter(compound == compound_name) %>%
    dplyr::arrange(score) %>%
    dplyr::mutate(
      carrier_label_wrapped = wrap_label(carrier_label, width = 26),
      carrier_label_wrapped = factor(carrier_label_wrapped, levels = carrier_label_wrapped)
    )
  
  ggplot(
    tmp,
    aes(x = carrier_label_wrapped, y = score, fill = carrier_label)
  ) +
    geom_col(
      width = 0.72,
      color = "white",
      linewidth = 0.6,
      alpha = 0.96
    ) +
    geom_text(
      aes(label = sprintf("%.2f", score)),
      hjust = -0.12,
      size = 4.0,
      fontface = "bold",
      color = "#333333"
    ) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = carrier_palette, guide = "none") +
    scale_y_continuous(
      limits = c(0, 1.12),
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
      title = "Carrier-backbone fit profile",
      subtitle = paste0("Compound: ", compound_name),
      x = NULL,
      y = "Carrier-fit score"
    ) +
    theme_acmpi(base_size = 13) +
    theme(
      axis.text.y = element_text(face = "bold", size = 10.8),
      plot.margin = margin(24, 105, 24, 28)
    )
}

############################################################
## 19. Figure 5: blueprint tile plot
############################################################

blueprint_long <- blueprint_table %>%
  dplyr::select(
    compound,
    Primary_Nanoformulation_Backbone,
    Secondary_Nanoformulation_Module,
    Conditional_Functional_Module,
    Ideal_Particle_Size,
    Ideal_PDI,
    Ideal_Zeta_Potential,
    Release_Strategy,
    Encapsulation_Driver,
    Co_loading_Need,
    Formulation_Priority_Class,
    Final_Development_Interpretation
  ) %>%
  tidyr::pivot_longer(
    cols = -compound,
    names_to = "Blueprint_Item",
    values_to = "Blueprint_Value"
  ) %>%
  dplyr::mutate(
    Blueprint_Item_Label = clean_label(Blueprint_Item),
    Blueprint_Item_Label = wrap_label(Blueprint_Item_Label, width = 18),
    Blueprint_Value_Label = wrap_label(Blueprint_Value, width = 40),
    section = dplyr::case_when(
      Blueprint_Item %in% c(
        "Primary_Nanoformulation_Backbone",
        "Secondary_Nanoformulation_Module",
        "Conditional_Functional_Module"
      ) ~ "Carrier architecture",
      Blueprint_Item %in% c(
        "Ideal_Particle_Size",
        "Ideal_PDI",
        "Ideal_Zeta_Potential"
      ) ~ "Material parameters",
      Blueprint_Item %in% c(
        "Formulation_Priority_Class",
        "Final_Development_Interpretation"
      ) ~ "Development priority",
      TRUE ~ "Design strategy"
    )
  )

blueprint_long$Blueprint_Item_Label <- factor(
  blueprint_long$Blueprint_Item_Label,
  levels = rev(unique(blueprint_long$Blueprint_Item_Label))
)

p_blueprint <- ggplot(
  blueprint_long,
  aes(x = compound, y = Blueprint_Item_Label, fill = section)
) +
  geom_tile(
    color = "white",
    linewidth = 0.8,
    width = 0.92,
    height = 0.86,
    alpha = 0.92
  ) +
  geom_text(
    aes(label = Blueprint_Value_Label),
    size = 3.15,
    color = "#2C2C2C",
    lineheight = 0.95
  ) +
  scale_fill_manual(values = section_palette, name = "Blueprint\nsection") +
  labs(
    title = "Ideal nanoformulation blueprint",
    subtitle = "Model-translated material parameters, formulation logic, and development priority",
    x = NULL,
    y = NULL
  ) +
  theme_acmpi(base_size = 13) +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 10.5),
    panel.border = element_blank(),
    plot.margin = margin(24, 110, 24, 44)
  )

############################################################
## 20. Figure 6: sensitivity interval plot
############################################################

sensitivity_plot_df <- sensitivity_table %>%
  dplyr::mutate(
    carrier_label_wrapped = wrap_label(carrier_label, width = 26),
    carrier_label_wrapped = factor(
      carrier_label_wrapped,
      levels = carrier_label_wrapped[order(baseline_score)]
    ),
    top_label = paste0("Top-rank: ", percent(top_rank_probability, accuracy = 1))
  )

p_sensitivity <- ggplot(
  sensitivity_plot_df,
  aes(y = carrier_label_wrapped)
) +
  geom_errorbarh(
    aes(xmin = sensitivity_p05, xmax = sensitivity_p95),
    height = 0.20,
    linewidth = 1.0,
    color = "#7A7A7A"
  ) +
  geom_point(
    aes(x = sensitivity_mean, fill = carrier_label),
    shape = 21,
    size = 5.0,
    color = "white",
    stroke = 0.7
  ) +
  geom_point(
    aes(x = baseline_score),
    shape = 23,
    size = 3.8,
    fill = "#333333",
    color = "white",
    stroke = 0.5
  ) +
  geom_text(
    aes(
      x = pmin(sensitivity_p95 + 0.04, 1.05),
      label = top_label
    ),
    hjust = 0,
    size = 3.6,
    fontface = "bold",
    color = "#333333"
  ) +
  scale_fill_manual(values = carrier_palette, guide = "none") +
  scale_x_continuous(
    limits = c(0, 1.18),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "Carrier-backbone sensitivity profile",
    subtitle = paste0(
      "Baseline score, sensitivity mean, 5-95% interval, and top-rank probability; ",
      N_SENSITIVITY,
      " iterations; +/-",
      WEIGHT_JITTER * 100,
      "% group-weight perturbation"
    ),
    x = "Carrier-fit score under weight perturbation",
    y = NULL
  ) +
  theme_acmpi(base_size = 13) +
  theme(
    axis.text.y = element_text(face = "bold", size = 10.8),
    plot.subtitle = element_text(size = 12.2, color = "#555555", lineheight = 1.15),
    plot.margin = margin(24, 135, 24, 28)
  )

############################################################
## 21. Figure 6b: functional module scores
############################################################

p_functional_module <- ggplot(
  functional_module_scores %>%
    dplyr::arrange(score) %>%
    dplyr::mutate(
      module_label_wrapped = wrap_label(module_label, width = 28),
      module_label_wrapped = factor(module_label_wrapped, levels = module_label_wrapped)
    ),
  aes(x = module_label_wrapped, y = score, fill = module_label)
) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.6,
    alpha = 0.96
  ) +
  geom_text(
    aes(label = paste0(sprintf("%.2f", score), " | ", recommendation)),
    hjust = -0.08,
    size = 3.7,
    fontface = "bold",
    color = "#333333"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = module_palette, guide = "none") +
  scale_y_continuous(
    limits = c(0, 1.18),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Functional module recommendation profile",
    subtitle = "Functional modules are scored separately from carrier backbones",
    x = NULL,
    y = "Functional module score"
  ) +
  theme_acmpi(base_size = 13) +
  theme(
    axis.text.y = element_text(face = "bold", size = 10.8),
    plot.margin = margin(24, 130, 24, 28)
  )

############################################################
## 22. Figure 7: ideal nanoformulation schematic
############################################################

make_schematic_plot <- function(blueprint_row) {
  
  comp <- blueprint_row$compound[1]
  primary <- blueprint_row$Primary_Nanoformulation_Backbone[1]
  secondary <- blueprint_row$Secondary_Nanoformulation_Module[1]
  module <- blueprint_row$Conditional_Functional_Module[1]
  size_text <- blueprint_row$Ideal_Particle_Size[1]
  zeta_text <- blueprint_row$Ideal_Zeta_Potential[1]
  release_text <- blueprint_row$Release_Strategy[1]
  driver_text <- blueprint_row$Encapsulation_Driver[1]
  priority_text <- blueprint_row$Final_Development_Interpretation[1]
  
  theta <- seq(0, 2 * pi, length.out = 250)
  
  make_circle <- function(x0, y0, r, layer) {
    data.frame(
      x = x0 + r * cos(theta),
      y = y0 + r * sin(theta),
      layer = layer
    )
  }
  
  circle_poly <- bind_rows(
    make_circle(0, 0, 1.60, "Functional surface module"),
    make_circle(0, 0, 1.23, "Secondary stabilization module"),
    make_circle(0, 0, 0.72, "Core carrier backbone")
  )
  
  layer_cols <- c(
    "Functional surface module" = "#CCEBC5",
    "Secondary stabilization module" = "#BEBADA",
    "Core carrier backbone" = "#8DD3C7"
  )
  
  anno <- data.frame(
    x = c(2.35, 2.35, 2.35, -2.25, -2.25),
    y = c(1.25, 0.35, -0.55, 0.75, -0.55),
    label = c(
      paste0("Primary backbone:\n", wrap_label(primary, 28)),
      paste0("Secondary module:\n", wrap_label(secondary, 28)),
      paste0("Functional module:\n", wrap_label(module, 28)),
      paste0("Size: ", size_text, "\nZeta: ", wrap_label(zeta_text, 24)),
      paste0("Release:\n", wrap_label(release_text, 28))
    ),
    stringsAsFactors = FALSE
  )
  
  ggplot() +
    geom_polygon(
      data = circle_poly,
      aes(x = x, y = y, fill = layer, group = layer),
      color = "white",
      linewidth = 1.1,
      alpha = 0.95
    ) +
    geom_text(
      aes(x = 0, y = 0.10),
      label = paste0(comp, "\ncore"),
      fontface = "bold",
      size = 4.4,
      color = "#2C2C2C",
      lineheight = 0.95
    ) +
    geom_segment(aes(x = 1.72, y = 1.00, xend = 2.08, yend = 1.18),
                 arrow = arrow(length = unit(0.16, "inches")),
                 linewidth = 0.65, color = "#555555") +
    geom_segment(aes(x = 1.42, y = 0.25, xend = 2.05, yend = 0.33),
                 arrow = arrow(length = unit(0.16, "inches")),
                 linewidth = 0.65, color = "#555555") +
    geom_segment(aes(x = 1.65, y = -0.60, xend = 2.05, yend = -0.55),
                 arrow = arrow(length = unit(0.16, "inches")),
                 linewidth = 0.65, color = "#555555") +
    geom_segment(aes(x = -1.58, y = 0.62, xend = -2.00, yend = 0.70),
                 arrow = arrow(length = unit(0.16, "inches")),
                 linewidth = 0.65, color = "#555555") +
    geom_segment(aes(x = -1.30, y = -0.45, xend = -1.98, yend = -0.52),
                 arrow = arrow(length = unit(0.16, "inches")),
                 linewidth = 0.65, color = "#555555") +
    geom_label(
      data = anno,
      aes(x = x, y = y, label = label),
      size = 3.55,
      lineheight = 0.96,
      fill = "white",
      label.size = 0.25,
      label.r = unit(0.15, "lines"),
      color = "#2C2C2C"
    ) +
    annotate(
      "text",
      x = 0,
      y = -2.03,
      label = paste0("Encapsulation driver: ", wrap_label(driver_text, 70)),
      size = 3.6,
      color = "#444444"
    ) +
    annotate(
      "text",
      x = 0,
      y = -2.28,
      label = paste0("Development interpretation: ", wrap_label(priority_text, 72)),
      size = 3.45,
      color = "#444444"
    ) +
    scale_fill_manual(values = layer_cols, name = "Schematic layer") +
    coord_equal(xlim = c(-3.25, 3.35), ylim = c(-2.55, 2.25), clip = "off") +
    labs(
      title = "Ideal nanoformulation schematic",
      subtitle = "Model-translated structural blueprint from compound-level features",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.02, color = "#2C2C2C"),
      plot.subtitle = element_text(size = 12.5, color = "#555555", hjust = 0.02),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      plot.margin = margin(24, 45, 24, 45)
    )
}

p_schematic <- make_schematic_plot(blueprint_table[1, , drop = FALSE])

############################################################
## 23. Save tables
############################################################

out_csv_stage2 <- file.path(table_dir, paste0(output_prefix, "_FullScoredData.csv"))
out_csv_classdef <- file.path(table_dir, paste0(output_prefix, "_NanoformulationClassDefinition.csv"))
out_csv_carrier <- file.path(table_dir, paste0(output_prefix, "_CarrierBackboneFitScores.csv"))
out_csv_contrib <- file.path(table_dir, paste0(output_prefix, "_CarrierFeatureContribution.csv"))
out_csv_module <- file.path(table_dir, paste0(output_prefix, "_FunctionalModuleScores.csv"))
out_csv_summary <- file.path(table_dir, paste0(output_prefix, "_TopCarrierSummary.csv"))
out_csv_blueprint <- file.path(table_dir, paste0(output_prefix, "_Blueprint.csv"))
out_csv_logic <- file.path(table_dir, paste0(output_prefix, "_DesignLogic.csv"))
out_csv_sensitivity <- file.path(table_dir, paste0(output_prefix, "_WeightSensitivity.csv"))
out_csv_carrier_qc <- file.path(table_dir, paste0(output_prefix, "_CarrierFeatureQC.csv"))
out_xlsx <- file.path(table_dir, paste0(output_prefix, "_Results.xlsx"))

write.csv(stage2, out_csv_stage2, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(class_definition_table, out_csv_classdef, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(carrier_scores, out_csv_carrier, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(carrier_feature_contribution, out_csv_contrib, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(functional_module_scores, out_csv_module, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(stage2_summary, out_csv_summary, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(blueprint_table, out_csv_blueprint, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(design_logic_table, out_csv_logic, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sensitivity_table, out_csv_sensitivity, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(carrier_feature_qc, out_csv_carrier_qc, row.names = FALSE, fileEncoding = "UTF-8")

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "ClassDefinition")
openxlsx::writeData(wb, "ClassDefinition", class_definition_table)

openxlsx::addWorksheet(wb, "TopCarrierSummary")
openxlsx::writeData(wb, "TopCarrierSummary", stage2_summary)

openxlsx::addWorksheet(wb, "CarrierBackboneFit")
openxlsx::writeData(wb, "CarrierBackboneFit", carrier_scores)

openxlsx::addWorksheet(wb, "CarrierContribution")
openxlsx::writeData(wb, "CarrierContribution", carrier_feature_contribution)

openxlsx::addWorksheet(wb, "FunctionalModules")
openxlsx::writeData(wb, "FunctionalModules", functional_module_scores)

openxlsx::addWorksheet(wb, "Blueprint")
openxlsx::writeData(wb, "Blueprint", blueprint_table)

openxlsx::addWorksheet(wb, "DesignLogic")
openxlsx::writeData(wb, "DesignLogic", design_logic_table)

openxlsx::addWorksheet(wb, "WeightSensitivity")
openxlsx::writeData(wb, "WeightSensitivity", sensitivity_table)

openxlsx::addWorksheet(wb, "CarrierFeatureQC")
openxlsx::writeData(wb, "CarrierFeatureQC", carrier_feature_qc)

openxlsx::addWorksheet(wb, "FullScoredData")
openxlsx::writeData(wb, "FullScoredData", stage2)

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

############################################################
## 24. Save figures
############################################################

fig4_png <- file.path(fig_png_dir, "Fig4_CarrierBackboneFit_Heatmap.png")
fig4_pdf <- file.path(fig_pdf_dir, "Fig4_CarrierBackboneFit_Heatmap.pdf")

fig5_png <- file.path(fig_png_dir, "Fig5_Ideal_Nanoformulation_Blueprint.png")
fig5_pdf <- file.path(fig_pdf_dir, "Fig5_Ideal_Nanoformulation_Blueprint.pdf")

fig6_png <- file.path(fig_png_dir, "Fig6_CarrierBackbone_Sensitivity_Profile.png")
fig6_pdf <- file.path(fig_pdf_dir, "Fig6_CarrierBackbone_Sensitivity_Profile.pdf")

fig6b_png <- file.path(fig_png_dir, "Fig6b_FunctionalModule_Profile.png")
fig6b_pdf <- file.path(fig_pdf_dir, "Fig6b_FunctionalModule_Profile.pdf")

fig7_png <- file.path(fig_png_dir, "Fig7_Ideal_Nanoformulation_Schematic.png")
fig7_pdf <- file.path(fig_pdf_dir, "Fig7_Ideal_Nanoformulation_Schematic.pdf")

save_plot_dual(
  plot = p_carrier_heatmap,
  png_path = fig4_png,
  pdf_path = fig4_pdf,
  width = 12.8,
  height = 7.0 + 0.25 * nrow(dat),
  dpi = FIG_DPI
)

save_plot_dual(
  plot = p_blueprint,
  png_path = fig5_png,
  pdf_path = fig5_pdf,
  width = 15.2,
  height = 10.4 + 0.15 * nrow(dat),
  dpi = FIG_DPI
)

save_plot_dual(
  plot = p_sensitivity,
  png_path = fig6_png,
  pdf_path = fig6_pdf,
  width = 13.8,
  height = 7.6,
  dpi = FIG_DPI
)

save_plot_dual(
  plot = p_functional_module,
  png_path = fig6b_png,
  pdf_path = fig6b_pdf,
  width = 13.2,
  height = 6.8,
  dpi = FIG_DPI
)

save_plot_dual(
  plot = p_schematic,
  png_path = fig7_png,
  pdf_path = fig7_pdf,
  width = 12.8,
  height = 8.6,
  dpi = FIG_DPI
)

for (compound_name in unique(carrier_scores$compound)) {
  
  p_bar <- make_carrier_bar_plot(compound_name)
  compound_id <- safe_name(compound_name)
  
  fig_bar_png <- file.path(
    fig_png_dir,
    paste0("Fig4b_CarrierBackboneFit_Barplot_", compound_id, ".png")
  )
  
  fig_bar_pdf <- file.path(
    fig_pdf_dir,
    paste0("Fig4b_CarrierBackboneFit_Barplot_", compound_id, ".pdf")
  )
  
  save_plot_dual(
    plot = p_bar,
    png_path = fig_bar_png,
    pdf_path = fig_bar_pdf,
    width = 12.8,
    height = 7.4,
    dpi = FIG_DPI
  )
}

############################################################
## 25. Print summary
############################################################

cat("\n============================================================\n")
cat("ACMPI-Nano Stage 2 v1.2: Ideal Nanoformulation Blueprint Model\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat(input_file, "\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Top carrier-backbone summary:\n")
print(stage2_summary, row.names = FALSE)

cat("\nCarrier-backbone fit scores:\n")
print(carrier_scores, row.names = FALSE)

cat("\nFunctional module scores:\n")
print(functional_module_scores, row.names = FALSE)

cat("\nBlueprint table:\n")
print(blueprint_table, row.names = FALSE)

cat("\nDesign logic table:\n")
print(design_logic_table, row.names = FALSE)

cat("\nCarrier feature contribution:\n")
print(carrier_feature_contribution, row.names = FALSE)

cat("\nWeight sensitivity summary:\n")
print(
  sensitivity_table[, c(
    "compound",
    "carrier_label",
    "baseline_score",
    "sensitivity_mean",
    "sensitivity_p05",
    "sensitivity_p95",
    "top_rank_probability"
  )],
  row.names = FALSE
)

cat("\nCarrier feature QC summary:\n")
print(
  carrier_feature_qc %>%
    dplyr::group_by(compound) %>%
    dplyr::summarise(
      n_features = dplyr::n(),
      n_missing = sum(is_missing),
      .groups = "drop"
    ),
  row.names = FALSE
)

cat("\nSaved tables:\n")
cat(out_csv_stage2, "\n")
cat(out_csv_classdef, "\n")
cat(out_csv_carrier, "\n")
cat(out_csv_contrib, "\n")
cat(out_csv_module, "\n")
cat(out_csv_summary, "\n")
cat(out_csv_blueprint, "\n")
cat(out_csv_logic, "\n")
cat(out_csv_sensitivity, "\n")
cat(out_csv_carrier_qc, "\n")
cat(out_xlsx, "\n")

cat("\nSaved figures in:\n")
cat(fig_png_dir, "\n")
cat(fig_pdf_dir, "\n")

cat("\nInterpretation guide:\n")
cat("- Carrier-backbone scores compare predefined primary nanoformulation backbones only.\n")
cat("- Functional module scores are evaluated separately and are not treated as competing carrier backbones.\n")
cat("- Multivalent_Assembly_Drive prevents small aromatic compounds from being over-classified as self-assembly candidates.\n")
cat("- Cyclodextrin_Cavity_Match and Small_Aromatic_Inclusion_Preference improve cyclodextrin inclusion recognition.\n")
cat("- Formulation_Priority_Index prevents structurally compatible but low-need compounds from being over-interpreted as high-priority nanoformulation candidates.\n")
cat("- CarrierFeatureContribution explains whether each carrier score is driven by primary, secondary, or modifier features.\n")
cat("- The schematic figure translates the blueprint into an interpretable ideal nanoformulation image.\n")
cat("- Safety and toxicity signals are interpreted as design-control requirements, not exclusion criteria.\n")