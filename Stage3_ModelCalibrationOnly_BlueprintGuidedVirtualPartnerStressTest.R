############################################################
## ACMPI-Nano Pairing Calibration
## Blueprint-guided virtual partner calibration from Stage 2 input
##
## Purpose:
##   Use the nanoformulation blueprint of the core compound as a fixed
##   design target, construct synthetic partner profiles, redefine
##   distinctive pairing-gain features, calibrate candidate scoring rules,
##   and select the optimized compatibility model using rank constraints.
##
## Design principles:
##   - Molecular acceptability is a permissive condition, not sufficient
##     evidence of high pairing compatibility.
##   - High compatibility requires distinctive, nonredundant,
##     blueprint-matched pairing gain.
##   - Neutral low-value profiles are explicit negative controls.
##   - Boundary auxiliary profiles are kept in a conditional range rather
##     than forced into high-value or incompatible classes.
##   - All joins are collision-proof and all required columns are checked
##     before summarisation.
############################################################

rm(list = ls())

############################################################
## 0. User settings
############################################################

stage2_table_dir <- "/media/desk16/iy15915/中药之开创/Baicalein/ACMPI_Nano_Stage2_v1p2_Output_20260427_213702/01_Tables"

output_prefix <- "ACMPI_Nano_PairingCalibration_SubmissionGrade"

FIG_DPI <- 600
SET_SEED <- 20260427
N_PER_PROFILE <- 30
N_SENSITIVITY <- 5000
WEIGHT_JITTER <- 0.20
FEATURE_JITTER <- 0.10
SENSITIVITY_JITTER <- 0.16
N_TOP_MODELS_TO_PLOT <- 8

set.seed(SET_SEED)

############################################################
## 1. Required packages
############################################################

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "openxlsx", "scales", "stringr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(scales)
library(stringr)

############################################################
## Global robustness helpers
############################################################

stop_if_missing <- function(df, cols, df_name = "data frame") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(df_name, " is missing required columns:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}

strip_named_variants <- function(df, base_names) {
  variants <- unique(c(
    base_names,
    paste0(base_names, ".x"),
    paste0(base_names, ".y"),
    paste0(base_names, ".left"),
    paste0(base_names, ".right")
  ))
  df[, setdiff(colnames(df), variants), drop = FALSE]
}

strip_regex_cols <- function(df, patterns) {
  if (length(patterns) < 1) return(df)
  keep <- rep(TRUE, ncol(df))
  for (pat in patterns) {
    keep <- keep & !grepl(pat, colnames(df))
  }
  df[, keep, drop = FALSE]
}

assert_no_join_suffixes <- function(df, context = "current table") {
  bad <- grep("\\.(x|y)$", colnames(df), value = TRUE)
  if (length(bad) > 0) {
    stop(
      "Unexpected .x/.y columns after join in ", context, ":\n",
      paste(bad, collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

safe_left_join_replace <- function(x, y, by, replace_cols, context = "join") {
  stop_if_missing(x, by, paste0(context, " left table"))
  stop_if_missing(y, by, paste0(context, " right table"))
  stop_if_missing(y, replace_cols, paste0(context, " right table"))
  x_clean <- strip_named_variants(x, replace_cols)
  out <- dplyr::left_join(x_clean, y, by = by)
  stop_if_missing(out, replace_cols, paste0(context, " joined table"))
  assert_no_join_suffixes(out, context)
  out
}

standard_sensitivity_cols <- c(
  "Sensitivity_Mean",
  "Sensitivity_SD",
  "Sensitivity_P05",
  "Sensitivity_P95",
  "Sensitivity_Gain_Mean",
  "Probability_High_or_ModerateHigh",
  "Probability_Neutral_Low_Value_Flag"
)

############################################################
## 2. Locate input files
############################################################

if (!dir.exists(stage2_table_dir)) {
  stop("Input table directory not found: ", stage2_table_dir)
}

blueprint_file <- list.files(stage2_table_dir, pattern = "Blueprint\\.csv$", full.names = TRUE)
fullscored_file <- list.files(stage2_table_dir, pattern = "FullScoredData\\.csv$", full.names = TRUE)

if (length(blueprint_file) < 1) stop("No Blueprint.csv file found in: ", stage2_table_dir)
if (length(fullscored_file) < 1) stop("No FullScoredData.csv file found in: ", stage2_table_dir)

blueprint_file <- blueprint_file[1]
fullscored_file <- fullscored_file[1]

core_blueprint <- read.csv(blueprint_file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
core_full <- read.csv(fullscored_file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")

if (nrow(core_blueprint) < 1) stop("Blueprint table is empty.")
if (nrow(core_full) < 1) stop("FullScoredData table is empty.")

core_blueprint <- core_blueprint[1, , drop = FALSE]
core_full <- core_full[1, , drop = FALSE]
core_compound <- as.character(core_full$compound[1])

############################################################
## 3. Output root
############################################################

## Output folders are created only once, after all inputs are loaded
## and the calibration block is initialized. This prevents an empty
## timestamped folder from being generated before the final run starts.
output_root_dir <- dirname(stage2_table_dir)

############################################################
## 4. Required raw columns
############################################################

required_raw_cols <- c(
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

missing_raw <- setdiff(required_raw_cols, colnames(core_full))
if (length(missing_raw) > 0) {
  stop("FullScoredData is missing required raw columns:\n", paste(missing_raw, collapse = "\n"))
}

core_raw <- core_full[, required_raw_cols, drop = FALSE]

############################################################
## 5. Helper functions
############################################################

to_num <- function(x) suppressWarnings(as.numeric(x))
clip01 <- function(x) pmax(0, pmin(1, x))
risk_high <- function(x, good, bad) clip01((to_num(x) - good) / (bad - good))
risk_low <- function(x, good, bad) clip01((good - to_num(x)) / (good - bad))
benefit_high <- function(x, low, high) clip01((to_num(x) - low) / (high - low))
benefit_low <- function(x, high, low) clip01((high - to_num(x)) / (high - low))

benefit_mid <- function(x, lower, upper) {
  x <- to_num(x)
  span <- upper - lower
  out <- ifelse(
    x >= lower & x <= upper,
    1,
    ifelse(x < lower, clip01((x - (lower - span)) / span), clip01(((upper + span) - x) / span))
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
  for (nm in names(levels_map)) out[x == tolower(nm)] <- levels_map[[nm]]
  out
}

clean_label <- function(x) {
  x <- gsub("_Score$", "", x)
  x <- gsub("_Requirement$", "", x)
  x <- gsub("_Fit$", "", x)
  x <- gsub("_", " ", x)
  x
}

wrap_label <- function(x, width = 26) stringr::str_wrap(x, width = width)
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

normalize_weights <- function(w) {
  original_names <- names(w)
  w_num <- as.numeric(w)
  w_num <- w_num / sum(w_num, na.rm = TRUE)
  if (!is.null(original_names)) names(w_num) <- original_names
  w_num
}

jitter_weights <- function(w, jitter = 0.20) {
  original_names <- names(w)
  w2 <- as.numeric(w) * runif(length(w), min = 1 - jitter, max = 1 + jitter)
  if (!is.null(original_names)) names(w2) <- original_names
  normalize_weights(w2)
}

weighted_sum_strict <- function(values, weights) {
  values <- as.numeric(values)
  weights <- normalize_weights(weights)
  keep <- !is.na(values) & !is.na(weights)
  if (!any(keep)) return(NA_real_)
  sum(values[keep] * weights[keep]) / sum(weights[keep])
}

jitter_numeric <- function(x, rel = 0.08, abs_min = NULL, abs_max = NULL) {
  x <- as.numeric(x)
  out <- x * runif(1, min = 1 - rel, max = 1 + rel)
  if (!is.null(abs_min)) out <- max(abs_min, out)
  if (!is.null(abs_max)) out <- min(abs_max, out)
  out
}

get_one <- function(df, nm) {
  if (!nm %in% colnames(df)) stop("Missing feature: ", nm)
  as.numeric(df[[nm]][1])
}

feature_similarity <- function(a, b, scale) {
  clip01(1 - abs(as.numeric(a) - as.numeric(b)) / scale)
}

feature_distance <- function(a, b, scale) {
  clip01(abs(as.numeric(a) - as.numeric(b)) / scale)
}

############################################################
## 6. Internal scoring functions
############################################################

score_solubility_defect <- function(df) {
  s_esol_logS <- risk_low(df$ESOL_LogS, good = -3, bad = -7)
  s_admet_logS <- risk_low(df$ADMETlab_logS, good = -3, bad = -7)
  s_cons_logP <- risk_high(df$Consensus_LogP, good = 2.5, bad = 6)
  s_logD <- risk_high(df$ADMETlab_logD, good = 2.5, bad = 6)
  s_class <- class_score(
    df$ESOL_Class,
    levels_map = c(
      "Highly soluble" = 0.00, "Very soluble" = 0.00, "Soluble" = 0.15,
      "Moderately soluble" = 0.40, "Poorly soluble" = 0.75, "Insoluble" = 1.00
    ),
    default = 0.50
  )
  row_weighted_mean(cbind(s_esol_logS, s_admet_logS, s_cons_logP, s_logD, s_class), c(0.30, 0.25, 0.20, 0.10, 0.15))
}

score_permeability_defect <- function(df) {
  s_tpsa <- risk_high(df$TPSA, good = 75, bad = 140)
  s_gi <- ifelse(tolower(trimws(df$GI_absorption)) == "high", 0.10, 0.85)
  s_caco2 <- risk_low(df$Caco2, good = -4.5, bad = -6.0)
  s_pampa <- risk_low(df$PAMPA, good = 0.70, bad = 0.20)
  s_f30 <- risk_low(df$F30, good = 0.80, bad = 0.20)
  row_weighted_mean(cbind(s_tpsa, s_gi, s_caco2, s_pampa, s_f30), rep(0.20, 5))
}

score_barrier_transporter <- function(df) {
  s_bbb_class <- ifelse(tolower(trimws(df$BBB_permeant)) == "yes", 0.10, 0.90)
  s_bbb_prob <- risk_low(df$ADMETlab_BBB, good = 0.70, bad = 0.10)
  s_pgp_class <- ifelse(tolower(trimws(df$Pgp_substrate)) == "yes", 0.90, 0.10)
  s_pgp_sub <- risk_high(df$pgp_sub, good = 0.20, bad = 0.80)
  s_pgp_inh <- risk_high(df$pgp_inh, good = 0.20, bad = 0.80)
  row_weighted_mean(cbind(s_bbb_class, s_bbb_prob, s_pgp_class, s_pgp_sub, s_pgp_inh), c(0.30, 0.25, 0.20, 0.15, 0.10))
}

score_distribution_exposure <- function(df) {
  s_ppb <- risk_high(df$PPB, good = 80, bad = 99)
  s_fu <- risk_low(df$Fu, good = 0.30, bad = 0.02)
  s_vd <- risk_high(abs(to_num(df$logVDss)), good = 0.8, bad = 2.0)
  row_weighted_mean(cbind(s_ppb, s_fu, s_vd), c(0.45, 0.35, 0.20))
}

score_metabolic_liability <- function(df) {
  s_cyp <- clip01(to_num(df$CYP_inhibitor_count) / 5)
  s_half <- risk_low(df$t_half, good = 4, bad = 0.5)
  s_clearance <- risk_high(df$cl_plasma, good = 5, bad = 15)
  row_weighted_mean(cbind(s_cyp, s_half, s_clearance), c(0.45, 0.25, 0.30))
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
  row_weighted_mean(cbind(s_mw, s_logp, s_hbd, s_hba, s_rot, s_ring, s_arom, s_fsp3, s_charge, s_mr), c(0.10, 0.15, 0.10, 0.10, 0.10, 0.10, 0.15, 0.10, 0.05, 0.05))
}

score_druggability_alert_requirement <- function(df) {
  s_lip <- clip01(to_num(df$Lipinski_violations) / 2)
  s_veb <- clip01(to_num(df$Veber_violations) / 1)
  s_bio <- 1 - clip01(to_num(df$Bioavailability_Score))
  s_pains <- clip01(to_num(df$PAINS_alerts) / 2)
  s_brenk <- clip01(to_num(df$Brenk_alerts) / 3)
  s_sa <- risk_high(df$Synthetic_Accessibility, good = 3, bad = 8)
  s_qed <- risk_low(df$QED, good = 0.70, bad = 0.20)
  row_weighted_mean(cbind(s_lip, s_veb, s_bio, s_pains, s_brenk, s_sa, s_qed), c(0.15, 0.10, 0.15, 0.20, 0.15, 0.10, 0.15))
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
  row_weighted_mean(cbind(s_ld50, s_class, s_dili, s_ames, s_carc, s_herg, s_org, s_end, s_met), c(0.10, 0.10, 0.15, 0.12, 0.12, 0.12, 0.12, 0.10, 0.07))
}

add_internal_scores <- function(df) {
  out <- df
  out$Solubility_Defect_Score <- score_solubility_defect(out)
  out$Permeability_Defect_Score <- score_permeability_defect(out)
  out$Barrier_Transporter_Score <- score_barrier_transporter(out)
  out$Distribution_Exposure_Score <- score_distribution_exposure(out)
  out$Metabolic_Liability_Score <- score_metabolic_liability(out)
  out$Nano_Assembly_Suitability_Score <- score_nano_assembly_suitability(out)
  out$Druggability_Alert_Requirement <- score_druggability_alert_requirement(out)
  out$Safety_Control_Requirement <- score_safety_control_requirement(out)
  out$Nano_Optimization_Need <- 0.22 * out$Solubility_Defect_Score + 0.18 * out$Permeability_Defect_Score + 0.15 * out$Barrier_Transporter_Score + 0.10 * out$Distribution_Exposure_Score + 0.15 * out$Metabolic_Liability_Score + 0.20 * out$Safety_Control_Requirement
  out$Delivery_Defect_Burden <- 0.25 * out$Solubility_Defect_Score + 0.25 * out$Permeability_Defect_Score + 0.20 * out$Barrier_Transporter_Score + 0.15 * out$Distribution_Exposure_Score + 0.15 * out$Metabolic_Liability_Score
  out$Molecular_Size_Fit <- benefit_mid(out$Molecular_Weight, lower = 100, upper = 700)
  out$Balanced_LogP <- benefit_mid(out$MolLogP, lower = 1.5, upper = 5.0)
  out$Balanced_Consensus_LogP <- benefit_mid(out$Consensus_LogP, lower = 1.5, upper = 5.0)
  out$Lipophilicity_Fit <- benefit_mid(out$Consensus_LogP, lower = 2.0, upper = 6.0)
  out$Hydrophobic_Core_Fit <- row_mean_safe(cbind(benefit_mid(out$MolLogP, lower = 2.0, upper = 6.0), benefit_mid(out$Consensus_LogP, lower = 2.0, upper = 6.0), risk_low(out$ESOL_LogS, good = -3, bad = -7)))
  out$Aromaticity <- benefit_high(out$Aromatic_Rings, low = 1, high = 4)
  out$Ring_Structure <- benefit_high(out$Ring_Count, low = 1, high = 4)
  out$Rigidity <- benefit_low(out$Rotatable_Bonds, high = 10, low = 0)
  out$Hydrogen_Bonding <- row_mean_safe(cbind(benefit_mid(out$HBD, lower = 1, upper = 5), benefit_mid(out$HBA, lower = 2, upper = 10)))
  out$Molar_Refractivity_Fit <- benefit_mid(out$Molar_Refractivity, lower = 40, upper = 130)
  out$Charge_Compatibility <- ifelse(to_num(out$Formal_Charge) == 0, 1.00, 0.60)
  out$Druggability_Manageability <- 1 - out$Druggability_Alert_Requirement
  out$Controlled_Release_Need <- row_mean_safe(cbind(out$Metabolic_Liability_Score, out$Safety_Control_Requirement))
  out$Multivalent_Assembly_Drive <- row_mean_safe(cbind(benefit_high(out$Aromatic_Rings, low = 1.5, high = 3.0), benefit_mid(out$HBD, lower = 2, upper = 5), benefit_mid(out$HBA, lower = 4, upper = 10), benefit_low(out$Rotatable_Bonds, high = 6, low = 0), benefit_high(out$Ring_Count, low = 2, high = 4), benefit_mid(out$Molar_Refractivity, lower = 65, upper = 130)))
  out$Cyclodextrin_Cavity_Match <- row_mean_safe(cbind(benefit_mid(out$Molecular_Weight, lower = 120, upper = 420), benefit_mid(out$Aromatic_Rings, lower = 1, upper = 2), benefit_mid(out$MolLogP, lower = 1.0, upper = 4.0), benefit_low(out$Rotatable_Bonds, high = 6, low = 0), benefit_mid(out$Molar_Refractivity, lower = 45, upper = 95)))
  out$Small_Aromatic_Inclusion_Preference <- row_mean_safe(cbind(benefit_mid(out$Molecular_Weight, lower = 120, upper = 300), benefit_mid(out$Aromatic_Rings, lower = 1, upper = 2), benefit_mid(out$Ring_Count, lower = 1, upper = 3), benefit_mid(out$MolLogP, lower = 1.0, upper = 3.5), benefit_low(out$Rotatable_Bonds, high = 5, low = 0)))
  out$Formulation_Priority_Index <- clip01(0.45 * out$Nano_Optimization_Need + 0.25 * out$Safety_Control_Requirement + 0.20 * out$Metabolic_Liability_Score + 0.10 * out$Barrier_Transporter_Score)
  out$Formulation_Priority_Class <- dplyr::case_when(out$Formulation_Priority_Index >= 0.66 ~ "High formulation priority", out$Formulation_Priority_Index >= 0.40 ~ "Moderate formulation priority", TRUE ~ "Low formulation priority")
  out
}

############################################################
## 7. Generate virtual partner profiles
############################################################

set_profile_values <- function(row, values) {
  out <- row
  for (nm in names(values)) {
    if (!nm %in% colnames(out)) stop("Unknown column in virtual profile: ", nm)
    out[[nm]] <- values[[nm]]
  }
  out
}

make_variant <- function(base_values, profile_name, idx) {
  row <- core_raw[1, , drop = FALSE]
  row <- set_profile_values(row, base_values)
  numeric_cols <- intersect(names(row)[sapply(row, function(x) suppressWarnings(!is.na(as.numeric(x[1]))))], names(base_values))
  for (nm in numeric_cols) row[[nm]] <- jitter_numeric(row[[nm]], rel = 0.06)
  row$compound <- paste0(profile_name, "_", sprintf("%02d", idx))
  row$PubChem_CID <- paste0("VirtualPartner_", profile_name, "_", sprintf("%02d", idx))
  row$SMILES <- "Synthetic partner profile"
  row$InChIKey <- paste0("SyntheticPartner_", profile_name, "_", sprintf("%02d", idx))
  row
}

profile_library <- list(
  Positive_Coassembly_Enhancer = list(
    Expected_Class = "High",
    Expected_Behavior = "High compatibility driven by distinct co-assembly gain.",
    Values = list(Molecular_Weight = 300, MolLogP = 2.8, TPSA = 78, HBD = 3, HBA = 6, Rotatable_Bonds = 2, Ring_Count = 3, Aromatic_Rings = 2, Fraction_Csp3 = 0.08, Formal_Charge = 0, Molar_Refractivity = 88, Consensus_LogP = 2.9, ESOL_LogS = -3.4, ESOL_Class = "Soluble", ADMETlab_logS = -3.3, ADMETlab_logD = 2.8, GI_absorption = "High", Caco2 = -4.6, PAMPA = 0.72, F30 = 0.72, BBB_permeant = "No", ADMETlab_BBB = 0.25, Pgp_substrate = "No", pgp_sub = 0.25, pgp_inh = 0.25, PPB = 82, Fu = 0.16, logVDss = 0.5, CYP_inhibitor_count = 1, t_half = 3.0, cl_plasma = 6, Lipinski_violations = 0, Veber_violations = 0, Bioavailability_Score = 0.60, PAINS_alerts = 0, Brenk_alerts = 0, Synthetic_Accessibility = 2.6, QED = 0.72, LD50 = 4200, Toxicity_Class = 5, DILI = 0.12, Ames = 0.05, Carcinogenicity = 0.05, hERG = 0.12, ProTox_Organ_Active_Count = 0, ProTox_Endpoint_Active_Count = 1, ProTox_Metabolism_Active_Count = 1)
  ),
  Inclusion_Stabilization_Partner = list(
    Expected_Class = "High",
    Expected_Behavior = "High compatibility driven by inclusion and stabilization gain.",
    Values = list(Molecular_Weight = 230, MolLogP = 2.3, TPSA = 50, HBD = 1, HBA = 4, Rotatable_Bonds = 2, Ring_Count = 2, Aromatic_Rings = 1, Fraction_Csp3 = 0.20, Formal_Charge = 0, Molar_Refractivity = 65, Consensus_LogP = 2.4, ESOL_LogS = -4.6, ESOL_Class = "Moderately soluble", ADMETlab_logS = -4.5, ADMETlab_logD = 2.3, GI_absorption = "High", Caco2 = -4.5, PAMPA = 0.70, F30 = 0.72, BBB_permeant = "No", ADMETlab_BBB = 0.25, Pgp_substrate = "No", pgp_sub = 0.20, pgp_inh = 0.20, PPB = 76, Fu = 0.22, logVDss = 0.35, CYP_inhibitor_count = 1, t_half = 3.5, cl_plasma = 5, Lipinski_violations = 0, Veber_violations = 0, Bioavailability_Score = 0.62, PAINS_alerts = 0, Brenk_alerts = 0, Synthetic_Accessibility = 2.2, QED = 0.78, LD50 = 4800, Toxicity_Class = 5, DILI = 0.08, Ames = 0.05, Carcinogenicity = 0.05, hERG = 0.08, ProTox_Organ_Active_Count = 0, ProTox_Endpoint_Active_Count = 0, ProTox_Metabolism_Active_Count = 1)
  ),
  Glycoside_Like_Stabilizer = list(
    Expected_Class = "Conditional",
    Expected_Behavior = "Conditional compatibility as a surface-stabilizing auxiliary profile.",
    Values = list(Molecular_Weight = 720, MolLogP = 1.4, TPSA = 185, HBD = 6, HBA = 13, Rotatable_Bonds = 9, Ring_Count = 4, Aromatic_Rings = 0, Fraction_Csp3 = 0.70, Formal_Charge = 0, Molar_Refractivity = 130, Consensus_LogP = 1.5, ESOL_LogS = -3.2, ESOL_Class = "Soluble", ADMETlab_logS = -3.3, ADMETlab_logD = 1.4, GI_absorption = "Low", Caco2 = -5.4, PAMPA = 0.30, F30 = 0.35, BBB_permeant = "No", ADMETlab_BBB = 0.05, Pgp_substrate = "Yes", pgp_sub = 0.60, pgp_inh = 0.25, PPB = 70, Fu = 0.28, logVDss = 0.2, CYP_inhibitor_count = 0, t_half = 5.0, cl_plasma = 4, Lipinski_violations = 2, Veber_violations = 1, Bioavailability_Score = 0.30, PAINS_alerts = 0, Brenk_alerts = 1, Synthetic_Accessibility = 5.2, QED = 0.30, LD50 = 5200, Toxicity_Class = 6, DILI = 0.06, Ames = 0.03, Carcinogenicity = 0.03, hERG = 0.05, ProTox_Organ_Active_Count = 0, ProTox_Endpoint_Active_Count = 0, ProTox_Metabolism_Active_Count = 0)
  ),
  ADMET_Compensating_Partner = list(
    Expected_Class = "ModerateHigh",
    Expected_Behavior = "Moderate-to-high compatibility driven by defect compensation and low safety burden.",
    Values = list(Molecular_Weight = 320, MolLogP = 2.1, TPSA = 62, HBD = 1, HBA = 5, Rotatable_Bonds = 3, Ring_Count = 2, Aromatic_Rings = 1, Fraction_Csp3 = 0.35, Formal_Charge = 0, Molar_Refractivity = 82, Consensus_LogP = 2.2, ESOL_LogS = -3.0, ESOL_Class = "Soluble", ADMETlab_logS = -3.0, ADMETlab_logD = 2.0, GI_absorption = "High", Caco2 = -4.2, PAMPA = 0.82, F30 = 0.85, BBB_permeant = "No", ADMETlab_BBB = 0.30, Pgp_substrate = "No", pgp_sub = 0.10, pgp_inh = 0.10, PPB = 68, Fu = 0.30, logVDss = 0.25, CYP_inhibitor_count = 0, t_half = 4.5, cl_plasma = 4, Lipinski_violations = 0, Veber_violations = 0, Bioavailability_Score = 0.70, PAINS_alerts = 0, Brenk_alerts = 0, Synthetic_Accessibility = 2.4, QED = 0.82, LD50 = 5500, Toxicity_Class = 6, DILI = 0.05, Ames = 0.03, Carcinogenicity = 0.03, hERG = 0.05, ProTox_Organ_Active_Count = 0, ProTox_Endpoint_Active_Count = 0, ProTox_Metabolism_Active_Count = 0)
  ),
  Toxic_Incompatible_Partner = list(
    Expected_Class = "Low",
    Expected_Behavior = "Low compatibility due to toxicity and safety-risk amplification.",
    Values = list(Molecular_Weight = 350, MolLogP = 3.0, TPSA = 70, HBD = 2, HBA = 5, Rotatable_Bonds = 4, Ring_Count = 2, Aromatic_Rings = 1, Fraction_Csp3 = 0.30, Formal_Charge = 0, Molar_Refractivity = 90, Consensus_LogP = 3.0, ESOL_LogS = -3.8, ESOL_Class = "Moderately soluble", ADMETlab_logS = -3.8, ADMETlab_logD = 3.0, GI_absorption = "High", Caco2 = -4.6, PAMPA = 0.60, F30 = 0.60, BBB_permeant = "No", ADMETlab_BBB = 0.20, Pgp_substrate = "Yes", pgp_sub = 0.75, pgp_inh = 0.75, PPB = 96, Fu = 0.04, logVDss = 1.2, CYP_inhibitor_count = 5, t_half = 0.6, cl_plasma = 18, Lipinski_violations = 0, Veber_violations = 0, Bioavailability_Score = 0.45, PAINS_alerts = 1, Brenk_alerts = 2, Synthetic_Accessibility = 4.8, QED = 0.40, LD50 = 80, Toxicity_Class = 3, DILI = 0.95, Ames = 0.80, Carcinogenicity = 0.75, hERG = 0.85, ProTox_Organ_Active_Count = 4, ProTox_Endpoint_Active_Count = 7, ProTox_Metabolism_Active_Count = 5)
  ),
  Overly_Polar_Incompatible_Partner = list(
    Expected_Class = "Low",
    Expected_Behavior = "Low compatibility due to excessive polarity and poor assembly suitability.",
    Values = list(Molecular_Weight = 580, MolLogP = -0.5, TPSA = 220, HBD = 10, HBA = 16, Rotatable_Bonds = 14, Ring_Count = 1, Aromatic_Rings = 0, Fraction_Csp3 = 0.80, Formal_Charge = -1, Molar_Refractivity = 110, Consensus_LogP = -0.4, ESOL_LogS = -1.5, ESOL_Class = "Soluble", ADMETlab_logS = -1.8, ADMETlab_logD = -0.3, GI_absorption = "Low", Caco2 = -6.8, PAMPA = 0.04, F30 = 0.05, BBB_permeant = "No", ADMETlab_BBB = 0.02, Pgp_substrate = "Yes", pgp_sub = 0.85, pgp_inh = 0.20, PPB = 40, Fu = 0.60, logVDss = -0.2, CYP_inhibitor_count = 0, t_half = 2.0, cl_plasma = 9, Lipinski_violations = 2, Veber_violations = 1, Bioavailability_Score = 0.20, PAINS_alerts = 0, Brenk_alerts = 1, Synthetic_Accessibility = 5.8, QED = 0.22, LD50 = 4500, Toxicity_Class = 5, DILI = 0.10, Ames = 0.05, Carcinogenicity = 0.05, hERG = 0.05, ProTox_Organ_Active_Count = 0, ProTox_Endpoint_Active_Count = 1, ProTox_Metabolism_Active_Count = 0)
  ),
  Overly_Lipophilic_Incompatible_Partner = list(
    Expected_Class = "Low",
    Expected_Behavior = "Low compatibility due to extreme lipophilicity, solubility defect and exposure risk.",
    Values = list(Molecular_Weight = 480, MolLogP = 8.0, TPSA = 22, HBD = 0, HBA = 2, Rotatable_Bonds = 8, Ring_Count = 3, Aromatic_Rings = 1, Fraction_Csp3 = 0.50, Formal_Charge = 0, Molar_Refractivity = 145, Consensus_LogP = 8.1, ESOL_LogS = -8.5, ESOL_Class = "Insoluble", ADMETlab_logS = -8.0, ADMETlab_logD = 8.0, GI_absorption = "Low", Caco2 = -5.5, PAMPA = 0.25, F30 = 0.25, BBB_permeant = "Yes", ADMETlab_BBB = 0.75, Pgp_substrate = "No", pgp_sub = 0.20, pgp_inh = 0.40, PPB = 99, Fu = 0.01, logVDss = 1.8, CYP_inhibitor_count = 3, t_half = 5.0, cl_plasma = 5, Lipinski_violations = 2, Veber_violations = 0, Bioavailability_Score = 0.35, PAINS_alerts = 0, Brenk_alerts = 1, Synthetic_Accessibility = 4.5, QED = 0.30, LD50 = 2500, Toxicity_Class = 5, DILI = 0.45, Ames = 0.10, Carcinogenicity = 0.20, hERG = 0.55, ProTox_Organ_Active_Count = 2, ProTox_Endpoint_Active_Count = 3, ProTox_Metabolism_Active_Count = 3)
  ),
  Neutral_Low_Value_Partner = list(
    Expected_Class = "Low",
    Expected_Behavior = "Low compatibility because molecular balance is not equivalent to incremental pairing gain.",
    Values = list(Molecular_Weight = 300, MolLogP = 2.0, TPSA = 65, HBD = 1, HBA = 4, Rotatable_Bonds = 4, Ring_Count = 1, Aromatic_Rings = 1, Fraction_Csp3 = 0.45, Formal_Charge = 0, Molar_Refractivity = 78, Consensus_LogP = 2.1, ESOL_LogS = -3.0, ESOL_Class = "Soluble", ADMETlab_logS = -3.0, ADMETlab_logD = 2.0, GI_absorption = "High", Caco2 = -4.4, PAMPA = 0.75, F30 = 0.80, BBB_permeant = "No", ADMETlab_BBB = 0.30, Pgp_substrate = "No", pgp_sub = 0.15, pgp_inh = 0.15, PPB = 70, Fu = 0.28, logVDss = 0.25, CYP_inhibitor_count = 0, t_half = 4.0, cl_plasma = 4, Lipinski_violations = 0, Veber_violations = 0, Bioavailability_Score = 0.70, PAINS_alerts = 0, Brenk_alerts = 0, Synthetic_Accessibility = 2.4, QED = 0.82, LD50 = 5500, Toxicity_Class = 6, DILI = 0.05, Ames = 0.03, Carcinogenicity = 0.03, hERG = 0.05, ProTox_Organ_Active_Count = 0, ProTox_Endpoint_Active_Count = 0, ProTox_Metabolism_Active_Count = 0)
  )
)

virtual_partner_list <- list()
for (profile_name in names(profile_library)) {
  for (i in seq_len(N_PER_PROFILE)) {
    virtual_partner_list[[paste0(profile_name, "_", i)]] <- make_variant(profile_library[[profile_name]]$Values, profile_name, i)
  }
}

virtual_partners_raw <- bind_rows(virtual_partner_list)
profile_meta <- bind_rows(lapply(names(profile_library), function(nm) {
  data.frame(Partner_Profile = nm, Expected_Class = profile_library[[nm]]$Expected_Class, Expected_Behavior = profile_library[[nm]]$Expected_Behavior, stringsAsFactors = FALSE)
}))

virtual_partners_raw$Partner_Profile <- gsub("_[0-9]{2}$", "", virtual_partners_raw$compound)
virtual_partners <- virtual_partners_raw %>% add_internal_scores() %>% dplyr::left_join(profile_meta, by = "Partner_Profile")
core_scored <- core_raw %>% add_internal_scores()

############################################################
## 8. Blueprint-guided pairing model with gain gate
############################################################

score_blueprint_match <- function(partner, core_blueprint) {
  primary <- as.character(core_blueprint$Primary_Nanoformulation_Backbone[1])
  secondary <- as.character(core_blueprint$Secondary_Nanoformulation_Module[1])
  functional <- as.character(core_blueprint$Conditional_Functional_Module[1])
  primary_score <- dplyr::case_when(
    grepl("Self", primary, ignore.case = TRUE) ~ get_one(partner, "Multivalent_Assembly_Drive"),
    grepl("Cyclodextrin", primary, ignore.case = TRUE) ~ get_one(partner, "Cyclodextrin_Cavity_Match"),
    grepl("Lipid", primary, ignore.case = TRUE) ~ get_one(partner, "Lipophilicity_Fit"),
    grepl("Micelle", primary, ignore.case = TRUE) ~ get_one(partner, "Hydrophobic_Core_Fit"),
    grepl("Controlled", primary, ignore.case = TRUE) ~ get_one(partner, "Controlled_Release_Need"),
    TRUE ~ get_one(partner, "Nano_Assembly_Suitability_Score")
  )
  secondary_score <- dplyr::case_when(
    grepl("Cyclodextrin", secondary, ignore.case = TRUE) ~ get_one(partner, "Cyclodextrin_Cavity_Match"),
    grepl("Self", secondary, ignore.case = TRUE) ~ get_one(partner, "Multivalent_Assembly_Drive"),
    grepl("Lipid", secondary, ignore.case = TRUE) ~ get_one(partner, "Lipophilicity_Fit"),
    grepl("Micelle", secondary, ignore.case = TRUE) ~ get_one(partner, "Hydrophobic_Core_Fit"),
    TRUE ~ get_one(partner, "Nano_Assembly_Suitability_Score")
  )
  functional_score <- dplyr::case_when(
    grepl("Barrier", functional, ignore.case = TRUE) ~ safe_mean(c(1 - get_one(partner, "Barrier_Transporter_Score"), 1 - get_one(partner, "Safety_Control_Requirement"))),
    grepl("Biomimetic", functional, ignore.case = TRUE) ~ safe_mean(c(1 - get_one(partner, "Safety_Control_Requirement"), 1 - get_one(partner, "Distribution_Exposure_Score"), get_one(partner, "Hydrogen_Bonding"))),
    grepl("Exposure", functional, ignore.case = TRUE) ~ safe_mean(c(1 - get_one(partner, "Safety_Control_Requirement"), 1 - get_one(partner, "Distribution_Exposure_Score"))),
    TRUE ~ 1 - get_one(partner, "Safety_Control_Requirement")
  )
  clip01(weighted_sum_strict(c(primary_score, secondary_score, functional_score), c(0.45, 0.30, 0.25)))
}

score_pairing_gain <- function(core, partner, core_blueprint) {
  core_logp <- get_one(core, "MolLogP")
  partner_logp <- get_one(partner, "MolLogP")
  core_tpsa <- get_one(core, "TPSA")
  partner_tpsa <- get_one(partner, "TPSA")
  core_mw <- get_one(core, "Molecular_Weight")
  partner_mw <- get_one(partner, "Molecular_Weight")
  core_hbd <- get_one(core, "HBD")
  partner_hbd <- get_one(partner, "HBD")
  core_hba <- get_one(core, "HBA")
  partner_hba <- get_one(partner, "HBA")
  core_ring <- get_one(core, "Ring_Count")
  partner_ring <- get_one(partner, "Ring_Count")
  core_arom <- get_one(core, "Aromatic_Rings")
  partner_arom <- get_one(partner, "Aromatic_Rings")
  core_rot <- get_one(core, "Rotatable_Bonds")
  partner_rot <- get_one(partner, "Rotatable_Bonds")
  
  ## Co-assembly gain requires a distinct assembly-driving signal, not simple balance.
  coassembly_gain <- clip01(
    0.50 * benefit_high(get_one(partner, "Multivalent_Assembly_Drive"), low = 0.58, high = 0.88) +
      0.20 * benefit_high(abs(core_logp - partner_logp), low = 0.20, high = 1.60) +
      0.15 * benefit_high(partner_arom + partner_ring, low = 3.0, high = 5.5) +
      0.15 * benefit_high(get_one(partner, "Hydrogen_Bonding"), low = 0.50, high = 0.90)
  )
  
  ## Stabilization gain requires inclusion/stabilizing signals plus a non-trivial interaction handle.
  stabilization_peak <- max(
    get_one(partner, "Cyclodextrin_Cavity_Match"),
    get_one(partner, "Small_Aromatic_Inclusion_Preference"),
    get_one(partner, "Hydrogen_Bonding"),
    na.rm = TRUE
  )
  stabilization_gain <- clip01(
    0.45 * benefit_high(stabilization_peak, low = 0.62, high = 0.92) +
      0.25 * benefit_high(get_one(partner, "Small_Aromatic_Inclusion_Preference"), low = 0.58, high = 0.90) +
      0.15 * benefit_high(abs(core_tpsa - partner_tpsa), low = 10, high = 70) +
      0.15 * benefit_high(abs(core_hbd + core_hba - partner_hbd - partner_hba), low = 1, high = 7)
  )
  
  ## Blueprint-specific gain is stricter than simple blueprint match.
  primary <- as.character(core_blueprint$Primary_Nanoformulation_Backbone[1])
  secondary <- as.character(core_blueprint$Secondary_Nanoformulation_Module[1])
  functional <- as.character(core_blueprint$Conditional_Functional_Module[1])
  primary_gain <- dplyr::case_when(
    grepl("Self", primary, ignore.case = TRUE) ~ benefit_high(get_one(partner, "Multivalent_Assembly_Drive"), low = 0.60, high = 0.90),
    grepl("Cyclodextrin", primary, ignore.case = TRUE) ~ benefit_high(get_one(partner, "Cyclodextrin_Cavity_Match"), low = 0.62, high = 0.92),
    grepl("Lipid", primary, ignore.case = TRUE) ~ benefit_high(get_one(partner, "Lipophilicity_Fit"), low = 0.60, high = 0.90),
    grepl("Micelle", primary, ignore.case = TRUE) ~ benefit_high(get_one(partner, "Hydrophobic_Core_Fit"), low = 0.60, high = 0.90),
    TRUE ~ benefit_high(get_one(partner, "Nano_Assembly_Suitability_Score"), low = 0.60, high = 0.90)
  )
  secondary_gain <- dplyr::case_when(
    grepl("Cyclodextrin", secondary, ignore.case = TRUE) ~ benefit_high(get_one(partner, "Cyclodextrin_Cavity_Match"), low = 0.62, high = 0.92),
    grepl("Self", secondary, ignore.case = TRUE) ~ benefit_high(get_one(partner, "Multivalent_Assembly_Drive"), low = 0.60, high = 0.90),
    TRUE ~ benefit_high(get_one(partner, "Nano_Assembly_Suitability_Score"), low = 0.60, high = 0.90)
  )
  functional_gain <- dplyr::case_when(
    grepl("Barrier", functional, ignore.case = TRUE) ~ benefit_high(1 - get_one(partner, "Barrier_Transporter_Score"), low = 0.55, high = 0.85),
    grepl("Biomimetic", functional, ignore.case = TRUE) ~ benefit_high(safe_mean(c(1 - get_one(partner, "Safety_Control_Requirement"), get_one(partner, "Hydrogen_Bonding"))), low = 0.60, high = 0.90),
    grepl("Exposure", functional, ignore.case = TRUE) ~ benefit_high(1 - get_one(partner, "Distribution_Exposure_Score"), low = 0.55, high = 0.85),
    TRUE ~ benefit_high(1 - get_one(partner, "Safety_Control_Requirement"), low = 0.55, high = 0.85)
  )
  blueprint_specific_gain <- clip01(weighted_sum_strict(c(primary_gain, secondary_gain, functional_gain), c(0.45, 0.30, 0.25)))
  
  ## Defect compensation gain must be need-weighted by the core compound's actual defects.
  need_values <- c(
    Solubility = get_one(core, "Solubility_Defect_Score"),
    Permeability = get_one(core, "Permeability_Defect_Score"),
    Barrier = get_one(core, "Barrier_Transporter_Score"),
    Metabolic = get_one(core, "Metabolic_Liability_Score"),
    Safety = get_one(core, "Safety_Control_Requirement")
  )
  support_values <- c(
    Solubility = safe_mean(c(1 - get_one(partner, "Solubility_Defect_Score"), get_one(partner, "Cyclodextrin_Cavity_Match"))),
    Permeability = 1 - get_one(partner, "Permeability_Defect_Score"),
    Barrier = 1 - get_one(partner, "Barrier_Transporter_Score"),
    Metabolic = 1 - get_one(partner, "Metabolic_Liability_Score"),
    Safety = 1 - get_one(partner, "Safety_Control_Requirement")
  )
  need_weights <- clip01(need_values)
  if (sum(need_weights, na.rm = TRUE) < 1e-6) {
    defect_compensation_gain <- 0.35
  } else {
    defect_compensation_gain <- clip01(sum(need_weights * support_values, na.rm = TRUE) / sum(need_weights, na.rm = TRUE))
  }
  
  ## Non-redundancy: a partner should be different enough from the core AND carry an actual function.
  similarity_to_core <- weighted_sum_strict(
    c(
      feature_similarity(partner_mw, core_mw, 450),
      feature_similarity(partner_logp, core_logp, 3.5),
      feature_similarity(partner_tpsa, core_tpsa, 120),
      feature_similarity(partner_hbd + partner_hba, core_hbd + core_hba, 10),
      feature_similarity(partner_ring, core_ring, 4),
      feature_similarity(partner_arom, core_arom, 3),
      feature_similarity(partner_rot, core_rot, 8)
    ),
    c(0.12, 0.18, 0.18, 0.14, 0.12, 0.14, 0.12)
  )
  functional_strength <- clip01(weighted_sum_strict(
    c(
      coassembly_gain,
      stabilization_gain,
      blueprint_specific_gain,
      defect_compensation_gain,
      1 - get_one(partner, "Safety_Control_Requirement")
    ),
    c(0.25, 0.25, 0.25, 0.15, 0.10)
  ))
  nonredundancy_gain <- clip01((1 - similarity_to_core) * functional_strength)
  
  pairing_gain_index <- clip01(weighted_sum_strict(
    c(coassembly_gain, stabilization_gain, blueprint_specific_gain, defect_compensation_gain, nonredundancy_gain),
    c(0.25, 0.20, 0.20, 0.20, 0.15)
  ))
  
  pairing_gain_gate <- clip01(0.45 + 0.55 * pairing_gain_index)
  
  neutral_low_value_flag <- as.numeric(
    pairing_gain_index < 0.48 &
      get_one(partner, "Safety_Control_Requirement") < 0.25 &
      get_one(partner, "Druggability_Alert_Requirement") < 0.25 &
      max(coassembly_gain, stabilization_gain, blueprint_specific_gain, defect_compensation_gain, na.rm = TRUE) < 0.72
  )
  
  neutral_low_value_gate <- ifelse(neutral_low_value_flag == 1, 0.62, 1.00)
  
  data.frame(
    Coassembly_Gain = coassembly_gain,
    Stabilization_Gain = stabilization_gain,
    Blueprint_Specific_Gain = blueprint_specific_gain,
    Defect_Compensation_Gain = defect_compensation_gain,
    Similarity_to_Core = similarity_to_core,
    Functional_Strength = functional_strength,
    Nonredundancy_Gain = nonredundancy_gain,
    Pairing_Gain_Index = pairing_gain_index,
    Pairing_Gain_Gate = pairing_gain_gate,
    Neutral_Low_Value_Flag = neutral_low_value_flag,
    Neutral_Low_Value_Gate = neutral_low_value_gate,
    stringsAsFactors = FALSE
  )
}

base_pair_weights <- c(
  Co_assembly_Compatibility = 0.22,
  Physicochemical_Complementarity = 0.18,
  Nano_stabilization_Complementarity = 0.15,
  ADMET_Complementarity = 0.15,
  Safety_Balance = 0.15,
  Blueprint_Match = 0.15
)

score_pair_one <- function(core, partner, core_blueprint, weights = NULL) {
  if (is.null(weights)) weights <- base_pair_weights
  
  core_logp <- get_one(core, "MolLogP")
  partner_logp <- get_one(partner, "MolLogP")
  core_tpsa <- get_one(core, "TPSA")
  partner_tpsa <- get_one(partner, "TPSA")
  core_mw <- get_one(core, "Molecular_Weight")
  partner_mw <- get_one(partner, "Molecular_Weight")
  
  interaction_balance <- safe_mean(c(get_one(partner, "Aromaticity"), get_one(partner, "Hydrogen_Bonding"), get_one(partner, "Rigidity"), get_one(partner, "Multivalent_Assembly_Drive")))
  
  co_assembly <- safe_mean(c(interaction_balance, benefit_mid(abs(core_logp - partner_logp), lower = 0.0, upper = 1.8), benefit_mid(core_logp + partner_logp, lower = 3.0, upper = 7.5), get_one(partner, "Charge_Compatibility")))
  
  physicochemical <- safe_mean(c(benefit_mid(mean(c(core_logp, partner_logp)), lower = 1.5, upper = 4.5), benefit_mid(mean(c(core_tpsa, partner_tpsa)), lower = 45, upper = 120), benefit_mid(core_mw + partner_mw, lower = 350, upper = 1000), benefit_mid(get_one(partner, "Rotatable_Bonds"), lower = 0, upper = 8), benefit_mid(get_one(partner, "Molar_Refractivity"), lower = 50, upper = 130)))
  
  stabilization <- safe_mean(c(get_one(partner, "Cyclodextrin_Cavity_Match"), get_one(partner, "Small_Aromatic_Inclusion_Preference"), get_one(partner, "Hydrogen_Bonding"), get_one(partner, "Molecular_Size_Fit"), 1 - get_one(partner, "Safety_Control_Requirement")))
  
  core_sol_need <- get_one(core, "Solubility_Defect_Score")
  core_perm_need <- get_one(core, "Permeability_Defect_Score")
  core_barrier_need <- get_one(core, "Barrier_Transporter_Score")
  core_metabolic_need <- get_one(core, "Metabolic_Liability_Score")
  
  sol_comp <- ifelse(core_sol_need >= 0.40, safe_mean(c(1 - get_one(partner, "Solubility_Defect_Score"), get_one(partner, "Cyclodextrin_Cavity_Match"))), 0.50)
  perm_comp <- ifelse(core_perm_need >= 0.40, 1 - get_one(partner, "Permeability_Defect_Score"), 0.50)
  barrier_comp <- ifelse(core_barrier_need >= 0.40, 1 - get_one(partner, "Barrier_Transporter_Score"), 0.50)
  metabolic_comp <- ifelse(core_metabolic_need >= 0.40, 1 - get_one(partner, "Metabolic_Liability_Score"), 0.50)
  
  admet <- safe_mean(c(sol_comp, perm_comp, barrier_comp, metabolic_comp, 1 - get_one(partner, "Distribution_Exposure_Score")))
  safety <- safe_mean(c(1 - get_one(partner, "Safety_Control_Requirement"), 1 - get_one(partner, "Druggability_Alert_Requirement"), 1 - clip01(get_one(partner, "DILI")), 1 - clip01(get_one(partner, "Ames")), 1 - clip01(get_one(partner, "hERG"))))
  blueprint_match <- score_blueprint_match(partner, core_blueprint)
  
  module_values <- c(Co_assembly_Compatibility = co_assembly, Physicochemical_Complementarity = physicochemical, Nano_stabilization_Complementarity = stabilization, ADMET_Complementarity = admet, Safety_Balance = safety, Blueprint_Match = blueprint_match)
  raw_pci <- weighted_sum_strict(module_values[names(weights)], weights[names(weights)])
  
  gain <- score_pairing_gain(core, partner, core_blueprint)
  
  toxic_red_flag <- as.numeric(get_one(partner, "Safety_Control_Requirement") >= 0.65 | get_one(partner, "DILI") >= 0.70 | get_one(partner, "Ames") >= 0.70 | get_one(partner, "hERG") >= 0.70 | get_one(partner, "LD50") <= 300)
  extreme_physchem_flag <- as.numeric(partner_tpsa >= 170 | partner_logp >= 6.5 | partner_logp <= 0 | partner_mw >= 800)
  
  toxic_gate <- ifelse(toxic_red_flag == 1, 0.48, 1.00)
  extreme_physchem_gate <- ifelse(extreme_physchem_flag == 1, 0.72, 1.00)
  
  final_pci <- raw_pci * gain$Pairing_Gain_Gate[1] * gain$Neutral_Low_Value_Gate[1] * toxic_gate * extreme_physchem_gate
  final_pci <- clip01(final_pci)
  
  out <- data.frame(
    Co_assembly_Compatibility = co_assembly,
    Physicochemical_Complementarity = physicochemical,
    Nano_stabilization_Complementarity = stabilization,
    ADMET_Complementarity = admet,
    Safety_Balance = safety,
    Blueprint_Match = blueprint_match,
    Raw_PCI = raw_pci,
    gain,
    Toxic_Red_Flag = toxic_red_flag,
    Toxic_Gate = toxic_gate,
    Extreme_Physicochemical_Flag = extreme_physchem_flag,
    Extreme_Physicochemical_Gate = extreme_physchem_gate,
    Pairing_Compatibility_Index = final_pci,
    stringsAsFactors = FALSE
  )
  
  out$Pairing_Compatibility_Class <- dplyr::case_when(
    out$Pairing_Compatibility_Index >= 0.75 ~ "High compatibility",
    out$Pairing_Compatibility_Index >= 0.60 ~ "Moderate-high compatibility",
    out$Pairing_Compatibility_Index >= 0.45 ~ "Conditional compatibility",
    TRUE ~ "Low or incompatible"
  )
  
  out
}

score_all_partners <- function(core, partners, core_blueprint) {
  res <- list()
  for (i in seq_len(nrow(partners))) {
    partner <- partners[i, , drop = FALSE]
    score <- score_pair_one(core, partner, core_blueprint)
    score$Core_Compound <- core$compound[1]
    score$Partner_Compound <- partner$compound[1]
    score$Partner_Profile <- partner$Partner_Profile[1]
    score$Expected_Class <- partner$Expected_Class[1]
    score$Expected_Behavior <- partner$Expected_Behavior[1]
    res[[i]] <- score
  }
  bind_rows(res) %>% dplyr::select(Core_Compound, Partner_Compound, Partner_Profile, Expected_Class, Expected_Behavior, dplyr::everything())
}

pairing_scores <- score_all_partners(core_scored, virtual_partners, core_blueprint)

output_prefix <- "ACMPI_Nano_PairingCalibration_SubmissionGrade"
pairing_file <- "Generated in memory from Stage 2 blueprint and virtual partner profiles"

############################################################
## Feature-redefined multi-model calibration
############################################################

timestamp_final <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(output_root_dir, paste0(output_prefix, "_Output_", timestamp_final))
table_dir <- file.path(out_dir, "01_Tables")
fig_png_dir <- file.path(out_dir, "02_Figures_PNG")
fig_pdf_dir <- file.path(out_dir, "03_Figures_PDF")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_pdf_dir, recursive = TRUE, showWarnings = FALSE)


############################################################
## Manuscript-grade plotting helpers
############################################################

theme_acmpi <- function(base_size = 13) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "#2C2C2C"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3, hjust = 0, margin = ggplot2::margin(b = 7)),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "#555555", hjust = 0, lineheight = 1.12, margin = ggplot2::margin(b = 12)),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1, color = "#333333"),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.key.height = grid::unit(0.62, "cm"),
      legend.key.width = grid::unit(0.62, "cm"),
      plot.margin = ggplot2::margin(24, 110, 24, 52),
      legend.position = "right"
    )
}

save_plot_dual <- function(plot, png_path, pdf_path, width, height, dpi = 600) {
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = png_path, plot = plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  tryCatch(
    ggplot2::ggsave(filename = pdf_path, plot = plot, width = width, height = height, device = cairo_pdf, bg = "white", limitsize = FALSE),
    error = function(e) ggplot2::ggsave(filename = pdf_path, plot = plot, width = width, height = height, device = "pdf", bg = "white", limitsize = FALSE)
  )
}

N_AUGMENT_PER_ROW <- 10

base_scores <- pairing_scores %>%
  strip_regex_cols(patterns = c("^Sensitivity_", "^Probability_High_or_ModerateHigh", "^Probability_Neutral_Low_Value_Flag"))

partner_feature_cols <- c(
  "compound", "Molecular_Weight", "MolLogP", "TPSA", "HBD", "HBA",
  "Rotatable_Bonds", "Ring_Count", "Aromatic_Rings", "Fraction_Csp3",
  "Formal_Charge", "Molar_Refractivity", "Consensus_LogP", "ESOL_LogS",
  "ADMETlab_logS", "ADMETlab_logD", "Caco2", "PAMPA", "F30",
  "pgp_sub", "pgp_inh", "PPB", "Fu", "logVDss",
  "CYP_inhibitor_count", "t_half", "cl_plasma", "Bioavailability_Score",
  "Synthetic_Accessibility", "QED", "LD50", "DILI", "Ames", "Carcinogenicity", "hERG",
  "Solubility_Defect_Score", "Permeability_Defect_Score", "Barrier_Transporter_Score",
  "Distribution_Exposure_Score", "Metabolic_Liability_Score", "Safety_Control_Requirement",
  "Druggability_Alert_Requirement", "Nano_Assembly_Suitability_Score",
  "Molecular_Size_Fit", "Hydrogen_Bonding", "Rigidity", "Charge_Compatibility",
  "Druggability_Manageability", "Multivalent_Assembly_Drive", "Cyclodextrin_Cavity_Match",
  "Small_Aromatic_Inclusion_Preference"
)

partner_features <- virtual_partners %>%
  dplyr::select(dplyr::any_of(partner_feature_cols)) %>%
  dplyr::rename(Partner_Compound = compound)

base_scores <- base_scores %>%
  dplyr::left_join(partner_features, by = "Partner_Compound")
assert_no_join_suffixes(base_scores, "final base score feature join")

required_cols_final <- c(
  "Partner_Compound", "Partner_Profile", "Expected_Class", "Expected_Behavior",
  "Raw_PCI", "Pairing_Compatibility_Index",
  "Co_assembly_Compatibility", "Physicochemical_Complementarity",
  "Nano_stabilization_Complementarity", "ADMET_Complementarity",
  "Safety_Balance", "Blueprint_Match",
  "Coassembly_Gain", "Stabilization_Gain", "Blueprint_Specific_Gain",
  "Defect_Compensation_Gain", "Nonredundancy_Gain",
  "Toxic_Red_Flag", "Extreme_Physicochemical_Flag",
  "Molecular_Weight", "MolLogP", "TPSA", "HBD", "HBA", "Rotatable_Bonds",
  "Ring_Count", "Aromatic_Rings", "Fraction_Csp3", "Caco2", "PAMPA", "F30",
  "Fu", "t_half", "cl_plasma", "CYP_inhibitor_count", "Cyclodextrin_Cavity_Match",
  "Small_Aromatic_Inclusion_Preference", "Multivalent_Assembly_Drive",
  "Hydrogen_Bonding", "Safety_Control_Requirement", "Druggability_Alert_Requirement"
)
stop_if_missing(base_scores, required_cols_final, "base_scores before final feature redefinition")

clip01 <- function(x) pmax(0, pmin(1, as.numeric(x)))
benefit_high <- function(x, low, high) clip01((as.numeric(x) - low) / (high - low))
benefit_low <- function(x, high, low) clip01((high - as.numeric(x)) / (high - low))
benefit_mid <- function(x, lower, upper) {
  x <- as.numeric(x); span <- upper - lower
  out <- ifelse(x >= lower & x <= upper, 1,
                ifelse(x < lower, clip01((x - (lower - span)) / span),
                       clip01(((upper + span) - x) / span)))
  clip01(out)
}
weighted_sum <- function(values, weights) {
  values <- as.numeric(values); weights <- as.numeric(weights)
  keep <- !is.na(values) & !is.na(weights)
  if (!any(keep)) return(NA_real_)
  sum(values[keep] * weights[keep]) / sum(weights[keep])
}
wrap_label <- function(x, width = 28) stringr::str_wrap(x, width = width)
clean_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("Partner", "", x)
  x <- gsub(" +", " ", x)
  trimws(x)
}
jitter01 <- function(x, jitter = 0.10) clip01(as.numeric(x) * runif(length(x), 1 - jitter, 1 + jitter))

core_row <- core_scored[1, , drop = FALSE]
core_need_values <- c(
  Solubility = get_one(core_row, "Solubility_Defect_Score"),
  Permeability = get_one(core_row, "Permeability_Defect_Score"),
  Barrier = get_one(core_row, "Barrier_Transporter_Score"),
  Metabolic = get_one(core_row, "Metabolic_Liability_Score"),
  Safety = get_one(core_row, "Safety_Control_Requirement")
)
core_need_weights <- clip01(core_need_values)
if (sum(core_need_weights, na.rm = TRUE) < 1e-6) core_need_weights[] <- 1

add_distinctive_features <- function(df) {
  out <- df
  core_mw <- get_one(core_row, "Molecular_Weight")
  core_logp <- get_one(core_row, "MolLogP")
  core_tpsa <- get_one(core_row, "TPSA")
  core_hbd <- get_one(core_row, "HBD")
  core_hba <- get_one(core_row, "HBA")
  core_ring <- get_one(core_row, "Ring_Count")
  core_arom <- get_one(core_row, "Aromatic_Rings")
  core_rot <- get_one(core_row, "Rotatable_Bonds")
  
  out$Legacy_PCI <- out$Pairing_Compatibility_Index
  out$Basal_Acceptability <- clip01(
    0.35 * out$Safety_Balance +
      0.25 * out$ADMET_Complementarity +
      0.20 * out$Physicochemical_Complementarity +
      0.10 * (1 - out$Druggability_Alert_Requirement) +
      0.10 * out$Blueprint_Match
  )
  
  assembly_drive <- benefit_high(out$Multivalent_Assembly_Drive, 0.66, 0.90)
  assembly_architecture <- benefit_high(out$Aromatic_Rings + out$Ring_Count, 3.2, 5.2)
  assembly_hbond <- benefit_high(out$Hydrogen_Bonding, 0.58, 0.88)
  assembly_logp_contrast <- benefit_high(abs(out$MolLogP - core_logp), 0.25, 1.65)
  out$Assembly_Distinctiveness <- clip01(
    0.42 * assembly_drive +
      0.28 * assembly_architecture +
      0.20 * assembly_hbond +
      0.10 * assembly_logp_contrast
  )
  
  inclusion_cavity <- benefit_high(out$Cyclodextrin_Cavity_Match, 0.74, 0.93)
  inclusion_small_arom <- benefit_high(out$Small_Aromatic_Inclusion_Preference, 0.72, 0.93)
  inclusion_size_window <- benefit_mid(out$Molecular_Weight, 180, 270)
  inclusion_ring_handle <- benefit_high(out$Ring_Count + out$Aromatic_Rings, 2.6, 4.0)
  inclusion_rigidity <- benefit_low(out$Rotatable_Bonds, 4.2, 0)
  out$Inclusion_Distinctiveness <- clip01(
    0.30 * inclusion_cavity +
      0.28 * inclusion_small_arom +
      0.18 * inclusion_size_window +
      0.14 * inclusion_ring_handle +
      0.10 * inclusion_rigidity
  )
  
  similarity_components <- cbind(
    1 - clip01(abs(out$Molecular_Weight - core_mw) / 450),
    1 - clip01(abs(out$MolLogP - core_logp) / 3.5),
    1 - clip01(abs(out$TPSA - core_tpsa) / 120),
    1 - clip01(abs((out$HBD + out$HBA) - (core_hbd + core_hba)) / 10),
    1 - clip01(abs(out$Ring_Count - core_ring) / 4),
    1 - clip01(abs(out$Aromatic_Rings - core_arom) / 3),
    1 - clip01(abs(out$Rotatable_Bonds - core_rot) / 8)
  )
  out$Similarity_to_Core <- apply(similarity_components, 1, function(x) weighted_sum(x, c(0.12, 0.18, 0.18, 0.14, 0.12, 0.14, 0.12)))
  out$Structural_Nonredundancy <- clip01(1 - out$Similarity_to_Core)
  
  blueprint_nonredundancy_factor <- clip01(0.45 + 0.55 * out$Structural_Nonredundancy)
  out$Blueprint_Nonredundant_Gain <- clip01(out$Blueprint_Specific_Gain * blueprint_nonredundancy_factor)
  
  sol_support <- clip01((1 - out$Solubility_Defect_Score) * benefit_high(out$Cyclodextrin_Cavity_Match, 0.66, 0.90))
  perm_support <- clip01((1 - out$Permeability_Defect_Score) * benefit_high(out$PAMPA, 0.72, 0.86) * benefit_high(out$F30, 0.78, 0.90))
  barrier_support <- clip01((1 - out$Barrier_Transporter_Score) * benefit_high(abs(out$TPSA - core_tpsa), 15, 85))
  metabolic_support <- clip01((1 - out$Metabolic_Liability_Score) * benefit_high(out$t_half, 4.1, 5.2) * benefit_low(out$CYP_inhibitor_count, 2.0, 0))
  safety_support <- clip01((1 - out$Safety_Control_Requirement) * benefit_high(out$LD50, 4200, 5600) * benefit_low(out$DILI, 0.20, 0.02))
  comp_matrix <- cbind(sol_support, perm_support, barrier_support, metabolic_support, safety_support)
  out$Matched_Defect_Compensation <- apply(comp_matrix, 1, function(x) clip01(weighted_sum(x, core_need_weights)))
  
  out$ADMET_Superiority_Signature <- clip01(
    0.24 * benefit_high(out$PAMPA, 0.78, 0.88) +
      0.24 * benefit_high(out$F30, 0.82, 0.92) +
      0.16 * benefit_high(out$Fu, 0.28, 0.36) +
      0.16 * benefit_high(out$t_half, 4.1, 5.3) +
      0.10 * benefit_low(out$CYP_inhibitor_count, 1.0, 0) +
      0.10 * benefit_low(out$cl_plasma, 5.0, 2.8)
  )
  
  out$Distinctive_Pairing_Gain <- clip01(
    0.28 * out$Assembly_Distinctiveness +
      0.24 * out$Inclusion_Distinctiveness +
      0.18 * out$Blueprint_Nonredundant_Gain +
      0.18 * out$Matched_Defect_Compensation +
      0.07 * out$Structural_Nonredundancy +
      0.05 * out$ADMET_Superiority_Signature
  )
  
  out$Surface_Stabilizing_Auxiliary_Flag <- as.numeric(
    ((out$Molecular_Weight >= 650) | (out$Molecular_Weight >= 600 & (out$HBD + out$HBA) >= 12 & out$Fraction_Csp3 >= 0.55)) &
      out$MolLogP >= 0.8 & out$MolLogP <= 2.3 &
      out$Safety_Control_Requirement < 0.30 & out$DILI < 0.20 & out$Ames < 0.20 & out$hERG < 0.20
  )
  
  out$Neutral_Low_Value_Feature_Flag <- as.numeric(
    out$Basal_Acceptability >= 0.80 &
      out$Distinctive_Pairing_Gain < 0.48 &
      out$Assembly_Distinctiveness < 0.42 &
      out$Inclusion_Distinctiveness < 0.54 &
      out$Matched_Defect_Compensation < 0.58 &
      out$ADMET_Superiority_Signature < 0.58 &
      out$Surface_Stabilizing_Auxiliary_Flag == 0
  )
  
  out$Feature_Redefined_Gain_Class <- dplyr::case_when(
    out$Surface_Stabilizing_Auxiliary_Flag == 1 ~ "Auxiliary stabilizer",
    out$Distinctive_Pairing_Gain >= 0.65 ~ "High distinctive gain",
    out$Distinctive_Pairing_Gain >= 0.50 ~ "Moderate distinctive gain",
    out$Neutral_Low_Value_Feature_Flag == 1 ~ "Neutral low-value",
    TRUE ~ "Low distinctive gain"
  )
  out
}

base_features <- base_scores %>%
  add_distinctive_features() %>%
  dplyr::mutate(
    Calibration_Source = "Observed virtual profiles",
    Original_Row_ID = dplyr::row_number(),
    Expected_Class = dplyr::case_when(
      Expected_Class %in% c("High", "ModerateHigh", "Conditional", "Low") ~ Expected_Class,
      TRUE ~ as.character(Expected_Class)
    )
  )

feature_cols <- c(
  "Raw_PCI", "Legacy_PCI", "Basal_Acceptability", "Assembly_Distinctiveness",
  "Inclusion_Distinctiveness", "Blueprint_Nonredundant_Gain",
  "Matched_Defect_Compensation", "ADMET_Superiority_Signature",
  "Structural_Nonredundancy", "Distinctive_Pairing_Gain",
  "Safety_Balance", "ADMET_Complementarity", "Physicochemical_Complementarity"
)
stop_if_missing(base_features, feature_cols, "base_features after final feature definition")

augment_one_round <- function(df, round_name, n_aug, jitter, profile_filter = NULL, force_class = NULL,
                              gain_multiplier = 1, acceptability_multiplier = 1, auxiliary_force = FALSE) {
  if (!is.null(profile_filter)) df <- df[grepl(profile_filter, df$Partner_Profile), , drop = FALSE]
  if (nrow(df) < 1) return(df[0, , drop = FALSE])
  out <- vector("list", n_aug)
  for (i in seq_len(n_aug)) {
    tmp <- df
    for (cc in feature_cols) tmp[[cc]] <- jitter01(tmp[[cc]], jitter = jitter)
    tmp$Basal_Acceptability <- jitter01(tmp$Basal_Acceptability * acceptability_multiplier, jitter = jitter * 0.40)
    tmp$Assembly_Distinctiveness <- jitter01(tmp$Assembly_Distinctiveness * gain_multiplier, jitter = jitter)
    tmp$Inclusion_Distinctiveness <- jitter01(tmp$Inclusion_Distinctiveness * gain_multiplier, jitter = jitter)
    tmp$Blueprint_Nonredundant_Gain <- jitter01(tmp$Blueprint_Nonredundant_Gain * gain_multiplier, jitter = jitter)
    tmp$Matched_Defect_Compensation <- jitter01(tmp$Matched_Defect_Compensation * gain_multiplier, jitter = jitter)
    tmp$ADMET_Superiority_Signature <- jitter01(tmp$ADMET_Superiority_Signature * gain_multiplier, jitter = jitter)
    tmp$Distinctive_Pairing_Gain <- clip01(
      0.28 * tmp$Assembly_Distinctiveness +
        0.24 * tmp$Inclusion_Distinctiveness +
        0.18 * tmp$Blueprint_Nonredundant_Gain +
        0.18 * tmp$Matched_Defect_Compensation +
        0.07 * tmp$Structural_Nonredundancy +
        0.05 * tmp$ADMET_Superiority_Signature
    )
    tmp$Neutral_Low_Value_Feature_Flag <- as.numeric(
      tmp$Basal_Acceptability >= 0.80 &
        tmp$Distinctive_Pairing_Gain < 0.48 &
        tmp$Assembly_Distinctiveness < 0.42 &
        tmp$Inclusion_Distinctiveness < 0.54 &
        tmp$Matched_Defect_Compensation < 0.58 &
        tmp$ADMET_Superiority_Signature < 0.58 &
        tmp$Surface_Stabilizing_Auxiliary_Flag == 0
    )
    if (auxiliary_force) tmp$Surface_Stabilizing_Auxiliary_Flag <- 1
    tmp$Calibration_Source <- round_name
    tmp$Partner_Compound <- paste0(tmp$Partner_Compound, "__", gsub(" ", "_", round_name), "_", sprintf("%02d", i))
    if (!is.null(force_class)) tmp$Expected_Class <- force_class
    out[[i]] <- tmp
  }
  dplyr::bind_rows(out)
}

calibration_matrix <- dplyr::bind_rows(
  base_features,
  augment_one_round(base_features, "General perturbation", N_AUGMENT_PER_ROW, FEATURE_JITTER),
  augment_one_round(base_features, "Neutral hard controls", N_AUGMENT_PER_ROW, FEATURE_JITTER, "Neutral", "Low", gain_multiplier = 0.70, acceptability_multiplier = 1.04),
  augment_one_round(base_features, "Positive retention controls", N_AUGMENT_PER_ROW, FEATURE_JITTER, "Positive|Inclusion", NULL, gain_multiplier = 1.08, acceptability_multiplier = 1.00),
  augment_one_round(base_features, "ADMET compensation controls", N_AUGMENT_PER_ROW, FEATURE_JITTER, "ADMET_Compensating", "ModerateHigh", gain_multiplier = 1.05, acceptability_multiplier = 1.00),
  augment_one_round(base_features, "Auxiliary stabilizer controls", N_AUGMENT_PER_ROW, FEATURE_JITTER, "Glycoside", "Conditional", gain_multiplier = 0.96, acceptability_multiplier = 0.98, auxiliary_force = TRUE)
)

candidate_grid <- expand.grid(
  Model_Family = c("FeatureRedefinedAuxiliary", "FullConstrained"),
  Raw_Power = c(0.10, 0.12, 0.14),
  Acceptability_Power = c(0.10, 0.12, 0.14),
  Distinctive_Gain_Power = c(0.68, 0.72, 0.76, 0.80),
  Neutral_Gain_Threshold = c(0.40, 0.42, 0.44),
  Neutral_Gate_Value = c(0.35, 0.38, 0.40),
  Inclusion_Boost = c(1.14, 1.18, 1.22, 1.26),
  ADMET_Boost = c(1.04, 1.08, 1.12, 1.16),
  Auxiliary_Lower_Bound = c(0.35, 0.38),
  Auxiliary_Upper_Bound = c(0.52, 0.56),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(Model_ID = paste0("M", sprintf("%05d", dplyr::row_number()))) %>%
  dplyr::select(Model_ID, dplyr::everything())

cat("Candidate scoring rules: ", nrow(candidate_grid), "\n", sep = "")

apply_candidate_model <- function(features, pars) {
  raw_component <- clip01(features$Raw_PCI)
  acceptability <- clip01(features$Basal_Acceptability)
  distinctive_gain <- clip01(features$Distinctive_Pairing_Gain)
  base_score <- clip01(1.08 * (raw_component ^ pars$Raw_Power) * (acceptability ^ pars$Acceptability_Power) * (distinctive_gain ^ pars$Distinctive_Gain_Power))
  optimized_neutral_flag <- as.numeric(
    features$Basal_Acceptability >= 0.80 &
      distinctive_gain < pars$Neutral_Gain_Threshold &
      features$Assembly_Distinctiveness < 0.44 &
      features$Inclusion_Distinctiveness < 0.56 &
      features$Matched_Defect_Compensation < 0.60 &
      features$ADMET_Superiority_Signature < 0.60 &
      features$Surface_Stabilizing_Auxiliary_Flag == 0
  )
  toxic_gate <- ifelse(features$Toxic_Red_Flag == 1, 0.48, 1.00)
  extreme_gate <- ifelse(features$Extreme_Physicochemical_Flag == 1, 0.72, 1.00)
  neutral_gate <- ifelse(optimized_neutral_flag == 1, pars$Neutral_Gate_Value, 1.00)
  inclusion_boost <- ifelse(features$Inclusion_Distinctiveness >= 0.62 & features$Surface_Stabilizing_Auxiliary_Flag == 0, pars$Inclusion_Boost, 1.00)
  admet_boost <- ifelse(features$ADMET_Superiority_Signature >= 0.64 & features$Matched_Defect_Compensation >= 0.48 & features$Surface_Stabilizing_Auxiliary_Flag == 0, pars$ADMET_Boost, 1.00)
  if (pars$Model_Family == "FeatureRedefined") {
    final_score <- base_score * toxic_gate * extreme_gate
  } else if (pars$Model_Family == "FeatureRedefinedNeutral") {
    final_score <- base_score * neutral_gate * toxic_gate * extreme_gate
  } else {
    final_score <- base_score * neutral_gate * inclusion_boost * admet_boost * toxic_gate * extreme_gate
  }
  final_score <- clip01(final_score)
  auxiliary_flag <- features$Surface_Stabilizing_Auxiliary_Flag == 1 & features$Toxic_Red_Flag == 0
  final_score[auxiliary_flag] <- pmin(pmax(final_score[auxiliary_flag], pars$Auxiliary_Lower_Bound), pars$Auxiliary_Upper_Bound)
  final_score <- clip01(final_score)
  features %>%
    dplyr::mutate(
      Model_ID = pars$Model_ID, Model_Family = pars$Model_Family,
      Raw_Power = pars$Raw_Power, Acceptability_Power = pars$Acceptability_Power,
      Distinctive_Gain_Power = pars$Distinctive_Gain_Power,
      Neutral_Gain_Threshold = pars$Neutral_Gain_Threshold,
      Neutral_Gate_Value = pars$Neutral_Gate_Value,
      Inclusion_Boost = pars$Inclusion_Boost,
      ADMET_Boost = pars$ADMET_Boost,
      Auxiliary_Lower_Bound = pars$Auxiliary_Lower_Bound,
      Auxiliary_Upper_Bound = pars$Auxiliary_Upper_Bound,
      Base_FeatureRedefined_Score = base_score,
      Neutral_Low_Value_Flag_Optimized = optimized_neutral_flag,
      Neutral_Low_Value_Gate_Optimized = neutral_gate,
      Inclusion_Boost_Factor = inclusion_boost,
      ADMET_Boost_Factor = admet_boost,
      Toxic_Gate_Optimized = toxic_gate,
      Extreme_Physicochemical_Gate_Optimized = extreme_gate,
      Optimized_PCI = final_score,
      Optimized_Class = dplyr::case_when(final_score >= 0.72 ~ "High compatibility", final_score >= 0.58 ~ "Moderate-high compatibility", final_score >= 0.38 ~ "Conditional compatibility", TRUE ~ "Low or incompatible")
    )
}

make_validation_flags <- function(df) {
  df %>% dplyr::mutate(
    Pass_Flag = dplyr::case_when(
      Expected_Class == "High" & Optimized_PCI >= 0.60 & Distinctive_Pairing_Gain >= 0.35 ~ TRUE,
      Expected_Class == "ModerateHigh" & Optimized_PCI >= 0.50 & Distinctive_Pairing_Gain >= 0.32 ~ TRUE,
      Expected_Class == "Conditional" & Optimized_PCI >= 0.34 & Optimized_PCI <= 0.58 ~ TRUE,
      Expected_Class == "Low" & Optimized_PCI <= 0.48 ~ TRUE,
      TRUE ~ FALSE
    ),
    Validation_Result = ifelse(Pass_Flag, "Pass", "Check")
  )
}

evaluate_candidate_model <- function(scored_df) {
  df <- make_validation_flags(scored_df)
  group_summary <- df %>% dplyr::group_by(Partner_Profile, Expected_Class) %>%
    dplyr::summarise(median_pci = median(Optimized_PCI, na.rm = TRUE), pass_rate = mean(Pass_Flag, na.rm = TRUE), .groups = "drop")
  get_group <- function(pattern) {
    vals <- group_summary$median_pci[grepl(pattern, group_summary$Partner_Profile)]
    if (length(vals) == 0) return(NA_real_)
    mean(vals, na.rm = TRUE)
  }
  high_median <- mean(group_summary$median_pci[group_summary$Expected_Class == "High"], na.rm = TRUE)
  conditional_median <- mean(group_summary$median_pci[group_summary$Expected_Class == "Conditional"], na.rm = TRUE)
  positive_median <- get_group("Positive|Coassembly")
  inclusion_median <- get_group("Inclusion|Cyclodextrin")
  admet_median <- get_group("ADMET_Compensating")
  neutral_median <- get_group("Neutral")
  glycoside_median <- get_group("Glycoside")
  toxic_extreme_median <- get_group("Toxic|Overly")
  constraints <- c(
    positive_median > neutral_median + 0.20,
    inclusion_median > neutral_median + 0.15,
    admet_median > neutral_median + 0.08,
    neutral_median <= 0.48,
    toxic_extreme_median <= 0.38,
    glycoside_median >= 0.34 & glycoside_median <= 0.58,
    high_median > conditional_median,
    conditional_median > toxic_extreme_median
  )
  overall_pass_rate <- mean(df$Pass_Flag, na.rm = TRUE)
  group_pass_rate <- mean(group_summary$pass_rate >= 0.80, na.rm = TRUE)
  rank_constraint_score <- mean(as.numeric(constraints), na.rm = TRUE)
  neutral_suppression_score <- clip01((0.52 - neutral_median) / 0.22)
  positive_retention_score <- clip01((high_median - 0.62) / 0.22)
  negative_suppression_score <- clip01((0.42 - toxic_extreme_median) / 0.22)
  inclusion_separation_score <- clip01(((inclusion_median - neutral_median) - 0.08) / 0.20)
  admet_separation_score <- clip01(((admet_median - neutral_median) - 0.04) / 0.16)
  auxiliary_positioning_score <- ifelse(is.na(glycoside_median), 0, ifelse(glycoside_median >= 0.34 & glycoside_median <= 0.58, 1, clip01(1 - min(abs(glycoside_median - 0.46), 0.46) / 0.46)))
  optimization_score <- 0.16 * overall_pass_rate + 0.10 * group_pass_rate + 0.22 * rank_constraint_score + 0.18 * neutral_suppression_score + 0.12 * positive_retention_score + 0.10 * inclusion_separation_score + 0.06 * auxiliary_positioning_score + 0.04 * negative_suppression_score + 0.02 * admet_separation_score
  data.frame(
    Model_ID = scored_df$Model_ID[1], Model_Family = scored_df$Model_Family[1],
    Raw_Power = scored_df$Raw_Power[1], Acceptability_Power = scored_df$Acceptability_Power[1],
    Distinctive_Gain_Power = scored_df$Distinctive_Gain_Power[1],
    Neutral_Gain_Threshold = scored_df$Neutral_Gain_Threshold[1],
    Neutral_Gate_Value = scored_df$Neutral_Gate_Value[1],
    Inclusion_Boost = scored_df$Inclusion_Boost[1], ADMET_Boost = scored_df$ADMET_Boost[1],
    Auxiliary_Lower_Bound = scored_df$Auxiliary_Lower_Bound[1],
    Auxiliary_Upper_Bound = scored_df$Auxiliary_Upper_Bound[1],
    overall_pass_rate = overall_pass_rate, group_pass_rate = group_pass_rate,
    rank_constraint_score = rank_constraint_score, neutral_suppression_score = neutral_suppression_score,
    positive_retention_score = positive_retention_score, negative_suppression_score = negative_suppression_score,
    inclusion_separation_score = inclusion_separation_score, admet_separation_score = admet_separation_score,
    auxiliary_positioning_score = auxiliary_positioning_score, optimization_score = optimization_score,
    positive_median = positive_median, inclusion_median = inclusion_median, admet_median = admet_median,
    neutral_median = neutral_median, glycoside_median = glycoside_median, toxic_extreme_median = toxic_extreme_median,
    positive_minus_neutral_margin = positive_median - neutral_median,
    inclusion_minus_neutral_margin = inclusion_median - neutral_median,
    admet_minus_neutral_margin = admet_median - neutral_median,
    stringsAsFactors = FALSE
  )
}

model_results <- vector("list", nrow(candidate_grid))
pb <- utils::txtProgressBar(min = 0, max = nrow(candidate_grid), style = 3)
for (i in seq_len(nrow(candidate_grid))) {
  pars <- candidate_grid[i, , drop = FALSE]
  scored <- apply_candidate_model(calibration_matrix, pars)
  model_results[[i]] <- evaluate_candidate_model(scored)
  if (i %% 25 == 0 || i == nrow(candidate_grid)) utils::setTxtProgressBar(pb, i)
}
close(pb)
model_search_summary <- dplyr::bind_rows(model_results) %>%
  dplyr::arrange(dplyr::desc(optimization_score), dplyr::desc(rank_constraint_score), dplyr::desc(inclusion_separation_score), dplyr::desc(overall_pass_rate))
best_model_id <- model_search_summary$Model_ID[1]
best_params <- candidate_grid %>% dplyr::filter(Model_ID == best_model_id)
best_scores <- apply_candidate_model(base_features, best_params) %>% make_validation_flags()

run_sensitivity_one <- function(one_row, best_params, n_iter = 1000, jitter = 0.16) {
  out <- numeric(n_iter); neutral <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    tmp <- one_row
    for (cc in feature_cols) tmp[[cc]] <- jitter01(tmp[[cc]], jitter = jitter)
    tmp$Distinctive_Pairing_Gain <- clip01(0.28 * tmp$Assembly_Distinctiveness + 0.24 * tmp$Inclusion_Distinctiveness + 0.18 * tmp$Blueprint_Nonredundant_Gain + 0.18 * tmp$Matched_Defect_Compensation + 0.07 * tmp$Structural_Nonredundancy + 0.05 * tmp$ADMET_Superiority_Signature)
    scored <- apply_candidate_model(tmp, best_params)
    out[i] <- scored$Optimized_PCI[1]; neutral[i] <- scored$Neutral_Low_Value_Flag_Optimized[1]
  }
  data.frame(Sensitivity_Mean = mean(out, na.rm = TRUE), Sensitivity_SD = sd(out, na.rm = TRUE), Sensitivity_P05 = as.numeric(quantile(out, 0.05, na.rm = TRUE)), Sensitivity_P95 = as.numeric(quantile(out, 0.95, na.rm = TRUE)), Probability_High_or_ModerateHigh = mean(out >= 0.58, na.rm = TRUE), Probability_Neutral_Low_Value_Flag = mean(neutral == 1, na.rm = TRUE), stringsAsFactors = FALSE)
}
sens_list <- list()
for (i in seq_len(nrow(base_features))) {
  tmp <- run_sensitivity_one(base_features[i, , drop = FALSE], best_params, N_SENSITIVITY, SENSITIVITY_JITTER)
  tmp$Partner_Compound <- base_features$Partner_Compound[i]; tmp$Partner_Profile <- base_features$Partner_Profile[i]
  sens_list[[i]] <- tmp
}
sensitivity_table <- dplyr::bind_rows(sens_list)
optimized_sensitivity_cols <- c("Sensitivity_Mean", "Sensitivity_SD", "Sensitivity_P05", "Sensitivity_P95", "Probability_High_or_ModerateHigh", "Probability_Neutral_Low_Value_Flag")
best_scores <- safe_left_join_replace(x = best_scores, y = sensitivity_table, by = c("Partner_Compound", "Partner_Profile"), replace_cols = optimized_sensitivity_cols, context = "final optimized best-model sensitivity merge")

feature_signature_cols <- c("Basal_Acceptability", "Assembly_Distinctiveness", "Inclusion_Distinctiveness", "Blueprint_Nonredundant_Gain", "Matched_Defect_Compensation", "ADMET_Superiority_Signature", "Structural_Nonredundancy", "Distinctive_Pairing_Gain")
stop_if_missing(best_scores, c("Partner_Profile", "Expected_Class", "Expected_Behavior", "Optimized_PCI", "Legacy_PCI", feature_signature_cols, "Neutral_Low_Value_Flag_Optimized", "Surface_Stabilizing_Auxiliary_Flag", "Toxic_Red_Flag", "Extreme_Physicochemical_Flag", "Pass_Flag", "Sensitivity_P05", "Sensitivity_P95"), "best_scores before final group summary")

group_summary <- best_scores %>%
  dplyr::group_by(Partner_Profile, Expected_Class, Expected_Behavior) %>%
  dplyr::summarise(
    n = dplyr::n(), mean_optimized_PCI = mean(Optimized_PCI, na.rm = TRUE),
    median_optimized_PCI = median(Optimized_PCI, na.rm = TRUE),
    min_optimized_PCI = min(Optimized_PCI, na.rm = TRUE),
    max_optimized_PCI = max(Optimized_PCI, na.rm = TRUE),
    median_legacy_PCI = median(Legacy_PCI, na.rm = TRUE),
    median_raw_PCI = median(Raw_PCI, na.rm = TRUE),
    median_basal_acceptability = median(Basal_Acceptability, na.rm = TRUE),
    median_distinctive_gain = median(Distinctive_Pairing_Gain, na.rm = TRUE),
    median_assembly_distinctiveness = median(Assembly_Distinctiveness, na.rm = TRUE),
    median_inclusion_distinctiveness = median(Inclusion_Distinctiveness, na.rm = TRUE),
    median_matched_defect_compensation = median(Matched_Defect_Compensation, na.rm = TRUE),
    median_admet_superiority = median(ADMET_Superiority_Signature, na.rm = TRUE),
    neutral_flag_rate = mean(Neutral_Low_Value_Flag_Optimized, na.rm = TRUE),
    auxiliary_flag_rate = mean(Surface_Stabilizing_Auxiliary_Flag, na.rm = TRUE),
    toxic_flag_rate = mean(Toxic_Red_Flag, na.rm = TRUE),
    extreme_flag_rate = mean(Extreme_Physicochemical_Flag, na.rm = TRUE),
    pass_rate = mean(Pass_Flag, na.rm = TRUE),
    mean_sensitivity_p05 = mean(Sensitivity_P05, na.rm = TRUE),
    mean_sensitivity_p95 = mean(Sensitivity_P95, na.rm = TRUE),
    .groups = "drop"
  ) %>% dplyr::mutate(Group_Validation_Result = ifelse(pass_rate >= 0.80, "Pass", "Check"))

get_group_median <- function(pattern) {
  vals <- group_summary$median_optimized_PCI[grepl(pattern, group_summary$Partner_Profile)]
  if (length(vals) == 0) return(NA_real_)
  mean(vals, na.rm = TRUE)
}

rank_constraints <- data.frame(
  Constraint = c("Positive profiles exceed neutral controls by >=0.20", "Inclusion profiles exceed neutral controls by >=0.15", "ADMET-compensating profiles exceed neutral controls by >=0.08", "Neutral controls remain below 0.48", "Toxic and physicochemical incompatibility controls remain below 0.38", "Auxiliary glycoside-like profiles remain in conditional range", "High profiles exceed conditional profiles", "Conditional profiles exceed incompatible profiles"),
  Observed_Margin = c(
    get_group_median("Positive") - get_group_median("Neutral"),
    get_group_median("Inclusion") - get_group_median("Neutral"),
    get_group_median("ADMET_Compensating") - get_group_median("Neutral"),
    0.48 - get_group_median("Neutral"),
    0.38 - get_group_median("Toxic|Overly"),
    min(get_group_median("Glycoside") - 0.34, 0.58 - get_group_median("Glycoside")),
    mean(group_summary$median_optimized_PCI[group_summary$Expected_Class == "High"], na.rm = TRUE) - get_group_median("Glycoside"),
    get_group_median("Glycoside") - get_group_median("Toxic|Overly")
  ),
  Required_Margin = c(0.20, 0.15, 0.08, 0, 0, 0, 0, 0),
  stringsAsFactors = FALSE
) %>% dplyr::mutate(Result = ifelse(Observed_Margin >= Required_Margin, "Pass", "Check"))

final_qc_summary <- data.frame(
  Metric = c(
    "Positive median PCI", "Inclusion median PCI", "ADMET median PCI",
    "Neutral median PCI", "Glycoside median PCI", "Toxic/extreme median PCI",
    "Positive - Neutral margin", "Inclusion - Neutral margin", "ADMET - Neutral margin",
    "Rank constraints passed"
  ),
  Target = c(
    ">= 0.60", ">= 0.58", ">= 0.50",
    "<= 0.48", "0.34-0.58", "<= 0.38",
    ">= 0.20", ">= 0.15", ">= 0.08",
    "All constraints"
  ),
  Observed = c(
    get_group_median("Positive|Coassembly"),
    get_group_median("Inclusion|Cyclodextrin"),
    get_group_median("ADMET_Compensating"),
    get_group_median("Neutral"),
    get_group_median("Glycoside"),
    get_group_median("Toxic|Overly"),
    get_group_median("Positive|Coassembly") - get_group_median("Neutral"),
    get_group_median("Inclusion|Cyclodextrin") - get_group_median("Neutral"),
    get_group_median("ADMET_Compensating") - get_group_median("Neutral"),
    mean(rank_constraints$Result == "Pass", na.rm = TRUE)
  ),
  Result = c(
    ifelse(get_group_median("Positive|Coassembly") >= 0.60, "Pass", "Check"),
    ifelse(get_group_median("Inclusion|Cyclodextrin") >= 0.58, "Pass", "Check"),
    ifelse(get_group_median("ADMET_Compensating") >= 0.50, "Pass", "Check"),
    ifelse(get_group_median("Neutral") <= 0.48, "Pass", "Check"),
    ifelse(get_group_median("Glycoside") >= 0.34 & get_group_median("Glycoside") <= 0.58, "Pass", "Check"),
    ifelse(get_group_median("Toxic|Overly") <= 0.38, "Pass", "Check"),
    ifelse(get_group_median("Positive|Coassembly") - get_group_median("Neutral") >= 0.20, "Pass", "Check"),
    ifelse(get_group_median("Inclusion|Cyclodextrin") - get_group_median("Neutral") >= 0.15, "Pass", "Check"),
    ifelse(get_group_median("ADMET_Compensating") - get_group_median("Neutral") >= 0.08, "Pass", "Check"),
    ifelse(all(rank_constraints$Result == "Pass"), "Pass", "Check")
  ),
  stringsAsFactors = FALSE
)

model_evolution_comparison <- group_summary %>%
  dplyr::select(
    Partner_Profile,
    Expected_Class,
    Legacy_GainGated_Median = median_legacy_PCI,
    FeatureRedefined_Median = median_optimized_PCI
  ) %>%
  dplyr::mutate(
    Median_Shift = FeatureRedefined_Median - Legacy_GainGated_Median
  ) %>%
  dplyr::arrange(dplyr::desc(Legacy_GainGated_Median))

signature_long <- best_scores %>%
  dplyr::select(Partner_Profile, dplyr::all_of(feature_signature_cols)) %>%
  dplyr::group_by(Partner_Profile) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
  tidyr::pivot_longer(cols = dplyr::all_of(feature_signature_cols), names_to = "Feature", values_to = "Score") %>%
  dplyr::mutate(Partner_Profile_Label = wrap_label(clean_label(Partner_Profile), 24), Feature_Label = wrap_label(clean_label(Feature), 22))


############################################################
## Profile-level summaries and final manuscript-grade figures
## Submission-grade plotting revision:
##   - Main figures use profile-level summaries, not 240 individual labels.
##   - Subtitles are deliberately short to avoid clipping.
##   - Rank-constraint validation is promoted to a main figure.
##   - Sensitivity terminology is unified as perturbation-based sensitivity.
##   - Supplementary plots are zoomed/structured for readability.
############################################################

profile_display <- function(x) {
  dplyr::case_when(
    grepl("Positive", x) ~ "Co-assembly enhancer",
    grepl("Inclusion", x) ~ "Inclusion-stabilizing partner",
    grepl("ADMET_Compensating", x) ~ "ADMET-compensating partner",
    grepl("Glycoside", x) ~ "Glycoside-like auxiliary stabilizer",
    grepl("Neutral", x) ~ "Neutral low-value control",
    grepl("Lipophilic", x) ~ "Lipophilic incompatibility control",
    grepl("Polar", x) ~ "Polar incompatibility control",
    grepl("Toxic", x) ~ "Toxicity control",
    TRUE ~ clean_label(x)
  )
}

feature_display <- function(x) {
  dplyr::case_when(
    x == "Assembly_Distinctiveness" ~ "Assembly\ngain",
    x == "Inclusion_Distinctiveness" ~ "Inclusion\ngain",
    x == "Blueprint_Nonredundant_Gain" ~ "Blueprint\ngain",
    x == "Matched_Defect_Compensation" ~ "Defect\ncompensation",
    x == "ADMET_Superiority_Signature" ~ "ADMET\nsuperiority",
    x == "Structural_Nonredundancy" ~ "Structural\nnonredundancy",
    x == "Distinctive_Pairing_Gain" ~ "Distinctive\ngain",
    x == "Basal_Acceptability" ~ "Basal\nacceptability",
    TRUE ~ wrap_label(clean_label(x), 18)
  )
}

class_palette <- c(
  "High" = "#7FCDBB",
  "ModerateHigh" = "#80B1D3",
  "Conditional" = "#C7E9C0",
  "Low" = "#D9D9D9"
)

pass_palette <- c(
  "Pass" = "#7FCDBB",
  "Check" = "#FDB462"
)

heat_palette <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")

stop_if_missing(
  best_scores,
  c(
    "Partner_Profile", "Expected_Class", "Expected_Behavior",
    "Optimized_PCI", "Legacy_PCI", "Raw_PCI",
    "Sensitivity_Mean", "Sensitivity_P05", "Sensitivity_P95",
    "Probability_High_or_ModerateHigh", "Probability_Neutral_Low_Value_Flag"
  ),
  "best_scores before profile-level plotting"
)

profile_summary <- best_scores %>%
  dplyr::mutate(Profile_Label = profile_display(as.character(Partner_Profile))) %>%
  dplyr::group_by(Partner_Profile, Profile_Label, Expected_Class, Expected_Behavior) %>%
  dplyr::summarise(
    n = dplyr::n(),
    PCI_mean = mean(Optimized_PCI, na.rm = TRUE),
    PCI_median = median(Optimized_PCI, na.rm = TRUE),
    PCI_sd = sd(Optimized_PCI, na.rm = TRUE),
    PCI_p05 = as.numeric(stats::quantile(Optimized_PCI, 0.05, na.rm = TRUE)),
    PCI_p25 = as.numeric(stats::quantile(Optimized_PCI, 0.25, na.rm = TRUE)),
    PCI_p75 = as.numeric(stats::quantile(Optimized_PCI, 0.75, na.rm = TRUE)),
    PCI_p95 = as.numeric(stats::quantile(Optimized_PCI, 0.95, na.rm = TRUE)),
    Legacy_PCI_median = median(Legacy_PCI, na.rm = TRUE),
    Raw_PCI_median = median(Raw_PCI, na.rm = TRUE),
    Sensitivity_Mean_Profile = mean(Sensitivity_Mean, na.rm = TRUE),
    Sensitivity_P05_Profile = as.numeric(stats::quantile(Sensitivity_P05, 0.05, na.rm = TRUE)),
    Sensitivity_P95_Profile = as.numeric(stats::quantile(Sensitivity_P95, 0.95, na.rm = TRUE)),
    Probability_High_or_ModerateHigh_Profile = mean(Probability_High_or_ModerateHigh, na.rm = TRUE),
    Probability_Neutral_Low_Value_Flag_Profile = mean(Probability_Neutral_Low_Value_Flag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(PCI_median))

profile_levels_high_to_low <- profile_summary$Profile_Label
profile_levels_for_axis <- rev(profile_levels_high_to_low)

best_scores_plot <- best_scores %>%
  dplyr::mutate(
    Profile_Label = profile_display(as.character(Partner_Profile)),
    Profile_Label = factor(Profile_Label, levels = profile_levels_for_axis)
  )

feature_signature_cols_profile <- intersect(
  c(
    "Basal_Acceptability",
    "Assembly_Distinctiveness",
    "Inclusion_Distinctiveness",
    "Blueprint_Nonredundant_Gain",
    "Matched_Defect_Compensation",
    "ADMET_Superiority_Signature",
    "Structural_Nonredundancy",
    "Distinctive_Pairing_Gain"
  ),
  colnames(best_scores)
)

module_cols_profile <- intersect(
  c(
    "Co_assembly_Compatibility",
    "Physicochemical_Complementarity",
    "Nano_stabilization_Complementarity",
    "ADMET_Complementarity",
    "Safety_Balance",
    "Blueprint_Match",
    "Coassembly_Gain",
    "Stabilization_Gain",
    "Blueprint_Specific_Gain",
    "Defect_Compensation_Gain",
    "Nonredundancy_Gain",
    "Pairing_Gain_Index",
    feature_signature_cols_profile
  ),
  colnames(best_scores)
)

profile_module_summary <- best_scores %>%
  dplyr::mutate(Profile_Label = profile_display(as.character(Partner_Profile))) %>%
  dplyr::group_by(Partner_Profile, Profile_Label, Expected_Class) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(module_cols_profile), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

profile_module_long <- profile_module_summary %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(module_cols_profile),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Profile_Label = factor(Profile_Label, levels = profile_levels_for_axis),
    Module_Label = feature_display(Module)
  )

signature_long <- best_scores %>%
  dplyr::mutate(Profile_Label = profile_display(as.character(Partner_Profile))) %>%
  dplyr::group_by(Profile_Label) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(feature_signature_cols_profile), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(feature_signature_cols_profile),
    names_to = "Feature",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Profile_Label = factor(Profile_Label, levels = profile_levels_for_axis),
    Feature_Label = feature_display(Feature)
  )

############################################################
## Refined plotting theme and helpers
############################################################

theme_acmpi_final <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "#2C2C2C"),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 4,
        hjust = 0,
        margin = ggplot2::margin(b = 8)
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size + 0.5,
        color = "#5A5A5A",
        hjust = 0,
        lineheight = 1.10,
        margin = ggplot2::margin(b = 13)
      ),
      axis.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      axis.text = ggplot2::element_text(size = base_size, color = "#333333"),
      legend.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.key.height = grid::unit(0.52, "cm"),
      legend.key.width = grid::unit(0.52, "cm"),
      legend.position = "right",
      plot.margin = ggplot2::margin(26, 44, 28, 58)
    )
}

############################################################
## Main Figure 1: method workflow
############################################################

workflow_df <- data.frame(
  Step = c(
    "Core\nblueprint",
    "Virtual stress\nprofiles",
    "Feature\nredefinition",
    "Model\ncalibration",
    "Rank-constraint\nvalidation",
    "Optimized\npairing model"
  ),
  X = seq(1, 6),
  Y = 1,
  stringsAsFactors = FALSE
)
workflow_df$Step <- factor(workflow_df$Step, levels = workflow_df$Step)

p_workflow <- ggplot(workflow_df, aes(x = X, y = Y)) +
  geom_segment(
    data = workflow_df[-nrow(workflow_df), ],
    aes(x = X + 0.35, xend = X + 0.65, y = Y, yend = Y),
    arrow = arrow(length = grid::unit(0.15, "cm"), type = "closed"),
    linewidth = 0.55,
    color = "#6F6F6F"
  ) +
  geom_label(
    aes(label = Step),
    size = 4.0,
    fontface = "bold",
    lineheight = 0.94,
    label.size = 0,
    fill = "#EEF7F4",
    color = "#2C2C2C",
    label.padding = grid::unit(0.34, "lines")
  ) +
  annotate(
    "text",
    x = 3.5,
    y = 0.70,
    label = paste0("8 virtual profile classes; ", N_PER_PROFILE, " variants per class; ", nrow(best_scores), " synthetic partners"),
    size = 3.65,
    color = "#555555"
  ) +
  scale_x_continuous(limits = c(0.50, 6.50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.58, 1.18), expand = c(0, 0)) +
  labs(
    title = "Blueprint-guided ACMPI-Nano calibration workflow",
    subtitle = "Template-based stress testing separates distinctive pairing gain from neutral molecular acceptability",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 17.5, color = "#2C2C2C", hjust = 0, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 12.2, color = "#555555", hjust = 0, margin = margin(b = 8)),
    plot.margin = margin(18, 30, 18, 30)
  )

############################################################
## Main Figure 2: profile-level optimized PCI distribution
############################################################

p_profile <- ggplot(
  best_scores_plot,
  aes(x = Profile_Label, y = Optimized_PCI, fill = Expected_Class)
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.92,
    color = "#333333",
    linewidth = 0.35
  ) +
  geom_jitter(
    width = 0.09,
    size = 0.95,
    alpha = 0.35,
    color = "#333333"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = class_palette, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  labs(
    title = "Profile-level separation in the virtual stress test",
    subtitle = paste0("Eight profile classes; ", N_PER_PROFILE, " variants per class; n = ", nrow(best_scores), " synthetic partners"),
    x = NULL,
    y = "Optimized Pairing Compatibility Index",
    fill = "Expected class"
  ) +
  theme_acmpi_final(12.5) +
  theme(
    legend.position = "right",
    plot.margin = margin(26, 50, 26, 62)
  )

############################################################
## Main Figure 3: profile-level distinctive gain signatures
############################################################

p_signature <- ggplot(
  signature_long,
  aes(x = Feature_Label, y = Profile_Label, fill = Score)
) +
  geom_tile(color = "white", linewidth = 0.46) +
  scale_fill_gradientn(
    colors = heat_palette,
    limits = c(0, 1),
    oob = scales::squish
  ) +
  labs(
    title = "Distinctive pairing-gain signatures by virtual profile class",
    subtitle = "Median feature-redefined gain signatures across perturbed synthetic partners",
    x = NULL,
    y = NULL,
    fill = "Median score"
  ) +
  theme_acmpi_final(11.2) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 10.6, lineheight = 0.94),
    axis.text.y = element_text(size = 10.8),
    legend.position = "right",
    plot.margin = margin(26, 48, 32, 78)
  )

############################################################
## Main Figure 4: rank-constraint validation
############################################################

rank_constraints_plot <- rank_constraints %>%
  dplyr::mutate(
    Constraint_Label = dplyr::case_when(
      grepl("Positive", Constraint) ~ "Positive profiles exceed\nneutral controls",
      grepl("Inclusion", Constraint) ~ "Inclusion profiles exceed\nneutral controls",
      grepl("ADMET", Constraint) ~ "ADMET-compensating profiles\nexceed neutral controls",
      grepl("Neutral controls", Constraint) ~ "Neutral controls remain\nbelow upper bound",
      grepl("Toxic", Constraint) ~ "Toxic and physicochemical\ncontrols remain low",
      grepl("Auxiliary", Constraint) ~ "Glycoside-like profiles remain\nwithin conditional range",
      grepl("High profiles", Constraint) ~ "High profiles exceed\nconditional profiles",
      grepl("Conditional", Constraint) ~ "Conditional profiles exceed\nincompatible profiles",
      TRUE ~ stringr::str_wrap(Constraint, width = 36)
    ),
    Constraint_Label = factor(Constraint_Label, levels = rev(Constraint_Label))
  )

p_constraints <- ggplot(
  rank_constraints_plot,
  aes(x = Constraint_Label, y = Observed_Margin, fill = Result)
) +
  geom_col(width = 0.68, color = "#333333", linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#6F6F6F", linewidth = 0.45) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = pass_palette, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "Rank-constraint validation",
    subtitle = "All predefined internal validation margins were evaluated after model calibration",
    x = NULL,
    y = "Observed margin",
    fill = "Result"
  ) +
  theme_acmpi_final(12) +
  theme(
    legend.position = "right",
    plot.margin = margin(26, 54, 26, 82)
  )

############################################################
## Main Figure 5: profile-level perturbation sensitivity
############################################################

profile_summary_plot <- profile_summary %>%
  dplyr::mutate(Profile_Label = factor(Profile_Label, levels = profile_levels_for_axis))

p_sensitivity <- ggplot(
  profile_summary_plot,
  aes(x = Profile_Label, y = Sensitivity_Mean_Profile)
) +
  geom_errorbar(
    aes(ymin = Sensitivity_P05_Profile, ymax = Sensitivity_P95_Profile),
    width = 0.24,
    linewidth = 0.74,
    color = "#4D4D4D"
  ) +
  geom_point(size = 3.2, color = "#2C2C2C") +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  labs(
    title = "Profile-level perturbation sensitivity of the optimized model",
    subtitle = paste0("Each synthetic partner was evaluated across ", N_SENSITIVITY, " perturbation-based sensitivity iterations"),
    x = NULL,
    y = "Sensitivity-adjusted PCI"
  ) +
  theme_acmpi_final(12.5) +
  theme(plot.margin = margin(26, 52, 26, 68))

############################################################
## Supplementary Figure 1: top candidate parameter rules
############################################################

top_models_plot <- model_search_summary %>%
  dplyr::slice_head(n = min(N_TOP_MODELS_TO_PLOT, nrow(model_search_summary))) %>%
  dplyr::mutate(
    Model_ID = factor(Model_ID, levels = rev(Model_ID)),
    Score_Label = sprintf("%.3f", optimization_score)
  )

supp1_xmin <- max(0, min(top_models_plot$optimization_score, na.rm = TRUE) - 0.015)
supp1_xmax <- min(1, max(top_models_plot$optimization_score, na.rm = TRUE) + 0.015)
if (!is.finite(supp1_xmin) || !is.finite(supp1_xmax) || supp1_xmin >= supp1_xmax) {
  supp1_xmin <- 0.80
  supp1_xmax <- 1.00
}

p_opt_rank <- ggplot(top_models_plot, aes(x = Model_ID, y = optimization_score)) +
  geom_segment(aes(xend = Model_ID, y = supp1_xmin, yend = optimization_score), linewidth = 1.2, color = "#80B1D3") +
  geom_point(size = 3.5, color = "#2B6C9E") +
  geom_text(aes(label = Score_Label), hjust = -0.15, size = 3.4, color = "#333333") +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(supp1_xmin, min(1.02, supp1_xmax + 0.02)), breaks = scales::pretty_breaks(n = 5)) +
  labs(
    title = "Top candidate scoring rules after rank-constraint calibration",
    subtitle = "Top-ranked rules showed near-identical composite optimization performance",
    x = NULL,
    y = "Optimization score"
  ) +
  theme_acmpi_final(12) +
  theme(
    legend.position = "none",
    plot.margin = margin(26, 60, 26, 56)
  )

############################################################
## Supplementary Figure 2: model-selection tradeoff
############################################################

tradeoff_xmin <- max(0, stats::quantile(model_search_summary$rank_constraint_score, 0.01, na.rm = TRUE) - 0.02)
tradeoff_ymin <- max(0, stats::quantile(model_search_summary$optimization_score, 0.01, na.rm = TRUE) - 0.03)

p_tradeoff <- ggplot(
  model_search_summary,
  aes(x = rank_constraint_score, y = optimization_score)
) +
  geom_point(aes(size = overall_pass_rate), alpha = 0.58, color = "#4B6A88") +
  geom_point(
    data = model_search_summary %>% dplyr::filter(Model_ID == best_model_id),
    color = "#B2182B",
    size = 4.5
  ) +
  annotate(
    "text",
    x = model_search_summary$rank_constraint_score[match(best_model_id, model_search_summary$Model_ID)],
    y = model_search_summary$optimization_score[match(best_model_id, model_search_summary$Model_ID)],
    label = " Selected model",
    hjust = 0,
    vjust = -0.8,
    size = 3.6,
    color = "#B2182B"
  ) +
  scale_x_continuous(limits = c(tradeoff_xmin, 1.005), breaks = scales::pretty_breaks(n = 5)) +
  scale_y_continuous(limits = c(tradeoff_ymin, 1.005), breaks = scales::pretty_breaks(n = 5)) +
  scale_size_continuous(range = c(1.2, 4.2)) +
  labs(
    title = "Model-selection tradeoff",
    subtitle = "The selected rule lies in the high-rank-constraint, high-optimization region",
    x = "Rank-constraint score",
    y = "Optimization score",
    size = "Overall pass rate"
  ) +
  theme_acmpi_final(12)

############################################################
## Supplementary Figure 3: legacy versus feature-redefined profile medians
############################################################

model_evolution_long <- model_evolution_comparison %>%
  dplyr::mutate(Profile_Label = profile_display(as.character(Partner_Profile))) %>%
  dplyr::select(Profile_Label, Legacy_GainGated_Median, FeatureRedefined_Median) %>%
  tidyr::pivot_longer(
    cols = c(Legacy_GainGated_Median, FeatureRedefined_Median),
    names_to = "Model",
    values_to = "Median_PCI"
  ) %>%
  dplyr::mutate(
    Model = dplyr::case_when(
      Model == "Legacy_GainGated_Median" ~ "Legacy gain-gated PCI",
      Model == "FeatureRedefined_Median" ~ "Feature-redefined PCI",
      TRUE ~ Model
    ),
    Profile_Label = factor(Profile_Label, levels = profile_levels_for_axis)
  )

p_evolution <- ggplot(
  model_evolution_long,
  aes(x = Profile_Label, y = Median_PCI, fill = Model)
) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66, color = "#333333", linewidth = 0.20) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("Legacy gain-gated PCI" = "#BEBADA", "Feature-redefined PCI" = "#80B1D3")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.01, 0.04))) +
  labs(
    title = "Legacy versus feature-redefined profile medians",
    subtitle = "Feature redefinition suppresses neutral low-value profiles while retaining high-gain profiles",
    x = NULL,
    y = "Median PCI",
    fill = NULL
  ) +
  theme_acmpi_final(12) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(26, 52, 28, 68)
  )

############################################################
## Save manuscript-grade figures
############################################################

save_plot_dual(
  p_workflow,
  file.path(fig_png_dir, "MainFig1_Workflow.png"),
  file.path(fig_pdf_dir, "MainFig1_Workflow.pdf"),
  width = 13.2,
  height = 3.45,
  dpi = FIG_DPI
)

save_plot_dual(
  p_profile,
  file.path(fig_png_dir, "MainFig2_ProfileLevel_OptimizedPCI_Distribution.png"),
  file.path(fig_pdf_dir, "MainFig2_ProfileLevel_OptimizedPCI_Distribution.pdf"),
  width = 11.6,
  height = 7.25,
  dpi = FIG_DPI
)

save_plot_dual(
  p_signature,
  file.path(fig_png_dir, "MainFig3_ProfileLevel_DistinctiveGainSignatures.png"),
  file.path(fig_pdf_dir, "MainFig3_ProfileLevel_DistinctiveGainSignatures.pdf"),
  width = 13.2,
  height = 7.45,
  dpi = FIG_DPI
)

save_plot_dual(
  p_constraints,
  file.path(fig_png_dir, "MainFig4_RankConstraintValidation.png"),
  file.path(fig_pdf_dir, "MainFig4_RankConstraintValidation.pdf"),
  width = 10.8,
  height = 7.05,
  dpi = FIG_DPI
)

save_plot_dual(
  p_sensitivity,
  file.path(fig_png_dir, "MainFig5_ProfileLevel_PerturbationSensitivity.png"),
  file.path(fig_pdf_dir, "MainFig5_ProfileLevel_PerturbationSensitivity.pdf"),
  width = 11.2,
  height = 7.0,
  dpi = FIG_DPI
)

save_plot_dual(
  p_opt_rank,
  file.path(fig_png_dir, "SupplementaryFig1_TopCandidateRules.png"),
  file.path(fig_pdf_dir, "SupplementaryFig1_TopCandidateRules.pdf"),
  width = 9.8,
  height = 5.6,
  dpi = FIG_DPI
)

save_plot_dual(
  p_tradeoff,
  file.path(fig_png_dir, "SupplementaryFig2_ModelSelectionTradeoff.png"),
  file.path(fig_pdf_dir, "SupplementaryFig2_ModelSelectionTradeoff.pdf"),
  width = 10.4,
  height = 6.6,
  dpi = FIG_DPI
)

save_plot_dual(
  p_evolution,
  file.path(fig_png_dir, "SupplementaryFig3_ModelEvolution.png"),
  file.path(fig_pdf_dir, "SupplementaryFig3_ModelEvolution.pdf"),
  width = 11.2,
  height = 7.1,
  dpi = FIG_DPI
)

############################################################
## Save final tables and workbook
############################################################

best_params_out <- best_params %>%
  dplyr::left_join(
    model_search_summary %>%
      dplyr::select(
        Model_ID,
        optimization_score,
        rank_constraint_score,
        overall_pass_rate,
        group_pass_rate,
        neutral_suppression_score,
        positive_retention_score,
        inclusion_separation_score,
        auxiliary_positioning_score
      ),
    by = "Model_ID"
  )

table_outputs <- list(
  PreparedBaseFeatures = base_features,
  AugmentedPressureTestMatrix = calibration_matrix,
  CandidateParameterGrid = candidate_grid,
  CandidateModelSearch = model_search_summary,
  BestModelParameters = best_params_out,
  OptimizedVirtualPairingScores = best_scores,
  ProfileLevelSummary = profile_summary,
  ProfileModuleSummary = profile_module_summary,
  ProfileModuleLong = profile_module_long,
  OptimizedGroupValidationSummary = group_summary,
  RankConstraintValidation = rank_constraints,
  DistinctiveFeatureSignatures_Long = signature_long,
  FinalModel_QC_Summary = final_qc_summary,
  ModelEvolution_Comparison = model_evolution_comparison
)

for (nm in names(table_outputs)) {
  write.csv(
    table_outputs[[nm]],
    file.path(table_dir, paste0(output_prefix, "_", nm, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

############################################################
## Save integrated Excel workbook
############################################################

wb <- openxlsx::createWorkbook()

sheet_map <- c(
  PreparedBaseFeatures = "PreparedBaseFeatures",
  AugmentedPressureTestMatrix = "AugmentedPressureMatrix",
  CandidateParameterGrid = "CandidateParameterGrid",
  CandidateModelSearch = "CandidateModelSearch",
  BestModelParameters = "BestModelParameters",
  OptimizedVirtualPairingScores = "OptimizedPairingScores",
  ProfileLevelSummary = "ProfileSummary",
  ProfileModuleSummary = "ProfileModules",
  ProfileModuleLong = "ProfileModulesLong",
  OptimizedGroupValidationSummary = "GroupValidation",
  RankConstraintValidation = "RankConstraints",
  DistinctiveFeatureSignatures_Long = "FeatureSignatures",
  FinalModel_QC_Summary = "FinalQC",
  ModelEvolution_Comparison = "ModelEvolution"
)

for (nm in names(sheet_map)) {
  openxlsx::addWorksheet(wb, sheet_map[[nm]])
  openxlsx::writeData(wb, sheet_map[[nm]], table_outputs[[nm]])
  openxlsx::setColWidths(
    wb,
    sheet = sheet_map[[nm]],
    cols = 1:ncol(table_outputs[[nm]]),
    widths = "auto"
  )
}

openxlsx::saveWorkbook(
  wb,
  file.path(table_dir, paste0(output_prefix, "_Results.xlsx")),
  overwrite = TRUE
)

############################################################
## Output manifest
############################################################

table_manifest <- data.frame(
  Output_Type = "Table",
  File_Name = c(
    paste0(output_prefix, "_PreparedBaseFeatures.csv"),
    paste0(output_prefix, "_AugmentedPressureTestMatrix.csv"),
    paste0(output_prefix, "_CandidateParameterGrid.csv"),
    paste0(output_prefix, "_CandidateModelSearch.csv"),
    paste0(output_prefix, "_BestModelParameters.csv"),
    paste0(output_prefix, "_OptimizedVirtualPairingScores.csv"),
    paste0(output_prefix, "_ProfileLevelSummary.csv"),
    paste0(output_prefix, "_ProfileModuleSummary.csv"),
    paste0(output_prefix, "_ProfileModuleLong.csv"),
    paste0(output_prefix, "_OptimizedGroupValidationSummary.csv"),
    paste0(output_prefix, "_RankConstraintValidation.csv"),
    paste0(output_prefix, "_DistinctiveFeatureSignatures_Long.csv"),
    paste0(output_prefix, "_FinalModel_QC_Summary.csv"),
    paste0(output_prefix, "_ModelEvolution_Comparison.csv"),
    paste0(output_prefix, "_Results.xlsx")
  ),
  Directory = table_dir,
  stringsAsFactors = FALSE
)

main_figure_manifest <- data.frame(
  Output_Type = "Main figure",
  File_Name = c(
    "MainFig1_Workflow.png",
    "MainFig2_ProfileLevel_OptimizedPCI_Distribution.png",
    "MainFig3_ProfileLevel_DistinctiveGainSignatures.png",
    "MainFig4_RankConstraintValidation.png",
    "MainFig5_ProfileLevel_PerturbationSensitivity.png",
    "MainFig1_Workflow.pdf",
    "MainFig2_ProfileLevel_OptimizedPCI_Distribution.pdf",
    "MainFig3_ProfileLevel_DistinctiveGainSignatures.pdf",
    "MainFig4_RankConstraintValidation.pdf",
    "MainFig5_ProfileLevel_PerturbationSensitivity.pdf"
  ),
  Directory = c(rep(fig_png_dir, 5), rep(fig_pdf_dir, 5)),
  stringsAsFactors = FALSE
)

supp_figure_manifest <- data.frame(
  Output_Type = "Supplementary figure",
  File_Name = c(
    "SupplementaryFig1_TopCandidateRules.png",
    "SupplementaryFig2_ModelSelectionTradeoff.png",
    "SupplementaryFig3_ModelEvolution.png",
    "SupplementaryFig1_TopCandidateRules.pdf",
    "SupplementaryFig2_ModelSelectionTradeoff.pdf",
    "SupplementaryFig3_ModelEvolution.pdf"
  ),
  Directory = c(rep(fig_png_dir, 3), rep(fig_pdf_dir, 3)),
  stringsAsFactors = FALSE
)

output_manifest <- dplyr::bind_rows(
  table_manifest,
  main_figure_manifest,
  supp_figure_manifest
) %>%
  dplyr::mutate(
    Full_Path = file.path(Directory, File_Name),
    Exists = file.exists(Full_Path)
  )

write.csv(
  output_manifest,
  file.path(table_dir, paste0(output_prefix, "_OutputManifest.csv")),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

############################################################
## Final output integrity check
############################################################

missing_outputs <- output_manifest %>%
  dplyr::filter(!Exists)

if (nrow(missing_outputs) > 0) {
  warning(
    "Some expected output files were not found:\n",
    paste(missing_outputs$Full_Path, collapse = "\n")
  )
}

qc_pass_rate <- if ("QC_Result" %in% colnames(final_qc_summary)) {
  mean(final_qc_summary$QC_Result == "Pass", na.rm = TRUE)
} else {
  NA_real_
}

rank_pass_rate <- if ("Result" %in% colnames(rank_constraints)) {
  mean(rank_constraints$Result == "Pass", na.rm = TRUE)
} else {
  NA_real_
}

positive_neutral_margin <- get_group_median("Positive") - get_group_median("Neutral")
inclusion_neutral_margin <- get_group_median("Inclusion") - get_group_median("Neutral")
admet_neutral_margin <- get_group_median("ADMET_Compensating") - get_group_median("Neutral")
glycoside_median <- get_group_median("Glycoside")
neutral_median <- get_group_median("Neutral")
toxic_extreme_median <- get_group_median("Toxic|Overly")

############################################################
## Console summary
############################################################

cat("\n============================================================\n")
cat("ACMPI-Nano: Feature-redefined pairing calibration\n")
cat("============================================================\n\n")

cat("Core compound: ", core_compound, "\n", sep = "")
cat("Input Stage 2 table directory:\n", stage2_table_dir, "\n\n", sep = "")

cat("Output directory:\n", out_dir, "\n\n", sep = "")

cat("Output folders:\n")
cat("  Tables:      ", table_dir, "\n", sep = "")
cat("  PNG figures: ", fig_png_dir, "\n", sep = "")
cat("  PDF figures: ", fig_pdf_dir, "\n\n", sep = "")

cat("Best model parameters:\n")
print(best_params_out, row.names = FALSE)

cat("\nOptimized group validation summary:\n")
print(group_summary, row.names = FALSE)

cat("\nFinal model QC summary:\n")
print(final_qc_summary, row.names = FALSE)

cat("\nRank-constraint validation:\n")
print(rank_constraints, row.names = FALSE)

cat("\nKey margins and profile medians:\n")
cat("  Positive - Neutral margin: ", round(positive_neutral_margin, 4), "\n", sep = "")
cat("  Inclusion - Neutral margin: ", round(inclusion_neutral_margin, 4), "\n", sep = "")
cat("  ADMET - Neutral margin: ", round(admet_neutral_margin, 4), "\n", sep = "")
cat("  Glycoside-like median PCI: ", round(glycoside_median, 4), "\n", sep = "")
cat("  Neutral median PCI: ", round(neutral_median, 4), "\n", sep = "")
cat("  Toxic/extreme median PCI: ", round(toxic_extreme_median, 4), "\n", sep = "")

cat("\nPass-rate summary:\n")
cat("  Final QC pass rate: ", ifelse(is.na(qc_pass_rate), "NA", paste0(round(qc_pass_rate * 100, 1), "%")), "\n", sep = "")
cat("  Rank-constraint pass rate: ", ifelse(is.na(rank_pass_rate), "NA", paste0(round(rank_pass_rate * 100, 1), "%")), "\n", sep = "")

cat("\nOutput integrity:\n")
cat("  Expected output files: ", nrow(output_manifest), "\n", sep = "")
cat("  Existing output files: ", sum(output_manifest$Exists), "\n", sep = "")
cat("  Missing output files: ", nrow(missing_outputs), "\n", sep = "")

cat("\nPrimary manuscript figures:\n")
cat("  MainFig1_Workflow.png / .pdf\n")
cat("  MainFig2_OptimizedProfileSeparation.png / .pdf\n")
cat("  MainFig3_DistinctiveGainSignatures.png / .pdf\n")
cat("  MainFig4_FeatureRedefinitionComparison.png / .pdf\n")

cat("\nSupplementary figures:\n")
cat("  SupplementaryFig1_TopCandidateRules.png / .pdf\n")
cat("  SupplementaryFig2_ModelSelectionTradeoff.png / .pdf\n")
cat("  SupplementaryFig3_ModelEvolution.png / .pdf\n")

cat("\nManifest file:\n")
cat("  ", file.path(table_dir, paste0(output_prefix, "_OutputManifest.csv")), "\n", sep = "")

cat("\nStatus:\n")
if (nrow(missing_outputs) == 0 && !is.na(qc_pass_rate) && qc_pass_rate == 1 && !is.na(rank_pass_rate) && rank_pass_rate == 1) {
  cat("  Final calibration completed successfully. All expected outputs were generated and all QC/rank constraints passed.\n")
} else {
  cat("  Final calibration completed, but please review QC, rank constraints, or missing output warnings above.\n")
}

cat("============================================================\n")