############################################################
## ACMPI-Nano Real 20-Compound Pairing Run
## Blueprint-guided real partner screening
##
## FULL PERFECT FIXED VERSION v4.0 | publication-safe no-clipping figures
##
## Purpose:
##   Use the Baicalein ideal nanoformulation blueprint as a fixed design target
##   and screen 20 real compounds using the unchanged ACMPI pairing algorithm.
##
## Important:
##   - The scoring algorithm is NOT changed.
##   - Virtual partner generation is removed.
##   - Real partners are loaded from MASTER_20Drugs_Model_Input.csv.
##   - Blueprint is loaded from:
##     ACMPI_Nano_Stage2_v1p2_IdealNanoformulationBlueprint_Blueprint.csv
##   - Strict QC prevents NA propagation and empty plots.
##   - Sensitivity weight names are preserved after jittering.
##   - Figure layout is upgraded for 20 compounds.
############################################################

rm(list = ls())

############################################################
## 0. User settings
############################################################

root_dir <- "/media/desk16/iy15915/中药之开创/20药物"

blueprint_file <- file.path(
  root_dir,
  "ACMPI_Nano_Stage2_v1p2_IdealNanoformulationBlueprint_Blueprint.csv"
)

real_partner_file <- file.path(
  root_dir,
  "MASTER_20Drugs_Model_Input.csv"
)

output_prefix <- "ACMPI_Nano_Real20Drugs_Pairing"

FIG_DPI <- 600
SET_SEED <- 20260429
N_SENSITIVITY <- 1000
WEIGHT_JITTER <- 0.20

CORE_COMPOUND <- "Calycosin"

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
  "ggrepel"
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
library(ggrepel)

############################################################
## 2. Input checks
############################################################

if (!dir.exists(root_dir)) {
  stop("Root directory not found: ", root_dir)
}

if (!file.exists(blueprint_file)) {
  stop("Blueprint file not found: ", blueprint_file)
}

if (!file.exists(real_partner_file)) {
  stop("Real partner file not found: ", real_partner_file)
}

############################################################
## 3. Helper functions
############################################################

read_csv_safe <- function(path) {
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8"),
    error = function(e1) {
      tryCatch(
        read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"),
        error = function(e2) {
          tryCatch(
            read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "GB18030"),
            error = function(e3) {
              read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
            }
          )
        }
      )
    }
  )
}

stop_if_missing <- function(df, cols, df_name = "data frame") {
  missing <- setdiff(cols, colnames(df))
  
  if (length(missing) > 0) {
    stop(
      df_name,
      " is missing required columns:\n",
      paste(missing, collapse = "\n"),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

to_num <- function(x) suppressWarnings(as.numeric(x))

clip01 <- function(x) {
  pmax(0, pmin(1, as.numeric(x)))
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

normalize_weights <- function(w) {
  original_names <- names(w)
  w_num <- suppressWarnings(as.numeric(w))
  
  if (length(w_num) == 0) {
    stop("Invalid weights: empty weight vector.", call. = FALSE)
  }
  
  if (any(is.na(w_num))) {
    stop("Invalid weights: NA detected.", call. = FALSE)
  }
  
  if (sum(w_num, na.rm = TRUE) <= 0) {
    stop("Invalid weights: non-positive sum.", call. = FALSE)
  }
  
  out <- w_num / sum(w_num)
  
  if (!is.null(original_names)) {
    names(out) <- original_names
  }
  
  out
}

jitter_weights <- function(w, jitter = 0.20) {
  original_names <- names(w)
  w_num <- suppressWarnings(as.numeric(w))
  
  if (length(w_num) == 0) {
    stop("Invalid weights: empty weight vector before jittering.", call. = FALSE)
  }
  
  if (any(is.na(w_num))) {
    stop("Invalid weights: NA detected before jittering.", call. = FALSE)
  }
  
  w2 <- w_num * runif(length(w_num), min = 1 - jitter, max = 1 + jitter)
  
  if (!is.null(original_names)) {
    names(w2) <- original_names
  }
  
  normalize_weights(w2)
}

weighted_sum_strict <- function(values, weights) {
  values_num <- suppressWarnings(as.numeric(values))
  weights_norm <- normalize_weights(weights)
  
  if (length(values_num) == 0 || length(weights_norm) == 0) {
    stop("weighted_sum_strict received empty values or weights.", call. = FALSE)
  }
  
  if (length(values_num) != length(weights_norm)) {
    stop(
      "weighted_sum_strict length mismatch: values = ",
      length(values_num),
      ", weights = ",
      length(weights_norm),
      call. = FALSE
    )
  }
  
  keep <- !is.na(values_num) & !is.na(weights_norm)
  
  if (!any(keep)) return(NA_real_)
  
  sum(values_num[keep] * weights_norm[keep]) / sum(weights_norm[keep])
}

get_one <- function(df, nm) {
  if (!nm %in% colnames(df)) {
    stop("Missing feature: ", nm, call. = FALSE)
  }
  as.numeric(df[[nm]][1])
}

feature_similarity <- function(a, b, scale) {
  clip01(1 - abs(as.numeric(a) - as.numeric(b)) / scale)
}

display_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x
}

short_class <- function(x) {
  dplyr::case_when(
    x == "High compatibility" ~ "High",
    x == "Moderate-high compatibility" ~ "Moderate-high",
    x == "Conditional compatibility" ~ "Conditional",
    x == "Low or incompatible" ~ "Low/incompatible",
    TRUE ~ x
  )
}

############################################################
## 4. Load blueprint and real 20-compound table
############################################################

core_blueprint <- read_csv_safe(blueprint_file)
real_raw <- read_csv_safe(real_partner_file)

if (nrow(core_blueprint) < 1) {
  stop("Blueprint table is empty.")
}

if (nrow(real_raw) < 1) {
  stop("MASTER_20Drugs_Model_Input.csv is empty.")
}

core_blueprint <- core_blueprint[1, , drop = FALSE]

############################################################
## 5. Required raw columns
############################################################

required_raw_cols <- c(
  "compound", "PubChem_CID", "SMILES",
  
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
  
  "Lipinski_violations", "Veber_violations",
  "Bioavailability_Score", "PAINS_alerts", "Brenk_alerts",
  "Synthetic_Accessibility", "QED",
  
  "LD50", "Toxicity_Class", "DILI", "Ames",
  "Carcinogenicity", "hERG",
  "ProTox_Organ_Active_Count",
  "ProTox_Endpoint_Active_Count",
  "ProTox_Metabolism_Active_Count"
)

stop_if_missing(real_raw, required_raw_cols, "MASTER_20Drugs_Model_Input.csv")

real_raw <- real_raw[, required_raw_cols, drop = FALSE]

############################################################
## 5.5 Strict input QC
############################################################

critical_numeric_cols <- c(
  "Molecular_Weight", "MolLogP", "TPSA", "HBD", "HBA",
  "Rotatable_Bonds", "Ring_Count", "Aromatic_Rings",
  "Fraction_Csp3", "Formal_Charge", "Molar_Refractivity",
  
  "Consensus_LogP", "ESOL_LogS", "ADMETlab_logS", "ADMETlab_logD",
  
  "Caco2", "PAMPA", "F30",
  
  "ADMETlab_BBB", "pgp_sub", "pgp_inh",
  
  "PPB", "Fu", "logVDss",
  
  "CYP_inhibitor_count", "t_half", "cl_plasma",
  
  "Lipinski_violations", "Veber_violations",
  "Bioavailability_Score", "PAINS_alerts", "Brenk_alerts",
  "Synthetic_Accessibility", "QED",
  
  "LD50", "Toxicity_Class", "DILI", "Ames",
  "Carcinogenicity", "hERG",
  "ProTox_Organ_Active_Count",
  "ProTox_Endpoint_Active_Count",
  "ProTox_Metabolism_Active_Count"
)

critical_character_cols <- c(
  "compound", "SMILES", "ESOL_Class", "GI_absorption",
  "BBB_permeant", "Pgp_substrate"
)

for (cc in critical_numeric_cols) {
  real_raw[[cc]] <- to_num(real_raw[[cc]])
}

############################################################
## Type-safe input QC
## Numeric and character fields are checked separately to avoid
## pivot_longer() type-combination errors.
############################################################

qc_numeric <- real_raw %>%
  dplyr::mutate(.row_id = dplyr::row_number()) %>%
  dplyr::select(.row_id, compound, dplyr::all_of(critical_numeric_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(critical_numeric_cols),
    names_to = "field",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    value_chr = as.character(value),
    is_missing = is.na(value) | value_chr == "" | value_chr == "NA"
  ) %>%
  dplyr::select(.row_id, compound, field, value_chr, is_missing)

qc_character_fields <- setdiff(critical_character_cols, "compound")

qc_character <- real_raw %>%
  dplyr::mutate(.row_id = dplyr::row_number()) %>%
  dplyr::select(.row_id, compound, dplyr::all_of(qc_character_fields)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(qc_character_fields),
    names_to = "field",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    value_chr = as.character(value),
    is_missing = is.na(value_chr) | value_chr == "" | value_chr == "NA"
  ) %>%
  dplyr::select(.row_id, compound, field, value_chr, is_missing)

qc_compound <- real_raw %>%
  dplyr::mutate(
    .row_id = dplyr::row_number(),
    field = "compound",
    value_chr = as.character(compound),
    is_missing = is.na(value_chr) | value_chr == "" | value_chr == "NA"
  ) %>%
  dplyr::select(.row_id, compound, field, value_chr, is_missing)

input_qc <- dplyr::bind_rows(qc_numeric, qc_character, qc_compound) %>%
  dplyr::filter(is_missing) %>%
  dplyr::select(.row_id, compound, field, value_chr)

if (nrow(input_qc) > 0) {
  message("\n[INPUT QC FAILED] Missing critical values detected:")
  print(input_qc)
  stop(
    "\nPlease fix MASTER_20Drugs_Model_Input.csv before running ACMPI pairing.\n",
    call. = FALSE
  )
}

if (any(duplicated(real_raw$compound))) {
  dup_names <- unique(real_raw$compound[duplicated(real_raw$compound)])
  stop(
    "Duplicate compound names detected:\n",
    paste(dup_names, collapse = "\n"),
    call. = FALSE
  )
}

if (nrow(real_raw) != 20) {
  warning(
    "Expected 20 compounds, but detected ", nrow(real_raw),
    ". The script will continue, but please verify the input table."
  )
}

cat("\n[INPUT QC OK] MASTER_20Drugs_Model_Input.csv has no missing critical fields.\n")
cat("[INPUT QC OK] Compounds detected: ", nrow(real_raw), "\n", sep = "")

############################################################
## 6. Internal scoring functions
##    Unchanged from pressure-test algorithm
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
    c(0.30, 0.25, 0.20, 0.10, 0.15)
  )
}

score_permeability_defect <- function(df) {
  s_tpsa <- risk_high(df$TPSA, good = 75, bad = 140)
  s_gi <- ifelse(tolower(trimws(df$GI_absorption)) == "high", 0.10, 0.85)
  s_caco2 <- risk_low(df$Caco2, good = -4.5, bad = -6.0)
  s_pampa <- risk_low(df$PAMPA, good = 0.70, bad = 0.20)
  s_f30 <- risk_low(df$F30, good = 0.80, bad = 0.20)
  
  row_weighted_mean(
    cbind(s_tpsa, s_gi, s_caco2, s_pampa, s_f30),
    rep(0.20, 5)
  )
}

score_barrier_transporter <- function(df) {
  s_bbb_class <- ifelse(tolower(trimws(df$BBB_permeant)) == "yes", 0.10, 0.90)
  s_bbb_prob <- risk_low(df$ADMETlab_BBB, good = 0.70, bad = 0.10)
  s_pgp_class <- ifelse(tolower(trimws(df$Pgp_substrate)) == "yes", 0.90, 0.10)
  s_pgp_sub <- risk_high(df$pgp_sub, good = 0.20, bad = 0.80)
  s_pgp_inh <- risk_high(df$pgp_inh, good = 0.20, bad = 0.80)
  
  row_weighted_mean(
    cbind(s_bbb_class, s_bbb_prob, s_pgp_class, s_pgp_sub, s_pgp_inh),
    c(0.30, 0.25, 0.20, 0.15, 0.10)
  )
}

score_distribution_exposure <- function(df) {
  s_ppb <- risk_high(df$PPB, good = 80, bad = 99)
  s_fu <- risk_low(df$Fu, good = 0.30, bad = 0.02)
  s_vd <- risk_high(abs(to_num(df$logVDss)), good = 0.8, bad = 2.0)
  
  row_weighted_mean(
    cbind(s_ppb, s_fu, s_vd),
    c(0.45, 0.35, 0.20)
  )
}

score_metabolic_liability <- function(df) {
  s_cyp <- clip01(to_num(df$CYP_inhibitor_count) / 5)
  s_half <- risk_low(df$t_half, good = 4, bad = 0.5)
  s_clearance <- risk_high(df$cl_plasma, good = 5, bad = 15)
  
  row_weighted_mean(
    cbind(s_cyp, s_half, s_clearance),
    c(0.45, 0.25, 0.30)
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
    c(0.10, 0.15, 0.10, 0.10, 0.10, 0.10, 0.15, 0.10, 0.05, 0.05)
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
    c(0.15, 0.10, 0.15, 0.20, 0.15, 0.10, 0.15)
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
    c(0.10, 0.10, 0.15, 0.12, 0.12, 0.12, 0.12, 0.10, 0.07)
  )
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
  
  out$Molecular_Size_Fit <- benefit_mid(out$Molecular_Weight, lower = 100, upper = 700)
  out$Balanced_LogP <- benefit_mid(out$MolLogP, lower = 1.5, upper = 5.0)
  out$Balanced_Consensus_LogP <- benefit_mid(out$Consensus_LogP, lower = 1.5, upper = 5.0)
  out$Lipophilicity_Fit <- benefit_mid(out$Consensus_LogP, lower = 2.0, upper = 6.0)
  
  out$Hydrophobic_Core_Fit <- row_mean_safe(cbind(
    benefit_mid(out$MolLogP, lower = 2.0, upper = 6.0),
    benefit_mid(out$Consensus_LogP, lower = 2.0, upper = 6.0),
    risk_low(out$ESOL_LogS, good = -3, bad = -7)
  ))
  
  out$Aromaticity <- benefit_high(out$Aromatic_Rings, low = 1, high = 4)
  out$Ring_Structure <- benefit_high(out$Ring_Count, low = 1, high = 4)
  out$Rigidity <- benefit_low(out$Rotatable_Bonds, high = 10, low = 0)
  
  out$Hydrogen_Bonding <- row_mean_safe(cbind(
    benefit_mid(out$HBD, lower = 1, upper = 5),
    benefit_mid(out$HBA, lower = 2, upper = 10)
  ))
  
  out$Molar_Refractivity_Fit <- benefit_mid(out$Molar_Refractivity, lower = 40, upper = 130)
  out$Charge_Compatibility <- ifelse(to_num(out$Formal_Charge) == 0, 1.00, 0.60)
  out$Druggability_Manageability <- 1 - out$Druggability_Alert_Requirement
  out$Controlled_Release_Need <- row_mean_safe(cbind(out$Metabolic_Liability_Score, out$Safety_Control_Requirement))
  
  out$Multivalent_Assembly_Drive <- row_mean_safe(cbind(
    benefit_high(out$Aromatic_Rings, low = 1.5, high = 3.0),
    benefit_mid(out$HBD, lower = 2, upper = 5),
    benefit_mid(out$HBA, lower = 4, upper = 10),
    benefit_low(out$Rotatable_Bonds, high = 6, low = 0),
    benefit_high(out$Ring_Count, low = 2, high = 4),
    benefit_mid(out$Molar_Refractivity, lower = 65, upper = 130)
  ))
  
  out$Cyclodextrin_Cavity_Match <- row_mean_safe(cbind(
    benefit_mid(out$Molecular_Weight, lower = 120, upper = 420),
    benefit_mid(out$Aromatic_Rings, lower = 1, upper = 2),
    benefit_mid(out$MolLogP, lower = 1.0, upper = 4.0),
    benefit_low(out$Rotatable_Bonds, high = 6, low = 0),
    benefit_mid(out$Molar_Refractivity, lower = 45, upper = 95)
  ))
  
  out$Small_Aromatic_Inclusion_Preference <- row_mean_safe(cbind(
    benefit_mid(out$Molecular_Weight, lower = 120, upper = 300),
    benefit_mid(out$Aromatic_Rings, lower = 1, upper = 2),
    benefit_mid(out$Ring_Count, lower = 1, upper = 3),
    benefit_mid(out$MolLogP, lower = 1.0, upper = 3.5),
    benefit_low(out$Rotatable_Bonds, high = 5, low = 0)
  ))
  
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
## 7. Core and real partner definition
############################################################

real_scored <- real_raw %>% add_internal_scores()

internal_score_cols <- c(
  "Solubility_Defect_Score",
  "Permeability_Defect_Score",
  "Barrier_Transporter_Score",
  "Distribution_Exposure_Score",
  "Metabolic_Liability_Score",
  "Nano_Assembly_Suitability_Score",
  "Druggability_Alert_Requirement",
  "Safety_Control_Requirement",
  "Nano_Optimization_Need",
  "Delivery_Defect_Burden",
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
  "Formulation_Priority_Index"
)

internal_na <- real_scored %>%
  dplyr::filter(dplyr::if_any(dplyr::all_of(internal_score_cols), is.na)) %>%
  dplyr::select(compound, dplyr::all_of(internal_score_cols))

if (nrow(internal_na) > 0) {
  message("\n[INTERNAL SCORE QC FAILED] NA values detected after internal scoring:")
  print(internal_na)
  stop(
    "\nPlease check input numeric fields before running pairing.\n",
    call. = FALSE
  )
}

if (!CORE_COMPOUND %in% real_scored$compound) {
  stop(
    "CORE_COMPOUND was not found in MASTER_20Drugs_Model_Input.csv: ",
    CORE_COMPOUND,
    "\nAvailable compounds:\n",
    paste(real_scored$compound, collapse = "\n"),
    call. = FALSE
  )
}

core_scored <- real_scored %>%
  dplyr::filter(compound == CORE_COMPOUND) %>%
  dplyr::slice(1)

real_partners <- real_scored %>%
  dplyr::filter(compound != CORE_COMPOUND)

real_partners$Partner_Profile <- "Real_Compound"
real_partners$Expected_Class <- "Real-world candidate"
real_partners$Expected_Behavior <- "Real compound screened against fixed Baicalein ideal nanoformulation blueprint."

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
    grepl("Barrier", functional, ignore.case = TRUE) ~ safe_mean(c(
      1 - get_one(partner, "Barrier_Transporter_Score"),
      1 - get_one(partner, "Safety_Control_Requirement")
    )),
    grepl("Biomimetic", functional, ignore.case = TRUE) ~ safe_mean(c(
      1 - get_one(partner, "Safety_Control_Requirement"),
      1 - get_one(partner, "Distribution_Exposure_Score"),
      get_one(partner, "Hydrogen_Bonding")
    )),
    grepl("Exposure", functional, ignore.case = TRUE) ~ safe_mean(c(
      1 - get_one(partner, "Safety_Control_Requirement"),
      1 - get_one(partner, "Distribution_Exposure_Score")
    )),
    TRUE ~ 1 - get_one(partner, "Safety_Control_Requirement")
  )
  
  clip01(
    weighted_sum_strict(
      c(primary_score, secondary_score, functional_score),
      c(0.45, 0.30, 0.25)
    )
  )
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
  
  coassembly_gain <- clip01(
    0.50 * benefit_high(get_one(partner, "Multivalent_Assembly_Drive"), low = 0.58, high = 0.88) +
      0.20 * benefit_high(abs(core_logp - partner_logp), low = 0.20, high = 1.60) +
      0.15 * benefit_high(partner_arom + partner_ring, low = 3.0, high = 5.5) +
      0.15 * benefit_high(get_one(partner, "Hydrogen_Bonding"), low = 0.50, high = 0.90)
  )
  
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
    grepl("Biomimetic", functional, ignore.case = TRUE) ~ benefit_high(safe_mean(c(
      1 - get_one(partner, "Safety_Control_Requirement"),
      get_one(partner, "Hydrogen_Bonding")
    )), low = 0.60, high = 0.90),
    grepl("Exposure", functional, ignore.case = TRUE) ~ benefit_high(1 - get_one(partner, "Distribution_Exposure_Score"), low = 0.55, high = 0.85),
    TRUE ~ benefit_high(1 - get_one(partner, "Safety_Control_Requirement"), low = 0.55, high = 0.85)
  )
  
  blueprint_specific_gain <- clip01(
    weighted_sum_strict(
      c(primary_gain, secondary_gain, functional_gain),
      c(0.45, 0.30, 0.25)
    )
  )
  
  need_values <- c(
    Solubility = get_one(core, "Solubility_Defect_Score"),
    Permeability = get_one(core, "Permeability_Defect_Score"),
    Barrier = get_one(core, "Barrier_Transporter_Score"),
    Metabolic = get_one(core, "Metabolic_Liability_Score"),
    Safety = get_one(core, "Safety_Control_Requirement")
  )
  
  support_values <- c(
    Solubility = safe_mean(c(
      1 - get_one(partner, "Solubility_Defect_Score"),
      get_one(partner, "Cyclodextrin_Cavity_Match")
    )),
    Permeability = 1 - get_one(partner, "Permeability_Defect_Score"),
    Barrier = 1 - get_one(partner, "Barrier_Transporter_Score"),
    Metabolic = 1 - get_one(partner, "Metabolic_Liability_Score"),
    Safety = 1 - get_one(partner, "Safety_Control_Requirement")
  )
  
  need_weights <- clip01(need_values)
  
  if (sum(need_weights, na.rm = TRUE) < 1e-6) {
    defect_compensation_gain <- 0.35
  } else {
    defect_compensation_gain <- clip01(
      sum(need_weights * support_values, na.rm = TRUE) /
        sum(need_weights, na.rm = TRUE)
    )
  }
  
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
  
  functional_strength <- clip01(
    weighted_sum_strict(
      c(
        coassembly_gain,
        stabilization_gain,
        blueprint_specific_gain,
        defect_compensation_gain,
        1 - get_one(partner, "Safety_Control_Requirement")
      ),
      c(0.25, 0.25, 0.25, 0.15, 0.10)
    )
  )
  
  nonredundancy_gain <- clip01((1 - similarity_to_core) * functional_strength)
  
  pairing_gain_index <- clip01(
    weighted_sum_strict(
      c(
        coassembly_gain,
        stabilization_gain,
        blueprint_specific_gain,
        defect_compensation_gain,
        nonredundancy_gain
      ),
      c(0.25, 0.20, 0.20, 0.20, 0.15)
    )
  )
  
  pairing_gain_gate <- clip01(0.45 + 0.55 * pairing_gain_index)
  
  neutral_low_value_flag <- as.numeric(
    pairing_gain_index < 0.48 &
      get_one(partner, "Safety_Control_Requirement") < 0.25 &
      get_one(partner, "Druggability_Alert_Requirement") < 0.25 &
      max(
        coassembly_gain,
        stabilization_gain,
        blueprint_specific_gain,
        defect_compensation_gain,
        na.rm = TRUE
      ) < 0.72
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
  if (is.null(weights)) {
    weights <- base_pair_weights
  }
  
  if (is.null(names(weights))) {
    names(weights) <- names(base_pair_weights)
  }
  
  weights <- normalize_weights(weights)
  
  missing_weight_names <- setdiff(names(base_pair_weights), names(weights))
  if (length(missing_weight_names) > 0) {
    stop(
      "Pairing weights are missing required names:
",
      paste(missing_weight_names, collapse = "
"),
      call. = FALSE
    )
  }
  
  weights <- weights[names(base_pair_weights)]
  
  core_logp <- get_one(core, "MolLogP")
  partner_logp <- get_one(partner, "MolLogP")
  core_tpsa <- get_one(core, "TPSA")
  partner_tpsa <- get_one(partner, "TPSA")
  core_mw <- get_one(core, "Molecular_Weight")
  partner_mw <- get_one(partner, "Molecular_Weight")
  
  interaction_balance <- safe_mean(c(
    get_one(partner, "Aromaticity"),
    get_one(partner, "Hydrogen_Bonding"),
    get_one(partner, "Rigidity"),
    get_one(partner, "Multivalent_Assembly_Drive")
  ))
  
  co_assembly <- safe_mean(c(
    interaction_balance,
    benefit_mid(abs(core_logp - partner_logp), lower = 0.0, upper = 1.8),
    benefit_mid(core_logp + partner_logp, lower = 3.0, upper = 7.5),
    get_one(partner, "Charge_Compatibility")
  ))
  
  physicochemical <- safe_mean(c(
    benefit_mid(mean(c(core_logp, partner_logp)), lower = 1.5, upper = 4.5),
    benefit_mid(mean(c(core_tpsa, partner_tpsa)), lower = 45, upper = 120),
    benefit_mid(core_mw + partner_mw, lower = 350, upper = 1000),
    benefit_mid(get_one(partner, "Rotatable_Bonds"), lower = 0, upper = 8),
    benefit_mid(get_one(partner, "Molar_Refractivity"), lower = 50, upper = 130)
  ))
  
  stabilization <- safe_mean(c(
    get_one(partner, "Cyclodextrin_Cavity_Match"),
    get_one(partner, "Small_Aromatic_Inclusion_Preference"),
    get_one(partner, "Hydrogen_Bonding"),
    get_one(partner, "Molecular_Size_Fit"),
    1 - get_one(partner, "Safety_Control_Requirement")
  ))
  
  core_sol_need <- get_one(core, "Solubility_Defect_Score")
  core_perm_need <- get_one(core, "Permeability_Defect_Score")
  core_barrier_need <- get_one(core, "Barrier_Transporter_Score")
  core_metabolic_need <- get_one(core, "Metabolic_Liability_Score")
  
  sol_comp <- ifelse(
    core_sol_need >= 0.40,
    safe_mean(c(
      1 - get_one(partner, "Solubility_Defect_Score"),
      get_one(partner, "Cyclodextrin_Cavity_Match")
    )),
    0.50
  )
  
  perm_comp <- ifelse(
    core_perm_need >= 0.40,
    1 - get_one(partner, "Permeability_Defect_Score"),
    0.50
  )
  
  barrier_comp <- ifelse(
    core_barrier_need >= 0.40,
    1 - get_one(partner, "Barrier_Transporter_Score"),
    0.50
  )
  
  metabolic_comp <- ifelse(
    core_metabolic_need >= 0.40,
    1 - get_one(partner, "Metabolic_Liability_Score"),
    0.50
  )
  
  admet <- safe_mean(c(
    sol_comp,
    perm_comp,
    barrier_comp,
    metabolic_comp,
    1 - get_one(partner, "Distribution_Exposure_Score")
  ))
  
  safety <- safe_mean(c(
    1 - get_one(partner, "Safety_Control_Requirement"),
    1 - get_one(partner, "Druggability_Alert_Requirement"),
    1 - clip01(get_one(partner, "DILI")),
    1 - clip01(get_one(partner, "Ames")),
    1 - clip01(get_one(partner, "hERG"))
  ))
  
  blueprint_match <- score_blueprint_match(partner, core_blueprint)
  
  module_values <- c(
    Co_assembly_Compatibility = co_assembly,
    Physicochemical_Complementarity = physicochemical,
    Nano_stabilization_Complementarity = stabilization,
    ADMET_Complementarity = admet,
    Safety_Balance = safety,
    Blueprint_Match = blueprint_match
  )
  
  raw_pci <- weighted_sum_strict(
    module_values[names(weights)],
    weights[names(weights)]
  )
  
  gain <- score_pairing_gain(core, partner, core_blueprint)
  
  toxic_red_flag <- as.numeric(
    get_one(partner, "Safety_Control_Requirement") >= 0.65 |
      get_one(partner, "DILI") >= 0.70 |
      get_one(partner, "Ames") >= 0.70 |
      get_one(partner, "hERG") >= 0.70 |
      get_one(partner, "LD50") <= 300
  )
  
  extreme_physchem_flag <- as.numeric(
    partner_tpsa >= 170 |
      partner_logp >= 6.5 |
      partner_logp <= 0 |
      partner_mw >= 800
  )
  
  toxic_gate <- ifelse(toxic_red_flag == 1, 0.48, 1.00)
  extreme_physchem_gate <- ifelse(extreme_physchem_flag == 1, 0.72, 1.00)
  
  final_pci <- raw_pci *
    gain$Pairing_Gain_Gate[1] *
    gain$Neutral_Low_Value_Gate[1] *
    toxic_gate *
    extreme_physchem_gate
  
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
  
  bind_rows(res) %>%
    dplyr::select(
      Core_Compound,
      Partner_Compound,
      Partner_Profile,
      Expected_Class,
      Expected_Behavior,
      dplyr::everything()
    )
}

############################################################
## 9. Real partner scoring
############################################################

pairing_scores <- score_all_partners(
  core = core_scored,
  partners = real_partners,
  core_blueprint = core_blueprint
)

pairing_scores <- pairing_scores %>%
  dplyr::arrange(dplyr::desc(Pairing_Compatibility_Index))

############################################################
## 10. Sensitivity analysis
############################################################

run_pairing_sensitivity_one <- function(core, partner, core_blueprint,
                                        n_iter = 1000,
                                        jitter = 0.20) {
  vals <- numeric(n_iter)
  
  for (i in seq_len(n_iter)) {
    w <- jitter_weights(base_pair_weights, jitter = jitter)
    
    tmp <- score_pair_one(
      core = core,
      partner = partner,
      core_blueprint = core_blueprint,
      weights = w
    )
    
    vals[i] <- tmp$Pairing_Compatibility_Index[1]
  }
  
  data.frame(
    Sensitivity_Mean = mean(vals, na.rm = TRUE),
    Sensitivity_SD = sd(vals, na.rm = TRUE),
    Sensitivity_P05 = as.numeric(quantile(vals, 0.05, na.rm = TRUE)),
    Sensitivity_P95 = as.numeric(quantile(vals, 0.95, na.rm = TRUE)),
    Probability_High_or_ModerateHigh = mean(vals >= 0.60, na.rm = TRUE),
    Probability_Conditional_or_Better = mean(vals >= 0.45, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

sensitivity_list <- list()

for (i in seq_len(nrow(real_partners))) {
  partner <- real_partners[i, , drop = FALSE]
  
  sens <- run_pairing_sensitivity_one(
    core = core_scored,
    partner = partner,
    core_blueprint = core_blueprint,
    n_iter = N_SENSITIVITY,
    jitter = WEIGHT_JITTER
  )
  
  sens$Partner_Compound <- partner$compound[1]
  sensitivity_list[[i]] <- sens
}

sensitivity_df <- bind_rows(sensitivity_list)

pairing_scores <- pairing_scores %>%
  dplyr::left_join(sensitivity_df, by = "Partner_Compound") %>%
  dplyr::arrange(dplyr::desc(Pairing_Compatibility_Index))

############################################################
## 10.5 Strict output QC before plotting
############################################################

result_qc <- pairing_scores %>%
  dplyr::filter(
    is.na(Pairing_Compatibility_Index) |
      is.na(Sensitivity_Mean) |
      is.na(Sensitivity_P05) |
      is.na(Sensitivity_P95)
  ) %>%
  dplyr::select(
    Partner_Compound,
    Pairing_Compatibility_Index,
    Sensitivity_Mean,
    Sensitivity_P05,
    Sensitivity_P95,
    Raw_PCI,
    Pairing_Gain_Index,
    Toxic_Red_Flag,
    Toxic_Gate,
    Extreme_Physicochemical_Flag,
    Extreme_Physicochemical_Gate
  )

if (nrow(result_qc) > 0) {
  message("\n[OUTPUT QC FAILED] NA values detected in final results:")
  print(result_qc)
  stop(
    "\nStop before plotting: PCI or sensitivity results contain NA.\n",
    "Please check MASTER_20Drugs_Model_Input.csv and rerun.\n",
    call. = FALSE
  )
}

cat("\n[OUTPUT QC OK] No NA values in PCI or sensitivity results.\n")

############################################################
## 11. Output folders
############################################################

timestamp_final <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_dir <- file.path(
  root_dir,
  paste0(output_prefix, "_Output_", timestamp_final)
)

table_dir <- file.path(out_dir, "01_Tables")
fig_png_dir <- file.path(out_dir, "02_Figures_PNG")
fig_pdf_dir <- file.path(out_dir, "03_Figures_PDF")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_pdf_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 12. Save tables
############################################################

out_score_csv <- file.path(table_dir, "ACMPI_Real20Drugs_PairingScores.csv")
out_score_xlsx <- file.path(table_dir, "ACMPI_Real20Drugs_PairingScores.xlsx")

write.csv(
  pairing_scores,
  out_score_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

openxlsx::write.xlsx(
  pairing_scores,
  out_score_xlsx,
  overwrite = TRUE
)

out_partner_csv <- file.path(table_dir, "ACMPI_Real20Drugs_InternalScoredPartners.csv")
out_core_csv <- file.path(table_dir, "ACMPI_Core_InternalScored.csv")

write.csv(
  real_partners,
  out_partner_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  core_scored,
  out_core_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

############################################################
## 13. Main plots
##     PERFECT main-figure system for 20-compound library
##
##     Figure strategy:
##       Fig1  Full-library final PCI ranking
##       Fig2  Top-5 module-level signature heatmap
##       Fig3  Top-5 weight-jitter robustness
##       Fig4  Raw-to-final gate-effect map
##       Fig5  Experimental decision landscape
##
##     Notes:
##       - The ACMPI scoring algorithm is not changed.
##       - Only figure layout, labeling, and visual encoding are upgraded.
##       - All figures use conservative, low-saturation publication palettes.
##       - Dynamic width/height and wrapped labels are used to prevent clipping.
############################################################

## ---------- 13.1 Figure helper functions ----------

safe_wrap <- function(x, width = 28) {
  vapply(
    as.character(x),
    function(z) stringr::str_wrap(z, width = width),
    character(1)
  )
}

fmt_num <- function(x, digits = 3) {
  formatC(as.numeric(x), format = "f", digits = digits)
}

clean_class <- function(x) {
  dplyr::case_when(
    x == "High compatibility" ~ "High",
    x == "Moderate-high compatibility" ~ "Moderate-high",
    x == "Conditional compatibility" ~ "Conditional",
    x == "Low or incompatible" ~ "Low / incompatible",
    TRUE ~ as.character(x)
  )
}

class_levels <- c("High", "Moderate-high", "Conditional", "Low / incompatible")

class_colors <- c(
  "High" = "#3F6C51",
  "Moderate-high" = "#4C78A8",
  "Conditional" = "#E6A65D",
  "Low / incompatible" = "#B9C3B0"
)

gate_colors <- c(
  "No penalty" = "#4C78A8",
  "Toxicity gate" = "#D88C63",
  "Physicochemical gate" = "#8AA6A3",
  "Dual gate" = "#7D6A91"
)

module_palette <- c(
  "0" = "#F7F8F3",
  "1" = "#315A7D"
)

theme_acmpi <- function(base_size = 11.5, legend_position = "bottom") {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "#2C2C2C"),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 2.2,
        hjust = 0,
        lineheight = 1.05,
        margin = ggplot2::margin(b = 5)
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 0.3,
        color = "#5A5A5A",
        hjust = 0,
        lineheight = 1.08,
        margin = ggplot2::margin(b = 9)
      ),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.8, color = "#333333"),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1.0),
      legend.text = ggplot2::element_text(size = base_size - 1.3),
      legend.position = legend_position,
      legend.box = "horizontal",
      legend.margin = ggplot2::margin(t = 2, r = 0, b = 0, l = 0),
      legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.margin = ggplot2::margin(24, 44, 26, 38)
    )
}

save_plot_dual <- function(plot, png_path, pdf_path, width, height, dpi = 600) {
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  
  tryCatch(
    ggplot2::ggsave(
      filename = pdf_path,
      plot = plot,
      width = width,
      height = height,
      device = cairo_pdf,
      bg = "white",
      limitsize = FALSE
    ),
    error = function(e) {
      ggplot2::ggsave(
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

pairing_scores_plot <- pairing_scores %>%
  dplyr::mutate(
    Partner_Compound_Display = display_name(Partner_Compound),
    Partner_Compound_Label = safe_wrap(Partner_Compound_Display, width = 23),
    Pairing_Compatibility_Class_Display = clean_class(Pairing_Compatibility_Class),
    Pairing_Compatibility_Class_Display = factor(
      Pairing_Compatibility_Class_Display,
      levels = class_levels
    ),
    Gate_Category = dplyr::case_when(
      Toxic_Red_Flag == 1 & Extreme_Physicochemical_Flag == 1 ~ "Dual gate",
      Toxic_Red_Flag == 1 & Extreme_Physicochemical_Flag == 0 ~ "Toxicity gate",
      Toxic_Red_Flag == 0 & Extreme_Physicochemical_Flag == 1 ~ "Physicochemical gate",
      TRUE ~ "No penalty"
    ),
    Gate_Category = factor(
      Gate_Category,
      levels = c("No penalty", "Toxicity gate", "Physicochemical gate", "Dual gate")
    ),
    PCI_Label = fmt_num(Pairing_Compatibility_Index, 3),
    Raw_PCI_Label = fmt_num(Raw_PCI, 3),
    Rank = dplyr::row_number()
  )

n_plot <- nrow(pairing_scores_plot)
top_n <- min(5, n_plot)

top_candidates <- pairing_scores_plot %>%
  dplyr::arrange(dplyr::desc(Pairing_Compatibility_Index)) %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::pull(Partner_Compound_Display)

priority_candidates <- pairing_scores_plot %>%
  dplyr::arrange(dplyr::desc(Pairing_Compatibility_Index)) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::pull(Partner_Compound_Display)

rank_width <- 8.2
rank_height <- max(7.4, 2.3 + 0.31 * n_plot)

heat_width <- 8.4
heat_height <- 4.9

sens_width <- 7.3
sens_height <- 4.6

scatter_width <- 8.8
scatter_height <- 6.2

decision_width <- 8.8
decision_height <- 6.2

## ---------- 13.2 Fig1. Full-library ranking ----------

plot_df <- pairing_scores_plot %>%
  dplyr::mutate(
    Partner_Compound_Label = factor(
      Partner_Compound_Label,
      levels = rev(Partner_Compound_Label)
    )
  )

p_rank <- ggplot(
  plot_df,
  aes(
    x = Partner_Compound_Label,
    y = Pairing_Compatibility_Index,
    fill = Pairing_Compatibility_Class_Display
  )
) +
  geom_vline(xintercept = 0, color = NA) +
  geom_hline(yintercept = c(0.45, 0.60, 0.75), linetype = "dashed", linewidth = 0.35, color = "#8C8C8C") +
  geom_col(width = 0.70, color = "#333333", linewidth = 0.20) +
  geom_text(
    aes(label = PCI_Label),
    hjust = -0.12,
    size = 3.2,
    color = "#333333"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 0.72),
    breaks = seq(0, 0.70, 0.10),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "ACMPI pairing compatibility",
    subtitle = paste0("Core: ", CORE_COMPOUND, " | Target: Baicalein-derived nanoformulation blueprint"),
    x = NULL,
    y = "Final pairing compatibility index",
    fill = "Class"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_acmpi(base_size = 11.3, legend_position = "bottom") +
  theme(
    axis.text.y = ggplot2::element_text(size = 9.5, lineheight = 0.94),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
    plot.margin = ggplot2::margin(20, 64, 22, 32)
  )

save_plot_dual(
  p_rank,
  file.path(fig_png_dir, "Fig1_ACMPI_FinalPCI_FullRanking.png"),
  file.path(fig_pdf_dir, "Fig1_ACMPI_FinalPCI_FullRanking.pdf"),
  width = rank_width,
  height = rank_height,
  dpi = FIG_DPI
)

## ---------- 13.3 Fig2. Top-5 module-level heatmap ----------

module_cols <- c(
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
  "Nonredundancy_Gain"
)

module_labels <- c(
  Co_assembly_Compatibility = "Co-assembly",
  Physicochemical_Complementarity = "Physicochemical",
  Nano_stabilization_Complementarity = "Nano-stabilization",
  ADMET_Complementarity = "ADMET",
  Safety_Balance = "Safety",
  Blueprint_Match = "Blueprint",
  Coassembly_Gain = "Coassembly gain",
  Stabilization_Gain = "Stabilization gain",
  Blueprint_Specific_Gain = "Blueprint gain",
  Defect_Compensation_Gain = "Defect compensation",
  Nonredundancy_Gain = "Nonredundancy"
)

module_group <- c(
  Co_assembly_Compatibility = "Compatibility modules",
  Physicochemical_Complementarity = "Compatibility modules",
  Nano_stabilization_Complementarity = "Compatibility modules",
  ADMET_Complementarity = "Compatibility modules",
  Safety_Balance = "Compatibility modules",
  Blueprint_Match = "Compatibility modules",
  Coassembly_Gain = "Gain modules",
  Stabilization_Gain = "Gain modules",
  Blueprint_Specific_Gain = "Gain modules",
  Defect_Compensation_Gain = "Gain modules",
  Nonredundancy_Gain = "Gain modules"
)

heat_df <- pairing_scores_plot %>%
  dplyr::arrange(dplyr::desc(Pairing_Compatibility_Index)) %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::select(Partner_Compound_Display, dplyr::all_of(module_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(module_cols),
    names_to = "Module_Raw",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Partner_Compound_Display = factor(
      Partner_Compound_Display,
      levels = rev(top_candidates)
    ),
    Module = factor(
      module_labels[Module_Raw],
      levels = module_labels[module_cols]
    ),
    Module_Group = factor(
      module_group[Module_Raw],
      levels = c("Compatibility modules", "Gain modules")
    ),
    Score_Label = fmt_num(Score, 2)
  )

p_heat <- ggplot(
  heat_df,
  aes(x = Module, y = Partner_Compound_Display, fill = Score)
) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = Score_Label), size = 2.7, color = "#2B2B2B") +
  facet_grid(. ~ Module_Group, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colors = c("#F7F8F3", "#D7E2E6", "#8FB3C5", "#315A7D"),
    values = scales::rescale(c(0, 0.45, 0.70, 1)),
    limits = c(0, 1),
    oob = scales::squish
  ) +
  labs(
    title = "Top-candidate module signatures",
    subtitle = "Top 5 partners ranked by final PCI; tile values show module scores",
    x = NULL,
    y = NULL,
    fill = "Score"
  ) +
  theme_acmpi(base_size = 10.7, legend_position = "right") +
  theme(
    strip.background = ggplot2::element_rect(fill = "#F1F1EC", color = NA),
    strip.text = ggplot2::element_text(face = "bold", size = 9.8, color = "#333333"),
    axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, vjust = 1, size = 8.5),
    axis.text.y = ggplot2::element_text(size = 9.7),
    panel.spacing.x = grid::unit(0.45, "lines"),
    plot.margin = ggplot2::margin(18, 24, 34, 28)
  )

save_plot_dual(
  p_heat,
  file.path(fig_png_dir, "Fig2_ACMPI_Top5_ModuleSignature_Heatmap.png"),
  file.path(fig_pdf_dir, "Fig2_ACMPI_Top5_ModuleSignature_Heatmap.pdf"),
  width = heat_width,
  height = heat_height,
  dpi = FIG_DPI
)

## ---------- 13.4 Fig3. Top-5 sensitivity robustness ----------

sens_plot_df <- pairing_scores_plot %>%
  dplyr::arrange(dplyr::desc(Pairing_Compatibility_Index)) %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::mutate(
    Partner_Compound_Display = factor(
      Partner_Compound_Display,
      levels = rev(top_candidates)
    )
  )

p_sens <- ggplot(sens_plot_df, aes(x = Partner_Compound_Display)) +
  geom_hline(yintercept = c(0.45, 0.60), linetype = "dashed", linewidth = 0.35, color = "#9A9A9A") +
  geom_errorbar(
    aes(ymin = Sensitivity_P05, ymax = Sensitivity_P95),
    width = 0.22,
    linewidth = 0.72,
    color = "#5A5A5A"
  ) +
  geom_point(aes(y = Sensitivity_Mean), size = 3.2, color = "#315A7D") +
  geom_point(aes(y = Pairing_Compatibility_Index), size = 2.9, shape = 23, fill = "#E6A65D", color = "#333333", stroke = 0.4) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(0.40, 0.67),
    breaks = seq(0.40, 0.66, 0.05),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = "Robustness of prioritized candidates",
    subtitle = paste0("Weight jitter: ±", WEIGHT_JITTER * 100, "% | N = ", N_SENSITIVITY, " iterations per partner"),
    x = NULL,
    y = "PCI under weight perturbation"
  ) +
  theme_acmpi(base_size = 11.0, legend_position = "none") +
  theme(
    axis.text.y = ggplot2::element_text(size = 10.2),
    plot.margin = ggplot2::margin(20, 28, 22, 32)
  )

save_plot_dual(
  p_sens,
  file.path(fig_png_dir, "Fig3_ACMPI_Top5_WeightSensitivity.png"),
  file.path(fig_pdf_dir, "Fig3_ACMPI_Top5_WeightSensitivity.pdf"),
  width = sens_width,
  height = sens_height,
  dpi = FIG_DPI
)

## ---------- 13.5 Fig4. Raw-to-final gate effect ----------

label_gate_compounds <- pairing_scores_plot %>%
  dplyr::filter(
    Partner_Compound_Display %in% unique(c(
      top_candidates,
      "Quercetin", "Kaempferol", "Luteolin", "Ginsenoside Rg3", "Glycyrrhizin", "Saikosaponin A"
    ))
  )

p_gate <- ggplot(
  pairing_scores_plot,
  aes(
    x = Raw_PCI,
    y = Pairing_Compatibility_Index,
    color = Gate_Category,
    shape = Gate_Category
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.55, color = "#7F7F7F") +
  geom_hline(yintercept = c(0.45, 0.60), linetype = "dashed", linewidth = 0.35, color = "#9A9A9A") +
  geom_point(size = 3.2, alpha = 0.92, stroke = 0.45) +
  ggrepel::geom_text_repel(
    data = label_gate_compounds,
    aes(label = safe_wrap(Partner_Compound_Display, width = 14)),
    size = 2.9,
    color = "#2F2F2F",
    min.segment.length = 0,
    segment.color = "#B0B0B0",
    segment.linewidth = 0.25,
    box.padding = 0.35,
    point.padding = 0.22,
    force = 1.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = gate_colors, drop = FALSE) +
  scale_shape_manual(values = c("No penalty" = 16, "Toxicity gate" = 17, "Physicochemical gate" = 15, "Dual gate" = 18), drop = FALSE) +
  scale_x_continuous(limits = c(0.45, 0.90), breaks = seq(0.45, 0.90, 0.10), expand = expansion(mult = c(0.02, 0.13))) +
  scale_y_continuous(limits = c(0.05, 0.68), breaks = seq(0.10, 0.65, 0.10), expand = expansion(mult = c(0.02, 0.06))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Gate-aware prioritization",
    subtitle = "Raw compatibility converted to gate-aware final PCI",
    x = "Raw PCI before gates",
    y = "Final PCI after gates",
    color = "Gate status",
    shape = "Gate status"
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_acmpi(base_size = 10.8, legend_position = "bottom") +
  theme(
    legend.box = "vertical",
    plot.margin = ggplot2::margin(22, 58, 28, 34)
  )

save_plot_dual(
  p_gate,
  file.path(fig_png_dir, "Fig4_ACMPI_RawToFinal_GateEffect.png"),
  file.path(fig_pdf_dir, "Fig4_ACMPI_RawToFinal_GateEffect.pdf"),
  width = scatter_width,
  height = scatter_height,
  dpi = FIG_DPI
)

## ---------- 13.6 Fig5. Decision landscape ----------

decision_df <- pairing_scores_plot %>%
  dplyr::mutate(
    Decision_Tier = dplyr::case_when(
      Partner_Compound_Display %in% priority_candidates[1] ~ "Primary candidate",
      Partner_Compound_Display %in% priority_candidates[-1] ~ "Validation candidate",
      Pairing_Compatibility_Index >= 0.45 ~ "Reserve candidate",
      TRUE ~ "Not prioritized"
    ),
    Decision_Tier = factor(
      Decision_Tier,
      levels = c("Primary candidate", "Validation candidate", "Reserve candidate", "Not prioritized")
    )
  )

decision_colors <- c(
  "Primary candidate" = "#3F6C51",
  "Validation candidate" = "#4C78A8",
  "Reserve candidate" = "#E6A65D",
  "Not prioritized" = "#B9C3B0"
)

label_decision <- decision_df %>%
  dplyr::filter(Decision_Tier %in% c("Primary candidate", "Validation candidate") |
                  Partner_Compound_Display %in% c("Glycyrrhetinic Acid", "Curcumin"))

p_decision <- ggplot(
  decision_df,
  aes(
    x = Pairing_Gain_Index,
    y = Pairing_Compatibility_Index,
    size = Raw_PCI,
    fill = Decision_Tier
  )
) +
  geom_hline(yintercept = c(0.45, 0.60), linetype = "dashed", linewidth = 0.35, color = "#9A9A9A") +
  geom_vline(xintercept = 0.50, linetype = "dashed", linewidth = 0.35, color = "#9A9A9A") +
  geom_point(shape = 21, color = "#333333", stroke = 0.35, alpha = 0.92) +
  ggrepel::geom_text_repel(
    data = label_decision,
    aes(label = safe_wrap(Partner_Compound_Display, width = 14)),
    size = 2.9,
    color = "#2F2F2F",
    min.segment.length = 0,
    segment.color = "#B0B0B0",
    segment.linewidth = 0.25,
    box.padding = 0.35,
    point.padding = 0.22,
    force = 1.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = decision_colors, drop = FALSE) +
  scale_size_continuous(range = c(2.5, 6.2), limits = range(decision_df$Raw_PCI, na.rm = TRUE)) +
  scale_x_continuous(limits = c(0.20, 0.62), breaks = seq(0.20, 0.60, 0.10), expand = expansion(mult = c(0.02, 0.16))) +
  scale_y_continuous(limits = c(0.05, 0.68), breaks = seq(0.10, 0.65, 0.10), expand = expansion(mult = c(0.02, 0.06))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Experimental decision landscape",
    subtitle = "Gain, final PCI, raw potential, and decision tiers",
    x = "Pairing gain index",
    y = "Final PCI",
    fill = "Decision tier",
    size = "Raw PCI"
  ) +
  guides(
    fill = guide_legend(nrow = 2, byrow = TRUE),
    size = guide_legend(nrow = 1)
  ) +
  theme_acmpi(base_size = 10.8, legend_position = "bottom") +
  theme(
    legend.box = "vertical",
    plot.margin = ggplot2::margin(22, 58, 28, 34)
  )

save_plot_dual(
  p_decision,
  file.path(fig_png_dir, "Fig5_ACMPI_DecisionLandscape.png"),
  file.path(fig_pdf_dir, "Fig5_ACMPI_DecisionLandscape.pdf"),
  width = decision_width,
  height = decision_height,
  dpi = FIG_DPI
)

## ---------- 13.7 Decision table for writing ----------

decision_table <- decision_df %>%
  dplyr::select(
    Rank,
    Partner_Compound,
    Pairing_Compatibility_Index,
    Pairing_Compatibility_Class,
    Raw_PCI,
    Pairing_Gain_Index,
    Sensitivity_Mean,
    Sensitivity_P05,
    Sensitivity_P95,
    Probability_High_or_ModerateHigh,
    Probability_Conditional_or_Better,
    Toxic_Red_Flag,
    Extreme_Physicochemical_Flag,
    Gate_Category,
    Decision_Tier
  ) %>%
  dplyr::arrange(Rank)

write.csv(
  decision_table,
  file.path(table_dir, "ACMPI_Real20Drugs_DecisionTable.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

openxlsx::write.xlsx(
  decision_table,
  file.path(table_dir, "ACMPI_Real20Drugs_DecisionTable.xlsx"),
  overwrite = TRUE
)


############################################################
## 14. Final summary
############################################################

cat("\n============================================================\n")
cat("ACMPI-Nano real 20-compound pairing run completed.\n")
cat("Root directory: ", root_dir, "\n", sep = "")
cat("Blueprint file: ", blueprint_file, "\n", sep = "")
cat("Real partner file: ", real_partner_file, "\n", sep = "")
cat("Core compound used: ", CORE_COMPOUND, "\n", sep = "")
cat("Total compounds in library: ", nrow(real_raw), "\n", sep = "")
cat("Real partners scored: ", nrow(real_partners), "\n", sep = "")
cat("Sensitivity iterations per compound: ", N_SENSITIVITY, "\n", sep = "")
cat("Weight jitter: ±", WEIGHT_JITTER * 100, "%\n", sep = "")

cat("\nTop-ranked candidates:\n")
print(
  pairing_scores[, c(
    "Partner_Compound",
    "Pairing_Compatibility_Index",
    "Pairing_Compatibility_Class",
    "Raw_PCI",
    "Pairing_Gain_Index",
    "Toxic_Red_Flag",
    "Extreme_Physicochemical_Flag",
    "Sensitivity_Mean",
    "Sensitivity_P05",
    "Sensitivity_P95"
  )],
  row.names = FALSE
)

cat("\nSaved outputs:\n")
cat(out_score_csv, "\n")
cat(out_score_xlsx, "\n")
cat(out_partner_csv, "\n")
cat(out_core_csv, "\n")
cat(fig_png_dir, "\n")
cat(fig_pdf_dir, "\n")
cat("============================================================\n")