############################################################
## ACMPI-Nano v1.0
## Stage 1: Core Compound Nanoformulation Value Model
## FULL FIXED VERSION
##
## Key fixes:
##   - Creates a timestamped output folder automatically
##   - Saves all tables and figures into separated subfolders
##   - Exports high-resolution PNG and PDF figures
##   - Uses top-journal pastel color palette
##   - Prevents text overflow by wrapping titles/subtitles
##   - Enlarges figure canvas and margins
##   - Generates per-compound radar and design-driver plots
##   - Generates one global nanoformulation value map
############################################################

rm(list = ls())

############################################################
## 0. User settings
############################################################

target_dir <- "/media/desk16/iy15915/中药之开创/Baicalein"
input_file <- file.path(target_dir, "Baicalein_Master_Model_Input.csv")

output_prefix <- "ACMPI_Nano_Stage1_CoreCompound_NanoformulationValue"

SAVE_OUTPUT <- TRUE
FIG_DPI <- 600

############################################################
## 1. Required packages
############################################################

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "openxlsx",
  "scales",
  "stringr"
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

############################################################
## 2. Output folders
############################################################

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_dir <- file.path(
  target_dir,
  paste0("ACMPI_Nano_Stage1_Output_", timestamp)
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
## 4. Expected official columns
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

class_score <- function(x, levels_map, default = 0.5) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(default, length(x))
  
  for (nm in names(levels_map)) {
    out[x == tolower(nm)] <- levels_map[[nm]]
  }
  
  out
}

score_to_class <- function(x, low_cut = 0.33, high_cut = 0.66) {
  cut(
    x,
    breaks = c(-Inf, low_cut, high_cut, Inf),
    labels = c("Low", "Moderate", "High")
  )
}

clean_module_label <- function(x) {
  x <- gsub("_Score$", "", x)
  x <- gsub("_Requirement$", "", x)
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

############################################################
## 6. Module scoring functions
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

############################################################
## 7. Stage 1 model
############################################################

run_stage1_nano_value_model <- function(df) {
  
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
  
  out$Design_Control_Requirement_Class <- dplyr::case_when(
    out$Safety_Control_Requirement >= 0.66 ~ "High",
    out$Safety_Control_Requirement >= 0.33 ~ "Moderate",
    TRUE ~ "Low"
  )
  
  out$Nano_Optimization_Need_Class <- as.character(score_to_class(out$Nano_Optimization_Need))
  out$Nano_Assembly_Suitability_Class <- as.character(score_to_class(out$Nano_Assembly_Suitability_Score))
  out$Druggability_Alert_Class <- as.character(score_to_class(out$Druggability_Alert_Requirement))
  
  bottleneck_mat <- out[, c(
    "Solubility_Defect_Score",
    "Permeability_Defect_Score",
    "Barrier_Transporter_Score",
    "Distribution_Exposure_Score",
    "Metabolic_Liability_Score",
    "Safety_Control_Requirement",
    "Druggability_Alert_Requirement"
  )]
  
  out$Main_Design_Bottleneck <- colnames(bottleneck_mat)[
    max.col(bottleneck_mat, ties.method = "first")
  ]
  
  driver_mat <- out[, c(
    "Nano_Optimization_Need",
    "Nano_Assembly_Suitability_Score",
    "Delivery_Defect_Burden",
    "Safety_Control_Requirement"
  )]
  
  out$Main_Nanoformulation_Value_Driver <- colnames(driver_mat)[
    max.col(driver_mat, ties.method = "first")
  ]
  
  advantage_mat <- data.frame(
    Assembly = out$Nano_Assembly_Suitability_Score,
    Aromaticity = benefit_high(out$Aromatic_Rings, low = 1, high = 4),
    Rigidity = benefit_low(out$Rotatable_Bonds, high = 10, low = 0),
    Hydrogen_Bonding = rowMeans(
      cbind(
        benefit_mid(out$HBD, lower = 1, upper = 5),
        benefit_mid(out$HBA, lower = 2, upper = 10)
      ),
      na.rm = TRUE
    ),
    Balanced_Lipophilicity = benefit_mid(out$MolLogP, lower = 1.5, upper = 5.0),
    stringsAsFactors = FALSE
  )
  
  out$Main_Structural_Advantage <- colnames(advantage_mat)[
    max.col(advantage_mat, ties.method = "first")
  ]
  
  out$Requires_Solubility_Enhancement <- out$Solubility_Defect_Score >= 0.50
  out$Requires_Permeability_Enhancement <- out$Permeability_Defect_Score >= 0.50
  out$Requires_Barrier_or_Transporter_Strategy <- out$Barrier_Transporter_Score >= 0.50
  out$Requires_Exposure_Optimization <- out$Distribution_Exposure_Score >= 0.50
  out$Requires_Metabolic_Control <- out$Metabolic_Liability_Score >= 0.50
  out$Requires_Safety_Control <- out$Safety_Control_Requirement >= 0.50
  out$Has_Strong_Nano_Assembly_Advantage <- out$Nano_Assembly_Suitability_Score >= 0.66
  
  out
}

############################################################
## 8. Run model
############################################################

stage1 <- run_stage1_nano_value_model(dat)

############################################################
## 9. Summary tables
############################################################

stage1_summary_cols <- c(
  "compound", "PubChem_CID", "InChIKey",
  
  "Solubility_Defect_Score",
  "Permeability_Defect_Score",
  "Barrier_Transporter_Score",
  "Distribution_Exposure_Score",
  "Metabolic_Liability_Score",
  "Nano_Assembly_Suitability_Score",
  "Druggability_Alert_Requirement",
  "Safety_Control_Requirement",
  
  "Delivery_Defect_Burden",
  "Nano_Optimization_Need",
  "Nanoformulation_Value_Index",
  "Nanoformulation_Development_Value",
  "Design_Control_Requirement_Class",
  
  "Nano_Optimization_Need_Class",
  "Nano_Assembly_Suitability_Class",
  "Druggability_Alert_Class",
  
  "Main_Nanoformulation_Value_Driver",
  "Main_Design_Bottleneck",
  "Main_Structural_Advantage",
  
  "Requires_Solubility_Enhancement",
  "Requires_Permeability_Enhancement",
  "Requires_Barrier_or_Transporter_Strategy",
  "Requires_Exposure_Optimization",
  "Requires_Metabolic_Control",
  "Requires_Safety_Control",
  "Has_Strong_Nano_Assembly_Advantage"
)

stage1_summary <- stage1[, stage1_summary_cols]

numeric_cols <- sapply(stage1_summary, is.numeric)
stage1_summary[, numeric_cols] <- lapply(stage1_summary[, numeric_cols], function(x) round(x, 4))

module_score_cols <- c(
  "Solubility_Defect_Score",
  "Permeability_Defect_Score",
  "Barrier_Transporter_Score",
  "Distribution_Exposure_Score",
  "Metabolic_Liability_Score",
  "Nano_Assembly_Suitability_Score",
  "Druggability_Alert_Requirement",
  "Safety_Control_Requirement",
  "Delivery_Defect_Burden",
  "Nano_Optimization_Need",
  "Nanoformulation_Value_Index"
)

module_table <- stage1 %>%
  dplyr::select(compound, dplyr::all_of(module_score_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(module_score_cols),
    names_to = "module",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    score = round(as.numeric(score), 4),
    module_label = clean_module_label(module)
  )

############################################################
## 10. Figure data
############################################################

driver_score_cols <- c(
  "Solubility_Defect_Score",
  "Permeability_Defect_Score",
  "Barrier_Transporter_Score",
  "Distribution_Exposure_Score",
  "Metabolic_Liability_Score",
  "Nano_Assembly_Suitability_Score",
  "Druggability_Alert_Requirement",
  "Safety_Control_Requirement"
)

driver_label_map <- c(
  "Solubility_Defect_Score" = "Solubility defect",
  "Permeability_Defect_Score" = "Permeability defect",
  "Barrier_Transporter_Score" = "Barrier / transporter",
  "Distribution_Exposure_Score" = "Distribution / exposure",
  "Metabolic_Liability_Score" = "Metabolic liability",
  "Nano_Assembly_Suitability_Score" = "Nano-assembly suitability",
  "Druggability_Alert_Requirement" = "Druggability alert",
  "Safety_Control_Requirement" = "Design-control requirement"
)

value_map_df <- data.frame(
  compound = stage1$compound,
  Nano_Optimization_Need = stage1$Nano_Optimization_Need,
  Nano_Assembly_Suitability_Score = stage1$Nano_Assembly_Suitability_Score,
  Safety_Control_Requirement = stage1$Safety_Control_Requirement,
  Nanoformulation_Value_Index = stage1$Nanoformulation_Value_Index,
  Nanoformulation_Development_Value = stage1$Nanoformulation_Development_Value,
  stringsAsFactors = FALSE
)

############################################################
## 11. Plot settings
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
      legend.key = element_blank(),
      plot.margin = margin(24, 36, 24, 36),
      legend.position = "right"
    )
}

pastel_palette <- c(
  "#8DD3C7",
  "#FDB462",
  "#BEBADA",
  "#FB8072",
  "#80B1D3",
  "#B3DE69",
  "#FCCDE5",
  "#BC80BD",
  "#CCEBC5",
  "#FFED6F",
  "#D9D9D9"
)

ndv_palette <- c(
  "High NDV" = "#8DD3C7",
  "Moderate NDV" = "#FDB462",
  "Low NDV" = "#BEBADA"
)

############################################################
## 12. Plot functions
############################################################

make_radar_plot <- function(one_row) {
  
  radar_modules <- c(
    "Solubility_Defect_Score",
    "Permeability_Defect_Score",
    "Barrier_Transporter_Score",
    "Distribution_Exposure_Score",
    "Metabolic_Liability_Score",
    "Nano_Assembly_Suitability_Score",
    "Druggability_Alert_Requirement",
    "Safety_Control_Requirement"
  )
  
  radar_df <- data.frame(
    compound = one_row$compound[1],
    module = radar_modules,
    score = as.numeric(one_row[1, radar_modules]),
    stringsAsFactors = FALSE
  )
  
  ## Short labels to avoid text overflow.
  radar_df$module_label <- c(
    "Solubility",
    "Permeability",
    "Barrier /\nTransporter",
    "Distribution /\nExposure",
    "Metabolism",
    "Nano-\nAssembly",
    "Drug\nAlerts",
    "Design\nControl"
  )
  
  radar_df$angle_id <- seq_len(nrow(radar_df))
  radar_df_closed <- rbind(radar_df, radar_df[1, ])
  
  ggplot(radar_df_closed, aes(x = angle_id, y = score)) +
    geom_polygon(
      fill = "#8DD3C7",
      color = "#2B8C7F",
      alpha = 0.36,
      linewidth = 1.15
    ) +
    geom_line(
      color = "#2B8C7F",
      linewidth = 1.25
    ) +
    geom_point(
      color = "#2B8C7F",
      fill = "white",
      size = 3.1,
      stroke = 1.1,
      shape = 21
    ) +
    geom_hline(
      yintercept = c(0.25, 0.50, 0.75, 1.00),
      color = "#D8D8D8",
      linewidth = 0.35
    ) +
    scale_x_continuous(
      breaks = radar_df$angle_id,
      labels = radar_df$module_label
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 1.00),
      labels = c("0", "0.25", "0.50", "0.75", "1.00")
    ) +
    coord_polar(start = -pi / 8, clip = "off") +
    labs(
      title = "Core compound diagnostic radar",
      subtitle = paste0("Compound: ", one_row$compound[1]),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.02),
      plot.subtitle = element_text(size = 13, color = "#555555", hjust = 0.02),
      axis.text.x = element_text(size = 9.5, face = "bold", color = "#333333", lineheight = 0.92),
      axis.text.y = element_text(size = 8.5, color = "#666666"),
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      plot.margin = margin(40, 80, 40, 80)
    )
}

make_driver_plot <- function(one_row) {
  
  driver_df <- data.frame(
    compound = one_row$compound[1],
    module = driver_score_cols,
    driver = unname(driver_label_map[driver_score_cols]),
    score = as.numeric(one_row[1, driver_score_cols]),
    stringsAsFactors = FALSE
  )
  
  driver_df$driver <- factor(
    driver_df$driver,
    levels = driver_df$driver[order(driver_df$score)]
  )
  
  main_driver <- wrap_label(clean_module_label(one_row$Main_Nanoformulation_Value_Driver[1]), width = 34)
  main_bottleneck <- wrap_label(clean_module_label(one_row$Main_Design_Bottleneck[1]), width = 34)
  
  subtitle_text <- paste0(
    "Main driver: ", main_driver,
    "\nMain bottleneck: ", main_bottleneck
  )
  
  ggplot(
    driver_df,
    aes(x = driver, y = score, fill = driver)
  ) +
    geom_col(
      width = 0.70,
      color = "white",
      linewidth = 0.5,
      alpha = 0.96
    ) +
    geom_text(
      aes(label = sprintf("%.2f", score)),
      hjust = -0.12,
      size = 4.0,
      color = "#333333",
      fontface = "bold"
    ) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = pastel_palette) +
    scale_y_continuous(
      limits = c(0, 1.15),
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0, 0.10))
    ) +
    labs(
      title = "Design-driver profile",
      subtitle = subtitle_text,
      x = NULL,
      y = "Module score"
    ) +
    theme_acmpi(base_size = 13) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 12.5, color = "#555555", lineheight = 1.15),
      plot.margin = margin(24, 95, 24, 24)
    )
}

############################################################
## 13. Figure 2: global value map
############################################################

p_value_map <- ggplot() +
  annotate(
    "rect",
    xmin = 0, xmax = 0.33, ymin = 0, ymax = 0.66,
    fill = "#F2F2F2",
    alpha = 0.55
  ) +
  annotate(
    "rect",
    xmin = 0.33, xmax = 1, ymin = 0.66, ymax = 1,
    fill = "#DDEFD8",
    alpha = 0.55
  ) +
  annotate(
    "rect",
    xmin = 0.33, xmax = 1, ymin = 0, ymax = 0.66,
    fill = "#FCE5CD",
    alpha = 0.45
  ) +
  annotate(
    "rect",
    xmin = 0, xmax = 0.33, ymin = 0.66, ymax = 1,
    fill = "#EADCF8",
    alpha = 0.45
  ) +
  geom_hline(
    yintercept = 0.66,
    linetype = "dashed",
    color = "#9E9E9E",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0.33,
    linetype = "dashed",
    color = "#9E9E9E",
    linewidth = 0.5
  ) +
  geom_point(
    data = value_map_df,
    aes(
      x = Nano_Optimization_Need,
      y = Nano_Assembly_Suitability_Score,
      size = Safety_Control_Requirement,
      fill = Nanoformulation_Development_Value
    ),
    shape = 21,
    color = "#444444",
    stroke = 0.8,
    alpha = 0.90
  ) +
  geom_text(
    data = value_map_df,
    aes(
      x = Nano_Optimization_Need,
      y = Nano_Assembly_Suitability_Score,
      label = compound
    ),
    nudge_y = 0.055,
    size = 4.0,
    fontface = "bold",
    color = "#333333",
    check_overlap = TRUE
  ) +
  annotate(
    "text",
    x = 0.67, y = 0.95,
    label = "High nanoformulation\nvalue zone",
    size = 3.9,
    fontface = "bold",
    color = "#406B45"
  ) +
  annotate(
    "text",
    x = 0.16, y = 0.91,
    label = "Assembly-rich /\nlow-need zone",
    size = 3.6,
    color = "#6B4F8A"
  ) +
  annotate(
    "text",
    x = 0.68, y = 0.12,
    label = "Need-driven /\ncarrier-assisted zone",
    size = 3.6,
    color = "#8A5A2B"
  ) +
  scale_fill_manual(values = ndv_palette, drop = FALSE) +
  scale_size_continuous(range = c(5, 12), limits = c(0, 1)) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.04, 0.12))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.04, 0.12))
  ) +
  labs(
    title = "Nanoformulation development value map",
    subtitle = "Positioning by nano-optimization need and nano-assembly suitability",
    x = "Nano-optimization need",
    y = "Nano-assembly suitability",
    fill = "Development value",
    size = "Design-control\nrequirement"
  ) +
  theme_acmpi(base_size = 13) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12.5, color = "#555555"),
    plot.margin = margin(24, 70, 24, 28)
  ) +
  coord_cartesian(clip = "off")

############################################################
## 14. Save tables
############################################################

out_csv_summary <- file.path(
  table_dir,
  paste0(output_prefix, "_Summary.csv")
)

out_csv_modules <- file.path(
  table_dir,
  paste0(output_prefix, "_ModuleScores.csv")
)

out_csv_full <- file.path(
  table_dir,
  paste0(output_prefix, "_FullScoredData.csv")
)

out_xlsx <- file.path(
  table_dir,
  paste0(output_prefix, "_Results.xlsx")
)

write.csv(
  stage1_summary,
  out_csv_summary,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  module_table,
  out_csv_modules,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  stage1,
  out_csv_full,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Stage1_Summary")
openxlsx::writeData(wb, "Stage1_Summary", stage1_summary)

openxlsx::addWorksheet(wb, "ModuleScores")
openxlsx::writeData(wb, "ModuleScores", module_table)

openxlsx::addWorksheet(wb, "FullScoredData")
openxlsx::writeData(wb, "FullScoredData", stage1)

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

############################################################
## 15. Save figures
############################################################

## Global value map
fig2_png <- file.path(fig_png_dir, "Fig2_Nanoformulation_Value_Map.png")
fig2_pdf <- file.path(fig_pdf_dir, "Fig2_Nanoformulation_Value_Map.pdf")

ggsave(
  filename = fig2_png,
  plot = p_value_map,
  width = 12.4,
  height = 8.0,
  dpi = FIG_DPI,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = fig2_pdf,
  plot = p_value_map,
  width = 12.4,
  height = 8.0,
  device = cairo_pdf,
  bg = "white",
  limitsize = FALSE
)

## Per-compound radar and design-driver plots
for (i in seq_len(nrow(stage1))) {
  
  one_row <- stage1[i, , drop = FALSE]
  compound_id <- safe_name(one_row$compound[1])
  
  p_radar <- make_radar_plot(one_row)
  p_driver <- make_driver_plot(one_row)
  
  fig1_png <- file.path(
    fig_png_dir,
    paste0("Fig1_CoreCompound_Diagnostic_Radar_", compound_id, ".png")
  )
  
  fig1_pdf <- file.path(
    fig_pdf_dir,
    paste0("Fig1_CoreCompound_Diagnostic_Radar_", compound_id, ".pdf")
  )
  
  fig3_png <- file.path(
    fig_png_dir,
    paste0("Fig3_Design_Driver_Profile_", compound_id, ".png")
  )
  
  fig3_pdf <- file.path(
    fig_pdf_dir,
    paste0("Fig3_Design_Driver_Profile_", compound_id, ".pdf")
  )
  
  ggsave(
    filename = fig1_png,
    plot = p_radar,
    width = 10.2,
    height = 10.2,
    dpi = FIG_DPI,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = fig1_pdf,
    plot = p_radar,
    width = 10.2,
    height = 10.2,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = fig3_png,
    plot = p_driver,
    width = 12.6,
    height = 7.8,
    dpi = FIG_DPI,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = fig3_pdf,
    plot = p_driver,
    width = 12.6,
    height = 7.8,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
}

############################################################
## 16. Print result summary
############################################################

cat("\n============================================================\n")
cat("ACMPI-Nano Stage 1: Core Compound Nanoformulation Value Model\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat(input_file, "\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Compact diagnostic summary:\n")
print(stage1_summary, row.names = FALSE)

cat("\nModule score table:\n")
print(module_table[, c("compound", "module", "score")], row.names = FALSE)

cat("\nSaved tables:\n")
cat(out_csv_summary, "\n")
cat(out_csv_modules, "\n")
cat(out_csv_full, "\n")
cat(out_xlsx, "\n")

cat("\nSaved figures in:\n")
cat(fig_png_dir, "\n")
cat(fig_pdf_dir, "\n")

cat("\nInterpretation guide:\n")
cat("- Nanoformulation_Value_Index indicates overall nanoformulation development value.\n")
cat("- Nanoformulation_Development_Value is classified as High, Moderate, or Low NDV.\n")
cat("- Safety_Control_Requirement is a design-control need, not an exclusion penalty.\n")
cat("- Main_Design_Bottleneck identifies the strongest design challenge.\n")
cat("- Main_Nanoformulation_Value_Driver identifies the dominant contributor to nanoformulation value.\n")
cat("- Main_Structural_Advantage identifies the strongest molecular feature supporting nanoformulation design.\n")