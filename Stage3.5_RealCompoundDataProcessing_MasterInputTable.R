############################################################
## Build master model input tables for 20 compounds
## FULL FIXED VERSION
##
## Purpose:
##   Merge RDKit, SwissADME, ADMETlab and ProTox outputs
##   for 20 real compounds.
##
## Key fixes:
##   - Root path updated to /media/desk16/iy15915/中药之开创/20药物
##   - Full 20-compound LD50 and Toxicity_Class table included
##   - Compatible with old and new RDKit column names
##   - One master file per compound
##   - One combined MASTER_20Drugs_Model_Input.csv/xlsx
##   - Strict missing-field QC
############################################################

rm(list = ls())

############################################################
## 0. Root directory
############################################################

root_dir <- "/media/desk16/iy15915/中药之开创/20药物"

if (!dir.exists(root_dir)) {
  stop("Root directory does not exist: ", root_dir)
}

############################################################
## 1. Full 20-compound toxicity parameters
## Units:
##   LD50 = mg/kg
##   Toxicity_Class = ProTox-style predicted toxicity class
############################################################

tox_manual <- data.frame(
  compound = c(
    "Apigenin",
    "Berberine",
    "Caffeic Acid",
    "Calycosin",
    "Chlorogenic Acid",
    "Cryptotanshinone",
    "Curcumin",
    "Emodin",
    "Ferulic Acid",
    "Formononetin",
    "Ginsenoside Rg3",
    "Glycyrrhetinic Acid",
    "Glycyrrhizin",
    "Kaempferol",
    "Liquiritigenin",
    "Luteolin",
    "Naringenin",
    "Quercetin",
    "Saikosaponin A",
    "Tanshinone IIA"
  ),
  LD50 = c(
    2500,
    200,
    2980,
    2500,
    5000,
    8000,
    2000,
    5000,
    1772,
    2500,
    4000,
    560,
    1750,
    3919,
    2000,
    3919,
    2000,
    159,
    2000,
    1230
  ),
  Toxicity_Class = c(
    5,
    3,
    5,
    5,
    5,
    6,
    4,
    5,
    4,
    5,
    5,
    4,
    4,
    5,
    4,
    5,
    4,
    3,
    4,
    4
  ),
  stringsAsFactors = FALSE
)

############################################################
## 2. Packages
############################################################

pkgs <- c("openxlsx")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
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

get1 <- function(dat, col, default = NA) {
  if (is.null(dat)) return(default)
  if (!col %in% colnames(dat)) return(default)
  val <- dat[[col]][1]
  if (length(val) == 0 || is.na(val) || identical(val, "")) return(default)
  val
}

get1_any <- function(dat, cols, default = NA) {
  if (is.null(dat)) return(default)
  for (cc in cols) {
    if (cc %in% colnames(dat)) {
      val <- dat[[cc]][1]
      if (length(val) > 0 && !is.na(val) && !identical(val, "")) {
        return(val)
      }
    }
  }
  return(default)
}

yes_to_one <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(x %in% c("yes", "active", "true", "1"), 1, 0)
}

active_count <- function(dat, class_exact) {
  if (is.null(dat)) return(NA_integer_)
  if (!all(c("Classification", "Prediction") %in% colnames(dat))) {
    return(NA_integer_)
  }
  idx <- trimws(dat$Classification) == class_exact
  sum(tolower(trimws(dat$Prediction[idx])) == "active", na.rm = TRUE)
}

clean_compound_name <- function(x) {
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x
}

safe_file_base <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

find_file <- function(folder, pattern_keywords) {
  files <- list.files(folder, recursive = FALSE, full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  files <- files[file.exists(files)]
  files <- files[!dir.exists(files)]
  if (length(files) == 0) return(NA_character_)
  file_names_low <- tolower(basename(files))
  idx <- rep(TRUE, length(files))
  for (kw in pattern_keywords) {
    idx <- idx & grepl(tolower(kw), file_names_low, fixed = TRUE)
  }
  matched <- files[idx]
  if (length(matched) == 0) return(NA_character_)
  matched[1]
}

############################################################
## 4. Official column order
############################################################

official_cols <- c(
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
  "Lipinski_violations", "Veber_violations", "Bioavailability_Score",
  "PAINS_alerts", "Brenk_alerts", "Synthetic_Accessibility", "QED",
  "LD50", "Toxicity_Class",
  "DILI", "Ames", "Carcinogenicity", "hERG",
  "ProTox_Organ_Active_Count",
  "ProTox_Endpoint_Active_Count",
  "ProTox_Metabolism_Active_Count"
)

############################################################
## 5. Process one compound folder
############################################################

process_one_compound <- function(compound_dir) {
  folder_name <- basename(compound_dir)
  compound_name <- clean_compound_name(folder_name)
  file_base <- safe_file_base(compound_name)
  
  cat("\n========================================\n")
  cat("Processing compound:", compound_name, "\n")
  cat("Folder:", compound_dir, "\n")
  
  file_rdkit  <- find_file(compound_dir, c("rdkit"))
  file_swiss  <- find_file(compound_dir, c("swiss"))
  file_admet  <- find_file(compound_dir, c("admetlab"))
  file_protox <- find_file(compound_dir, c("protox"))
  
  required_files <- c(file_rdkit, file_swiss, file_admet, file_protox)
  
  if (any(is.na(required_files))) {
    warning("Missing one or more required files for: ", compound_name)
    return(list(
      master = NULL,
      log = data.frame(
        compound = compound_name,
        Status = "Missing_required_file",
        RDKit = file_rdkit,
        SwissADME = file_swiss,
        ADMETlab = file_admet,
        ProTox = file_protox,
        Output_CSV = NA,
        Output_XLSX = NA,
        Missing_Field_Count = NA,
        stringsAsFactors = FALSE
      ),
      missing = data.frame(
        compound = compound_name,
        field = "Required_file",
        value = NA,
        is_missing = TRUE,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  rdkit  <- read_csv_safe(file_rdkit)
  swiss  <- read_csv_safe(file_swiss)
  admet  <- read_csv_safe(file_admet)
  protox <- read_csv_safe(file_protox)
  
  swiss_cyp_cols <- c(
    "CYP1A2 inhibitor",
    "CYP2C19 inhibitor",
    "CYP2C9 inhibitor",
    "CYP2D6 inhibitor",
    "CYP3A4 inhibitor"
  )
  
  cyp_inhibitor_count <- sum(
    sapply(swiss_cyp_cols, function(cc) yes_to_one(get1(swiss, cc, "No"))),
    na.rm = TRUE
  )
  
  protox_organ_active_count <- active_count(protox, "Organ toxicity")
  protox_endpoint_active_count <- active_count(protox, "Toxicity end points")
  protox_metabolism_active_count <- active_count(protox, "Metabolism")
  
  tox_row <- tox_manual[tox_manual$compound == compound_name, , drop = FALSE]
  if (nrow(tox_row) == 0) {
    manual_ld50 <- NA
    manual_toxicity_class <- NA
  } else {
    manual_ld50 <- tox_row$LD50[1]
    manual_toxicity_class <- tox_row$Toxicity_Class[1]
  }
  
  master <- data.frame(
    compound = compound_name,
    PubChem_CID = get1_any(rdkit, c("PubChem_CID", "PubChem CID", "CID")),
    SMILES = get1_any(rdkit, c("SMILES", "Input_SMILES", "Canonical_SMILES", "Canonical SMILES", "smiles"), default = NA),
    Molecular_Weight = get1_any(rdkit, c("Molecular_Weight", "Molecular Weight", "MW")),
    MolLogP = get1_any(rdkit, c("MolLogP", "LogP", "Mol LogP")),
    TPSA = get1_any(rdkit, c("TPSA", "Topological_Polar_Surface_Area")),
    HBD = get1_any(rdkit, c("HBD", "NumHDonors")),
    HBA = get1_any(rdkit, c("HBA", "NumHAcceptors")),
    Rotatable_Bonds = get1_any(rdkit, c("Rotatable_Bonds", "Rotatable Bonds", "NumRotatableBonds")),
    Ring_Count = get1_any(rdkit, c("Ring_Count", "Ring Count")),
    Aromatic_Rings = get1_any(rdkit, c("Aromatic_Rings", "Aromatic Rings")),
    Fraction_Csp3 = get1_any(rdkit, c("Fraction_Csp3", "FractionCSP3")),
    Formal_Charge = get1_any(rdkit, c("Formal_Charge", "Formal Charge")),
    Molar_Refractivity = get1_any(rdkit, c("Molar_Refractivity", "Molar Refractivity", "MolMR")),
    Consensus_LogP = get1(swiss, "Consensus Log P"),
    ESOL_LogS = get1(swiss, "ESOL Log S"),
    ESOL_Class = get1(swiss, "ESOL Class"),
    ADMETlab_logS = get1(admet, "logS"),
    ADMETlab_logD = get1(admet, "logD"),
    GI_absorption = get1(swiss, "GI absorption"),
    Caco2 = get1(admet, "caco2"),
    PAMPA = get1(admet, "PAMPA"),
    F30 = get1(admet, "f30"),
    BBB_permeant = get1(swiss, "BBB permeant"),
    ADMETlab_BBB = get1(admet, "BBB"),
    Pgp_substrate = get1(swiss, "Pgp substrate"),
    pgp_sub = get1(admet, "pgp_sub"),
    pgp_inh = get1(admet, "pgp_inh"),
    PPB = get1(admet, "PPB"),
    Fu = get1(admet, "Fu"),
    logVDss = get1(admet, "logVDss"),
    CYP_inhibitor_count = cyp_inhibitor_count,
    t_half = get1_any(admet, c("t0.5", "t_half", "half_life")),
    cl_plasma = get1_any(admet, c("cl-plasma", "cl_plasma", "CLplasma")),
    Lipinski_violations = get1_any(swiss, c("Lipinski #violations", "Lipinski_Violations", "Lipinski violations")),
    Veber_violations = get1_any(swiss, c("Veber #violations", "Veber_Violations", "Veber violations")),
    Bioavailability_Score = get1(swiss, "Bioavailability Score"),
    PAINS_alerts = get1(swiss, "PAINS #alerts"),
    Brenk_alerts = get1(swiss, "Brenk #alerts"),
    Synthetic_Accessibility = get1(swiss, "Synthetic Accessibility"),
    QED = get1(admet, "QED"),
    LD50 = manual_ld50,
    Toxicity_Class = manual_toxicity_class,
    DILI = get1(admet, "DILI"),
    Ames = get1(admet, "Ames"),
    Carcinogenicity = get1(admet, "Carcinogenicity"),
    hERG = get1(admet, "hERG"),
    ProTox_Organ_Active_Count = protox_organ_active_count,
    ProTox_Endpoint_Active_Count = protox_endpoint_active_count,
    ProTox_Metabolism_Active_Count = protox_metabolism_active_count,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  master <- master[, official_cols]
  
  missing_check <- data.frame(
    compound = compound_name,
    field = colnames(master),
    value = as.character(master[1, ]),
    is_missing = is.na(as.character(master[1, ])) | as.character(master[1, ]) == "",
    stringsAsFactors = FALSE
  )
  
  missing_rows <- missing_check[missing_check$is_missing, ]
  
  out_csv <- file.path(compound_dir, paste0(file_base, "_Master_Model_Input.csv"))
  out_xlsx <- file.path(compound_dir, paste0(file_base, "_Master_Model_Input.xlsx"))
  out_missing <- file.path(compound_dir, paste0(file_base, "_Missing_Check.csv"))
  
  write.csv(master, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
  openxlsx::write.xlsx(master, out_xlsx, overwrite = TRUE)
  write.csv(missing_check, out_missing, row.names = FALSE, fileEncoding = "UTF-8")
  
  cat("[OK] Saved:", out_csv, "\n")
  cat("[OK] Saved:", out_xlsx, "\n")
  
  if (nrow(missing_rows) == 0) {
    cat("[OK] No missing fields.\n")
  } else {
    cat("[WARNING] Missing fields:", nrow(missing_rows), "\n")
  }
  
  log <- data.frame(
    compound = compound_name,
    Status = "OK",
    RDKit = file_rdkit,
    SwissADME = file_swiss,
    ADMETlab = file_admet,
    ProTox = file_protox,
    Output_CSV = out_csv,
    Output_XLSX = out_xlsx,
    Missing_Field_Count = nrow(missing_rows),
    stringsAsFactors = FALSE
  )
  
  return(list(master = master, log = log, missing = missing_check))
}

############################################################
## 6. Batch run
############################################################

compound_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
compound_dirs <- compound_dirs[basename(compound_dirs) %in% tox_manual$compound]

if (length(compound_dirs) == 0) {
  stop("No target compound folders found under: ", root_dir)
}

compound_dirs <- compound_dirs[order(match(basename(compound_dirs), tox_manual$compound))]

all_master <- list()
all_log <- list()
all_missing <- list()

for (d in compound_dirs) {
  res <- process_one_compound(d)
  all_log[[basename(d)]] <- res$log
  if (!is.null(res$master)) all_master[[basename(d)]] <- res$master
  if (!is.null(res$missing)) all_missing[[basename(d)]] <- res$missing
}

############################################################
## 7. Combine outputs
############################################################

if (length(all_master) > 0) {
  master_all <- do.call(rbind, all_master)
  rownames(master_all) <- NULL
} else {
  master_all <- data.frame()
}

log_all <- do.call(rbind, all_log)
missing_all <- do.call(rbind, all_missing)
rownames(log_all) <- NULL
rownames(missing_all) <- NULL

############################################################
## 8. Final QC
############################################################

expected_compounds <- tox_manual$compound
merged_compounds <- master_all$compound
missing_compounds <- setdiff(expected_compounds, merged_compounds)
extra_compounds <- setdiff(merged_compounds, expected_compounds)

if (length(missing_compounds) > 0) {
  warning("Some expected compounds were not merged:\n", paste(missing_compounds, collapse = "\n"))
}

if (length(extra_compounds) > 0) {
  warning("Unexpected compounds detected in merged table:\n", paste(extra_compounds, collapse = "\n"))
}

critical_fields <- official_cols
critical_missing_long <- data.frame()

if (nrow(master_all) > 0) {
  for (cc in critical_fields) {
    idx <- is.na(master_all[[cc]]) | as.character(master_all[[cc]]) == ""
    if (any(idx)) {
      critical_missing_long <- rbind(
        critical_missing_long,
        data.frame(
          compound = master_all$compound[idx],
          field = cc,
          value = as.character(master_all[[cc]][idx]),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

############################################################
## 9. Save combined outputs
############################################################

out_master_csv <- file.path(root_dir, "MASTER_20Drugs_Model_Input.csv")
out_master_xlsx <- file.path(root_dir, "MASTER_20Drugs_Model_Input.xlsx")
out_log_csv <- file.path(root_dir, "MASTER_20Drugs_Run_Log.csv")
out_missing_csv <- file.path(root_dir, "MASTER_20Drugs_Missing_Check.csv")
out_critical_missing_csv <- file.path(root_dir, "MASTER_20Drugs_Critical_Missing_Check.csv")

if (nrow(master_all) > 0) {
  write.csv(master_all, out_master_csv, row.names = FALSE, fileEncoding = "UTF-8")
  openxlsx::write.xlsx(master_all, out_master_xlsx, overwrite = TRUE)
} else {
  warning("No compound was successfully merged. MASTER table was not created.")
}

write.csv(log_all, out_log_csv, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(missing_all, out_missing_csv, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(critical_missing_long, out_critical_missing_csv, row.names = FALSE, fileEncoding = "UTF-8")

############################################################
## 10. Print final summary
############################################################

cat("\n============================================================\n")
cat("Batch processing completed.\n")
cat("Root directory:", root_dir, "\n")
cat("Target compound folders:", length(compound_dirs), "\n")
cat("Successfully merged compounds:", nrow(master_all), "\n")

cat("\nExpected compounds:\n")
print(expected_compounds)

cat("\nMerged compounds:\n")
print(master_all$compound)

cat("\nSaved combined files:\n")
if (nrow(master_all) > 0) {
  cat(out_master_csv, "\n")
  cat(out_master_xlsx, "\n")
}
cat(out_log_csv, "\n")
cat(out_missing_csv, "\n")
cat(out_critical_missing_csv, "\n")

if (nrow(master_all) > 0) {
  cat("\nCompound-level LD50 and toxicity class:\n")
  print(master_all[, c("compound", "LD50", "Toxicity_Class")], row.names = FALSE)
}

cat("\nRun log:\n")
print(log_all[, c("compound", "Status", "Missing_Field_Count")], row.names = FALSE)

cat("\nMissing-required-file compounds, if any:\n")
print(
  log_all[
    log_all$Status != "OK",
    c("compound", "Status", "RDKit", "SwissADME", "ADMETlab", "ProTox")
  ],
  row.names = FALSE
)

if (nrow(critical_missing_long) == 0) {
  cat("\n[QC OK] No critical missing fields in MASTER_20Drugs_Model_Input.\n")
} else {
  cat("\n[QC WARNING] Critical missing fields remain. See:\n")
  cat(out_critical_missing_csv, "\n")
  print(critical_missing_long)
}

cat("\nAll done.\n")
cat("============================================================\n")
