############################################################
## Build master model input table for Baicalein
## Fixed version:
##   - No PubChem file required
##   - No ProTox oral toxicity txt required
##   - InChIKey, LD50, Toxicity_Class are filled from verified values
############################################################

target_dir <- "/media/desk16/iy15915/中药之开创/Baicalein"

############################################################
## Verified manually curated identity / acute toxicity values
############################################################

verified_inchikey <- "FXNFHKRTJBSTCS-UHFFFAOYSA-N"
verified_ld50 <- 3919
verified_toxicity_class <- 5

############################################################
## Packages
############################################################

pkgs <- c("openxlsx")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

############################################################
## Input files
############################################################

file_rdkit  <- file.path(target_dir, "RDKit.csv")
file_swiss  <- file.path(target_dir, "SwissADME.csv")
file_admet  <- file.path(target_dir, "ADMETlab 3.0.csv")
file_protox <- file.path(target_dir, "ProTox 3.0.csv")

required_files <- c(file_rdkit, file_swiss, file_admet, file_protox)

for (x in required_files) {
  if (!file.exists(x)) {
    stop("Missing required file: ", x)
  }
}

############################################################
## Helper functions
############################################################

read_csv_safe <- function(path) {
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8"),
    error = function(e1) {
      tryCatch(
        read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"),
        error = function(e2) {
          read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
        }
      )
    }
  )
}

get1 <- function(dat, col, default = NA) {
  if (!col %in% colnames(dat)) return(default)
  val <- dat[[col]][1]
  if (length(val) == 0 || is.na(val) || identical(val, "")) return(default)
  val
}

yes_to_one <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(x %in% c("yes", "active", "true", "1"), 1, 0)
}

active_count <- function(dat, class_exact) {
  if (!all(c("Classification", "Prediction") %in% colnames(dat))) {
    return(NA_integer_)
  }
  
  idx <- trimws(dat$Classification) == class_exact
  
  sum(
    tolower(trimws(dat$Prediction[idx])) == "active",
    na.rm = TRUE
  )
}

############################################################
## Read data
############################################################

rdkit  <- read_csv_safe(file_rdkit)
swiss  <- read_csv_safe(file_swiss)
admet  <- read_csv_safe(file_admet)
protox <- read_csv_safe(file_protox)

############################################################
## Derived variables
############################################################

## SwissADME CYP inhibitor count
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

## ProTox active counts
protox_organ_active_count <- active_count(protox, "Organ toxicity")
protox_endpoint_active_count <- active_count(protox, "Toxicity end points")
protox_metabolism_active_count <- active_count(protox, "Metabolism")

############################################################
## Build master model input table
############################################################

master <- data.frame(
  ## Identity fields
  compound = get1(rdkit, "compound", "Baicalein"),
  PubChem_CID = get1(rdkit, "PubChem_CID"),
  SMILES = get1(rdkit, "SMILES"),
  InChIKey = verified_inchikey,
  
  ## Molecular structure
  Molecular_Weight = get1(rdkit, "Molecular_Weight"),
  MolLogP = get1(rdkit, "MolLogP"),
  TPSA = get1(rdkit, "TPSA"),
  HBD = get1(rdkit, "HBD"),
  HBA = get1(rdkit, "HBA"),
  Rotatable_Bonds = get1(rdkit, "Rotatable_Bonds"),
  Ring_Count = get1(rdkit, "Ring_Count"),
  Aromatic_Rings = get1(rdkit, "Aromatic_Rings"),
  Fraction_Csp3 = get1(rdkit, "Fraction_Csp3"),
  Formal_Charge = get1(rdkit, "Formal_Charge"),
  Molar_Refractivity = get1(rdkit, "Molar_Refractivity"),
  
  ## Solubility / lipophilicity
  Consensus_LogP = get1(swiss, "Consensus Log P"),
  ESOL_LogS = get1(swiss, "ESOL Log S"),
  ESOL_Class = get1(swiss, "ESOL Class"),
  ADMETlab_logS = get1(admet, "logS"),
  ADMETlab_logD = get1(admet, "logD"),
  
  ## Absorption / permeability
  GI_absorption = get1(swiss, "GI absorption"),
  Caco2 = get1(admet, "caco2"),
  PAMPA = get1(admet, "PAMPA"),
  F30 = get1(admet, "f30"),
  
  ## Barrier / transporter
  BBB_permeant = get1(swiss, "BBB permeant"),
  ADMETlab_BBB = get1(admet, "BBB"),
  Pgp_substrate = get1(swiss, "Pgp substrate"),
  pgp_sub = get1(admet, "pgp_sub"),
  pgp_inh = get1(admet, "pgp_inh"),
  
  ## Distribution / exposure
  PPB = get1(admet, "PPB"),
  Fu = get1(admet, "Fu"),
  logVDss = get1(admet, "logVDss"),
  
  ## Metabolism / clearance
  CYP_inhibitor_count = cyp_inhibitor_count,
  t_half = get1(admet, "t0.5"),
  cl_plasma = get1(admet, "cl-plasma"),
  
  ## Druggability / alerts
  Lipinski_violations = get1(swiss, "Lipinski #violations"),
  Veber_violations = get1(swiss, "Veber #violations"),
  Bioavailability_Score = get1(swiss, "Bioavailability Score"),
  PAINS_alerts = get1(swiss, "PAINS #alerts"),
  Brenk_alerts = get1(swiss, "Brenk #alerts"),
  Synthetic_Accessibility = get1(swiss, "Synthetic Accessibility"),
  QED = get1(admet, "QED"),
  
  ## Toxicity / safety
  LD50 = verified_ld50,
  Toxicity_Class = verified_toxicity_class,
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

############################################################
## Official column order
############################################################

official_cols <- c(
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
  
  "LD50", "Toxicity_Class", "DILI", "Ames", "Carcinogenicity", "hERG",
  "ProTox_Organ_Active_Count",
  "ProTox_Endpoint_Active_Count",
  "ProTox_Metabolism_Active_Count"
)

master <- master[, official_cols]

############################################################
## Missing value check
############################################################

missing_check <- data.frame(
  field = colnames(master),
  value = as.character(master[1, ]),
  is_missing = is.na(as.character(master[1, ])) | as.character(master[1, ]) == "",
  stringsAsFactors = FALSE
)

############################################################
## Print result
############################################################

cat("\nMaster model input table:\n")
print(master, row.names = FALSE)

cat("\nNumber of fields:", ncol(master), "\n")
cat("Identity fields: 4\n")
cat("Model scoring fields:", ncol(master) - 4, "\n")

cat("\nMissing fields:\n")
missing_rows <- missing_check[missing_check$is_missing, ]

if (nrow(missing_rows) == 0) {
  cat("No missing fields.\n")
} else {
  print(missing_rows, row.names = FALSE)
}

############################################################
## Save outputs
############################################################

out_csv <- file.path(target_dir, "Baicalein_Master_Model_Input.csv")
out_xlsx <- file.path(target_dir, "Baicalein_Master_Model_Input.xlsx")

write.csv(
  master,
  out_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

openxlsx::write.xlsx(
  master,
  out_xlsx,
  overwrite = TRUE
)

cat("\nSaved files:\n")
cat(out_csv, "\n")
cat(out_xlsx, "\n")