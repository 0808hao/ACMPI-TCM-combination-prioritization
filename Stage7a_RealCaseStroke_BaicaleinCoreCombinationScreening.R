############################################################
## Real six-drug ACMPI-NPS combination analysis
## Core drug: Baicalein
## FINAL no-overflow version
############################################################

rm(list = ls())

############################################################
## 0. Settings
############################################################
ROOT_DIR <- "/media/desk16/iy15915/中药之开创/脑卒中疾病基因"
DRUG_ROOT <- "/media/desk16/iy15915/中药之开创/真实药物验证"

STRING_SCORE_MIN <- 700
CORE_DRUG <- "Baicalein"
SET_SEED <- 20260504
set.seed(SET_SEED)

B_PERM_COMBINATION <- 500
DEGREE_BINS <- 10
DRIVER_EPS <- 0.005
NEAR_DISTANCE_CUTOFF <- 1
FAR_DISTANCE_CUTOFF <- 3

DRUGS_KEEP <- c(
  "Baicalein",
  "Caffeic Acid",
  "Curcumin",
  "Ferulic Acid",
  "Glycyrrhetinic Acid",
  "Liquiritigenin"
)

############################################################
## 1. Packages
############################################################
pkgs <- c("data.table", "ggplot2", "scales")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(data.table)
library(ggplot2)
library(scales)

############################################################
## 2. Locate latest stroke disease module
############################################################
find_latest_disease_dir <- function(root_dir, score_min = 700) {
  dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  pattern <- paste0("ANPS_StrokeDisease_STRING", score_min, "_")
  dirs <- dirs[grepl(pattern, basename(dirs), fixed = TRUE)]
  
  req <- c(
    "01_StrokeDisease_GeneMapping_STRING.tsv",
    "02_STRING700_NodeDistanceToStrokeDisease.tsv",
    "04_StrokeDisease_PPI_QC.tsv"
  )
  
  valid <- dirs[vapply(dirs, function(d) {
    all(file.exists(file.path(d, req)))
  }, logical(1))]
  
  if (length(valid) == 0) {
    stop("No valid stroke disease ANPS folder found.", call. = FALSE)
  }
  
  valid[which.max(file.info(valid)$mtime)]
}

DISEASE_DIR <- find_latest_disease_dir(ROOT_DIR, STRING_SCORE_MIN)

NODE_DISTANCE_FILE <- file.path(
  DISEASE_DIR,
  "02_STRING700_NodeDistanceToStrokeDisease.tsv"
)

OUT_DIR <- file.path(
  DRUG_ROOT,
  paste0(
    "RealDrugCombination_BaicaleinCore_Final_NoOverflow_STRING",
    STRING_SCORE_MIN,
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

FIG_DIR <- file.path(OUT_DIR, "Figures_NoOverflow")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\nDisease module folder:\n", DISEASE_DIR, "\n", sep = "")
cat("Output folder:\n", OUT_DIR, "\n", sep = "")

############################################################
## 3. Read stroke node-distance table
############################################################
node_dt <- fread(NODE_DISTANCE_FILE, sep = "\t", data.table = TRUE)

required_cols <- c(
  "string_id",
  "preferred_name",
  "distance_to_stroke_disease",
  "degree_STRING700",
  "is_stroke_disease_node"
)

miss <- setdiff(required_cols, colnames(node_dt))
if (length(miss) > 0) {
  stop("Missing columns in node distance file:\n", paste(miss, collapse = "\n"))
}

node_dt[, string_id := as.character(string_id)]
node_dt[, preferred_name := as.character(preferred_name)]
node_dt[, distance_to_stroke_disease := as.numeric(distance_to_stroke_disease)]
node_dt[, degree_STRING700 := as.numeric(degree_STRING700)]

if (!is.logical(node_dt$is_stroke_disease_node)) {
  node_dt[, is_stroke_disease_node := as.logical(is_stroke_disease_node)]
}

node_dt <- node_dt[
  !is.na(string_id) &
    !is.na(distance_to_stroke_disease) &
    !is.na(degree_STRING700)
]

node_dt[, degree_rank := frank(degree_STRING700, ties.method = "average")]
node_dt[, degree_bin := ceiling(degree_rank / .N * DEGREE_BINS)]
node_dt[degree_bin < 1, degree_bin := 1]
node_dt[degree_bin > DEGREE_BINS, degree_bin := DEGREE_BINS]

background_pool <- node_dt[
  is_stroke_disease_node == FALSE &
    !is.na(distance_to_stroke_disease) &
    !is.na(degree_STRING700)
]

all_background_ids <- background_pool$string_id
degree_bin_pools <- split(background_pool$string_id, background_pool$degree_bin)

dist_lookup <- node_dt$distance_to_stroke_disease
names(dist_lookup) <- node_dt$string_id

degree_lookup <- node_dt$degree_STRING700
names(degree_lookup) <- node_dt$string_id

bin_lookup <- node_dt$degree_bin
names(bin_lookup) <- node_dt$string_id

gene_lookup <- node_dt$preferred_name
names(gene_lookup) <- node_dt$string_id

disease_flag_lookup <- node_dt$is_stroke_disease_node
names(disease_flag_lookup) <- node_dt$string_id

############################################################
## 4. Read six real drug STRING mapping files
############################################################
drug_dirs <- list.dirs(DRUG_ROOT, recursive = FALSE, full.names = TRUE)
drug_dirs <- drug_dirs[basename(drug_dirs) %in% DRUGS_KEEP]

find_latest_drug_anps <- function(drug_dir) {
  dirs <- list.dirs(drug_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl("ANPS_.*_STRING700_", basename(dirs))]
  dirs <- dirs[file.exists(file.path(dirs, "01_GeneMapping_STRING.tsv"))]
  if (length(dirs) == 0) stop("No ANPS folder found in: ", drug_dir)
  dirs[which.max(file.info(dirs)$mtime)]
}

drug_map_list <- list()

for (d in drug_dirs) {
  drug_name <- basename(d)
  anps_dir <- find_latest_drug_anps(d)
  map_file <- file.path(anps_dir, "01_GeneMapping_STRING.tsv")
  
  x <- fread(map_file, sep = "\t", data.table = TRUE)
  
  x <- x[
    Mapped_to_STRING == TRUE &
      Present_in_STRING700_PPI == TRUE &
      !is.na(string_id)
  ]
  
  drug_map_list[[drug_name]] <- data.table(
    Drug = drug_name,
    Input_Gene = x$Input_Gene,
    string_id = as.character(x$string_id),
    preferred_name = as.character(x$preferred_name)
  )
  
  cat(drug_name, "targets in STRING700:", length(unique(x$string_id)), "\n")
}

drug_targets <- unique(rbindlist(drug_map_list, use.names = TRUE, fill = TRUE))
drug_names <- DRUGS_KEEP[DRUGS_KEEP %in% unique(drug_targets$Drug)]

if (!CORE_DRUG %in% drug_names) {
  stop("Core drug not found: ", CORE_DRUG)
}

partner_drugs <- setdiff(drug_names, CORE_DRUG)

############################################################
## 5. Helper functions
############################################################
mean_distance <- function(ids) {
  ids <- unique(ids)
  d <- as.numeric(dist_lookup[ids])
  d <- d[!is.na(d)]
  if (length(d) == 0) return(NA_real_)
  mean(d)
}

median_distance <- function(ids) {
  ids <- unique(ids)
  d <- as.numeric(dist_lookup[ids])
  d <- d[!is.na(d)]
  if (length(d) == 0) return(NA_real_)
  median(d)
}

get_drug_ids <- function(drugs) {
  drug_targets[Drug %in% drugs, unique(string_id)]
}

sample_safe <- function(x, n, replace = FALSE) {
  x <- unique(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) stop("Empty sampling pool.")
  if (length(x) < n && !replace) sample(x, n, replace = TRUE) else sample(x, n, replace = replace)
}

complete_to_n <- function(ids, n, fill_pool) {
  ids <- unique(ids)
  ids <- ids[!is.na(ids)]
  if (length(ids) >= n) return(ids[seq_len(n)])
  
  pool <- setdiff(unique(fill_pool), ids)
  extra <- sample_safe(pool, n - length(ids), replace = length(pool) < (n - length(ids)))
  out <- unique(c(ids, extra))
  
  if (length(out) < n) {
    out <- c(out, sample_safe(out, n - length(out), replace = TRUE))
  }
  
  out[seq_len(n)]
}

sample_degree_matched <- function(reference_ids, n_out) {
  reference_ids <- unique(reference_ids)
  reference_ids <- reference_ids[reference_ids %in% node_dt$string_id]
  
  ref_bins <- bin_lookup[reference_ids]
  ref_bins <- ref_bins[!is.na(ref_bins)]
  
  if (length(ref_bins) == 0) {
    stop("No degree bins found for reference_ids.", call. = FALSE)
  }
  
  sampled_bins <- sample(ref_bins, n_out, replace = TRUE)
  out <- character(n_out)
  
  for (i in seq_along(sampled_bins)) {
    b <- as.character(sampled_bins[i])
    pool <- degree_bin_pools[[b]]
    if (is.null(pool) || length(pool) == 0) pool <- all_background_ids
    out[i] <- sample(pool, 1)
  }
  
  complete_to_n(out, n_out, all_background_ids)
}

permute_combination_distance <- function(ids, b_perm = 500) {
  ids <- unique(ids)
  vals <- numeric(b_perm)
  
  for (b in seq_len(b_perm)) {
    sampled <- sample_degree_matched(ids, length(ids))
    vals[b] <- mean_distance(sampled)
  }
  
  vals <- vals[!is.na(vals)]
  
  list(
    random_mean = mean(vals),
    random_sd = sd(vals),
    random_values = vals
  )
}

classify_distance <- function(d) {
  ifelse(
    is.na(d), "Unknown",
    ifelse(
      d <= NEAR_DISTANCE_CUTOFF,
      "Disease_close",
      ifelse(d >= FAR_DISTANCE_CUTOFF, "Disease_far", "Intermediate")
    )
  )
}

############################################################
## 6. Build all Baicalein-core combinations
############################################################
combo_list <- list()
idx <- 1

for (k in 0:length(partner_drugs)) {
  cmbs <- combn(partner_drugs, k, simplify = FALSE)
  
  for (cmb in cmbs) {
    drugs <- c(CORE_DRUG, cmb)
    
    combo_list[[idx]] <- data.table(
      Combination_ID = sprintf("COMB_%03d", idx),
      Combination = paste(drugs, collapse = " + "),
      Included_Drugs = paste(drugs, collapse = ";"),
      N_Drugs = length(drugs)
    )
    
    idx <- idx + 1
  }
}

combination_def <- rbindlist(combo_list)

############################################################
## 7. Combination scoring + added target contribution
############################################################
core_ids <- get_drug_ids(CORE_DRUG)
core_distance <- mean_distance(core_ids)

core_d <- as.numeric(dist_lookup[core_ids])
core_close_n <- sum(core_d <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE)
core_far_n <- sum(core_d >= FAR_DISTANCE_CUTOFF, na.rm = TRUE)

combination_results <- rbindlist(lapply(seq_len(nrow(combination_def)), function(i) {
  
  row <- combination_def[i]
  drugs <- unlist(strsplit(row$Included_Drugs, ";", fixed = TRUE))
  ids <- unique(get_drug_ids(drugs))
  added_ids <- setdiff(ids, core_ids)
  
  d <- as.numeric(dist_lookup[ids])
  d <- d[!is.na(d)]
  
  ad <- as.numeric(dist_lookup[added_ids])
  ad <- ad[!is.na(ad)]
  
  combo_close_n <- sum(d <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE)
  combo_far_n <- sum(d >= FAR_DISTANCE_CUTOFF, na.rm = TRUE)
  
  added_close_n <- sum(ad <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE)
  added_far_n <- sum(ad >= FAR_DISTANCE_CUTOFF, na.rm = TRUE)
  added_intermediate_n <- sum(ad > NEAR_DISTANCE_CUTOFF & ad < FAR_DISTANCE_CUTOFF, na.rm = TRUE)
  
  data.table(
    Combination_ID = row$Combination_ID,
    Combination = row$Combination,
    Included_Drugs = row$Included_Drugs,
    N_Drugs = row$N_Drugs,
    Target_N = length(ids),
    
    Mean_Distance = mean_distance(ids),
    Median_Distance = median_distance(ids),
    Min_Distance = min(d, na.rm = TRUE),
    Max_Distance = max(d, na.rm = TRUE),
    
    Combo_Close_Target_N = combo_close_n,
    Combo_Far_Target_N = combo_far_n,
    Fraction_Near_Disease = mean(d <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE),
    Fraction_Far_From_Disease = mean(d >= FAR_DISTANCE_CUTOFF, na.rm = TRUE),
    
    Core_Target_N = length(core_ids),
    Core_Mean_Distance = core_distance,
    Core_Close_Target_N = core_close_n,
    Core_Far_Target_N = core_far_n,
    
    Network_Benefit_vs_Baicalein = core_distance - mean_distance(ids),
    
    Added_Target_N = length(added_ids),
    Added_Close_Target_N = added_close_n,
    Added_Intermediate_Target_N = added_intermediate_n,
    Added_Far_Target_N = added_far_n,
    Added_Mean_Distance = ifelse(length(ad) > 0, mean(ad, na.rm = TRUE), NA_real_),
    Added_Median_Distance = ifelse(length(ad) > 0, median(ad, na.rm = TRUE), NA_real_),
    Added_Close_Gain_vs_Core = combo_close_n - core_close_n,
    Added_Far_Gain_vs_Core = combo_far_n - core_far_n,
    Added_Net_CloseMinusFar = added_close_n - added_far_n,
    Effective_Target_Ratio = ifelse(length(added_ids) > 0, added_close_n / length(added_ids), NA_real_),
    Signal_to_Noise = added_close_n / (added_far_n + 1)
  )
}))

############################################################
## 8. Degree-matched permutation
############################################################
cat("\nRunning degree-matched permutation for real combinations...\n")

perm_results <- list()

for (i in seq_len(nrow(combination_results))) {
  
  cat("Permutation:", i, "/", nrow(combination_results), "\n")
  
  row <- combination_results[i]
  drugs <- unlist(strsplit(row$Included_Drugs, ";", fixed = TRUE))
  ids <- get_drug_ids(drugs)
  
  pr <- permute_combination_distance(ids, B_PERM_COMBINATION)
  
  z <- ifelse(
    is.na(pr$random_sd) || pr$random_sd == 0,
    NA_real_,
    (row$Mean_Distance - pr$random_mean) / pr$random_sd
  )
  
  p_close <- mean(pr$random_values <= row$Mean_Distance, na.rm = TRUE)
  
  perm_results[[i]] <- data.table(
    Combination_ID = row$Combination_ID,
    Combination = row$Combination,
    B_Permutation = B_PERM_COMBINATION,
    Observed_Distance = row$Mean_Distance,
    Random_Mean_Distance = pr$random_mean,
    Random_SD_Distance = pr$random_sd,
    Z_Score = z,
    Empirical_P_Closer = p_close
  )
}

permutation_results <- rbindlist(perm_results)

combination_results <- merge(
  combination_results,
  permutation_results[, .(
    Combination_ID,
    Random_Mean_Distance,
    Random_SD_Distance,
    Z_Score,
    Empirical_P_Closer
  )],
  by = "Combination_ID",
  all.x = TRUE
)

############################################################
## 9. Composite ranking score
############################################################
scale01_good_low <- function(x) {
  if (all(is.na(x)) || max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(rep(0.5, length(x)))
  1 - (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

scale01_good_high <- function(x) {
  if (all(is.na(x)) || max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(rep(0.5, length(x)))
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

combination_results[, Score_Distance := scale01_good_low(Mean_Distance)]
combination_results[, Score_Z := scale01_good_low(Z_Score)]
combination_results[, Score_AddedClose := scale01_good_high(Added_Close_Target_N)]
combination_results[, Score_NetAdded := scale01_good_high(Added_Net_CloseMinusFar)]
combination_results[, Score_AddedMeanDistance := scale01_good_low(Added_Mean_Distance)]
combination_results[is.na(Score_AddedMeanDistance), Score_AddedMeanDistance := Score_Distance]
combination_results[, Score_EffectiveTargetRatio := scale01_good_high(Effective_Target_Ratio)]
combination_results[is.na(Score_EffectiveTargetRatio), Score_EffectiveTargetRatio := Score_Distance]

combination_results[, ACMPI_NPS_CompositeScore :=
                      0.30 * Score_Distance +
                      0.20 * Score_Z +
                      0.18 * Score_AddedClose +
                      0.14 * Score_NetAdded +
                      0.10 * Score_AddedMeanDistance +
                      0.08 * Score_EffectiveTargetRatio]

setorder(
  combination_results,
  -ACMPI_NPS_CompositeScore,
  Mean_Distance,
  -Added_Close_Target_N,
  Added_Far_Target_N
)

############################################################
## 10. Incremental drug contribution
############################################################
incremental_results <- list()
idx <- 1

for (i in seq_len(nrow(combination_results))) {
  
  row <- combination_results[i]
  drugs <- unlist(strsplit(row$Included_Drugs, ";", fixed = TRUE))
  
  if (length(drugs) <= 1) next
  
  full_ids <- get_drug_ids(drugs)
  full_d <- mean_distance(full_ids)
  
  for (drug in setdiff(drugs, CORE_DRUG)) {
    base_drugs <- setdiff(drugs, drug)
    base_ids <- get_drug_ids(base_drugs)
    added_ids <- setdiff(get_drug_ids(drug), base_ids)
    
    base_d <- mean_distance(base_ids)
    ad <- as.numeric(dist_lookup[added_ids])
    ad <- ad[!is.na(ad)]
    
    incremental_results[[idx]] <- data.table(
      Combination_ID = row$Combination_ID,
      Combination = row$Combination,
      Removed_or_Added_Drug = drug,
      Base_Combination = paste(base_drugs, collapse = " + "),
      Base_Distance = base_d,
      Full_Distance = full_d,
      Incremental_Benefit = base_d - full_d,
      Newly_Added_Target_N = length(unique(added_ids)),
      Newly_Added_Close_Target_N = sum(ad <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE),
      Newly_Added_Intermediate_Target_N = sum(ad > NEAR_DISTANCE_CUTOFF & ad < FAR_DISTANCE_CUTOFF, na.rm = TRUE),
      Newly_Added_Far_Target_N = sum(ad >= FAR_DISTANCE_CUTOFF, na.rm = TRUE),
      Newly_Added_Mean_Distance = ifelse(length(ad) > 0, mean(ad), NA_real_),
      Newly_Added_Net_CloseMinusFar =
        sum(ad <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE) -
        sum(ad >= FAR_DISTANCE_CUTOFF, na.rm = TRUE),
      Newly_Effective_Target_Ratio =
        ifelse(length(unique(added_ids)) > 0,
               sum(ad <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE) / length(unique(added_ids)),
               NA_real_),
      Newly_Signal_to_Noise =
        sum(ad <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE) /
        (sum(ad >= FAR_DISTANCE_CUTOFF, na.rm = TRUE) + 1)
    )
    
    idx <- idx + 1
  }
}

incremental_results <- rbindlist(incremental_results)

incremental_summary <- incremental_results[
  ,
  .(
    N_Combination_Context = .N,
    Median_Incremental_Benefit = median(Incremental_Benefit, na.rm = TRUE),
    Mean_Incremental_Benefit = mean(Incremental_Benefit, na.rm = TRUE),
    Median_Newly_Added_Close_Target_N = median(Newly_Added_Close_Target_N, na.rm = TRUE),
    Median_Newly_Added_Far_Target_N = median(Newly_Added_Far_Target_N, na.rm = TRUE),
    Median_Newly_Added_Net_CloseMinusFar = median(Newly_Added_Net_CloseMinusFar, na.rm = TRUE),
    Median_Newly_Added_Mean_Distance = median(Newly_Added_Mean_Distance, na.rm = TRUE),
    Median_Newly_Effective_Target_Ratio = median(Newly_Effective_Target_Ratio, na.rm = TRUE),
    Median_Newly_Signal_to_Noise = median(Newly_Signal_to_Noise, na.rm = TRUE),
    Positive_Addition_Count = sum(Incremental_Benefit > DRIVER_EPS, na.rm = TRUE),
    Negative_Addition_Count = sum(Incremental_Benefit < -DRIVER_EPS, na.rm = TRUE)
  ),
  by = Removed_or_Added_Drug
][order(-Median_Incremental_Benefit, -Median_Newly_Signal_to_Noise)]

############################################################
## 11. Gene-level driver analysis
############################################################
driver_list <- list()

for (i in seq_len(nrow(combination_results))) {
  
  row <- combination_results[i]
  drugs <- unlist(strsplit(row$Included_Drugs, ";", fixed = TRUE))
  ids <- unique(get_drug_ids(drugs))
  
  d <- as.numeric(dist_lookup[ids])
  keep <- !is.na(d)
  ids <- ids[keep]
  d <- d[keep]
  
  if (length(ids) <= 2) next
  
  n <- length(ids)
  full_distance <- mean(d)
  distance_without <- (sum(d) - d) / (n - 1)
  driver_score <- distance_without - full_distance
  
  tmp <- data.table(
    Combination_ID = row$Combination_ID,
    Combination = row$Combination,
    Included_Drugs = row$Included_Drugs,
    string_id = ids,
    preferred_name = gene_lookup[ids],
    distance_to_stroke_disease = d,
    degree_STRING700 = degree_lookup[ids],
    is_stroke_disease_node = disease_flag_lookup[ids],
    Full_Combination_Distance = full_distance,
    Distance_Without_Gene = distance_without,
    DriverScore = driver_score
  )
  
  tmp[, Target_Distance_Class := classify_distance(distance_to_stroke_disease)]
  
  src <- drug_targets[
    string_id %in% ids,
    .(
      Source_Drugs = paste(sort(unique(Drug)), collapse = ";"),
      Source_Input_Genes = paste(sort(unique(Input_Gene)), collapse = ";")
    ),
    by = string_id
  ]
  
  tmp <- merge(tmp, src, by = "string_id", all.x = TRUE)
  driver_list[[i]] <- tmp
}

driver_results <- rbindlist(driver_list, use.names = TRUE, fill = TRUE)

driver_results[, Driver_Class := fifelse(
  DriverScore > DRIVER_EPS,
  "Beneficial_driver",
  fifelse(DriverScore < -DRIVER_EPS, "Dilution_gene", "Neutral_or_redundant")
)]

top_driver_results <- driver_results[
  order(Combination_ID, -DriverScore)
][
  ,
  head(.SD, 30),
  by = Combination_ID
]

top_dilution_results <- driver_results[
  order(Combination_ID, DriverScore)
][
  ,
  head(.SD, 30),
  by = Combination_ID
]

best_by_distance <- combination_results[order(Mean_Distance)][1]
best_by_composite <- combination_results[order(-ACMPI_NPS_CompositeScore)][1]

best_driver_genes <- driver_results[
  Combination_ID == best_by_composite$Combination_ID
][order(-DriverScore)]

############################################################
## 12. Added target detail table
############################################################
added_target_detail <- list()

for (i in seq_len(nrow(combination_results))) {
  
  row <- combination_results[i]
  drugs <- unlist(strsplit(row$Included_Drugs, ";", fixed = TRUE))
  ids <- unique(get_drug_ids(drugs))
  added_ids <- setdiff(ids, core_ids)
  
  if (length(added_ids) == 0) next
  
  tmp <- data.table(
    Combination_ID = row$Combination_ID,
    Combination = row$Combination,
    string_id = added_ids,
    preferred_name = gene_lookup[added_ids],
    distance_to_stroke_disease = as.numeric(dist_lookup[added_ids]),
    degree_STRING700 = as.numeric(degree_lookup[added_ids]),
    is_stroke_disease_node = as.logical(disease_flag_lookup[added_ids])
  )
  
  tmp[, Target_Distance_Class := classify_distance(distance_to_stroke_disease)]
  
  src <- drug_targets[
    string_id %in% added_ids,
    .(
      Source_Drugs = paste(sort(unique(Drug)), collapse = ";"),
      Source_Input_Genes = paste(sort(unique(Input_Gene)), collapse = ";")
    ),
    by = string_id
  ]
  
  tmp <- merge(tmp, src, by = "string_id", all.x = TRUE)
  added_target_detail[[i]] <- tmp
}

added_target_detail <- rbindlist(added_target_detail, use.names = TRUE, fill = TRUE)

############################################################
## 13. Save outputs
############################################################
out_comb <- file.path(OUT_DIR, "01_RealDrug_CombinationRanking_AddedTarget.tsv")
out_perm <- file.path(OUT_DIR, "02_RealDrug_CombinationPermutation.tsv")
out_inc <- file.path(OUT_DIR, "03_RealDrug_IncrementalDrugBenefit.tsv")
out_inc_sum <- file.path(OUT_DIR, "04_RealDrug_IncrementalDrugSummary.tsv")
out_driver <- file.path(OUT_DIR, "05_RealDrug_GeneDriverScores.tsv")
out_top_driver <- file.path(OUT_DIR, "06_RealDrug_TopBeneficialDrivers.tsv")
out_top_dilution <- file.path(OUT_DIR, "07_RealDrug_TopDilutionGenes.tsv")
out_best <- file.path(OUT_DIR, "08_RealDrug_BestCompositeCombination_DriverGenes.tsv")
out_added_detail <- file.path(OUT_DIR, "09_RealDrug_AddedTargetDetail.tsv")

fwrite(combination_results, out_comb, sep = "\t")
fwrite(permutation_results, out_perm, sep = "\t")
fwrite(incremental_results, out_inc, sep = "\t")
fwrite(incremental_summary, out_inc_sum, sep = "\t")
fwrite(driver_results, out_driver, sep = "\t")
fwrite(top_driver_results, out_top_driver, sep = "\t")
fwrite(top_dilution_results, out_top_dilution, sep = "\t")
fwrite(best_driver_genes, out_best, sep = "\t")
fwrite(added_target_detail, out_added_detail, sep = "\t")

############################################################
## 14. No-overflow figures
############################################################
FIG_DIR <- file.path(OUT_DIR, "Figures_NoOverflow")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

combination_results <- combination_results[order(-ACMPI_NPS_CompositeScore)]
combination_results[, Combo_Short := paste0("C", seq_len(.N))]

combo_key <- combination_results[, .(
  Combo_Short,
  Combination_ID,
  Combination,
  Included_Drugs,
  N_Drugs,
  Target_N,
  Mean_Distance,
  ACMPI_NPS_CompositeScore,
  Added_Close_Target_N,
  Added_Far_Target_N,
  Added_Net_CloseMinusFar,
  Effective_Target_Ratio,
  Signal_to_Noise,
  Z_Score
)]

out_key <- file.path(OUT_DIR, "10_CombinationShortLabel_Key.tsv")
fwrite(combo_key, out_key, sep = "\t")

caption_file <- file.path(OUT_DIR, "11_FigureCaptions.txt")
writeLines(c(
  "Fig1: Top 20 combinations ranked by mean distance. Lower values indicate closer proximity to the stroke disease module.",
  "Fig2: Top 20 combinations ranked by ACMPI-NPS composite score.",
  "Fig3: Added disease-close versus disease-far target balance relative to Baicalein.",
  "Fig4: Incremental benefit of each added drug.",
  "Fig5: Multi-metric heatmap of all combinations.",
  "Fig6: Top driver genes in the best composite combination."
), caption_file)

PAL_N <- c(
  "1" = "#BFDCC6",
  "2" = "#C7D8F2",
  "3" = "#F5C7CF",
  "4" = "#F7DFB2",
  "5" = "#D7C7F2",
  "6" = "#BFE5EA"
)

PAL_CLASS <- c(
  "Disease_close" = "#F2A7A3",
  "Intermediate" = "#CDB8E9",
  "Disease_far" = "#AFCBEF",
  "Unknown" = "#D9D9D9"
)

theme_clean <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "grey20"),
      axis.text.y = element_text(size = base_size),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(12, 50, 12, 12)
    )
}

save_both <- function(p, name, width, height) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")),
         p, width = width, height = height, dpi = 600,
         bg = "white", limitsize = FALSE)
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")),
         p, width = width, height = height,
         bg = "white", limitsize = FALSE)
}

fig1_dt <- combination_results[order(Mean_Distance)][1:min(.N, 20)]
fig1_dt[, Combo_Short := factor(Combo_Short, levels = rev(Combo_Short))]
fig1_dt[, N_Drugs_Factor := factor(N_Drugs)]

p1 <- ggplot(fig1_dt, aes(Combo_Short, Mean_Distance, fill = N_Drugs_Factor)) +
  geom_col(width = 0.65, color = "grey35", linewidth = 0.18) +
  geom_text(aes(label = sprintf("%.3f", Mean_Distance)), hjust = -0.08, size = 3.1, color = "grey25") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = PAL_N, name = "Drug number") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Top combinations by mean distance", x = NULL, y = "Mean distance") +
  theme_clean(12)

save_both(p1, "Fig1_Top20_MeanDistance", 8.6, 6.2)

fig2_dt <- combination_results[order(-ACMPI_NPS_CompositeScore)][1:min(.N, 20)]
fig2_dt[, Combo_Short := factor(Combo_Short, levels = rev(Combo_Short))]
fig2_dt[, N_Drugs_Factor := factor(N_Drugs)]

p2 <- ggplot(fig2_dt, aes(Combo_Short, ACMPI_NPS_CompositeScore, fill = N_Drugs_Factor)) +
  geom_col(width = 0.65, color = "grey35", linewidth = 0.18) +
  geom_text(aes(label = sprintf("%.3f", ACMPI_NPS_CompositeScore)), hjust = -0.08, size = 3.1, color = "grey25") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = PAL_N, name = "Drug number") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Top combinations by composite score", x = NULL, y = "Composite score") +
  theme_clean(12)

save_both(p2, "Fig2_Top20_CompositeScore", 8.6, 6.2)

fig3_dt <- combination_results[order(-Added_Net_CloseMinusFar)][1:min(.N, 20)]
fig3_dt[, Combo_Short := factor(Combo_Short, levels = rev(Combo_Short))]

fig3_long <- melt(
  fig3_dt,
  id.vars = "Combo_Short",
  measure.vars = c("Added_Close_Target_N", "Added_Far_Target_N"),
  variable.name = "Target_Class",
  value.name = "Target_N"
)

fig3_long[, Target_Class := fifelse(
  Target_Class == "Added_Close_Target_N",
  "Close targets",
  "Far targets"
)]
fig3_long[Target_Class == "Far targets", Target_N := -Target_N]

p3 <- ggplot(fig3_long, aes(Combo_Short, Target_N, fill = Target_Class)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.15) +
  coord_flip() +
  scale_fill_manual(
    values = c("Close targets" = "#F2A7A3", "Far targets" = "#AFCBEF"),
    name = NULL
  ) +
  labs(title = "Added-target signal and dilution", x = NULL, y = "Added target count") +
  theme_clean(12)

save_both(p3, "Fig3_AddedTargetBalance", 8.6, 6.2)

inc_order <- incremental_summary[order(-Median_Incremental_Benefit)]$Removed_or_Added_Drug
incremental_results[, Removed_or_Added_Drug := factor(Removed_or_Added_Drug, levels = inc_order)]

p4 <- ggplot(
  incremental_results,
  aes(Removed_or_Added_Drug, Incremental_Benefit, fill = Removed_or_Added_Drug)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.86, color = "grey25") +
  geom_jitter(width = 0.10, size = 1.5, alpha = 0.45, color = "grey30") +
  scale_fill_manual(values = hue_pal(l = 82, c = 45)(length(inc_order))) +
  labs(title = "Incremental benefit of added drugs", x = NULL, y = "Incremental benefit") +
  theme_clean(12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "none"
  )

save_both(p4, "Fig4_IncrementalDrugBenefit", 8.8, 5.4)

heat_dt <- copy(combination_results)[order(-ACMPI_NPS_CompositeScore)]
heat_dt[, Combo_Short := factor(Combo_Short, levels = rev(Combo_Short))]

heat_long <- heat_dt[, .(
  Combo_Short,
  Composite = ACMPI_NPS_CompositeScore,
  Distance = -Mean_Distance,
  CloseTargets = Added_Close_Target_N,
  NetCloseFar = Added_Net_CloseMinusFar,
  EffectiveRatio = Effective_Target_Ratio,
  SignalNoise = Signal_to_Noise,
  Zscore = -Z_Score
)]

heat_long <- melt(
  heat_long,
  id.vars = "Combo_Short",
  variable.name = "Metric",
  value.name = "Value"
)

heat_long[, Scaled_Value := as.numeric(scale(Value)), by = Metric]

p5 <- ggplot(heat_long, aes(Metric, Combo_Short, fill = Scaled_Value)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#AFCBEF",
    mid = "white",
    high = "#F2A7A3",
    midpoint = 0,
    name = "Scaled\nscore"
  ) +
  labs(title = "Multi-metric combination landscape", x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 7),
    legend.position = "right",
    plot.margin = margin(12, 22, 12, 12)
  )

save_both(p5, "Fig5_AllCombination_Heatmap", 9.2, 8.8)

best_plot <- best_driver_genes[order(-DriverScore)][1:min(.N, 20)]
best_plot[, preferred_name := factor(preferred_name, levels = rev(preferred_name))]

p6 <- ggplot(best_plot, aes(preferred_name, DriverScore, fill = Target_Distance_Class)) +
  geom_col(width = 0.65, color = "grey35", linewidth = 0.18) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = PAL_CLASS, name = "Target class") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Top driver genes", x = NULL, y = "DriverScore") +
  theme_clean(12)

save_both(p6, "Fig6_TopDriverGenes", 8.6, 6.2)

############################################################
## 15. Final report
############################################################
cat("\n============================================================\n")
cat("Real six-drug ACMPI-NPS final no-overflow analysis completed.\n")
cat("============================================================\n")

cat("\nBest by raw mean distance:\n")
print(best_by_distance)

cat("\nBest by composite score:\n")
print(best_by_composite)

cat("\nIncremental drug summary:\n")
print(incremental_summary)

cat("\nSaved files:\n")
cat(out_comb, "\n")
cat(out_perm, "\n")
cat(out_inc, "\n")
cat(out_inc_sum, "\n")
cat(out_driver, "\n")
cat(out_top_driver, "\n")
cat(out_top_dilution, "\n")
cat(out_best, "\n")
cat(out_added_detail, "\n")
cat(out_key, "\n")
cat(caption_file, "\n")

cat("\nFigures saved to:\n")
cat(FIG_DIR, "\n")
cat("============================================================\n")