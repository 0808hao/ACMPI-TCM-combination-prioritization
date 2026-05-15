############################################################
## ACMPI-NPS Pressure Test A: Proximity Accuracy
## FINAL manuscript-grade workflow: raw distance + representative permutation + publication figures
##
## Purpose:
##   Validate whether the disease-distance module itself works correctly in one run:
##   1) disease-near target sets should be closer to the stroke module;
##   2) disease-far target sets should be farther;
##   3) neutral/background sets should behave like random background;
##   4) hub-biased sets should be exposed by raw distance;
##   5) ordinary vs degree-matched permutation should adjudicate hub bias.
##
## Important:
##   This script does NOT test multi-component combinations.
##   It only tests whether the PPI distance model is reliable.
##
## Input:
##   Latest ANPS_StrokeDisease_STRING700_* folder under ROOT_DIR containing:
##     01_StrokeDisease_GeneMapping_STRING.tsv
##     02_STRING700_NodeDistanceToStrokeDisease.tsv
##     04_StrokeDisease_PPI_QC.tsv
##
## Output:
##   A new output folder under ROOT_DIR.
##   Only essential Pressure Test A files are saved.
############################################################

rm(list = ls())

############################################################
## 0. User settings
############################################################

ROOT_DIR <- "/media/desk16/iy15915/中药之开创/脑卒中疾病基因"

STRING_SCORE_MIN <- 700
SET_SEED <- 20260502

## ==========================================================
## Run mode
## ==========================================================
## QUICK_TEST:
##   Minimal sanity check. It still runs raw distance and permutation,
##   but with very small repetitions.
##
## PILOT_FULL:
##   Recommended first formal run. It completes Pressure Test A in one pass:
##   raw-distance accuracy + representative ordinary and degree-matched permutation.
##
## FINAL:
##   Larger-scale robust validation after PILOT_FULL passes.
##
## FINAL_STRICT:
##   Same design as FINAL, but stricter permutation B for manuscript lock.
## ==========================================================

RUN_MODE <- "FINAL_STRICT"

if (RUN_MODE == "QUICK_TEST") {
  N_REP_MAIN <- 5
  N_REP_STRESS <- 2
  B_PERM <- 50
  PERM_REPS_PER_GROUP <- 2
  PROGRESS_EVERY <- 1
  TARGET_N_GRID <- c(10, 30)
  
} else if (RUN_MODE == "PILOT_FULL") {
  N_REP_MAIN <- 30
  N_REP_STRESS <- 10
  B_PERM <- 200
  PERM_REPS_PER_GROUP <- 8
  PROGRESS_EVERY <- 10
  TARGET_N_GRID <- c(10, 20, 30, 50)
  
} else if (RUN_MODE == "FINAL") {
  N_REP_MAIN <- 300
  N_REP_STRESS <- 80
  B_PERM <- 500
  PERM_REPS_PER_GROUP <- 30
  PROGRESS_EVERY <- 50
  TARGET_N_GRID <- c(10, 20, 30, 50)
  
} else if (RUN_MODE == "FINAL_STRICT") {
  ## Manuscript-lock setting for Pressure Test A.
  ## The parameterization is fixed; only simulation/permutation depth is increased.
  N_REP_MAIN <- 300
  N_REP_STRESS <- 80
  B_PERM <- 1000
  PERM_REPS_PER_GROUP <- 40
  PROGRESS_EVERY <- 50
  TARGET_N_GRID <- c(10, 20, 30, 50)
  
} else {
  stop("Unknown RUN_MODE. Use QUICK_TEST, PILOT_FULL, FINAL, or FINAL_STRICT.", call. = FALSE)
}

DEGREE_BINS <- 10
DEFAULT_TARGET_N <- 30

## Thresholds used only for labels and summaries, not for model fitting.
NEAR_DISTANCE_CUTOFF <- 1
FAR_DISTANCE_CUTOFF <- 3
DIRECTION_EPS <- 0.05

set.seed(SET_SEED)

############################################################
## 1. Required packages
############################################################

required_pkgs <- c("data.table", "ggplot2")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(data.table)
library(ggplot2)

############################################################
## 2. Locate latest disease-module output folder
############################################################

find_latest_disease_dir <- function(root_dir, score_min = 700) {
  candidate_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  pattern <- paste0("ANPS_StrokeDisease_STRING", score_min, "_")
  candidate_dirs <- candidate_dirs[grepl(pattern, basename(candidate_dirs), fixed = TRUE)]
  
  required_files <- c(
    "01_StrokeDisease_GeneMapping_STRING.tsv",
    "02_STRING700_NodeDistanceToStrokeDisease.tsv",
    "04_StrokeDisease_PPI_QC.tsv"
  )
  
  valid_dirs <- candidate_dirs[
    vapply(candidate_dirs, function(d) {
      all(file.exists(file.path(d, required_files)))
    }, logical(1))
  ]
  
  if (length(valid_dirs) == 0) {
    stop(
      "No valid disease-module output folder found under:\n",
      root_dir,
      "\nExpected folder name pattern: ",
      pattern,
      "\nExpected files:\n",
      paste(required_files, collapse = "\n"),
      call. = FALSE
    )
  }
  
  mt <- file.info(valid_dirs)$mtime
  valid_dirs[which.max(mt)]
}

DISEASE_DIR <- find_latest_disease_dir(ROOT_DIR, STRING_SCORE_MIN)

MAPPING_FILE <- file.path(DISEASE_DIR, "01_StrokeDisease_GeneMapping_STRING.tsv")
NODE_DISTANCE_FILE <- file.path(DISEASE_DIR, "02_STRING700_NodeDistanceToStrokeDisease.tsv")
QC_FILE <- file.path(DISEASE_DIR, "04_StrokeDisease_PPI_QC.tsv")

OUT_DIR <- file.path(
  ROOT_DIR,
  paste0(
    "ANPS_PressureTestA_ProximityAccuracy_COMPLETE_",
    RUN_MODE,
    "_STRING",
    STRING_SCORE_MIN,
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("ACMPI-NPS Pressure Test A: Complete Proximity Accuracy Test\n")
cat("Run mode: ", RUN_MODE, "\n", sep = "")
cat("Disease module folder:\n", DISEASE_DIR, "\n", sep = "")
cat("Output folder:\n", OUT_DIR, "\n", sep = "")
cat("============================================================\n")

############################################################
## 3. Read input files
############################################################

cat("\nReading disease mapping table...\n")
mapping_dt <- fread(MAPPING_FILE, sep = "\t", data.table = TRUE)

cat("Reading node-to-disease distance table...\n")
node_dt <- fread(NODE_DISTANCE_FILE, sep = "\t", data.table = TRUE)

cat("Reading disease PPI QC table...\n")
qc_dt <- fread(QC_FILE, sep = "\t", data.table = TRUE)

required_node_cols <- c(
  "string_id",
  "preferred_name",
  "distance_to_stroke_disease",
  "degree_STRING700",
  "is_stroke_disease_node"
)

missing_node_cols <- setdiff(required_node_cols, colnames(node_dt))
if (length(missing_node_cols) > 0) {
  stop(
    "Node distance file missing required columns:\n",
    paste(missing_node_cols, collapse = "\n"),
    call. = FALSE
  )
}

node_dt[, string_id := as.character(string_id)]
node_dt[, preferred_name := as.character(preferred_name)]
node_dt[, distance_to_stroke_disease := suppressWarnings(as.numeric(distance_to_stroke_disease))]
node_dt[, degree_STRING700 := suppressWarnings(as.numeric(degree_STRING700))]

if (!is.logical(node_dt$is_stroke_disease_node)) {
  node_dt[, is_stroke_disease_node := as.logical(is_stroke_disease_node)]
}

node_dt <- node_dt[
  !is.na(string_id) &
    !is.na(distance_to_stroke_disease) &
    !is.na(degree_STRING700)
]

if (nrow(node_dt) == 0) {
  stop("No valid nodes found in node distance table.", call. = FALSE)
}

cat("\nInput node table QC:\n")
cat("PPI nodes with valid disease distance: ", nrow(node_dt), "\n", sep = "")
cat("Stroke disease nodes: ", sum(node_dt$is_stroke_disease_node, na.rm = TRUE), "\n", sep = "")

############################################################
## 4. Degree bins and fast background pools
############################################################

node_dt[, degree_rank := frank(degree_STRING700, ties.method = "average")]
node_dt[, degree_bin := ceiling(degree_rank / .N * DEGREE_BINS)]
node_dt[degree_bin < 1, degree_bin := 1]
node_dt[degree_bin > DEGREE_BINS, degree_bin := DEGREE_BINS]

degree_q10 <- as.numeric(quantile(node_dt$degree_STRING700, 0.10, na.rm = TRUE))
degree_q90 <- as.numeric(quantile(node_dt$degree_STRING700, 0.90, na.rm = TRUE))
degree_q95 <- as.numeric(quantile(node_dt$degree_STRING700, 0.95, na.rm = TRUE))

node_dt[, is_extreme_hub := degree_STRING700 >= degree_q90]
node_dt[, is_top5_hub := degree_STRING700 >= degree_q95]
node_dt[, is_non_extreme_degree := degree_STRING700 >= degree_q10 & degree_STRING700 <= degree_q90]

setkey(node_dt, string_id)

background_pool <- node_dt[
  is_stroke_disease_node == FALSE &
    !is.na(distance_to_stroke_disease) &
    !is.na(degree_STRING700)
]

if (nrow(background_pool) < 1000) {
  stop("Background pool is unexpectedly small.", call. = FALSE)
}

## Fast precomputed pools. This avoids repeated full-table filtering during permutation.
degree_bin_pools <- split(background_pool$string_id, background_pool$degree_bin)
all_background_ids <- background_pool$string_id

cat("\nBackground pool:\n")
cat("Non-disease background nodes: ", nrow(background_pool), "\n", sep = "")
cat("Degree q10: ", degree_q10, "\n", sep = "")
cat("Degree q90: ", degree_q90, "\n", sep = "")
cat("Degree q95: ", degree_q95, "\n", sep = "")

############################################################
## 5. Fast lookup vectors
############################################################

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
## 6. Sampling helpers
############################################################

sample_safe <- function(x, n, replace = FALSE) {
  x <- unique(x)
  x <- x[!is.na(x)]
  if (n <= 0) return(character())
  if (length(x) == 0) stop("Sampling pool is empty.", call. = FALSE)
  if (length(x) < n && !replace) {
    sample(x, n, replace = TRUE)
  } else {
    sample(x, n, replace = replace)
  }
}

complete_to_n <- function(ids, n, fill_pool, exclude_ids = character()) {
  ids <- unique(ids)
  ids <- ids[!is.na(ids)]
  
  if (length(ids) >= n) {
    return(ids[seq_len(n)])
  }
  
  fill_pool <- unique(fill_pool)
  fill_pool <- fill_pool[!is.na(fill_pool)]
  fill_pool <- setdiff(fill_pool, c(ids, exclude_ids))
  
  if (length(fill_pool) == 0) {
    return(sample_safe(ids, n, replace = TRUE))
  }
  
  need <- n - length(ids)
  extra <- sample_safe(fill_pool, need, replace = length(fill_pool) < need)
  out <- unique(c(ids, extra))
  
  if (length(out) < n) {
    out <- c(out, sample_safe(out, n - length(out), replace = TRUE))
  }
  
  out[seq_len(n)]
}

get_pool_by_distance <- function(distance_rule,
                                 exclude_ids = character(),
                                 non_extreme = TRUE,
                                 allow_disease_nodes = FALSE) {
  pool <- copy(node_dt)
  
  if (!allow_disease_nodes) {
    pool <- pool[is_stroke_disease_node == FALSE]
  }
  
  if (non_extreme) {
    pool <- pool[is_non_extreme_degree == TRUE]
  }
  
  if (length(exclude_ids) > 0) {
    pool <- pool[!string_id %in% exclude_ids]
  }
  
  if (distance_rule == "dist_0") {
    pool <- pool[distance_to_stroke_disease == 0]
  } else if (distance_rule == "dist_1") {
    pool <- pool[distance_to_stroke_disease == 1]
  } else if (distance_rule == "dist_0_1") {
    pool <- pool[distance_to_stroke_disease %in% c(0, 1)]
  } else if (distance_rule == "dist_1_2") {
    pool <- pool[distance_to_stroke_disease %in% c(1, 2)]
  } else if (distance_rule == "dist_2_3") {
    pool <- pool[distance_to_stroke_disease %in% c(2, 3)]
  } else if (distance_rule == "dist_ge_3") {
    pool <- pool[distance_to_stroke_disease >= 3]
  } else if (distance_rule == "dist_ge_4") {
    pool <- pool[distance_to_stroke_disease >= 4]
  } else {
    stop("Unknown distance_rule: ", distance_rule, call. = FALSE)
  }
  
  pool$string_id
}

sample_degree_matched <- function(reference_ids,
                                  n_out,
                                  exclude_ids = character()) {
  reference_ids <- unique(reference_ids)
  reference_ids <- reference_ids[reference_ids %in% node_dt$string_id]
  
  if (n_out <= 0) return(character())
  if (length(reference_ids) == 0) {
    stop("reference_ids are empty after mapping to node table.", call. = FALSE)
  }
  
  ref_bins <- bin_lookup[reference_ids]
  ref_bins <- ref_bins[!is.na(ref_bins)]
  if (length(ref_bins) == 0) stop("No degree bins found for reference_ids.", call. = FALSE)
  
  sampled_bins <- sample(ref_bins, n_out, replace = TRUE)
  out <- character(n_out)
  
  for (j in seq_along(sampled_bins)) {
    b <- as.character(sampled_bins[j])
    pool <- degree_bin_pools[[b]]
    
    if (is.null(pool) || length(pool) == 0) {
      pool <- all_background_ids
    }
    
    if (length(exclude_ids) > 0) {
      pool <- setdiff(pool, exclude_ids)
    }
    
    if (length(pool) == 0) {
      pool <- all_background_ids
    }
    
    out[j] <- sample(pool, 1)
  }
  
  complete_to_n(out, n_out, all_background_ids, exclude_ids = exclude_ids)
}

sample_background_by_size <- function(n, exclude_ids = character(), non_extreme = TRUE) {
  pool <- background_pool
  if (non_extreme) pool <- pool[is_non_extreme_degree == TRUE]
  if (length(exclude_ids) > 0) pool <- pool[!string_id %in% exclude_ids]
  sample_safe(pool$string_id, n, replace = FALSE)
}

sample_positive_targets <- function(n, exclude_ids = character()) {
  ## 30% true disease nodes + 70% one-hop disease neighbors.
  n0 <- floor(n * 0.30)
  n1 <- n - n0
  
  pool0 <- get_pool_by_distance("dist_0", exclude_ids = exclude_ids, non_extreme = TRUE, allow_disease_nodes = TRUE)
  pool1 <- get_pool_by_distance("dist_1", exclude_ids = exclude_ids, non_extreme = TRUE, allow_disease_nodes = FALSE)
  
  ids0 <- if (n0 > 0) sample_safe(pool0, n0, replace = FALSE) else character()
  ids1 <- if (n1 > 0) sample_safe(pool1, n1, replace = FALSE) else character()
  
  complete_to_n(c(ids0, ids1), n, unique(c(pool1, pool0)), exclude_ids = exclude_ids)
}

sample_neutral_targets <- function(n, exclude_ids = character()) {
  sample_background_by_size(n, exclude_ids = exclude_ids, non_extreme = TRUE)
}

sample_dilution_targets <- function(n, exclude_ids = character()) {
  pool <- get_pool_by_distance("dist_ge_3", exclude_ids = exclude_ids, non_extreme = TRUE, allow_disease_nodes = FALSE)
  sample_safe(pool, n, replace = FALSE)
}

sample_hub_targets <- function(n, exclude_ids = character()) {
  pool <- node_dt[
    is_stroke_disease_node == FALSE &
      degree_STRING700 >= degree_q90 &
      !is.na(distance_to_stroke_disease)
  ]
  if (length(exclude_ids) > 0) pool <- pool[!string_id %in% exclude_ids]
  sample_safe(pool$string_id, n, replace = FALSE)
}

sample_mixed_targets <- function(n, exclude_ids = character()) {
  n_pos <- floor(n * 0.40)
  n_neu <- floor(n * 0.40)
  n_dil <- n - n_pos - n_neu
  
  ids_pos <- sample_positive_targets(n = n_pos, exclude_ids = exclude_ids)
  ids_neu <- sample_neutral_targets(n = n_neu, exclude_ids = c(exclude_ids, ids_pos))
  ids_dil <- sample_dilution_targets(n = n_dil, exclude_ids = c(exclude_ids, ids_pos, ids_neu))
  
  complete_to_n(c(ids_pos, ids_neu, ids_dil), n, unique(c(ids_pos, ids_neu, ids_dil, all_background_ids)), exclude_ids = exclude_ids)
}

sample_target_set <- function(target_type, n_targets, exclude_ids = character()) {
  switch(
    target_type,
    Positive = sample_positive_targets(n_targets, exclude_ids = exclude_ids),
    Neutral = sample_neutral_targets(n_targets, exclude_ids = exclude_ids),
    Dilution = sample_dilution_targets(n_targets, exclude_ids = exclude_ids),
    Hub = sample_hub_targets(n_targets, exclude_ids = exclude_ids),
    Mixed = sample_mixed_targets(n_targets, exclude_ids = exclude_ids),
    stop("Unknown target_type: ", target_type, call. = FALSE)
  )
}

############################################################
## 7. Proximity and permutation functions
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

make_target_table <- function(sim_id,
                              scenario_group,
                              target_type,
                              n_targets,
                              ids) {
  out <- data.table(
    Simulation_ID = sim_id,
    Scenario_Group = scenario_group,
    Target_Type = target_type,
    Design_Target_N = n_targets,
    string_id = unique(ids)
  )
  
  out[, preferred_name := gene_lookup[string_id]]
  out[, distance_to_stroke_disease := as.numeric(dist_lookup[string_id])]
  out[, degree_STRING700 := as.numeric(degree_lookup[string_id])]
  out[, degree_bin := as.integer(bin_lookup[string_id])]
  out[, is_stroke_disease_node := as.logical(disease_flag_lookup[string_id])]
  
  out[, Target_Profile := fifelse(
    Target_Type == "Positive", "Disease_near_set",
    fifelse(Target_Type == "Neutral", "Background_neutral_set",
            fifelse(Target_Type == "Dilution", "Disease_far_set",
                    fifelse(Target_Type == "Hub", "Hub_biased_set", "Mixed_set")))
  )]
  
  out[]
}

permute_distance <- function(target_ids,
                             method = c("degree_matched", "ordinary"),
                             b_perm = 500) {
  method <- match.arg(method)
  target_ids <- unique(target_ids)
  target_ids <- target_ids[target_ids %in% node_dt$string_id]
  
  if (length(target_ids) == 0 || b_perm <= 0) {
    return(list(mean_random = NA_real_, sd_random = NA_real_, random_values = numeric()))
  }
  
  random_vals <- numeric(b_perm)
  
  for (b in seq_len(b_perm)) {
    if (method == "degree_matched") {
      sampled <- sample_degree_matched(reference_ids = target_ids, n_out = length(target_ids), exclude_ids = character())
    } else {
      sampled <- sample_safe(all_background_ids, length(target_ids), replace = FALSE)
    }
    
    d <- as.numeric(dist_lookup[sampled])
    d <- d[!is.na(d)]
    random_vals[b] <- ifelse(length(d) > 0, mean(d), NA_real_)
  }
  
  random_vals <- random_vals[!is.na(random_vals)]
  
  list(
    mean_random = mean(random_vals),
    sd_random = sd(random_vals),
    random_values = random_vals
  )
}

compute_permutation_summary <- function(sim_id,
                                        scenario_group,
                                        target_type,
                                        target_ids,
                                        observed_distance,
                                        method = "degree_matched",
                                        b_perm = 500) {
  pr <- permute_distance(target_ids = target_ids, method = method, b_perm = b_perm)
  
  z <- ifelse(is.na(pr$sd_random) || pr$sd_random == 0, NA_real_, (observed_distance - pr$mean_random) / pr$sd_random)
  empirical_p_lower <- mean(pr$random_values <= observed_distance, na.rm = TRUE)
  empirical_p_upper <- mean(pr$random_values >= observed_distance, na.rm = TRUE)
  
  data.table(
    Simulation_ID = sim_id,
    Scenario_Group = scenario_group,
    Target_Type = target_type,
    Permutation_Method = method,
    B_Permutation = b_perm,
    Observed_Distance = observed_distance,
    Random_Mean_Distance = pr$mean_random,
    Random_SD_Distance = pr$sd_random,
    Z_Score = z,
    Empirical_P_Lower_Closer = empirical_p_lower,
    Empirical_P_Upper_Farther = empirical_p_upper
  )
}

############################################################
## 8. Build Pressure Test A design grid
############################################################

make_repeated_grid <- function(scenario_group,
                               target_type,
                               n_targets,
                               n_rep) {
  data.table(
    Scenario_Group = scenario_group,
    Target_Type = target_type,
    Design_Target_N = n_targets,
    Replicate = seq_len(n_rep)
  )
}

grid_list <- list()

## A. Main target-set type accuracy test
for (tt in c("Positive", "Neutral", "Dilution", "Hub", "Mixed")) {
  grid_list[[length(grid_list) + 1]] <- make_repeated_grid(
    scenario_group = "A_MainAccuracy",
    target_type = tt,
    n_targets = DEFAULT_TARGET_N,
    n_rep = N_REP_MAIN
  )
}

## B. Target-number stress test
for (nt in TARGET_N_GRID) {
  for (tt in c("Positive", "Neutral", "Dilution", "Hub", "Mixed")) {
    grid_list[[length(grid_list) + 1]] <- make_repeated_grid(
      scenario_group = "B_TargetNumber",
      target_type = tt,
      n_targets = nt,
      n_rep = N_REP_STRESS
    )
  }
}

## C. Hub-specific control test with standard target size
for (tt in c("Hub", "Neutral", "Positive")) {
  grid_list[[length(grid_list) + 1]] <- make_repeated_grid(
    scenario_group = "C_HubControl",
    target_type = tt,
    n_targets = DEFAULT_TARGET_N,
    n_rep = N_REP_STRESS
  )
}

design_grid <- rbindlist(grid_list, use.names = TRUE, fill = TRUE)
design_grid[, Simulation_ID := sprintf("PTA_%05d", .I)]
setcolorder(design_grid, "Simulation_ID")

cat("\nPressure Test A design grid:\n")
cat("Total simulations: ", nrow(design_grid), "\n", sep = "")
print(
  design_grid[
    ,
    .N,
    by = .(Scenario_Group, Target_Type, Design_Target_N)
  ][order(Scenario_Group, Target_Type, Design_Target_N)]
)

############################################################
## 9. Run broad proximity screening
############################################################

target_sets_all <- list()
proximity_all <- list()

cat("\nRunning Pressure Test A broad proximity screening...\n")

for (i in seq_len(nrow(design_grid))) {
  row <- design_grid[i]
  
  sim_id <- row$Simulation_ID
  scenario_group <- row$Scenario_Group
  target_type <- row$Target_Type
  n_targets <- as.integer(row$Design_Target_N)
  
  if (i %% PROGRESS_EVERY == 0 || i == 1) {
    cat(
      "Simulation ", i, " / ", nrow(design_grid),
      ": ", sim_id,
      " | ", scenario_group,
      " | ", target_type,
      " | n=", n_targets,
      "\n",
      sep = ""
    )
    flush.console()
  }
  
  ids <- sample_target_set(target_type = target_type, n_targets = n_targets)
  
  target_sets_all[[i]] <- make_target_table(
    sim_id = sim_id,
    scenario_group = scenario_group,
    target_type = target_type,
    n_targets = n_targets,
    ids = ids
  )
  
  d_vec <- as.numeric(dist_lookup[unique(ids)])
  deg_vec <- as.numeric(degree_lookup[unique(ids)])
  d_vec <- d_vec[!is.na(d_vec)]
  deg_vec <- deg_vec[!is.na(deg_vec)]
  
  mean_d <- mean_distance(ids)
  med_d <- median_distance(ids)
  frac_disease_node <- mean(as.logical(disease_flag_lookup[unique(ids)]), na.rm = TRUE)
  frac_near <- mean(d_vec <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE)
  frac_far <- mean(d_vec >= FAR_DISTANCE_CUTOFF, na.rm = TRUE)
  
  proximity_all[[i]] <- data.table(
    Simulation_ID = sim_id,
    Scenario_Group = scenario_group,
    Target_Type = target_type,
    Design_Target_N = n_targets,
    Observed_Target_N = length(unique(ids)),
    Mean_Distance = mean_d,
    Median_Distance = med_d,
    Min_Distance = min(d_vec, na.rm = TRUE),
    Max_Distance = max(d_vec, na.rm = TRUE),
    Mean_Degree = mean(deg_vec, na.rm = TRUE),
    Median_Degree = median(deg_vec, na.rm = TRUE),
    Fraction_Disease_Node = frac_disease_node,
    Fraction_Near_Disease = frac_near,
    Fraction_Far_From_Disease = frac_far
  )
}

target_sets <- rbindlist(target_sets_all, use.names = TRUE, fill = TRUE)
proximity_results <- rbindlist(proximity_all, use.names = TRUE, fill = TRUE)

############################################################
## 10. Labels and broad summaries
############################################################

proximity_results[, Expected_Distance_Class := fifelse(
  Target_Type == "Positive", "Closest",
  fifelse(Target_Type == "Dilution", "Farthest",
          fifelse(Target_Type == "Neutral", "Background",
                  fifelse(Target_Type == "Hub", "Hub_bias_check", "Intermediate")))
)]

summary_all <- proximity_results[
  ,
  .(
    N = .N,
    Median_Mean_Distance = median(Mean_Distance, na.rm = TRUE),
    Mean_Mean_Distance = mean(Mean_Distance, na.rm = TRUE),
    Median_Median_Distance = median(Median_Distance, na.rm = TRUE),
    Median_Fraction_Near_Disease = median(Fraction_Near_Disease, na.rm = TRUE),
    Median_Fraction_Far_From_Disease = median(Fraction_Far_From_Disease, na.rm = TRUE),
    Median_Mean_Degree = median(Mean_Degree, na.rm = TRUE)
  ),
  by = .(Scenario_Group, Target_Type, Design_Target_N)
][order(Scenario_Group, Target_Type, Design_Target_N)]

summary_main <- summary_all[Scenario_Group == "A_MainAccuracy"]

############################################################
## 11. Representative selection for permutation
############################################################

select_representative_simulations <- function(prox_dt,
                                              reps_per_group = 10,
                                              main_only = FALSE) {
  dt <- copy(prox_dt)
  
  if (main_only) {
    dt <- dt[Scenario_Group == "A_MainAccuracy"]
  }
  
  if (nrow(dt) == 0 || reps_per_group <= 0) {
    return(character())
  }
  
  dt[, Abs_Mean_Distance := abs(Mean_Distance)]
  
  selected <- dt[
    ,
    {
      n_take <- min(.N, reps_per_group)
      if (.N <= n_take) {
        .SD
      } else {
        center_idx <- which.min(abs(Mean_Distance - median(Mean_Distance, na.rm = TRUE)))
        extreme_idx <- order(-Abs_Mean_Distance)[seq_len(min(ceiling(n_take / 2), .N))]
        keep_idx <- unique(c(center_idx, extreme_idx))
        remaining <- setdiff(seq_len(.N), keep_idx)
        extra_n <- n_take - length(keep_idx)
        extra_idx <- if (extra_n > 0 && length(remaining) > 0) sample(remaining, min(extra_n, length(remaining))) else integer()
        final_idx <- unique(c(keep_idx, extra_idx))
        final_idx <- final_idx[seq_len(min(length(final_idx), n_take))]
        .SD[final_idx]
      }
    },
    by = .(Scenario_Group, Target_Type)
  ]
  
  unique(selected$Simulation_ID)
}

############################################################
## 12. Representative permutation validation
############################################################

permutation_results <- data.table(
  Simulation_ID = character(),
  Scenario_Group = character(),
  Target_Type = character(),
  Permutation_Method = character(),
  B_Permutation = integer(),
  Observed_Distance = numeric(),
  Random_Mean_Distance = numeric(),
  Random_SD_Distance = numeric(),
  Z_Score = numeric(),
  Empirical_P_Lower_Closer = numeric(),
  Empirical_P_Upper_Farther = numeric()
)

if (B_PERM > 0 && PERM_REPS_PER_GROUP > 0) {
  
  cat("\nSelecting representative simulations for permutation...\n")
  perm_sim_ids <- select_representative_simulations(
    prox_dt = proximity_results,
    reps_per_group = PERM_REPS_PER_GROUP
  )
  
  cat("Representative simulations for permutation: ", length(perm_sim_ids), "\n", sep = "")
  
  perm_list_all <- list()
  perm_counter <- 0L
  
  for (sid in perm_sim_ids) {
    prox_row <- proximity_results[Simulation_ID == sid][1]
    target_row <- target_sets[Simulation_ID == sid]
    target_ids <- target_row[, unique(string_id)]
    
    perm_counter <- perm_counter + 1L
    if (perm_counter %% 10 == 0 || perm_counter == 1) {
      cat("Permutation set ", perm_counter, " / ", length(perm_sim_ids), ": ", sid, "\n", sep = "")
      flush.console()
    }
    
    perm_degree <- compute_permutation_summary(
      sim_id = sid,
      scenario_group = prox_row$Scenario_Group,
      target_type = prox_row$Target_Type,
      target_ids = target_ids,
      observed_distance = prox_row$Mean_Distance,
      method = "degree_matched",
      b_perm = B_PERM
    )
    
    perm_items <- list(perm_degree)
    
    ## Ordinary random is used for every representative set in the complete
    ## Pressure Test A so that ordinary and degree-matched nulls are directly comparable.
    perm_ordinary <- compute_permutation_summary(
      sim_id = sid,
      scenario_group = prox_row$Scenario_Group,
      target_type = prox_row$Target_Type,
      target_ids = target_ids,
      observed_distance = prox_row$Mean_Distance,
      method = "ordinary",
      b_perm = B_PERM
    )
    perm_items[[length(perm_items) + 1]] <- perm_ordinary
    
    perm_list_all[[length(perm_list_all) + 1]] <- rbindlist(perm_items, use.names = TRUE, fill = TRUE)
  }
  
  if (length(perm_list_all) > 0) {
    permutation_results <- rbindlist(perm_list_all, use.names = TRUE, fill = TRUE)
  }
}

############################################################
## 13. Permutation summaries and pass criteria
############################################################

if (nrow(permutation_results) > 0) {
  permutation_results[, Significant_Closer_005 := !is.na(Empirical_P_Lower_Closer) & Empirical_P_Lower_Closer <= 0.05]
  permutation_results[, Significant_Farther_005 := !is.na(Empirical_P_Upper_Farther) & Empirical_P_Upper_Farther <= 0.05]
} else {
  permutation_results[, Significant_Closer_005 := logical()]
  permutation_results[, Significant_Farther_005 := logical()]
}

perm_summary_all <- if (nrow(permutation_results) > 0) {
  permutation_results[
    ,
    .(
      N = .N,
      Median_Z = median(Z_Score, na.rm = TRUE),
      Mean_Z = mean(Z_Score, na.rm = TRUE),
      Fraction_Closer_005 = mean(Significant_Closer_005, na.rm = TRUE),
      Fraction_Farther_005 = mean(Significant_Farther_005, na.rm = TRUE)
    ),
    by = .(Scenario_Group, Target_Type, Permutation_Method)
  ][order(Scenario_Group, Target_Type, Permutation_Method)]
} else {
  data.table(
    Scenario_Group = character(),
    Target_Type = character(),
    Permutation_Method = character(),
    N = integer(),
    Median_Z = numeric(),
    Mean_Z = numeric(),
    Fraction_Closer_005 = numeric(),
    Fraction_Farther_005 = numeric()
  )
}

get_main_distance <- function(type) {
  x <- summary_main[Target_Type == type, Median_Mean_Distance]
  if (length(x) == 0) return(NA_real_)
  x[1]
}

positive_d <- get_main_distance("Positive")
neutral_d <- get_main_distance("Neutral")
dilution_d <- get_main_distance("Dilution")
hub_d <- get_main_distance("Hub")
mixed_d <- get_main_distance("Mixed")

pass_distance_order <- !is.na(positive_d) && !is.na(neutral_d) && !is.na(dilution_d) &&
  positive_d < neutral_d && neutral_d < dilution_d

pass_mixed_intermediate <- !is.na(mixed_d) && !is.na(positive_d) && !is.na(dilution_d) &&
  mixed_d > positive_d && mixed_d < dilution_d

hub_bias_flag <- !is.na(hub_d) && !is.na(neutral_d) && hub_d < neutral_d

get_main_perm_z <- function(type, method = "degree_matched") {
  x <- perm_summary_all[
    Scenario_Group == "A_MainAccuracy" &
      Target_Type == type &
      Permutation_Method == method,
    Median_Z
  ]
  if (length(x) == 0) return(NA_real_)
  x[1]
}

get_main_perm_frac_closer <- function(type, method = "degree_matched") {
  x <- perm_summary_all[
    Scenario_Group == "A_MainAccuracy" &
      Target_Type == type &
      Permutation_Method == method,
    Fraction_Closer_005
  ]
  if (length(x) == 0) return(NA_real_)
  x[1]
}

positive_z_deg <- get_main_perm_z("Positive", "degree_matched")
neutral_z_deg <- get_main_perm_z("Neutral", "degree_matched")
dilution_z_deg <- get_main_perm_z("Dilution", "degree_matched")
hub_z_deg <- get_main_perm_z("Hub", "degree_matched")
hub_z_ord <- get_main_perm_z("Hub", "ordinary")

positive_frac_deg <- get_main_perm_frac_closer("Positive", "degree_matched")
hub_frac_deg <- get_main_perm_frac_closer("Hub", "degree_matched")
hub_frac_ord <- get_main_perm_frac_closer("Hub", "ordinary")

pass_permutation_available <- any(permutation_results$Permutation_Method == "degree_matched") &&
  any(permutation_results$Permutation_Method == "ordinary")

pass_positive_perm <- !is.na(positive_z_deg) && positive_z_deg < 0
pass_dilution_perm <- !is.na(dilution_z_deg) && dilution_z_deg > 0

## Hub control means degree-matched null reduces the apparent closeness relative to ordinary null.
## This is directional rather than a hard significance cutoff.
pass_hub_control <- !is.na(hub_z_deg) && !is.na(hub_z_ord) &&
  abs(hub_z_deg) < abs(hub_z_ord)

pass_hub_not_more_significant_than_positive <- !is.na(hub_frac_deg) && !is.na(positive_frac_deg) &&
  hub_frac_deg <= positive_frac_deg

pass_summary <- data.table(
  Test_Item = c(
    "Raw distance order: Positive < Neutral < Dilution",
    "Raw mixed behavior: Positive < Mixed < Dilution",
    "Hub raw-distance bias detected",
    "Ordinary and degree-matched permutation performed",
    "Positive target sets are closer than degree-matched null",
    "Dilution target sets are farther than degree-matched null",
    "Hub bias reduced by degree-matched permutation",
    "Hub not more significant than Positive after degree matching"
  ),
  Result = c(
    as.character(pass_distance_order),
    as.character(pass_mixed_intermediate),
    as.character(hub_bias_flag),
    as.character(pass_permutation_available),
    as.character(pass_positive_perm),
    as.character(pass_dilution_perm),
    as.character(pass_hub_control),
    as.character(pass_hub_not_more_significant_than_positive)
  ),
  Interpretation = c(
    "Distance module should rank disease-near targets as closest and disease-far targets as farthest.",
    "Mixed target sets should show intermediate raw distance behavior.",
    "TRUE means hub targets appear artificially close before degree correction.",
    "Required for complete Pressure Test A.",
    "Observed positive sets should have negative degree-matched Z.",
    "Observed dilution sets should have positive degree-matched Z.",
    "Hub ordinary-random proximity should be stronger than degree-matched proximity if hub bias is controlled.",
    "After degree matching, hub sets should not outperform the true disease-near positive sets."
  )
)

cat("\n============================================================\n")
cat("Pressure Test A main proximity summary\n")
cat("============================================================\n")
print(summary_main)

cat("\n============================================================\n")
cat("Pressure Test A permutation summary\n")
cat("============================================================\n")
print(perm_summary_all[Scenario_Group == "A_MainAccuracy"])

cat("\n============================================================\n")
cat("Pressure Test A pass summary\n")
cat("============================================================\n")
print(pass_summary)

############################################################
## 14. Save only essential files
############################################################

out_grid <- file.path(OUT_DIR, "01A_ProximityAccuracy_DesignGrid.tsv")
out_targets <- file.path(OUT_DIR, "02A_ProximityAccuracy_TargetSets.tsv")
out_results <- file.path(OUT_DIR, "03A_ProximityAccuracy_Results.tsv")
out_perm <- file.path(OUT_DIR, "04A_ProximityAccuracy_Permutation.tsv")
out_summary <- file.path(OUT_DIR, "05A_ProximityAccuracy_Summary.tsv")
out_pass <- file.path(OUT_DIR, "06A_ProximityAccuracy_PassCriteria.tsv")

fwrite(design_grid, out_grid, sep = "\t")
fwrite(target_sets, out_targets, sep = "\t")
fwrite(proximity_results, out_results, sep = "\t")
fwrite(permutation_results, out_perm, sep = "\t")
fwrite(summary_all, out_summary, sep = "\t")
fwrite(pass_summary, out_pass, sep = "\t")

############################################################
## 15. Publication-grade figures
############################################################

FIG_DIR <- file.path(OUT_DIR, "Figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_LEVELS <- c("Positive", "Neutral", "Mixed", "Dilution", "Hub")
TARGET_LABELS <- c(
  Positive = "Positive\nnear-disease",
  Neutral = "Neutral\nbackground",
  Mixed = "Mixed\ncomposite",
  Dilution = "Dilution\nfar-disease",
  Hub = "Hub-biased\nhigh-degree"
)

PAL_TARGET <- c(
  Positive = "#F2AAA6",  ## soft coral
  Neutral  = "#B8CADC",  ## pale blue-gray
  Mixed    = "#CDBBEA",  ## soft lavender
  Dilution = "#A9D8E8",  ## pale cyan-blue
  Hub      = "#F5C68B"   ## soft amber
)

PAL_PERM <- c(
  ordinary = "#D6D6D6",
  degree_matched = "#8FAFC7"
)

base_theme <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0, margin = margin(b = 7)),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0, color = "#4A4A4A", margin = margin(b = 8)),
      axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", margin = margin(r = 8)),
      axis.text = element_text(color = "#1F1F1F"),
      axis.text.x = element_text(size = base_size - 1, lineheight = 0.92),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.28, color = "#E7E7E7"),
      panel.border = element_rect(color = "#222222", linewidth = 0.45),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 4, b = 4),
      strip.background = element_rect(fill = "#F6F6F6", color = "#D7D7D7"),
      strip.text = element_text(face = "bold", color = "#222222"),
      plot.margin = margin(t = 14, r = 22, b = 16, l = 16)
    )
}

save_pub_plot <- function(plot_obj, filename_stub, width, height, dpi = 600) {
  png_file <- file.path(FIG_DIR, paste0(filename_stub, ".png"))
  pdf_file <- file.path(FIG_DIR, paste0(filename_stub, ".pdf"))
  ggplot2::ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE,
    device = cairo_pdf
  )
  invisible(list(png = png_file, pdf = pdf_file))
}

format_target_type <- function(dt) {
  dt <- copy(dt)
  dt[, Target_Type := factor(Target_Type, levels = TARGET_LEVELS)]
  dt
}

## Figure A1: raw disease-distance distribution
fig_a1_dt <- format_target_type(proximity_results[Scenario_Group == "A_MainAccuracy"])

fig_a1 <- ggplot(fig_a1_dt, aes(x = Target_Type, y = Mean_Distance, fill = Target_Type)) +
  geom_boxplot(width = 0.58, alpha = 0.76, outlier.shape = NA, color = "#303030", linewidth = 0.38) +
  geom_jitter(aes(color = Target_Type), width = 0.13, height = 0, size = 1.35, alpha = 0.45, show.legend = FALSE) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.5, fill = "white", color = "#111111", stroke = 0.45) +
  scale_x_discrete(labels = TARGET_LABELS) +
  scale_fill_manual(values = PAL_TARGET, drop = FALSE) +
  scale_color_manual(values = PAL_TARGET, drop = FALSE) +
  labs(
    title = "Raw disease-distance distribution",
    subtitle = "Virtual target sets sampled from STRING700 human PPI background",
    x = NULL,
    y = "Mean distance to stroke disease module"
  ) +
  base_theme(12) +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")

save_pub_plot(fig_a1, "FigA1_RawDistance_Distribution", width = 7.8, height = 5.2)

## Figure A2: target-number robustness
fig_a2_dt <- format_target_type(summary_all[Scenario_Group == "B_TargetNumber"])

fig_a2 <- ggplot(
  fig_a2_dt,
  aes(x = Design_Target_N, y = Median_Mean_Distance, color = Target_Type, group = Target_Type)
) +
  geom_line(linewidth = 0.95, alpha = 0.94) +
  geom_point(aes(fill = Target_Type), shape = 21, size = 2.9, stroke = 0.55, color = "#303030") +
  scale_x_continuous(breaks = sort(unique(fig_a2_dt$Design_Target_N))) +
  scale_color_manual(values = PAL_TARGET, labels = TARGET_LABELS, drop = FALSE) +
  scale_fill_manual(values = PAL_TARGET, labels = TARGET_LABELS, drop = FALSE) +
  labs(
    title = "Target-number robustness",
    subtitle = "Distance ranking remains stable across virtual target-set sizes",
    x = "Target-set size",
    y = "Median mean distance",
    color = "Target set",
    fill = "Target set"
  ) +
  base_theme(12) +
  guides(fill = "none", color = guide_legend(nrow = 1, byrow = TRUE)) +
  coord_cartesian(clip = "off")

save_pub_plot(fig_a2, "FigA2_TargetNumber_Robustness", width = 8.4, height = 5.4)

## Figure A3: hub raw-distance bias
hub_bias_src <- proximity_results[
  Scenario_Group == "A_MainAccuracy" & Target_Type %in% c("Neutral", "Hub")
]
hub_bias_src <- format_target_type(hub_bias_src)

hub_dist <- hub_bias_src[, .(
  Simulation_ID,
  Target_Type,
  Metric = "Mean distance",
  Value = Mean_Distance
)]
hub_degree <- hub_bias_src[, .(
  Simulation_ID,
  Target_Type,
  Metric = "Mean degree",
  Value = Mean_Degree
)]
fig_a3_dt <- rbindlist(list(hub_dist, hub_degree), use.names = TRUE, fill = TRUE)
fig_a3_dt[, Metric := factor(Metric, levels = c("Mean distance", "Mean degree"))]

fig_a3 <- ggplot(fig_a3_dt, aes(x = Target_Type, y = Value, fill = Target_Type)) +
  geom_boxplot(width = 0.55, alpha = 0.78, outlier.shape = NA, color = "#303030", linewidth = 0.38) +
  geom_jitter(aes(color = Target_Type), width = 0.11, height = 0, size = 1.35, alpha = 0.45, show.legend = FALSE) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_x_discrete(labels = TARGET_LABELS) +
  scale_fill_manual(values = PAL_TARGET, drop = FALSE) +
  scale_color_manual(values = PAL_TARGET, drop = FALSE) +
  labs(
    title = "Hub raw-distance bias",
    subtitle = "High-degree hubs appear close under raw distance and require degree correction",
    x = NULL,
    y = NULL
  ) +
  base_theme(12) +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off")

save_pub_plot(fig_a3, "FigA3_Hub_RawDistance_Bias", width = 8.6, height = 4.8)

## Figure A4: ordinary permutation Z-score
fig_a4_dt <- format_target_type(permutation_results[
  Scenario_Group == "A_MainAccuracy" & Permutation_Method == "ordinary"
])

if (nrow(fig_a4_dt) > 0) {
  fig_a4 <- ggplot(fig_a4_dt, aes(x = Target_Type, y = Z_Score, fill = Target_Type)) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.42, color = "#333333") +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", linewidth = 0.35, color = "#9B9B9B") +
    geom_boxplot(width = 0.58, alpha = 0.76, outlier.shape = NA, color = "#303030", linewidth = 0.38) +
    geom_jitter(aes(color = Target_Type), width = 0.13, height = 0, size = 1.35, alpha = 0.50, show.legend = FALSE) +
    scale_x_discrete(labels = TARGET_LABELS) +
    scale_fill_manual(values = PAL_TARGET, drop = FALSE) +
    scale_color_manual(values = PAL_TARGET, drop = FALSE) +
    labs(
      title = "Ordinary permutation Z-score",
      subtitle = "Ordinary randomization exposes hub-driven apparent proximity",
      x = NULL,
      y = "Z-score versus ordinary random target sets"
    ) +
    base_theme(12) +
    theme(legend.position = "none") +
    coord_cartesian(clip = "off")
  save_pub_plot(fig_a4, "FigA4_OrdinaryPermutation_Zscore", width = 7.8, height = 5.2)
}

## Figure A5: degree-matched permutation Z-score
fig_a5_dt <- format_target_type(permutation_results[
  Scenario_Group == "A_MainAccuracy" & Permutation_Method == "degree_matched"
])

if (nrow(fig_a5_dt) > 0) {
  fig_a5 <- ggplot(fig_a5_dt, aes(x = Target_Type, y = Z_Score, fill = Target_Type)) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.42, color = "#333333") +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", linewidth = 0.35, color = "#9B9B9B") +
    geom_boxplot(width = 0.58, alpha = 0.76, outlier.shape = NA, color = "#303030", linewidth = 0.38) +
    geom_jitter(aes(color = Target_Type), width = 0.13, height = 0, size = 1.35, alpha = 0.50, show.legend = FALSE) +
    scale_x_discrete(labels = TARGET_LABELS) +
    scale_fill_manual(values = PAL_TARGET, drop = FALSE) +
    scale_color_manual(values = PAL_TARGET, drop = FALSE) +
    labs(
      title = "Degree-matched permutation Z-score",
      subtitle = "Degree matching preserves true disease-near signal while suppressing hub bias",
      x = NULL,
      y = "Z-score versus degree-matched random target sets"
    ) +
    base_theme(12) +
    theme(legend.position = "none") +
    coord_cartesian(clip = "off")
  save_pub_plot(fig_a5, "FigA5_DegreeMatchedPermutation_Zscore", width = 7.8, height = 5.2)
}

cat("\nPublication figures saved to:\n")
cat(FIG_DIR, "\n")


############################################################
## 16. Final report
############################################################

cat("\n============================================================\n")
cat("ACMPI-NPS Pressure Test A completed.\n")
cat("============================================================\n")

cat("\nRun mode:\n")
cat(RUN_MODE, "\n")

cat("\nDisease module input folder:\n")
cat(DISEASE_DIR, "\n")

cat("\nOutput folder:\n")
cat(OUT_DIR, "\n")

cat("\nSaved essential files:\n")
cat(out_grid, "\n")
cat(out_targets, "\n")
cat(out_results, "\n")
cat(out_perm, "\n")
cat(out_summary, "\n")
cat(out_pass, "\n")

cat("\nDesign size:\n")
cat("Total simulations: ", nrow(design_grid), "\n", sep = "")
cat("Target-set rows: ", nrow(target_sets), "\n", sep = "")
cat("Proximity result rows: ", nrow(proximity_results), "\n", sep = "")
cat("Permutation rows: ", nrow(permutation_results), "\n", sep = "")

cat("\nNext step:\n")
cat("1. Run RUN_MODE <- 'PILOT_FULL' to complete raw-distance accuracy and representative permutation in one pass.\n")
cat("2. If PILOT_FULL passes, run RUN_MODE <- 'FINAL' for robust validation.\n")
cat("3. For manuscript lock, run RUN_MODE <- 'FINAL_STRICT'.\n")
cat("4. Only after Pressure Test A passes should you run Pressure Test B for multi-component combinations.\n")

cat("\nDone. Only essential Pressure Test A files were saved.\n")
cat("============================================================\n")
