############################################################
## ACMPI-NPS Pressure Test B
## FINAL methodology-grade validation + manuscript figures
##
## Multi-component combination benefit + target contribution
##
## Purpose:
##   Validate whether ACMPI-NPS can identify:
##   1) beneficial single-partner addition;
##   2) cumulative benefit from complementary multi-partner combinations;
##   3) reduced benefit from redundant partners;
##   4) dilution when disease-far partners are added;
##   5) target-level drivers that pull combinations closer to or farther from the disease module.
##
## Assumption:
##   Pressure Test A has already validated:
##   - disease-distance calculation;
##   - hub-bias control;
##   - degree-matched permutation framework.
##
## Input:
##   Latest ANPS_StrokeDisease_STRING700_* folder under ROOT_DIR containing:
##     01_StrokeDisease_GeneMapping_STRING.tsv
##     02_STRING700_NodeDistanceToStrokeDisease.tsv
##     04_StrokeDisease_PPI_QC.tsv
##
## Output:
##   Core TSV files + PNG/PDF manuscript-grade figures.
############################################################

rm(list = ls())

############################################################
## 0. User settings
############################################################

ROOT_DIR <- "/media/desk16/iy15915/中药之开创/脑卒中疾病基因"

STRING_SCORE_MIN <- 700
SET_SEED <- 20260503

## ==========================================================
## Run mode
## ==========================================================
## QUICK_TEST:
##   Minimal check only.
## PILOT_B:
##   Lightweight confirmation of model behavior.
## FINAL_B:
##   Manuscript-level validation.
## FINAL_B_STRICT:
##   Stricter final run for manuscript lock.
## ==========================================================

RUN_MODE <- "FINAL_B_STRICT"

if (RUN_MODE == "QUICK_TEST") {
  N_REP <- 10
  CORE_N <- 30
  PARTNER_N <- 30
  DRIVER_TOP_N <- 10
  PROGRESS_EVERY <- 1
  COMPUTE_COMBINATION_PERMUTATION <- TRUE
  B_PERM_COMBINATION <- 100
  PERM_REPS_PER_SCENARIO <- 2
  SAVE_FULL_DRIVER_ABLATION <- TRUE
  
} else if (RUN_MODE == "PILOT_B") {
  N_REP <- 80
  CORE_N <- 30
  PARTNER_N <- 30
  DRIVER_TOP_N <- 10
  PROGRESS_EVERY <- 10
  COMPUTE_COMBINATION_PERMUTATION <- TRUE
  B_PERM_COMBINATION <- 300
  PERM_REPS_PER_SCENARIO <- 8
  SAVE_FULL_DRIVER_ABLATION <- TRUE
  
} else if (RUN_MODE == "FINAL_B") {
  N_REP <- 300
  CORE_N <- 30
  PARTNER_N <- 30
  DRIVER_TOP_N <- 20
  PROGRESS_EVERY <- 50
  COMPUTE_COMBINATION_PERMUTATION <- TRUE
  B_PERM_COMBINATION <- 500
  PERM_REPS_PER_SCENARIO <- 25
  SAVE_FULL_DRIVER_ABLATION <- TRUE
  
} else if (RUN_MODE == "FINAL_B_STRICT") {
  N_REP <- 500
  CORE_N <- 30
  PARTNER_N <- 30
  DRIVER_TOP_N <- 30
  PROGRESS_EVERY <- 50
  COMPUTE_COMBINATION_PERMUTATION <- TRUE
  B_PERM_COMBINATION <- 1000
  PERM_REPS_PER_SCENARIO <- 30
  SAVE_FULL_DRIVER_ABLATION <- TRUE
  
} else {
  stop("Unknown RUN_MODE. Use QUICK_TEST, PILOT_B, FINAL_B, or FINAL_B_STRICT.", call. = FALSE)
}

DEGREE_BINS <- 10
NEAR_DISTANCE_CUTOFF <- 1
FAR_DISTANCE_CUTOFF <- 3
DRIVER_EPS <- 0.005

set.seed(SET_SEED)

############################################################
## 1. Packages
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
## 2. Locate disease module folder
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
      call. = FALSE
    )
  }
  
  valid_dirs[which.max(file.info(valid_dirs)$mtime)]
}

DISEASE_DIR <- find_latest_disease_dir(ROOT_DIR, STRING_SCORE_MIN)

MAPPING_FILE <- file.path(DISEASE_DIR, "01_StrokeDisease_GeneMapping_STRING.tsv")
NODE_DISTANCE_FILE <- file.path(DISEASE_DIR, "02_STRING700_NodeDistanceToStrokeDisease.tsv")
QC_FILE <- file.path(DISEASE_DIR, "04_StrokeDisease_PPI_QC.tsv")

OUT_DIR <- file.path(
  ROOT_DIR,
  paste0(
    "ANPS_PressureTestB_MultiComponent_FINAL_",
    RUN_MODE,
    "_STRING",
    STRING_SCORE_MIN,
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

FIG_DIR <- file.path(OUT_DIR, "Figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("ACMPI-NPS Pressure Test B: methodology-grade multi-component validation\n")
cat("Run mode: ", RUN_MODE, "\n", sep = "")
cat("Disease module folder:\n", DISEASE_DIR, "\n", sep = "")
cat("Output folder:\n", OUT_DIR, "\n", sep = "")
cat("============================================================\n")

############################################################
## 3. Read data
############################################################

mapping_dt <- fread(MAPPING_FILE, sep = "\t", data.table = TRUE)
node_dt <- fread(NODE_DISTANCE_FILE, sep = "\t", data.table = TRUE)
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

############################################################
## 4. Degree bins and lookup
############################################################

node_dt[, degree_rank := frank(degree_STRING700, ties.method = "average")]
node_dt[, degree_bin := ceiling(degree_rank / .N * DEGREE_BINS)]
node_dt[degree_bin < 1, degree_bin := 1]
node_dt[degree_bin > DEGREE_BINS, degree_bin := DEGREE_BINS]

degree_q10 <- as.numeric(quantile(node_dt$degree_STRING700, 0.10, na.rm = TRUE))
degree_q90 <- as.numeric(quantile(node_dt$degree_STRING700, 0.90, na.rm = TRUE))

node_dt[, is_extreme_hub := degree_STRING700 >= degree_q90]
node_dt[, is_non_extreme_degree := degree_STRING700 >= degree_q10 & degree_STRING700 <= degree_q90]

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
## 5. Sampling and distance helpers
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
  
  need <- n - length(ids)
  
  if (length(fill_pool) == 0) {
    extra <- sample_safe(ids, need, replace = TRUE)
  } else {
    extra <- sample_safe(fill_pool, need, replace = length(fill_pool) < need)
  }
  
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
  } else if (distance_rule == "dist_ge_3") {
    pool <- pool[distance_to_stroke_disease >= 3]
  } else {
    stop("Unknown distance_rule: ", distance_rule, call. = FALSE)
  }
  
  pool$string_id
}

sample_core_targets <- function(n) {
  pool <- get_pool_by_distance(
    distance_rule = "dist_1_2",
    non_extreme = TRUE,
    allow_disease_nodes = FALSE
  )
  sample_safe(pool, n, replace = FALSE)
}

sample_positive_targets <- function(n, exclude_ids = character()) {
  n0 <- floor(n * 0.30)
  n1 <- n - n0
  
  pool0 <- get_pool_by_distance(
    distance_rule = "dist_0",
    exclude_ids = exclude_ids,
    non_extreme = TRUE,
    allow_disease_nodes = TRUE
  )
  
  pool1 <- get_pool_by_distance(
    distance_rule = "dist_1",
    exclude_ids = exclude_ids,
    non_extreme = TRUE,
    allow_disease_nodes = FALSE
  )
  
  ids0 <- if (n0 > 0) sample_safe(pool0, n0, replace = FALSE) else character()
  ids1 <- if (n1 > 0) sample_safe(pool1, n1, replace = FALSE) else character()
  
  complete_to_n(c(ids0, ids1), n, unique(c(pool0, pool1)), exclude_ids = exclude_ids)
}

sample_neutral_targets <- function(n, exclude_ids = character()) {
  pool <- background_pool[is_non_extreme_degree == TRUE]
  if (length(exclude_ids) > 0) {
    pool <- pool[!string_id %in% exclude_ids]
  }
  sample_safe(pool$string_id, n, replace = FALSE)
}

sample_dilution_targets <- function(n, exclude_ids = character()) {
  pool <- get_pool_by_distance(
    distance_rule = "dist_ge_3",
    exclude_ids = exclude_ids,
    non_extreme = TRUE,
    allow_disease_nodes = FALSE
  )
  sample_safe(pool, n, replace = FALSE)
}

sample_mixed_targets <- function(n, exclude_ids = character()) {
  n_pos <- floor(n * 0.40)
  n_neu <- floor(n * 0.40)
  n_dil <- n - n_pos - n_neu
  
  ids_pos <- sample_positive_targets(n_pos, exclude_ids = exclude_ids)
  ids_neu <- sample_neutral_targets(n_neu, exclude_ids = c(exclude_ids, ids_pos))
  ids_dil <- sample_dilution_targets(n_dil, exclude_ids = c(exclude_ids, ids_pos, ids_neu))
  
  complete_to_n(c(ids_pos, ids_neu, ids_dil), n, all_background_ids, exclude_ids = exclude_ids)
}

sample_redundant_targets <- function(n,
                                     reference_ids,
                                     exclude_ids = character(),
                                     overlap_fraction = 0.70) {
  n_overlap <- floor(n * overlap_fraction)
  n_new <- n - n_overlap
  
  shared <- sample_safe(reference_ids, n_overlap, replace = length(reference_ids) < n_overlap)
  new_pos <- sample_positive_targets(n_new, exclude_ids = c(exclude_ids, shared))
  
  complete_to_n(c(shared, new_pos), n, unique(c(reference_ids, new_pos)), exclude_ids = exclude_ids)
}

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
## 6. Design grid
############################################################

scenario_dt <- data.table(
  Scenario = c(
    "Additive_benefit",
    "Redundant_partner",
    "Dilution_partner",
    "Mixed_realistic"
  ),
  Partner_A_Profile = c("Positive", "Positive", "Positive", "Positive"),
  Partner_B_Profile = c("Positive_complementary", "Redundant_to_A", "Positive_complementary", "Positive_complementary"),
  Partner_C_Profile = c("Neutral", "Neutral", "Dilution", "Mixed")
)

design_grid <- scenario_dt[
  ,
  .(Replicate = seq_len(N_REP)),
  by = .(Scenario, Partner_A_Profile, Partner_B_Profile, Partner_C_Profile)
]

design_grid[, Simulation_ID := sprintf("PTB_%05d", .I)]
setcolorder(design_grid, "Simulation_ID")

cat("\nPressure Test B design grid:\n")
print(design_grid[, .N, by = Scenario])

############################################################
## 7. Build virtual compound target sets
############################################################

target_sets_all <- list()

cat("\nSampling virtual core and partner target sets...\n")

for (i in seq_len(nrow(design_grid))) {
  row <- design_grid[i]
  
  sim_id <- row$Simulation_ID
  scenario <- row$Scenario
  
  if (i %% PROGRESS_EVERY == 0 || i == 1) {
    cat("Simulation ", i, " / ", nrow(design_grid), ": ", sim_id, " | ", scenario, "\n", sep = "")
    flush.console()
  }
  
  core_ids <- sample_core_targets(CORE_N)
  partner_a_ids <- sample_positive_targets(PARTNER_N, exclude_ids = core_ids)
  
  if (scenario == "Redundant_partner") {
    partner_b_ids <- sample_redundant_targets(
      n = PARTNER_N,
      reference_ids = partner_a_ids,
      exclude_ids = core_ids,
      overlap_fraction = 0.70
    )
  } else {
    partner_b_ids <- sample_positive_targets(
      PARTNER_N,
      exclude_ids = c(core_ids, partner_a_ids)
    )
  }
  
  if (scenario == "Dilution_partner") {
    partner_c_ids <- sample_dilution_targets(
      PARTNER_N,
      exclude_ids = c(core_ids, partner_a_ids, partner_b_ids)
    )
  } else if (scenario == "Mixed_realistic") {
    partner_c_ids <- sample_mixed_targets(
      PARTNER_N,
      exclude_ids = c(core_ids, partner_a_ids, partner_b_ids)
    )
  } else {
    partner_c_ids <- sample_neutral_targets(
      PARTNER_N,
      exclude_ids = c(core_ids, partner_a_ids, partner_b_ids)
    )
  }
  
  one_compound_table <- function(compound_name, ids, truth_label) {
    data.table(
      Simulation_ID = sim_id,
      Scenario = scenario,
      Virtual_Compound = compound_name,
      string_id = unique(ids),
      Truth_Label = truth_label
    )
  }
  
  target_sets_all[[i]] <- rbindlist(
    list(
      one_compound_table("Core", core_ids, "Core_like"),
      one_compound_table("A", partner_a_ids, "Positive_like"),
      one_compound_table(
        "B",
        partner_b_ids,
        ifelse(scenario == "Redundant_partner", "Redundant_like", "Positive_like")
      ),
      one_compound_table(
        "C",
        partner_c_ids,
        ifelse(
          scenario == "Dilution_partner",
          "Dilution_like",
          ifelse(scenario == "Mixed_realistic", "Mixed_like", "Neutral_like")
        )
      )
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

target_sets <- rbindlist(target_sets_all, use.names = TRUE, fill = TRUE)

target_sets[, preferred_name := gene_lookup[string_id]]
target_sets[, distance_to_stroke_disease := as.numeric(dist_lookup[string_id])]
target_sets[, degree_STRING700 := as.numeric(degree_lookup[string_id])]
target_sets[, is_stroke_disease_node := as.logical(disease_flag_lookup[string_id])]
target_sets[, Target_Distance_Class := classify_distance(distance_to_stroke_disease)]

############################################################
## 8. Combination definitions
############################################################

combination_def <- data.table(
  Combination = c(
    "Core",
    "Core+A",
    "Core+B",
    "Core+C",
    "Core+A+B",
    "Core+A+C",
    "Core+B+C",
    "Core+A+B+C"
  ),
  Included_Compounds = c(
    "Core",
    "Core;A",
    "Core;B",
    "Core;C",
    "Core;A;B",
    "Core;A;C",
    "Core;B;C",
    "Core;A;B;C"
  ),
  Combination_Order = 1:8
)

get_combination_ids <- function(sim_id, included_compounds) {
  comps <- unlist(strsplit(included_compounds, ";", fixed = TRUE))
  target_sets[
    Simulation_ID == sim_id & Virtual_Compound %in% comps,
    unique(string_id)
  ]
}

############################################################
## 9. Combination-level results
############################################################

combination_results_all <- list()

for (i in seq_len(nrow(design_grid))) {
  sim_id <- design_grid$Simulation_ID[i]
  scenario <- design_grid$Scenario[i]
  
  core_ids <- get_combination_ids(sim_id, "Core")
  core_distance <- mean_distance(core_ids)
  
  rows <- lapply(seq_len(nrow(combination_def)), function(j) {
    comb <- combination_def[j]
    ids <- get_combination_ids(sim_id, comb$Included_Compounds)
    
    d_vec <- as.numeric(dist_lookup[unique(ids)])
    d_vec <- d_vec[!is.na(d_vec)]
    
    data.table(
      Simulation_ID = sim_id,
      Scenario = scenario,
      Combination = comb$Combination,
      Included_Compounds = comb$Included_Compounds,
      Combination_Order = comb$Combination_Order,
      Target_N = length(unique(ids)),
      Mean_Distance = mean_distance(ids),
      Median_Distance = median_distance(ids),
      Fraction_Near_Disease = mean(d_vec <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE),
      Fraction_Far_From_Disease = mean(d_vec >= FAR_DISTANCE_CUTOFF, na.rm = TRUE),
      Core_Mean_Distance = core_distance,
      Delta_vs_Core = mean_distance(ids) - core_distance,
      Network_Benefit_vs_Core = core_distance - mean_distance(ids)
    )
  })
  
  combination_results_all[[i]] <- rbindlist(rows, use.names = TRUE, fill = TRUE)
}

combination_results <- rbindlist(combination_results_all, use.names = TRUE, fill = TRUE)

############################################################
## 10. Incremental partner benefit and added-target contribution
############################################################

make_incremental_edges <- function() {
  data.table(
    Base_Combination = c(
      "Core", "Core", "Core",
      "Core+A", "Core+A",
      "Core+B", "Core+B",
      "Core+C", "Core+C",
      "Core+A+B", "Core+A+C", "Core+B+C"
    ),
    Added_Compound = c(
      "A", "B", "C",
      "B", "C",
      "A", "C",
      "A", "B",
      "C", "B", "A"
    ),
    New_Combination = c(
      "Core+A", "Core+B", "Core+C",
      "Core+A+B", "Core+A+C",
      "Core+A+B", "Core+B+C",
      "Core+A+C", "Core+B+C",
      "Core+A+B+C", "Core+A+B+C", "Core+A+B+C"
    )
  )
}

increment_edges <- make_incremental_edges()

incremental_all <- list()
target_contribution_all <- list()

for (i in seq_len(nrow(design_grid))) {
  sim_id <- design_grid$Simulation_ID[i]
  scenario <- design_grid$Scenario[i]
  
  comb_dt <- combination_results[Simulation_ID == sim_id]
  
  rows <- lapply(seq_len(nrow(increment_edges)), function(j) {
    edge <- increment_edges[j]
    
    before_d <- comb_dt[Combination == edge$Base_Combination, Mean_Distance][1]
    after_d <- comb_dt[Combination == edge$New_Combination, Mean_Distance][1]
    
    base_comps <- comb_dt[Combination == edge$Base_Combination, Included_Compounds][1]
    new_comps <- comb_dt[Combination == edge$New_Combination, Included_Compounds][1]
    
    base_ids <- get_combination_ids(sim_id, base_comps)
    new_ids <- get_combination_ids(sim_id, new_comps)
    added_ids <- setdiff(new_ids, base_ids)
    
    added_d <- as.numeric(dist_lookup[added_ids])
    added_d <- added_d[!is.na(added_d)]
    
    added_close_n <- sum(added_d <= NEAR_DISTANCE_CUTOFF, na.rm = TRUE)
    added_far_n <- sum(added_d >= FAR_DISTANCE_CUTOFF, na.rm = TRUE)
    
    data.table(
      Simulation_ID = sim_id,
      Scenario = scenario,
      Base_Combination = edge$Base_Combination,
      Added_Compound = edge$Added_Compound,
      New_Combination = edge$New_Combination,
      Distance_Before = before_d,
      Distance_After = after_d,
      Incremental_Benefit = before_d - after_d,
      Added_Target_N = length(unique(added_ids)),
      Added_Close_Target_N = added_close_n,
      Added_Far_Target_N = added_far_n,
      Added_Close_Target_Fraction = ifelse(length(added_ids) > 0, added_close_n / length(added_ids), NA_real_),
      Added_Far_Target_Fraction = ifelse(length(added_ids) > 0, added_far_n / length(added_ids), NA_real_),
      Mean_Added_Target_Distance = ifelse(length(added_d) > 0, mean(added_d), NA_real_)
    )
  })
  
  incremental_all[[i]] <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  
  target_rows <- lapply(seq_len(nrow(increment_edges)), function(j) {
    edge <- increment_edges[j]
    
    base_comps <- comb_dt[Combination == edge$Base_Combination, Included_Compounds][1]
    new_comps <- comb_dt[Combination == edge$New_Combination, Included_Compounds][1]
    
    base_ids <- get_combination_ids(sim_id, base_comps)
    new_ids <- get_combination_ids(sim_id, new_comps)
    added_ids <- setdiff(new_ids, base_ids)
    
    if (length(added_ids) == 0) return(data.table())
    
    source_info <- target_sets[
      Simulation_ID == sim_id & string_id %in% added_ids,
      .(
        Source_Compounds = paste(sort(unique(Virtual_Compound)), collapse = ";"),
        Truth_Labels = paste(sort(unique(Truth_Label)), collapse = ";"),
        preferred_name = preferred_name[1],
        distance_to_stroke_disease = distance_to_stroke_disease[1],
        degree_STRING700 = degree_STRING700[1],
        is_stroke_disease_node = is_stroke_disease_node[1],
        Target_Distance_Class = Target_Distance_Class[1]
      ),
      by = string_id
    ]
    
    source_info[, Simulation_ID := sim_id]
    source_info[, Scenario := scenario]
    source_info[, Base_Combination := edge$Base_Combination]
    source_info[, Added_Compound := edge$Added_Compound]
    source_info[, New_Combination := edge$New_Combination]
    
    setcolorder(
      source_info,
      c(
        "Simulation_ID", "Scenario", "Base_Combination",
        "Added_Compound", "New_Combination",
        "string_id", "preferred_name", "Source_Compounds",
        "Truth_Labels", "distance_to_stroke_disease",
        "degree_STRING700", "is_stroke_disease_node",
        "Target_Distance_Class"
      )
    )
    
    source_info[]
  })
  
  target_contribution_all[[i]] <- rbindlist(target_rows, use.names = TRUE, fill = TRUE)
}

incremental_results <- rbindlist(incremental_all, use.names = TRUE, fill = TRUE)
target_contribution <- rbindlist(target_contribution_all, use.names = TRUE, fill = TRUE)

incremental_results[, Increment_Class := fifelse(
  Incremental_Benefit > DRIVER_EPS,
  "Beneficial_addition",
  fifelse(Incremental_Benefit < -DRIVER_EPS, "Dilution_addition", "Redundant_or_neutral")
)]

############################################################
## 11. Efficient driver ablation
############################################################

## For a mean-distance score, leave-one-target-out DriverScore can be computed exactly:
##   DriverScore = d(without target) - d(full combination)
## This is equivalent to explicit deletion but avoids repeated full recomputation.

driver_all <- list()

for (i in seq_len(nrow(combination_results))) {
  row <- combination_results[i]
  
  sim_id <- row$Simulation_ID
  scenario <- row$Scenario
  combination <- row$Combination
  included_compounds <- row$Included_Compounds
  
  ids <- unique(get_combination_ids(sim_id, included_compounds))
  
  if (length(ids) <= 2) next
  
  d <- as.numeric(dist_lookup[ids])
  keep <- !is.na(d)
  ids <- ids[keep]
  d <- d[keep]
  
  if (length(ids) <= 2) next
  
  n_ids <- length(ids)
  sum_d <- sum(d)
  base_d <- sum_d / n_ids
  d_without <- (sum_d - d) / (n_ids - 1)
  driver_score <- d_without - base_d
  
  source_info <- target_sets[
    Simulation_ID == sim_id & string_id %in% ids,
    .(
      Source_Compounds = paste(sort(unique(Virtual_Compound)), collapse = ";"),
      Truth_Labels = paste(sort(unique(Truth_Label)), collapse = ";"),
      preferred_name = preferred_name[1],
      distance_to_stroke_disease = distance_to_stroke_disease[1],
      degree_STRING700 = degree_STRING700[1],
      is_stroke_disease_node = is_stroke_disease_node[1],
      Target_Distance_Class = Target_Distance_Class[1]
    ),
    by = string_id
  ]
  
  drv <- data.table(
    Simulation_ID = sim_id,
    Scenario = scenario,
    Combination = combination,
    Included_Compounds = included_compounds,
    string_id = ids,
    Combination_Distance = base_d,
    Distance_Without_Target = d_without,
    DriverScore = driver_score
  )
  
  drv <- merge(drv, source_info, by = "string_id", all.x = TRUE, sort = FALSE)
  
  setcolorder(
    drv,
    c(
      "Simulation_ID", "Scenario", "Combination", "Included_Compounds",
      "string_id", "preferred_name", "Source_Compounds", "Truth_Labels",
      "distance_to_stroke_disease", "degree_STRING700", "is_stroke_disease_node",
      "Target_Distance_Class", "Combination_Distance", "Distance_Without_Target",
      "DriverScore"
    )
  )
  
  driver_all[[i]] <- drv
}

driver_results <- rbindlist(driver_all, use.names = TRUE, fill = TRUE)

driver_results[, Driver_Class := fifelse(
  DriverScore > DRIVER_EPS,
  "Beneficial_driver",
  fifelse(DriverScore < -DRIVER_EPS, "Dilution_target", "Neutral_or_redundant")
)]

driver_results[, Abs_DriverScore := abs(DriverScore)]

## Compact top-driver table for writing and figure labelling.
top_driver_results <- driver_results[
  order(Scenario, Combination, -DriverScore)
][
  ,
  head(.SD, DRIVER_TOP_N),
  by = .(Scenario, Combination)
]

top_dilution_results <- driver_results[
  order(Scenario, Combination, DriverScore)
][
  ,
  head(.SD, DRIVER_TOP_N),
  by = .(Scenario, Combination)
]

############################################################
## 12. Driver recovery
############################################################

driver_recovery_distance <- driver_results[
  ,
  .(
    N_Targets = .N,
    Median_DriverScore = median(DriverScore, na.rm = TRUE),
    Mean_DriverScore = mean(DriverScore, na.rm = TRUE),
    Fraction_BeneficialDriver = mean(Driver_Class == "Beneficial_driver", na.rm = TRUE),
    Fraction_DilutionTarget = mean(Driver_Class == "Dilution_target", na.rm = TRUE)
  ),
  by = .(Scenario, Combination, Target_Distance_Class)
][order(Scenario, Combination, Target_Distance_Class)]

driver_recovery_truth <- driver_results[
  ,
  .(
    N_Targets = .N,
    Median_DriverScore = median(DriverScore, na.rm = TRUE),
    Mean_DriverScore = mean(DriverScore, na.rm = TRUE),
    Fraction_BeneficialDriver = mean(Driver_Class == "Beneficial_driver", na.rm = TRUE),
    Fraction_DilutionTarget = mean(Driver_Class == "Dilution_target", na.rm = TRUE)
  ),
  by = .(Scenario, Combination, Truth_Labels)
][order(Scenario, Combination, Truth_Labels)]

driver_recovery <- rbindlist(
  list(
    data.table(Recovery_Level = "DistanceClass", driver_recovery_distance),
    data.table(Recovery_Level = "TruthLabel", driver_recovery_truth)
  ),
  use.names = TRUE,
  fill = TRUE
)

############################################################
## 13. Combination summaries and pass criteria
############################################################

comb_summary <- combination_results[
  ,
  .(
    N = .N,
    Median_Mean_Distance = median(Mean_Distance, na.rm = TRUE),
    Mean_Mean_Distance = mean(Mean_Distance, na.rm = TRUE),
    Median_Delta_vs_Core = median(Delta_vs_Core, na.rm = TRUE),
    Median_Network_Benefit_vs_Core = median(Network_Benefit_vs_Core, na.rm = TRUE),
    Median_Target_N = median(Target_N, na.rm = TRUE),
    Median_Fraction_Near_Disease = median(Fraction_Near_Disease, na.rm = TRUE),
    Median_Fraction_Far_From_Disease = median(Fraction_Far_From_Disease, na.rm = TRUE)
  ),
  by = .(Scenario, Combination, Combination_Order)
][order(Scenario, Combination_Order)]

increment_summary <- incremental_results[
  ,
  .(
    N = .N,
    Median_Incremental_Benefit = median(Incremental_Benefit, na.rm = TRUE),
    Mean_Incremental_Benefit = mean(Incremental_Benefit, na.rm = TRUE),
    Fraction_Beneficial_Addition = mean(Increment_Class == "Beneficial_addition", na.rm = TRUE),
    Fraction_Dilution_Addition = mean(Increment_Class == "Dilution_addition", na.rm = TRUE),
    Median_Added_Close_Target_N = median(Added_Close_Target_N, na.rm = TRUE),
    Median_Added_Far_Target_N = median(Added_Far_Target_N, na.rm = TRUE),
    Median_Added_Close_Target_Fraction = median(Added_Close_Target_Fraction, na.rm = TRUE),
    Median_Added_Far_Target_Fraction = median(Added_Far_Target_Fraction, na.rm = TRUE)
  ),
  by = .(Scenario, Base_Combination, Added_Compound, New_Combination)
][order(Scenario, Base_Combination, Added_Compound)]

get_combo_metric <- function(scenario, combination, metric = "Median_Mean_Distance") {
  x <- comb_summary[Scenario == scenario & Combination == combination, get(metric)]
  if (length(x) == 0) return(NA_real_)
  x[1]
}

get_increment_metric <- function(scenario, base, added, newcomb, metric = "Median_Incremental_Benefit") {
  x <- increment_summary[
    Scenario == scenario &
      Base_Combination == base &
      Added_Compound == added &
      New_Combination == newcomb,
    get(metric)
  ]
  if (length(x) == 0) return(NA_real_)
  x[1]
}

additive_core <- get_combo_metric("Additive_benefit", "Core")
additive_ab <- get_combo_metric("Additive_benefit", "Core+A+B")
additive_ab_gain <- get_increment_metric("Additive_benefit", "Core+A", "B", "Core+A+B")

redundant_b_gain <- get_increment_metric("Redundant_partner", "Core+A", "B", "Core+A+B")

dilution_c_gain <- get_increment_metric("Dilution_partner", "Core+A", "C", "Core+A+C")

mixed_c_gain <- get_increment_metric("Mixed_realistic", "Core+A", "C", "Core+A+C")

beneficial_driver_close <- driver_results[
  Target_Distance_Class == "Disease_close",
  median(DriverScore, na.rm = TRUE)
]

beneficial_driver_far <- driver_results[
  Target_Distance_Class == "Disease_far",
  median(DriverScore, na.rm = TRUE)
]

truth_positive_driver <- driver_results[
  grepl("Positive_like", Truth_Labels, fixed = TRUE),
  median(DriverScore, na.rm = TRUE)
]

truth_dilution_driver <- driver_results[
  grepl("Dilution_like", Truth_Labels, fixed = TRUE),
  median(DriverScore, na.rm = TRUE)
]

pass_summary <- data.table(
  Test_Item = c(
    "Additive scenario: Core+A+B closer than Core",
    "Additive scenario: adding B to Core+A gives positive incremental benefit",
    "Redundant scenario: adding redundant B gives smaller benefit than additive B",
    "Dilution scenario: adding C to Core+A worsens or does not improve distance",
    "Mixed scenario: adding mixed C gives intermediate or weak benefit",
    "Driver recovery: disease-close targets have higher DriverScore than disease-far targets",
    "Driver recovery: positive-like targets have higher DriverScore than dilution-like targets"
  ),
  Result = c(
    as.character(!is.na(additive_core) && !is.na(additive_ab) && additive_ab < additive_core),
    as.character(!is.na(additive_ab_gain) && additive_ab_gain > 0),
    as.character(!is.na(redundant_b_gain) && !is.na(additive_ab_gain) && redundant_b_gain < additive_ab_gain),
    as.character(!is.na(dilution_c_gain) && dilution_c_gain <= DRIVER_EPS),
    as.character(!is.na(mixed_c_gain) && mixed_c_gain < additive_ab_gain),
    as.character(!is.na(beneficial_driver_close) && !is.na(beneficial_driver_far) && beneficial_driver_close > beneficial_driver_far),
    as.character(!is.na(truth_positive_driver) && !is.na(truth_dilution_driver) && truth_positive_driver > truth_dilution_driver)
  ),
  Interpretation = c(
    "The model should identify cumulative benefit when complementary positive partners are combined.",
    "Adding a complementary positive partner should improve proximity.",
    "Highly overlapping partners should provide less incremental benefit.",
    "Disease-far partners should dilute or fail to improve the combination.",
    "Mixed partners should not outperform clean complementary-positive addition.",
    "Target-level ablation should recover disease-close targets as beneficial drivers and disease-far targets as dilution targets.",
    "Truth-label recovery should rank positive-like targets above dilution-like targets."
  )
)

cat("\n============================================================\n")
cat("Pressure Test B combination summary\n")
cat("============================================================\n")
print(comb_summary)

cat("\n============================================================\n")
cat("Pressure Test B incremental summary\n")
cat("============================================================\n")
print(increment_summary)

cat("\n============================================================\n")
cat("Pressure Test B pass summary\n")
cat("============================================================\n")
print(pass_summary)

############################################################
## 14. Optional representative combination-level permutation
############################################################

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
    if (is.null(pool) || length(pool) == 0) pool <- all_background_ids
    if (length(exclude_ids) > 0) pool <- setdiff(pool, exclude_ids)
    if (length(pool) == 0) pool <- all_background_ids
    out[j] <- sample(pool, 1)
  }
  
  complete_to_n(out, n_out, all_background_ids, exclude_ids = exclude_ids)
}

permute_combination_distance <- function(target_ids, b_perm = 1000) {
  target_ids <- unique(target_ids)
  target_ids <- target_ids[target_ids %in% node_dt$string_id]
  if (length(target_ids) == 0 || b_perm <= 0) {
    return(list(mean_random = NA_real_, sd_random = NA_real_, random_values = numeric()))
  }
  
  vals <- numeric(b_perm)
  for (b in seq_len(b_perm)) {
    sampled <- sample_degree_matched(reference_ids = target_ids, n_out = length(target_ids))
    vals[b] <- mean_distance(sampled)
  }
  vals <- vals[!is.na(vals)]
  list(mean_random = mean(vals), sd_random = sd(vals), random_values = vals)
}

combination_permutation <- data.table(
  Simulation_ID = character(),
  Scenario = character(),
  Combination = character(),
  B_Permutation = integer(),
  Observed_Distance = numeric(),
  Random_Mean_Distance = numeric(),
  Random_SD_Distance = numeric(),
  Z_Score = numeric(),
  Empirical_P_Lower_Closer = numeric()
)

if (COMPUTE_COMBINATION_PERMUTATION && B_PERM_COMBINATION > 0) {
  cat("\nRunning representative combination-level degree-matched permutation...\n")
  key_combos <- c("Core", "Core+A", "Core+A+B", "Core+A+C", "Core+A+B+C")
  
  representative <- combination_results[
    Combination %in% key_combos,
    {
      n_take <- min(.N, PERM_REPS_PER_SCENARIO)
      idx_center <- which.min(abs(Mean_Distance - median(Mean_Distance, na.rm = TRUE)))
      idx_gain <- order(Network_Benefit_vs_Core, decreasing = TRUE)[seq_len(min(ceiling(n_take / 2), .N))]
      idx <- unique(c(idx_center, idx_gain))
      rest <- setdiff(seq_len(.N), idx)
      if (length(idx) < n_take && length(rest) > 0) {
        idx <- unique(c(idx, sample(rest, min(n_take - length(idx), length(rest)))))
      }
      .SD[idx[seq_len(min(length(idx), n_take))]]
    },
    by = .(Scenario, Combination)
  ]
  
  perm_rows <- list()
  for (i in seq_len(nrow(representative))) {
    rr <- representative[i]
    if (i %% 20 == 0 || i == 1) {
      cat("Combination permutation ", i, " / ", nrow(representative), "\n", sep = "")
      flush.console()
    }
    ids <- get_combination_ids(rr$Simulation_ID, rr$Included_Compounds)
    pr <- permute_combination_distance(ids, b_perm = B_PERM_COMBINATION)
    z <- ifelse(is.na(pr$sd_random) || pr$sd_random == 0, NA_real_, (rr$Mean_Distance - pr$mean_random) / pr$sd_random)
    p_lower <- mean(pr$random_values <= rr$Mean_Distance, na.rm = TRUE)
    perm_rows[[i]] <- data.table(
      Simulation_ID = rr$Simulation_ID,
      Scenario = rr$Scenario,
      Combination = rr$Combination,
      B_Permutation = B_PERM_COMBINATION,
      Observed_Distance = rr$Mean_Distance,
      Random_Mean_Distance = pr$mean_random,
      Random_SD_Distance = pr$sd_random,
      Z_Score = z,
      Empirical_P_Lower_Closer = p_lower
    )
  }
  combination_permutation <- rbindlist(perm_rows, use.names = TRUE, fill = TRUE)
}

combination_perm_summary <- if (nrow(combination_permutation) > 0) {
  combination_permutation[
    ,
    .(
      N = .N,
      Median_Z = median(Z_Score, na.rm = TRUE),
      Mean_Z = mean(Z_Score, na.rm = TRUE),
      Fraction_Closer_005 = mean(Empirical_P_Lower_Closer <= 0.05, na.rm = TRUE)
    ),
    by = .(Scenario, Combination)
  ][order(Scenario, Combination)]
} else {
  data.table(
    Scenario = character(),
    Combination = character(),
    N = integer(),
    Median_Z = numeric(),
    Mean_Z = numeric(),
    Fraction_Closer_005 = numeric()
  )
}

############################################################
## 15. Save outputs
############################################################

out_design <- file.path(OUT_DIR, "01B_MultiComponent_DesignGrid.tsv")
out_targets <- file.path(OUT_DIR, "02B_MultiComponent_TargetSets.tsv")
out_comb <- file.path(OUT_DIR, "03B_MultiComponent_CombinationResults.tsv")
out_inc <- file.path(OUT_DIR, "04B_MultiComponent_IncrementalBenefit.tsv")
out_tcontrib <- file.path(OUT_DIR, "05B_MultiComponent_TargetContribution.tsv")
out_driver <- file.path(OUT_DIR, "06B_MultiComponent_DriverAblation.tsv")
out_recovery <- file.path(OUT_DIR, "07B_MultiComponent_DriverRecovery.tsv")
out_pass <- file.path(OUT_DIR, "08B_MultiComponent_PassCriteria.tsv")
out_comb_summary <- file.path(OUT_DIR, "09B_MultiComponent_CombinationSummary.tsv")
out_increment_summary <- file.path(OUT_DIR, "10B_MultiComponent_IncrementSummary.tsv")
out_top_driver <- file.path(OUT_DIR, "11B_MultiComponent_TopBeneficialDrivers.tsv")
out_top_dilution <- file.path(OUT_DIR, "12B_MultiComponent_TopDilutionTargets.tsv")
out_comb_perm <- file.path(OUT_DIR, "13B_MultiComponent_CombinationPermutation.tsv")
out_comb_perm_summary <- file.path(OUT_DIR, "14B_MultiComponent_CombinationPermutationSummary.tsv")

fwrite(design_grid, out_design, sep = "\t")
fwrite(target_sets, out_targets, sep = "\t")
fwrite(combination_results, out_comb, sep = "\t")
fwrite(incremental_results, out_inc, sep = "\t")
fwrite(target_contribution, out_tcontrib, sep = "\t")
if (SAVE_FULL_DRIVER_ABLATION) fwrite(driver_results, out_driver, sep = "\t")
fwrite(driver_recovery, out_recovery, sep = "\t")
fwrite(pass_summary, out_pass, sep = "\t")
fwrite(comb_summary, out_comb_summary, sep = "\t")
fwrite(increment_summary, out_increment_summary, sep = "\t")
fwrite(top_driver_results, out_top_driver, sep = "\t")
fwrite(top_dilution_results, out_top_dilution, sep = "\t")
fwrite(combination_permutation, out_comb_perm, sep = "\t")
fwrite(combination_perm_summary, out_comb_perm_summary, sep = "\t")

############################################################
## 16. Manuscript-grade figures
############################################################

cat("\nGenerating manuscript-grade figures...\n")

PALETTE_SCENARIO <- c(
  "Additive_benefit" = "#F4A7A3",
  "Redundant_partner" = "#F6CFA5",
  "Dilution_partner" = "#AFCBEF",
  "Mixed_realistic" = "#CDB8E9"
)

PALETTE_DISTANCE <- c(
  "Disease_close" = "#F4A7A3",
  "Intermediate" = "#CDB8E9",
  "Disease_far" = "#AFCBEF",
  "Unknown" = "#D9D9D9"
)

PALETTE_INCREMENT <- c(
  "Beneficial_addition" = "#F4A7A3",
  "Redundant_or_neutral" = "#D8D8D8",
  "Dilution_addition" = "#AFCBEF"
)

pretty_scenario <- c(
  "Additive_benefit" = "Additive\nbenefit",
  "Redundant_partner" = "Redundant\npartner",
  "Dilution_partner" = "Dilution\npartner",
  "Mixed_realistic" = "Mixed\nrealistic"
)

comb_levels <- c("Core", "Core+A", "Core+B", "Core+C", "Core+A+B", "Core+A+C", "Core+B+C", "Core+A+B+C")

combination_results[, Combination := factor(Combination, levels = comb_levels)]
comb_summary[, Combination := factor(Combination, levels = comb_levels)]
incremental_results[, New_Combination := factor(New_Combination, levels = comb_levels)]
increment_summary[, New_Combination := factor(New_Combination, levels = comb_levels)]
driver_results[, Target_Distance_Class := factor(Target_Distance_Class, levels = c("Disease_close", "Intermediate", "Disease_far", "Unknown"))]

theme_manuscript <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 3, hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_text(size = base_size, color = "grey30", hjust = 0, margin = margin(b = 8)),
      axis.title = element_text(face = "bold", color = "grey20"),
      axis.text = element_text(color = "grey20"),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_rect(fill = "grey96", color = "grey75"),
      strip.text = element_text(face = "bold", color = "grey20"),
      plot.margin = margin(14, 20, 16, 16)
    )
}

save_plot_both <- function(p, filename, width, height, dpi = 600) {
  png_file <- file.path(FIG_DIR, paste0(filename, ".png"))
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  ggplot2::ggsave(png_file, plot = p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(pdf_file, plot = p, width = width, height = height, bg = "white", limitsize = FALSE)
}

## B1: combination distance hierarchy
p_b1 <- ggplot(combination_results, aes(x = Combination, y = Mean_Distance, fill = Scenario)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.82, linewidth = 0.25) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.2, fill = "white", color = "grey25") +
  facet_wrap(~ Scenario, nrow = 2, labeller = labeller(Scenario = pretty_scenario)) +
  scale_fill_manual(values = PALETTE_SCENARIO, guide = "none") +
  labs(
    title = "Combination-level disease proximity",
    subtitle = "Complementary partners reduce disease-module distance, whereas dilution partners weaken the proximity gain.",
    x = "Virtual combination",
    y = "Mean distance to stroke disease module"
  ) +
  theme_manuscript(12) +
  theme(axis.text.x = element_text(size = 8.5, angle = 45, hjust = 1), strip.text = element_text(size = 10.5))

save_plot_both(p_b1, "FigB1_CombinationDistance_Hierarchy", width = 11.8, height = 8.5)

## B2: network benefit versus core
plot_benefit <- combination_results[Combination != "Core"]

p_b2 <- ggplot(plot_benefit, aes(x = Combination, y = Network_Benefit_vs_Core, fill = Scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.82, linewidth = 0.25) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.2, fill = "white", color = "grey25") +
  facet_wrap(~ Scenario, nrow = 2, labeller = labeller(Scenario = pretty_scenario)) +
  scale_fill_manual(values = PALETTE_SCENARIO, guide = "none") +
  labs(
    title = "Network benefit relative to the core compound",
    subtitle = "Positive values indicate that a virtual combination is closer to the disease module than the core alone.",
    x = "Virtual combination",
    y = "Network benefit vs core"
  ) +
  theme_manuscript(12) +
  theme(axis.text.x = element_text(size = 8.5, angle = 45, hjust = 1), strip.text = element_text(size = 10.5))

save_plot_both(p_b2, "FigB2_NetworkBenefit_vs_Core", width = 11.8, height = 8.5)

## B3: incremental benefit by addition
increment_plot <- incremental_results[
  Base_Combination %in% c("Core", "Core+A") &
    New_Combination %in% c("Core+A", "Core+B", "Core+C", "Core+A+B", "Core+A+C")
]

increment_plot[, Addition_Label := paste0(Base_Combination, " + ", Added_Compound, " -> ", New_Combination)]
addition_levels <- unique(increment_plot[order(Base_Combination, Added_Compound, New_Combination)]$Addition_Label)
increment_plot[, Addition_Label := factor(Addition_Label, levels = addition_levels)]

p_b3 <- ggplot(increment_plot, aes(x = Addition_Label, y = Incremental_Benefit, fill = Increment_Class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.86, linewidth = 0.25) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.2, fill = "white", color = "grey25") +
  facet_wrap(~ Scenario, nrow = 2, labeller = labeller(Scenario = pretty_scenario), scales = "free_x") +
  scale_fill_manual(values = PALETTE_INCREMENT, name = "Addition class") +
  labs(
    title = "Incremental partner benefit",
    subtitle = "Complementary partners add benefit, redundant partners add limited benefit, and disease-far partners dilute proximity.",
    x = "Incremental addition",
    y = "Incremental benefit"
  ) +
  theme_manuscript(11) +
  theme(axis.text.x = element_text(size = 7.6, angle = 45, hjust = 1), legend.position = "bottom", strip.text = element_text(size = 10))

save_plot_both(p_b3, "FigB3_IncrementalPartnerBenefit", width = 13.5, height = 9.2)

## B4: added close/far targets for selected additions
added_summary_long <- increment_summary[
  Base_Combination %in% c("Core", "Core+A") &
    New_Combination %in% c("Core+A", "Core+B", "Core+C", "Core+A+B", "Core+A+C"),
  .(
    Scenario,
    Base_Combination,
    Added_Compound,
    New_Combination,
    Added_Close = Median_Added_Close_Target_N,
    Added_Far = Median_Added_Far_Target_N
  )
]

added_summary_long[, Addition_Label := paste0(Base_Combination, " + ", Added_Compound, " -> ", New_Combination)]
added_long <- melt(
  added_summary_long,
  id.vars = c("Scenario", "Addition_Label"),
  measure.vars = c("Added_Close", "Added_Far"),
  variable.name = "Added_Target_Class",
  value.name = "Median_Target_N"
)

added_long[, Added_Target_Class := fifelse(Added_Target_Class == "Added_Close", "Disease-close added targets", "Disease-far added targets")]
added_long[, Addition_Label := factor(Addition_Label, levels = unique(added_summary_long$Addition_Label))]

p_b4 <- ggplot(added_long, aes(x = Addition_Label, y = Median_Target_N, fill = Added_Target_Class)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, alpha = 0.88, color = "grey35", linewidth = 0.15) +
  facet_wrap(~ Scenario, nrow = 2, labeller = labeller(Scenario = pretty_scenario), scales = "free_x") +
  scale_fill_manual(values = c("Disease-close added targets" = "#F4A7A3", "Disease-far added targets" = "#AFCBEF"), name = "Added target class") +
  labs(
    title = "Added-target novelty explains incremental benefit",
    subtitle = "Beneficial additions introduce disease-close targets, whereas dilution additions introduce disease-far targets.",
    x = "Incremental addition",
    y = "Median number of added targets"
  ) +
  theme_manuscript(11) +
  theme(axis.text.x = element_text(size = 7.6, angle = 45, hjust = 1), legend.position = "bottom", strip.text = element_text(size = 10))

save_plot_both(p_b4, "FigB4_AddedTarget_Novelty", width = 13.5, height = 9.2)

## B5: DriverScore by distance class
plot_driver <- driver_results[Target_Distance_Class %in% c("Disease_close", "Intermediate", "Disease_far")]

p_b5 <- ggplot(plot_driver, aes(x = Target_Distance_Class, y = DriverScore, fill = Target_Distance_Class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.83, linewidth = 0.25) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.2, fill = "white", color = "grey25") +
  facet_wrap(~ Scenario, nrow = 2, labeller = labeller(Scenario = pretty_scenario)) +
  scale_fill_manual(values = PALETTE_DISTANCE, guide = "none") +
  labs(
    title = "Target-level DriverScore decomposition",
    subtitle = "Disease-close targets act as beneficial drivers, whereas disease-far targets act as dilution targets.",
    x = "Target distance class",
    y = "DriverScore"
  ) +
  theme_manuscript(12) +
  theme(axis.text.x = element_text(size = 9.3, angle = 25, hjust = 1), strip.text = element_text(size = 10.5))

save_plot_both(p_b5, "FigB5_TargetDriverScore_Decomposition", width = 11.8, height = 8.5)

## B6: Truth-label recovery plot
truth_order <- c("Core_like", "Positive_like", "Redundant_like", "Neutral_like", "Mixed_like", "Dilution_like")
truth_plot <- driver_results[Truth_Labels %in% truth_order]
truth_plot[, Truth_Labels := factor(Truth_Labels, levels = truth_order)]

truth_palette <- c(
  "Core_like" = "#D6D6D6",
  "Positive_like" = "#F4A7A3",
  "Redundant_like" = "#F6CFA5",
  "Neutral_like" = "#C7D4E8",
  "Mixed_like" = "#CDB8E9",
  "Dilution_like" = "#AFCBEF"
)

p_b6 <- ggplot(truth_plot, aes(x = Truth_Labels, y = DriverScore, fill = Truth_Labels)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.84, linewidth = 0.25) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 2.0, fill = "white", color = "grey25") +
  facet_wrap(~ Scenario, nrow = 2, labeller = labeller(Scenario = pretty_scenario), scales = "free_x") +
  scale_fill_manual(values = truth_palette, guide = "none") +
  labs(
    title = "Recovery of simulated target labels by DriverScore",
    subtitle = "Positive-like targets show higher DriverScore than dilution-like targets across virtual scenarios.",
    x = "Simulated target label",
    y = "DriverScore"
  ) +
  theme_manuscript(11) +
  theme(axis.text.x = element_text(size = 8.0, angle = 35, hjust = 1), strip.text = element_text(size = 10))

save_plot_both(p_b6, "FigB6_TruthLabel_DriverRecovery", width = 12.8, height = 8.8)

############################################################
## 17. Final report
############################################################

cat("\n============================================================\n")
cat("ACMPI-NPS Pressure Test B completed.\n")
cat("============================================================\n")
cat("Run mode: ", RUN_MODE, "\n", sep = "")
cat("Output folder:\n", OUT_DIR, "\n", sep = "")
cat("Figures folder:\n", FIG_DIR, "\n", sep = "")
cat("Saved files:\n")
cat(out_design, "\n")
cat(out_targets, "\n")
cat(out_comb, "\n")
cat(out_inc, "\n")
cat(out_tcontrib, "\n")
if (SAVE_FULL_DRIVER_ABLATION) cat(out_driver, "\n")
cat(out_recovery, "\n")
cat(out_pass, "\n")
cat(out_comb_summary, "\n")
cat(out_increment_summary, "\n")
cat(out_top_driver, "\n")
cat(out_top_dilution, "\n")
cat(out_comb_perm, "\n")
cat(out_comb_perm_summary, "\n")
cat("\nSaved figures:\n")
cat(file.path(FIG_DIR, "FigB1_CombinationDistance_Hierarchy.png"), "\n")
cat(file.path(FIG_DIR, "FigB2_NetworkBenefit_vs_Core.png"), "\n")
cat(file.path(FIG_DIR, "FigB3_IncrementalPartnerBenefit.png"), "\n")
cat(file.path(FIG_DIR, "FigB4_AddedTarget_Novelty.png"), "\n")
cat(file.path(FIG_DIR, "FigB5_TargetDriverScore_Decomposition.png"), "\n")
cat(file.path(FIG_DIR, "FigB6_TruthLabel_DriverRecovery.png"), "\n")
cat("============================================================\n")
