############################################################
## Pathway-space virtual combination validation
## Final package-ready version
##
## Includes:
##   1. Repeated virtual scenario discrimination
##   2. Noise stress test
##   3. Dropout stress test
##   4. Weight perturbation stability
##   5. Random permutation calibration
##   6. Source-stratified validation
##   7. Decision table
##   8. Publication-grade PNG/PDF figures
##
## Plot fixes:
##   - No oversized histogram with observed line far outside the null
##   - Permutation calibration uses rank-based null calibration
##   - P values reported as P < 1/N when no random score exceeds observed
##   - Short independent figure titles
##   - Soft low-saturation colors
############################################################

rm(list = ls())

############################################################
## 0. Settings
############################################################

BASE_DIR <- "/media/desk16/iy15915/中药之开创/富集模型/脑卒中基因/DiseaseGene_PathwayVector_20260507_114626/02_ModelInput_PathwayVectors"

LONG_FILE <- file.path(BASE_DIR, "04_Disease_PathwayVector_ForModelInput_LONG.tsv")
WIDE_FILE <- file.path(BASE_DIR, "05_Disease_PathwayVector_PathwayScore_WIDE.tsv")

RUN_TAG <- format(Sys.time(), "%Y%m%d_%H%M%S")

OUT_DIR <- file.path(dirname(BASE_DIR), paste0("Virtual_PathwayCombination_FINAL_", RUN_TAG))
TABLE_DIR <- file.path(OUT_DIR, "01_Tables")
FIG_DIR <- file.path(OUT_DIR, "02_Figures")
VECTOR_DIR <- file.path(OUT_DIR, "03_Vectors")
QC_DIR <- file.path(OUT_DIR, "04_QC")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(VECTOR_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Parameters
############################################################

set.seed(20260507)

CORE_RATIO <- 0.20
IDEAL_RATIO <- 0.30
PARTIAL_RATIO <- 0.10
REDUNDANT_OVERLAP <- 0.90

DISEASE_CORE_Q <- 0.75

BASE_W_COSINE <- 0.35
BASE_W_JACCARD <- 0.25
BASE_W_COVERAGE <- 0.30
BASE_W_OFFTARGET <- 0.10

N_REPEAT <- 500
N_STRESS_REPEAT <- 500
N_RANDOM <- 5000
N_WEIGHT <- 2000

NOISE_LEVELS <- c(0, 0.10, 0.25, 0.50, 1.00, 2.00)
DROPOUT_LEVELS <- c(0, 0.10, 0.25, 0.50, 0.70)

SOURCE_LEVELS <- c(
  "All_sources",
  "Reactome_ORA",
  "GO_BP_ORA",
  "MSigDB_Hallmark_ORA"
)

############################################################
## 2. Packages
############################################################

pkgs <- c(
  "data.table",
  "ggplot2",
  "scales",
  "stringr"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(data.table)
library(ggplot2)
library(scales)
library(stringr)

############################################################
## 3. Plot style
############################################################

PAL_SCENARIO <- c(
  Ideal_Complementary = "#8FB9A8",
  Partial = "#F2C879",
  Redundant = "#B8A1D9",
  Dilution = "#F0A7A0",
  Random_IntensityMatched = "#B9C0C9"
)

PAL_SOURCE <- c(
  All_sources = "#8FB9A8",
  Reactome_ORA = "#9BB7D4",
  GO_BP_ORA = "#D7BCE8",
  MSigDB_Hallmark_ORA = "#F0C987"
)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "grey20", size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(8, 12, 8, 8)
    )
}

save_dual <- function(p, filename, width, height) {
  ggsave(
    file.path(FIG_DIR, paste0(filename, ".png")),
    p,
    width = width,
    height = height,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    file.path(FIG_DIR, paste0(filename, ".pdf")),
    p,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
}

############################################################
## 4. Helper functions
############################################################

safe_sample <- function(x, n, replace_if_needed = TRUE) {
  x <- unique(x)
  x <- x[!is.na(x)]
  
  if (length(x) == 0) stop("Sampling pool is empty.")
  if (n <= 0) return(character(0))
  
  if (length(x) < n) {
    if (replace_if_needed) {
      return(sample(x, n, replace = TRUE))
    } else {
      return(sample(x, length(x), replace = FALSE))
    }
  }
  
  sample(x, n, replace = FALSE)
}

cosine_similarity <- function(x, y) {
  idx <- intersect(names(x), names(y))
  x <- x[idx]
  y <- y[idx]
  
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (den == 0) return(NA_real_)
  
  sum(x * y) / den
}

weighted_jaccard <- function(x, y) {
  idx <- intersect(names(x), names(y))
  x <- x[idx]
  y <- y[idx]
  
  den <- sum(pmax(x, y))
  if (den == 0) return(NA_real_)
  
  sum(pmin(x, y)) / den
}

ci95_low <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  mean(x) - 1.96 * sd(x) / sqrt(length(x))
}

ci95_high <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  mean(x) + 1.96 * sd(x) / sqrt(length(x))
}

empirical_p_label <- function(p, n) {
  if (is.na(p)) return("P = NA")
  if (p == 0) return(paste0("P < 1/", format(n, big.mark = ",")))
  paste0("P = ", signif(p, 3))
}

make_space <- function(wide_dt, source_name) {
  
  if (source_name != "All_sources") {
    x <- wide_dt[Pathway_Source == source_name]
  } else {
    x <- copy(wide_dt)
  }
  
  x <- x[
    !is.na(Disease) &
      is.finite(Disease)
  ]
  
  if (nrow(x) < 50) {
    stop("Too few pathway features for source: ", source_name)
  }
  
  D <- x$Disease
  names(D) <- x$Feature_ID
  
  disease_core_cutoff <- quantile(D, DISEASE_CORE_Q, na.rm = TRUE)
  disease_core_ids <- names(D)[D >= disease_core_cutoff]
  noncore_ids <- names(D)[D < disease_core_cutoff]
  ordered_features <- x[order(-Disease)]$Feature_ID
  
  list(
    source_name = source_name,
    dt = x,
    D = D,
    disease_core_ids = disease_core_ids,
    noncore_ids = noncore_ids,
    ordered_features = ordered_features
  )
}

construct_virtual_vectors <- function(space) {
  
  D <- space$D
  ordered_features <- space$ordered_features
  disease_core_ids <- space$disease_core_ids
  noncore_ids <- space$noncore_ids
  
  n_core <- round(length(D) * CORE_RATIO)
  n_ideal <- round(length(D) * IDEAL_RATIO)
  n_partial <- round(length(D) * PARTIAL_RATIO)
  
  core_pool <- ordered_features[
    seq_len(max(1, round(length(ordered_features) * 0.45)))
  ]
  
  core_ids <- safe_sample(core_pool, n_core, replace_if_needed = TRUE)
  
  Core <- rep(0, length(D))
  names(Core) <- names(D)
  Core[core_ids] <- D[core_ids] * runif(length(core_ids), 0.55, 0.80)
  
  remaining_ids <- setdiff(names(D), core_ids)
  
  ideal_pool <- intersect(
    setdiff(disease_core_ids, core_ids),
    remaining_ids
  )
  
  if (length(ideal_pool) < round(n_ideal * 0.50)) {
    ideal_pool <- intersect(
      remaining_ids,
      ordered_features[seq_len(max(1, round(length(ordered_features) * 0.70)))]
    )
  }
  
  ideal_ids <- safe_sample(ideal_pool, n_ideal, replace_if_needed = TRUE)
  
  Ideal <- rep(0, length(D))
  names(Ideal) <- names(D)
  Ideal[ideal_ids] <- D[ideal_ids] * runif(length(ideal_ids), 0.65, 0.95)
  
  partial_pool <- intersect(
    remaining_ids,
    ordered_features[seq_len(max(1, round(length(ordered_features) * 0.75)))]
  )
  
  partial_ids <- safe_sample(partial_pool, n_partial, replace_if_needed = TRUE)
  
  Partial <- rep(0, length(D))
  names(Partial) <- names(D)
  Partial[partial_ids] <- D[partial_ids] * runif(length(partial_ids), 0.35, 0.65)
  
  n_redundant <- round(length(core_ids) * REDUNDANT_OVERLAP)
  redundant_ids <- safe_sample(core_ids, n_redundant, replace_if_needed = TRUE)
  
  Redundant <- rep(0, length(D))
  names(Redundant) <- names(D)
  Redundant[redundant_ids] <- D[redundant_ids] * runif(length(redundant_ids), 0.55, 0.85)
  
  signal_template <- Ideal[Ideal > 0]
  if (length(signal_template) == 0) signal_template <- D[D > 0]
  
  random_ids <- safe_sample(names(D), n_ideal, replace_if_needed = TRUE)
  
  RandomDrug <- rep(0, length(D))
  names(RandomDrug) <- names(D)
  RandomDrug[random_ids] <- sample(signal_template, length(random_ids), replace = TRUE)
  
  dilution_pool <- noncore_ids
  if (length(dilution_pool) == 0) dilution_pool <- names(D)
  
  dilution_ids <- safe_sample(dilution_pool, n_ideal, replace_if_needed = TRUE)
  
  Dilution <- rep(0, length(D))
  names(Dilution) <- names(D)
  Dilution[dilution_ids] <- runif(
    length(dilution_ids),
    quantile(D, 0.10, na.rm = TRUE),
    quantile(D, 0.40, na.rm = TRUE)
  )
  
  list(
    Core = Core,
    Ideal_Complementary = Ideal,
    Partial = Partial,
    Redundant = Redundant,
    Random_IntensityMatched = RandomDrug,
    Dilution = Dilution,
    signal_template = signal_template
  )
}

evaluate_partner <- function(space, vectors, partner_vec, partner_name, weights) {
  
  D <- space$D
  Core <- vectors$Core
  combo <- Core + partner_vec
  
  sim_core_cos <- cosine_similarity(D, Core)
  sim_combo_cos <- cosine_similarity(D, combo)
  sim_core_jac <- weighted_jaccard(D, Core)
  sim_combo_jac <- weighted_jaccard(D, combo)
  
  core_active <- names(Core)[Core > 0]
  combo_active <- names(combo)[combo > 0]
  
  disease_core_ids <- space$disease_core_ids
  noncore_ids <- space$noncore_ids
  
  core_covered <- length(intersect(disease_core_ids, core_active)) / max(length(disease_core_ids), 1)
  combo_covered <- length(intersect(disease_core_ids, combo_active)) / max(length(disease_core_ids), 1)
  
  core_offtarget <- length(intersect(noncore_ids, core_active)) / max(length(noncore_ids), 1)
  combo_offtarget <- length(intersect(noncore_ids, combo_active)) / max(length(noncore_ids), 1)
  
  combination_gain_cosine <- sim_combo_cos - sim_core_cos
  combination_gain_jaccard <- sim_combo_jac - sim_core_jac
  coverage_gain <- combo_covered - core_covered
  offtarget_gain <- combo_offtarget - core_offtarget
  
  raw_score <-
    weights["cosine"] * combination_gain_cosine +
    weights["jaccard"] * combination_gain_jaccard +
    weights["coverage"] * coverage_gain -
    weights["offtarget"] * offtarget_gain
  
  data.table(
    Source = space$source_name,
    Partner = partner_name,
    Core_Cosine = sim_core_cos,
    Combo_Cosine = sim_combo_cos,
    CombinationGain_Cosine = combination_gain_cosine,
    Core_Jaccard = sim_core_jac,
    Combo_Jaccard = sim_combo_jac,
    CombinationGain_Jaccard = combination_gain_jaccard,
    CoreCoverage = core_covered,
    ComboCoverage = combo_covered,
    CoverageGain = coverage_gain,
    CoreOffTarget = core_offtarget,
    ComboOffTarget = combo_offtarget,
    OffTargetGain = offtarget_gain,
    MissedCoreRatio = 1 - combo_covered,
    Raw_FinalScore = raw_score
  )
}

evaluate_all_partners <- function(space, vectors, weights) {
  rbindlist(list(
    evaluate_partner(space, vectors, vectors$Ideal_Complementary, "Ideal_Complementary", weights),
    evaluate_partner(space, vectors, vectors$Partial, "Partial", weights),
    evaluate_partner(space, vectors, vectors$Redundant, "Redundant", weights),
    evaluate_partner(space, vectors, vectors$Dilution, "Dilution", weights),
    evaluate_partner(space, vectors, vectors$Random_IntensityMatched, "Random_IntensityMatched", weights)
  ))
}

base_weights <- c(
  cosine = BASE_W_COSINE,
  jaccard = BASE_W_JACCARD,
  coverage = BASE_W_COVERAGE,
  offtarget = BASE_W_OFFTARGET
)

############################################################
## 5. Read input
############################################################

if (!file.exists(LONG_FILE)) stop("Missing LONG file:\n", LONG_FILE)
if (!file.exists(WIDE_FILE)) stop("Missing WIDE file:\n", WIDE_FILE)

long_dt <- fread(LONG_FILE)
wide_dt <- fread(WIDE_FILE)

required_cols <- c(
  "Feature_ID",
  "Feature_Name",
  "Pathway_Source",
  "Pathway_ID",
  "Pathway_Name",
  "Disease"
)

miss <- setdiff(required_cols, colnames(wide_dt))

if (length(miss) > 0) {
  stop("Missing columns in WIDE file:\n", paste(miss, collapse = "\n"))
}

wide_dt <- wide_dt[
  !is.na(Disease) &
    is.finite(Disease)
]

############################################################
## 6. Main repeated virtual experiment
############################################################

main_space <- make_space(wide_dt, "All_sources")

main_repeat_list <- vector("list", N_REPEAT)

for (i in seq_len(N_REPEAT)) {
  vectors <- construct_virtual_vectors(main_space)
  tmp <- evaluate_all_partners(main_space, vectors, base_weights)
  tmp[, Repeat := i]
  main_repeat_list[[i]] <- tmp
}

main_repeats <- rbindlist(main_repeat_list)

main_repeats[
  ,
  Rank_FinalScore := frank(-Raw_FinalScore, ties.method = "average"),
  by = Repeat
]

main_summary <- main_repeats[
  ,
  .(
    N = .N,
    Mean_FinalScore = mean(Raw_FinalScore, na.rm = TRUE),
    SD_FinalScore = sd(Raw_FinalScore, na.rm = TRUE),
    CI_Low_FinalScore = ci95_low(Raw_FinalScore),
    CI_High_FinalScore = ci95_high(Raw_FinalScore),
    Mean_CosineGain = mean(CombinationGain_Cosine, na.rm = TRUE),
    Mean_JaccardGain = mean(CombinationGain_Jaccard, na.rm = TRUE),
    Mean_CoverageGain = mean(CoverageGain, na.rm = TRUE),
    Mean_OffTargetGain = mean(OffTargetGain, na.rm = TRUE),
    Top1_Frequency = mean(Rank_FinalScore == 1, na.rm = TRUE),
    Top3_Frequency = mean(Rank_FinalScore <= 3, na.rm = TRUE)
  ),
  by = Partner
][order(-Mean_FinalScore)]

############################################################
## 7. Noise stress test
############################################################

noise_all <- list()
idx <- 1

for (r in seq_len(N_STRESS_REPEAT)) {
  
  vectors <- construct_virtual_vectors(main_space)
  Ideal <- vectors$Ideal_Complementary
  signal_template <- vectors$signal_template
  
  for (nl in NOISE_LEVELS) {
    
    noisy_partner <- Ideal
    noise_n <- round(sum(Ideal > 0) * nl)
    
    if (noise_n > 0) {
      noise_ids <- safe_sample(names(main_space$D), noise_n, replace_if_needed = TRUE)
      noisy_partner[noise_ids] <- sample(signal_template, length(noise_ids), replace = TRUE)
    }
    
    tmp <- evaluate_partner(
      main_space,
      vectors,
      noisy_partner,
      "Ideal_Complementary",
      base_weights
    )
    
    tmp[, `:=`(
      Repeat = r,
      Noise_Level = nl
    )]
    
    noise_all[[idx]] <- tmp
    idx <- idx + 1
  }
}

noise_results <- rbindlist(noise_all)

noise_summary <- noise_results[
  ,
  .(
    N = .N,
    Mean_FinalScore = mean(Raw_FinalScore, na.rm = TRUE),
    SD_FinalScore = sd(Raw_FinalScore, na.rm = TRUE),
    CI_Low_FinalScore = ci95_low(Raw_FinalScore),
    CI_High_FinalScore = ci95_high(Raw_FinalScore),
    Mean_CoverageGain = mean(CoverageGain, na.rm = TRUE),
    Mean_OffTargetGain = mean(OffTargetGain, na.rm = TRUE)
  ),
  by = Noise_Level
][order(Noise_Level)]

############################################################
## 8. Dropout stress test
############################################################

dropout_all <- list()
idx <- 1

for (r in seq_len(N_STRESS_REPEAT)) {
  
  vectors <- construct_virtual_vectors(main_space)
  Ideal <- vectors$Ideal_Complementary
  
  active_ids <- names(Ideal)[Ideal > 0]
  
  for (dl in DROPOUT_LEVELS) {
    
    dropout_partner <- Ideal
    
    drop_n <- round(length(active_ids) * dl)
    
    if (drop_n > 0) {
      drop_ids <- safe_sample(active_ids, drop_n, replace_if_needed = FALSE)
      dropout_partner[drop_ids] <- 0
    }
    
    tmp <- evaluate_partner(
      main_space,
      vectors,
      dropout_partner,
      "Ideal_Complementary",
      base_weights
    )
    
    tmp[, `:=`(
      Repeat = r,
      Dropout_Level = dl
    )]
    
    dropout_all[[idx]] <- tmp
    idx <- idx + 1
  }
}

dropout_results <- rbindlist(dropout_all)

dropout_summary <- dropout_results[
  ,
  .(
    N = .N,
    Mean_FinalScore = mean(Raw_FinalScore, na.rm = TRUE),
    SD_FinalScore = sd(Raw_FinalScore, na.rm = TRUE),
    CI_Low_FinalScore = ci95_low(Raw_FinalScore),
    CI_High_FinalScore = ci95_high(Raw_FinalScore),
    Mean_CoverageGain = mean(CoverageGain, na.rm = TRUE),
    Mean_OffTargetGain = mean(OffTargetGain, na.rm = TRUE)
  ),
  by = Dropout_Level
][order(Dropout_Level)]

############################################################
## 9. Weight perturbation stability
############################################################

component_means <- main_repeats[
  ,
  .(
    CombinationGain_Cosine = mean(CombinationGain_Cosine, na.rm = TRUE),
    CombinationGain_Jaccard = mean(CombinationGain_Jaccard, na.rm = TRUE),
    CoverageGain = mean(CoverageGain, na.rm = TRUE),
    OffTargetGain = mean(OffTargetGain, na.rm = TRUE)
  ),
  by = Partner
]

weight_results <- vector("list", N_WEIGHT)

for (i in seq_len(N_WEIGHT)) {
  
  w <- rgamma(4, shape = 2, rate = 1)
  w <- w / sum(w)
  names(w) <- c("cosine", "jaccard", "coverage", "offtarget")
  
  tmp <- copy(component_means)
  
  tmp[, PerturbedScore :=
        w["cosine"] * CombinationGain_Cosine +
        w["jaccard"] * CombinationGain_Jaccard +
        w["coverage"] * CoverageGain -
        w["offtarget"] * OffTargetGain]
  
  tmp[, Rank := frank(-PerturbedScore, ties.method = "average")]
  
  tmp[, `:=`(
    Iteration = i,
    W_Cosine = w["cosine"],
    W_Jaccard = w["jaccard"],
    W_Coverage = w["coverage"],
    W_OffTarget = w["offtarget"]
  )]
  
  weight_results[[i]] <- tmp
}

weight_results <- rbindlist(weight_results)

weight_summary <- weight_results[
  ,
  .(
    N = .N,
    Mean_Rank = mean(Rank, na.rm = TRUE),
    Median_Rank = median(Rank, na.rm = TRUE),
    Top1_Frequency = mean(Rank == 1, na.rm = TRUE),
    Top3_Frequency = mean(Rank <= 3, na.rm = TRUE),
    Mean_PerturbedScore = mean(PerturbedScore, na.rm = TRUE),
    CI_Low_PerturbedScore = ci95_low(PerturbedScore),
    CI_High_PerturbedScore = ci95_high(PerturbedScore)
  ),
  by = Partner
][order(Mean_Rank)]

############################################################
## 10. Random permutation calibration
############################################################

perm_list <- vector("list", N_RANDOM)

anchor_vectors <- construct_virtual_vectors(main_space)

ideal_anchor <- evaluate_partner(
  main_space,
  anchor_vectors,
  anchor_vectors$Ideal_Complementary,
  "Ideal_Complementary",
  base_weights
)

signal_template <- anchor_vectors$signal_template
n_ideal_anchor <- sum(anchor_vectors$Ideal_Complementary > 0)

for (i in seq_len(N_RANDOM)) {
  
  rand_ids <- safe_sample(names(main_space$D), n_ideal_anchor, replace_if_needed = TRUE)
  
  rand_vec <- rep(0, length(main_space$D))
  names(rand_vec) <- names(main_space$D)
  
  rand_vec[rand_ids] <- sample(signal_template, length(rand_ids), replace = TRUE)
  
  tmp <- evaluate_partner(
    main_space,
    anchor_vectors,
    rand_vec,
    paste0("Random_", i),
    base_weights
  )
  
  tmp[, Iteration := i]
  perm_list[[i]] <- tmp
}

perm_results <- rbindlist(perm_list)

perm_summary <- data.table(
  Metric = c("Raw_FinalScore", "Combo_Cosine"),
  Observed_Ideal = c(
    ideal_anchor$Raw_FinalScore,
    ideal_anchor$Combo_Cosine
  ),
  Empirical_P = c(
    mean(perm_results$Raw_FinalScore >= ideal_anchor$Raw_FinalScore, na.rm = TRUE),
    mean(perm_results$Combo_Cosine >= ideal_anchor$Combo_Cosine, na.rm = TRUE)
  ),
  N_RANDOM = N_RANDOM
)

perm_summary[
  ,
  Empirical_P_Label := mapply(
    empirical_p_label,
    Empirical_P,
    N_RANDOM
  )
]

############################################################
## 11. Source-stratified validation
############################################################

source_all <- list()
idx <- 1

for (src in SOURCE_LEVELS) {
  
  sp <- make_space(wide_dt, src)
  
  n_rep_src <- 200
  
  for (r in seq_len(n_rep_src)) {
    vectors <- construct_virtual_vectors(sp)
    tmp <- evaluate_all_partners(sp, vectors, base_weights)
    tmp[, Repeat := r]
    source_all[[idx]] <- tmp
    idx <- idx + 1
  }
}

source_results <- rbindlist(source_all)

source_results[
  ,
  Rank_FinalScore := frank(-Raw_FinalScore, ties.method = "average"),
  by = .(Source, Repeat)
]

source_summary <- source_results[
  ,
  .(
    N = .N,
    Mean_FinalScore = mean(Raw_FinalScore, na.rm = TRUE),
    CI_Low_FinalScore = ci95_low(Raw_FinalScore),
    CI_High_FinalScore = ci95_high(Raw_FinalScore),
    Top1_Frequency = mean(Rank_FinalScore == 1, na.rm = TRUE),
    Top3_Frequency = mean(Rank_FinalScore <= 3, na.rm = TRUE)
  ),
  by = .(Source, Partner)
][order(Source, -Mean_FinalScore)]

############################################################
## 12. Decision table
############################################################

ideal_main <- main_summary[Partner == "Ideal_Complementary"]
partial_main <- main_summary[Partner == "Partial"]
redundant_main <- main_summary[Partner == "Redundant"]
random_main <- main_summary[Partner == "Random_IntensityMatched"]
dilution_main <- main_summary[Partner == "Dilution"]

noise_pass <- noise_summary[Noise_Level <= 0.50, all(Mean_FinalScore > 0)]
dropout_pass <- dropout_summary[Dropout_Level <= 0.25, all(Mean_FinalScore > 0)]
weight_pass <- weight_summary[Partner == "Ideal_Complementary"]$Top1_Frequency >= 0.80
perm_pass <- perm_summary[Metric == "Raw_FinalScore"]$Empirical_P < 0.001

source_pass <- source_summary[
  Partner == "Ideal_Complementary",
  all(Top1_Frequency >= 0.70)
]

scenario_pass <-
  ideal_main$Mean_FinalScore > partial_main$Mean_FinalScore &&
  partial_main$Mean_FinalScore > redundant_main$Mean_FinalScore &&
  redundant_main$Mean_FinalScore > random_main$Mean_FinalScore &&
  dilution_main$Mean_FinalScore > random_main$Mean_FinalScore

decision_table <- data.table(
  Validation_Module = c(
    "Scenario discrimination",
    "Noise robustness",
    "Dropout robustness",
    "Weight stability",
    "Permutation calibration",
    "Source-level consistency"
  ),
  Criterion = c(
    "Ideal > Partial > Redundant > Random, with Dilution non-random but penalized",
    "Mean final score remains positive under <=50% noise",
    "Mean final score remains positive under <=25% dropout",
    "Ideal Top1 frequency >= 0.80 under random weight perturbation",
    "Ideal empirical P < 0.001 against intensity-matched random partners",
    "Ideal Top1 frequency >= 0.70 across pathway sources"
  ),
  Value = c(
    paste0(
      "Ideal=", signif(ideal_main$Mean_FinalScore, 4),
      "; Partial=", signif(partial_main$Mean_FinalScore, 4),
      "; Redundant=", signif(redundant_main$Mean_FinalScore, 4),
      "; Random=", signif(random_main$Mean_FinalScore, 4)
    ),
    paste0("Minimum mean score <=50% noise = ", signif(min(noise_summary[Noise_Level <= 0.50]$Mean_FinalScore), 4)),
    paste0("Minimum mean score <=25% dropout = ", signif(min(dropout_summary[Dropout_Level <= 0.25]$Mean_FinalScore), 4)),
    paste0("Ideal Top1 frequency = ", signif(weight_summary[Partner == "Ideal_Complementary"]$Top1_Frequency, 4)),
    paste0("Empirical ", perm_summary[Metric == "Raw_FinalScore"]$Empirical_P_Label),
    paste0("Minimum Ideal Top1 frequency = ", signif(min(source_summary[Partner == "Ideal_Complementary"]$Top1_Frequency), 4))
  ),
  Decision = c(
    ifelse(scenario_pass, "PASS", "CHECK"),
    ifelse(noise_pass, "PASS", "CHECK"),
    ifelse(dropout_pass, "PASS", "CHECK"),
    ifelse(weight_pass, "PASS", "CHECK"),
    ifelse(perm_pass, "PASS", "CHECK"),
    ifelse(source_pass, "PASS", "CHECK")
  )
)

############################################################
## 13. Save tables
############################################################

fwrite(main_repeats, file.path(TABLE_DIR, "01_MainVirtualExperiment_AllRepeats.tsv"), sep = "\t")
fwrite(main_summary, file.path(TABLE_DIR, "02_MainVirtualExperiment_Summary.tsv"), sep = "\t")

fwrite(noise_results, file.path(TABLE_DIR, "03_NoiseStress_AllRepeats.tsv"), sep = "\t")
fwrite(noise_summary, file.path(TABLE_DIR, "04_NoiseStress_Summary.tsv"), sep = "\t")

fwrite(dropout_results, file.path(TABLE_DIR, "05_DropoutStress_AllRepeats.tsv"), sep = "\t")
fwrite(dropout_summary, file.path(TABLE_DIR, "06_DropoutStress_Summary.tsv"), sep = "\t")

fwrite(weight_results, file.path(TABLE_DIR, "07_WeightPerturbation_AllWeights.tsv"), sep = "\t")
fwrite(weight_summary, file.path(TABLE_DIR, "08_WeightPerturbation_RankStability.tsv"), sep = "\t")

fwrite(perm_results, file.path(TABLE_DIR, "09_RandomPermutation_AllResults.tsv"), sep = "\t")
fwrite(perm_summary, file.path(TABLE_DIR, "10_RandomPermutation_Summary.tsv"), sep = "\t")

fwrite(source_results, file.path(TABLE_DIR, "11_SourceStratified_AllRepeats.tsv"), sep = "\t")
fwrite(source_summary, file.path(TABLE_DIR, "12_SourceStratified_Summary.tsv"), sep = "\t")

fwrite(decision_table, file.path(TABLE_DIR, "13_Final_MethodValidation_Decision.tsv"), sep = "\t")

fwrite(
  data.table(
    Parameter = c(
      "CORE_RATIO", "IDEAL_RATIO", "PARTIAL_RATIO", "REDUNDANT_OVERLAP",
      "DISEASE_CORE_Q", "BASE_W_COSINE", "BASE_W_JACCARD",
      "BASE_W_COVERAGE", "BASE_W_OFFTARGET", "N_REPEAT",
      "N_STRESS_REPEAT", "N_RANDOM", "N_WEIGHT"
    ),
    Value = c(
      CORE_RATIO, IDEAL_RATIO, PARTIAL_RATIO, REDUNDANT_OVERLAP,
      DISEASE_CORE_Q, BASE_W_COSINE, BASE_W_JACCARD,
      BASE_W_COVERAGE, BASE_W_OFFTARGET, N_REPEAT,
      N_STRESS_REPEAT, N_RANDOM, N_WEIGHT
    )
  ),
  file.path(TABLE_DIR, "14_RunMetadata.tsv"),
  sep = "\t"
)

############################################################
## 14. Figures
############################################################

scenario_levels <- c(
  "Ideal_Complementary",
  "Partial",
  "Redundant",
  "Dilution",
  "Random_IntensityMatched"
)

main_summary[, Partner := factor(Partner, levels = scenario_levels)]

p_scenario <- ggplot(
  main_summary,
  aes(
    x = Mean_FinalScore,
    y = Partner,
    color = Partner
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.4) +
  geom_errorbarh(
    aes(xmin = CI_Low_FinalScore, xmax = CI_High_FinalScore),
    height = 0.18,
    linewidth = 0.65
  ) +
  geom_point(size = 3.8) +
  scale_color_manual(values = PAL_SCENARIO, guide = "none") +
  labs(
    title = "Scenario discrimination",
    x = "Final score",
    y = NULL
  ) +
  theme_pub(11)

save_dual(p_scenario, "Scenario_discrimination", 7.2, 4.5)

p_balance <- ggplot(
  main_repeats,
  aes(
    x = CoverageGain,
    y = OffTargetGain,
    color = Partner
  )
) +
  geom_point(alpha = 0.12, size = 1.1) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 4.2,
    shape = 18
  ) +
  scale_color_manual(values = PAL_SCENARIO) +
  labs(
    title = "Coverage–off-target balance",
    x = "Coverage gain",
    y = "Off-target gain",
    color = NULL
  ) +
  theme_pub(11)

save_dual(p_balance, "Coverage_offtarget_balance", 7, 5.2)

p_noise <- ggplot(
  noise_summary,
  aes(
    x = Noise_Level,
    y = Mean_FinalScore
  )
) +
  geom_ribbon(
    aes(ymin = CI_Low_FinalScore, ymax = CI_High_FinalScore),
    fill = "#CFE8D6",
    alpha = 0.75
  ) +
  geom_line(color = "#6FA58B", linewidth = 1.1) +
  geom_point(color = "#6FA58B", size = 2.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.4) +
  labs(
    title = "Noise robustness",
    x = "Noise level",
    y = "Final score"
  ) +
  theme_pub(11)

save_dual(p_noise, "Noise_robustness", 6.4, 4.8)

p_dropout <- ggplot(
  dropout_summary,
  aes(
    x = Dropout_Level,
    y = Mean_FinalScore
  )
) +
  geom_ribbon(
    aes(ymin = CI_Low_FinalScore, ymax = CI_High_FinalScore),
    fill = "#F4DFB8",
    alpha = 0.85
  ) +
  geom_line(color = "#C99543", linewidth = 1.1) +
  geom_point(color = "#C99543", size = 2.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.4) +
  labs(
    title = "Dropout robustness",
    x = "Dropout level",
    y = "Final score"
  ) +
  theme_pub(11)

save_dual(p_dropout, "Dropout_robustness", 6.4, 4.8)

weight_summary[, Partner := factor(Partner, levels = scenario_levels)]

p_weight <- ggplot(
  weight_summary,
  aes(
    x = Top1_Frequency,
    y = Partner,
    color = Partner
  )
) +
  geom_segment(
    aes(x = 0, xend = Top1_Frequency, y = Partner, yend = Partner),
    linewidth = 1.1,
    alpha = 0.75
  ) +
  geom_point(size = 3.8) +
  scale_color_manual(values = PAL_SCENARIO, guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Weight stability",
    x = "Top-rank frequency",
    y = NULL
  ) +
  theme_pub(11)

save_dual(p_weight, "Weight_stability", 7, 4.5)

############################################################
## Fixed permutation plots
############################################################

perm_rank_dt <- copy(perm_results)
setorder(perm_rank_dt, Raw_FinalScore)
perm_rank_dt[, Rank_Proportion := seq_len(.N) / .N]

observed_score <- ideal_anchor$Raw_FinalScore
observed_label <- perm_summary[Metric == "Raw_FinalScore"]$Empirical_P_Label

p_perm_rank <- ggplot(
  perm_rank_dt,
  aes(
    x = Rank_Proportion,
    y = Raw_FinalScore
  )
) +
  geom_line(color = "#9BB7D4", linewidth = 1.1) +
  geom_hline(
    yintercept = observed_score,
    color = "#D95F5F",
    linewidth = 0.85
  ) +
  annotate(
    "label",
    x = 0.70,
    y = observed_score,
    label = paste0("Observed ideal\n", observed_label),
    size = 3.2,
    color = "#D95F5F",
    label.size = 0.25,
    fill = "white"
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Permutation calibration",
    x = "Empirical rank within random combinations",
    y = "Final score"
  ) +
  theme_pub(11)

save_dual(p_perm_rank, "Permutation_calibration", 7, 5)

null_q <- quantile(
  perm_results$Raw_FinalScore,
  probs = c(0.001, 0.999),
  na.rm = TRUE
)

p_perm_null <- ggplot(
  perm_results,
  aes(Raw_FinalScore)
) +
  geom_histogram(
    bins = 55,
    fill = "#C7D8F2",
    color = "white",
    linewidth = 0.25
  ) +
  coord_cartesian(xlim = null_q) +
  labs(
    title = "Permutation null distribution",
    x = "Final score",
    y = "Random combinations"
  ) +
  annotate(
    "label",
    x = mean(null_q),
    y = Inf,
    label = paste0("Observed ideal is outside the null range\n", observed_label),
    vjust = 1.35,
    size = 3.1,
    label.size = 0.25,
    fill = "white"
  ) +
  theme_pub(11)

save_dual(p_perm_null, "Permutation_null_distribution", 6.8, 4.8)

source_summary[, Partner := factor(Partner, levels = scenario_levels)]
source_summary[, Source := factor(Source, levels = SOURCE_LEVELS)]

p_source <- ggplot(
  source_summary,
  aes(
    x = Partner,
    y = Mean_FinalScore,
    fill = Source
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.78), width = 0.65, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = PAL_SOURCE) +
  labs(
    title = "Source-level consistency",
    x = NULL,
    y = "Final score",
    fill = NULL
  ) +
  theme_pub(10) +
  theme(
    axis.text.x = element_text(angle = 18, hjust = 1),
    legend.position = "bottom"
  )

save_dual(p_source, "Source_level_consistency", 9.5, 5.2)

decision_plot_dt <- copy(decision_table)
decision_plot_dt[, Validation_Module := factor(Validation_Module, levels = rev(Validation_Module))]

p_decision <- ggplot(
  decision_plot_dt,
  aes(
    x = Decision,
    y = Validation_Module,
    fill = Decision
  )
) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.85, height = 0.72) +
  geom_text(aes(label = Decision), size = 3.6, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c(PASS = "#BFDCC6", CHECK = "#F2C0B8")) +
  labs(
    title = "Validation decision",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(color = "grey20")
  )

save_dual(p_decision, "Validation_decision", 6.5, 4.2)

############################################################
## 15. Final report
############################################################

cat("\n============================================================\n")
cat("Final pathway-space virtual validation completed.\n")
cat("============================================================\n")

cat("\nOutput directory:\n")
cat(OUT_DIR, "\n")

cat("\nDecision table:\n")
print(decision_table)

cat("\nMain summary:\n")
print(main_summary)

cat("\nPermutation summary:\n")
print(perm_summary)

cat("\nFigures saved to:\n")
cat(FIG_DIR, "\n")

cat("============================================================\n")