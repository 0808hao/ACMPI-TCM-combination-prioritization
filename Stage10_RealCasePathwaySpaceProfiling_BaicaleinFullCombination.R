############################################################
## Real disease-specific baicalein-based full-combination analysis
##
## Purpose:
##   Each disease is analysed independently.
##   Baicalein is fixed as the anchor compound.
##   The other five compounds are evaluated as all non-empty subsets:
##     C(5,1) + C(5,2) + C(5,3) + C(5,4) + C(5,5) = 31 combinations.
##
## Key rule:
##   Input directory names are kept identical to the local folders.
##   Output labels, tables, figure titles, and file names are in English.
##
## Critical:
##   Algorithm parameters and scoring weights are not changed.
############################################################

rm(list = ls())

############################################################
## 0. Settings
############################################################

ROOT_DIR <- "/media/desk16/iy15915/中药之开创/富集模型"

RUN_TAG <- format(Sys.time(), "%Y%m%d_%H%M%S")

OUT_ROOT <- file.path(
  ROOT_DIR,
  paste0("Real_Baicalein_FullCombination_OutputEnglish_", RUN_TAG)
)

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

DISEASE_CONFIG <- data.table::data.table(
  Disease_Label = c(
    "Stroke",
    "Rheumatoid_Arthritis",
    "Diabetes_Mellitus",
    "Osteoarthritis"
  ),
  Disease_DirName = c(
    "脑卒中基因",
    "类风湿关节炎",
    "糖尿病",
    "骨关节炎"
  )
)

DRUG_ROOT_DIRNAME <- "药物基因"

ANCHOR_DRUG <- "Baicalein"

PARTNER_DRUGS <- c(
  "Caffeic Acid",
  "Curcumin",
  "Ferulic Acid",
  "Glycyrrhetinic Acid",
  "Liquiritigenin"
)

ALL_DRUGS <- c(ANCHOR_DRUG, PARTNER_DRUGS)

############################################################
## 1. Parameters
##    These values are kept unchanged from the original algorithm.
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

base_weights <- c(
  cosine = BASE_W_COSINE,
  jaccard = BASE_W_JACCARD,
  coverage = BASE_W_COVERAGE,
  offtarget = BASE_W_OFFTARGET
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
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(data.table)
library(ggplot2)
library(scales)
library(stringr)

############################################################
## 3. Plot style and utility functions
############################################################

PAL_SIZE <- c(
  "1" = "#8FB9A8",
  "2" = "#F2C879",
  "3" = "#B8A1D9",
  "4" = "#F0A7A0",
  "5" = "#9BB7D4"
)

PAL_SOURCE <- c(
  All_sources = "#8FB9A8",
  Reactome_ORA = "#9BB7D4",
  GO_BP_ORA = "#D7BCE8",
  MSigDB_Hallmark_ORA = "#F0C987"
)

PAL_DECOMP <- c(
  Weighted_CosineGain = "#9BB7D4",
  Weighted_JaccardGain = "#D7BCE8",
  Weighted_CoverageGain = "#8FB9A8",
  Weighted_OffTargetPenalty = "#F0A7A0"
)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 1,
        hjust = 0,
        color = "grey15",
        lineheight = 0.95,
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        size = base_size - 1,
        hjust = 0,
        color = "grey30",
        lineheight = 0.95,
        margin = margin(b = 8)
      ),
      axis.title = element_text(
        face = "bold",
        size = base_size,
        color = "grey20"
      ),
      axis.text = element_text(
        color = "grey20",
        size = base_size - 1
      ),
      legend.title = element_text(
        face = "bold",
        size = base_size - 1
      ),
      legend.text = element_text(
        size = base_size - 1,
        color = "grey20"
      ),
      legend.position = "right",
      panel.grid.major.y = element_line(
        color = "grey92",
        linewidth = 0.25
      ),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(12, 24, 16, 16)
    )
}

save_dual <- function(p, filename, width, height, fig_dir) {
  ggsave(
    file.path(fig_dir, paste0(filename, ".png")),
    p,
    width = width,
    height = height,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    file.path(fig_dir, paste0(filename, ".pdf")),
    p,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
}

safe_name <- function(x) {
  gsub("[/\\\\:*?\"<>| ]+", "_", x)
}

wrap_text <- function(x, width = 42) {
  stringr::str_wrap(x, width = width)
}

make_plot_title <- function(title, width = 72) {
  stringr::str_wrap(title, width = width)
}

make_dynamic_height <- function(n_labels, base = 4.8, per_label = 0.34, min_height = 5.8, max_height = 13.5) {
  pmin(max_height, pmax(min_height, base + per_label * n_labels))
}

make_dynamic_width <- function(labels, base = 8.8, per_char = 0.045, min_width = 9.8, max_width = 13.8) {
  max_chars <- max(nchar(as.character(labels)), na.rm = TRUE)
  pmin(max_width, pmax(min_width, base + per_char * max_chars))
}

get_disease_dirname <- function(disease_label) {
  out <- DISEASE_CONFIG[Disease_Label == disease_label, Disease_DirName]
  if (length(out) != 1 || is.na(out)) {
    stop("Cannot map disease label to input directory: ", disease_label)
  }
  out
}

############################################################
## 4. Core helper functions
##    The scoring logic is kept unchanged.
############################################################

cosine_similarity <- function(x, y) {
  idx <- intersect(names(x), names(y))
  x <- x[idx]
  y <- y[idx]
  
  den <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (den == 0) {
    return(NA_real_)
  }
  
  sum(x * y) / den
}

weighted_jaccard <- function(x, y) {
  idx <- intersect(names(x), names(y))
  x <- x[idx]
  y <- y[idx]
  
  den <- sum(pmax(x, y))
  if (den == 0) {
    return(NA_real_)
  }
  
  sum(pmin(x, y)) / den
}

make_space <- function(wide_dt, source_name, disease_label) {
  if (source_name != "All_sources") {
    x <- wide_dt[Pathway_Source == source_name]
  } else {
    x <- copy(wide_dt)
  }
  
  x <- x[!is.na(Disease) & is.finite(Disease)]
  
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
    disease_label = disease_label,
    dt = x,
    D = D,
    disease_core_ids = disease_core_ids,
    noncore_ids = noncore_ids,
    ordered_features = ordered_features
  )
}

make_partner_subsets <- function(partner_drugs) {
  out <- list()
  idx <- 1
  
  for (k in seq_along(partner_drugs)) {
    cmb <- combn(partner_drugs, k, simplify = FALSE)
    for (x in cmb) {
      out[[idx]] <- x
      idx <- idx + 1
    }
  }
  
  out
}

evaluate_real_combination <- function(space, core_vec, partner_set_vec, partner_set_name, partner_set_size, weights) {
  D <- space$D
  
  Core <- core_vec[names(D)]
  PartnerSet <- partner_set_vec[names(D)]
  
  Core[is.na(Core)] <- 0
  PartnerSet[is.na(PartnerSet)] <- 0
  
  names(Core) <- names(D)
  names(PartnerSet) <- names(D)
  
  combo <- Core + PartnerSet
  
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
  
  weighted_cosine <- weights["cosine"] * combination_gain_cosine
  weighted_jaccard_value <- weights["jaccard"] * combination_gain_jaccard
  weighted_coverage <- weights["coverage"] * coverage_gain
  weighted_offtarget_penalty <- -weights["offtarget"] * offtarget_gain
  
  raw_score <-
    weighted_cosine +
    weighted_jaccard_value +
    weighted_coverage +
    weighted_offtarget_penalty
  
  data.table(
    Disease = space$disease_label,
    Source = space$source_name,
    Anchor = ANCHOR_DRUG,
    Partner_Set = partner_set_name,
    Partner_Set_Size = partner_set_size,
    Combination = paste(ANCHOR_DRUG, partner_set_name, sep = " + "),
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
    Weighted_CosineGain = weighted_cosine,
    Weighted_JaccardGain = weighted_jaccard_value,
    Weighted_CoverageGain = weighted_coverage,
    Weighted_OffTargetPenalty = weighted_offtarget_penalty,
    Raw_FinalScore = raw_score,
    Active_Core_N = length(core_active),
    Active_Combination_N = length(combo_active),
    Disease_Core_N = length(disease_core_ids),
    Disease_NonCore_N = length(noncore_ids)
  )
}

############################################################
## 5. Input discovery
############################################################

find_latest_vector_dir <- function(base_dir) {
  if (!dir.exists(base_dir)) {
    return(NA_character_)
  }
  
  candidate_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  candidate_dirs <- candidate_dirs[
    grepl("PathwayVector", basename(candidate_dirs), ignore.case = TRUE)
  ]
  
  if (length(candidate_dirs) == 0) {
    return(NA_character_)
  }
  
  info <- file.info(candidate_dirs)
  candidate_dirs[order(info$mtime, decreasing = TRUE)][1]
}

locate_entity_input <- function(entity_type, entity_label, entity_dirname = entity_label) {
  base_dir <- if (entity_type == "Disease") {
    file.path(ROOT_DIR, entity_dirname)
  } else {
    file.path(ROOT_DIR, DRUG_ROOT_DIRNAME, entity_dirname)
  }
  
  latest_dir <- find_latest_vector_dir(base_dir)
  
  if (is.na(latest_dir)) {
    stop(
      "No PathwayVector folder found for entity: ",
      entity_label,
      " under: ",
      base_dir
    )
  }
  
  model_dir <- file.path(latest_dir, "02_ModelInput_PathwayVectors")
  
  long_candidates <- c(
    file.path(model_dir, "04_Entity_PathwayVector_ForModelInput_LONG.tsv"),
    file.path(model_dir, "04_Disease_PathwayVector_ForModelInput_LONG.tsv")
  )
  
  wide_candidates <- c(
    file.path(model_dir, "05_Entity_PathwayVector_PathwayScore_WIDE.tsv"),
    file.path(model_dir, "05_Disease_PathwayVector_PathwayScore_WIDE.tsv")
  )
  
  long_file <- long_candidates[file.exists(long_candidates)][1]
  wide_file <- wide_candidates[file.exists(wide_candidates)][1]
  
  if (is.na(long_file) || !file.exists(long_file)) {
    stop("LONG input file not found for entity: ", entity_label)
  }
  
  if (is.na(wide_file) || !file.exists(wide_file)) {
    stop("WIDE input file not found for entity: ", entity_label)
  }
  
  data.table(
    Entity_Type = entity_type,
    Entity_Label = entity_label,
    Entity_DirName = entity_dirname,
    Base_Dir = base_dir,
    Latest_PathwayVector_Dir = latest_dir,
    ModelInput_Dir = model_dir,
    LONG_File = long_file,
    WIDE_File = wide_file
  )
}

input_map <- rbindlist(
  c(
    lapply(seq_len(nrow(DISEASE_CONFIG)), function(i) {
      locate_entity_input(
        entity_type = "Disease",
        entity_label = DISEASE_CONFIG$Disease_Label[i],
        entity_dirname = DISEASE_CONFIG$Disease_DirName[i]
      )
    }),
    lapply(ALL_DRUGS, function(x) {
      locate_entity_input(
        entity_type = "Drug",
        entity_label = x,
        entity_dirname = x
      )
    })
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  input_map,
  file.path(OUT_ROOT, "00_Input_File_Map.tsv"),
  sep = "\t"
)

############################################################
## 6. Read pathway vectors
############################################################

read_entity_wide <- function(file, entity_label, entity_dirname = entity_label) {
  x <- fread(file)
  
  required_meta <- c(
    "Feature_ID",
    "Feature_Name",
    "Pathway_Source",
    "Pathway_ID",
    "Pathway_Name"
  )
  
  miss <- setdiff(required_meta, colnames(x))
  
  if (length(miss) > 0) {
    stop("Missing columns in ", file, ":\n", paste(miss, collapse = "\n"))
  }
  
  value_candidates <- unique(c(
    entity_label,
    entity_dirname,
    "Disease",
    "Value"
  ))
  
  value_col <- value_candidates[value_candidates %in% colnames(x)][1]
  
  if (is.na(value_col)) {
    non_meta <- setdiff(colnames(x), required_meta)
    numeric_cols <- non_meta[sapply(x[, ..non_meta], is.numeric)]
    
    if (length(numeric_cols) == 0) {
      stop("Cannot identify numeric pathway-score column in: ", file)
    }
    
    value_col <- numeric_cols[1]
  }
  
  out <- x[, .(
    Feature_ID,
    Feature_Name,
    Pathway_Source,
    Pathway_ID,
    Pathway_Name,
    Value = as.numeric(get(value_col))
  )]
  
  out <- out[!is.na(Value) & is.finite(Value)]
  out <- unique(out, by = "Feature_ID")
  out
}

disease_vectors <- list()
drug_vectors <- list()

for (i in seq_len(nrow(DISEASE_CONFIG))) {
  disease_label <- DISEASE_CONFIG$Disease_Label[i]
  disease_dirname <- DISEASE_CONFIG$Disease_DirName[i]
  file <- input_map[Entity_Type == "Disease" & Entity_Label == disease_label]$WIDE_File
  disease_vectors[[disease_label]] <- read_entity_wide(file, disease_label, disease_dirname)
}

for (drug in ALL_DRUGS) {
  file <- input_map[Entity_Type == "Drug" & Entity_Label == drug]$WIDE_File
  drug_vectors[[drug]] <- read_entity_wide(file, drug, drug)
}

############################################################
## 7. Vector alignment
############################################################

build_disease_drug_space <- function(disease_label) {
  disease_dt <- copy(disease_vectors[[disease_label]])
  setnames(disease_dt, "Value", "Disease")
  
  out <- disease_dt
  
  for (drug in ALL_DRUGS) {
    drug_dt <- copy(drug_vectors[[drug]])[, .(Feature_ID, Drug_Value = Value)]
    setnames(drug_dt, "Drug_Value", safe_name(drug))
    out <- merge(out, drug_dt, by = "Feature_ID", all.x = TRUE)
  }
  
  drug_cols <- safe_name(ALL_DRUGS)
  
  for (cc in drug_cols) {
    out[is.na(get(cc)), (cc) := 0]
  }
  
  out
}

extract_vector <- function(aligned_dt, drug_name) {
  v <- aligned_dt[[safe_name(drug_name)]]
  names(v) <- aligned_dt$Feature_ID
  v[is.na(v)] <- 0
  v
}

make_partner_set_vector <- function(aligned_dt, partner_set) {
  if (length(partner_set) == 0) {
    v <- rep(0, nrow(aligned_dt))
    names(v) <- aligned_dt$Feature_ID
    return(v)
  }
  
  cols <- safe_name(partner_set)
  mat <- as.matrix(aligned_dt[, ..cols])
  mat[is.na(mat)] <- 0
  
  v <- rowSums(mat)
  names(v) <- aligned_dt$Feature_ID
  v
}

############################################################
## 8. Main disease-specific full-combination analysis
############################################################

partner_sets <- make_partner_subsets(PARTNER_DRUGS)

all_disease_results <- list()
all_source_results <- list()
all_best_summary <- list()
all_added_pathways <- list()
all_driver_pathways <- list()

for (disease_label in DISEASE_CONFIG$Disease_Label) {
  disease_out <- file.path(OUT_ROOT, safe_name(disease_label))
  table_dir <- file.path(disease_out, "01_Tables")
  fig_dir <- file.path(disease_out, "02_Figures")
  vector_dir <- file.path(disease_out, "03_Vectors")
  qc_dir <- file.path(disease_out, "04_QC")
  
  dir.create(disease_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(vector_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  
  aligned_dt <- build_disease_drug_space(disease_label)
  
  fwrite(
    aligned_dt,
    file.path(vector_dir, "00_Aligned_DiseaseDrug_PathwaySpace.tsv"),
    sep = "\t"
  )
  
  core_vec_all <- extract_vector(aligned_dt, ANCHOR_DRUG)
  main_space <- make_space(aligned_dt, "All_sources", disease_label)
  
  ############################################################
  ## 8.1 All-source full 31-combination ranking
  ############################################################
  
  main_results <- rbindlist(
    lapply(partner_sets, function(pset) {
      partner_set_name <- paste(pset, collapse = " + ")
      partner_set_vec <- make_partner_set_vector(aligned_dt, pset)
      
      evaluate_real_combination(
        space = main_space,
        core_vec = core_vec_all,
        partner_set_vec = partner_set_vec,
        partner_set_name = partner_set_name,
        partner_set_size = length(pset),
        weights = base_weights
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
  
  main_results[, Rank_FinalScore := frank(-Raw_FinalScore, ties.method = "average")]
  
  main_results[
    ,
    Rank_Within_Size := frank(-Raw_FinalScore, ties.method = "average"),
    by = Partner_Set_Size
  ]
  
  setorder(main_results, Rank_FinalScore, -Raw_FinalScore)
  
  ############################################################
  ## 8.2 Source-stratified full-combination scores
  ############################################################
  
  source_results <- list()
  
  for (src in SOURCE_LEVELS) {
    tmp <- tryCatch({
      sp <- make_space(aligned_dt, src, disease_label)
      core_vec_src <- extract_vector(sp$dt, ANCHOR_DRUG)
      
      rbindlist(
        lapply(partner_sets, function(pset) {
          partner_set_name <- paste(pset, collapse = " + ")
          partner_set_vec <- make_partner_set_vector(sp$dt, pset)
          
          evaluate_real_combination(
            space = sp,
            core_vec = core_vec_src,
            partner_set_vec = partner_set_vec,
            partner_set_name = partner_set_name,
            partner_set_size = length(pset),
            weights = base_weights
          )
        }),
        use.names = TRUE,
        fill = TRUE
      )
    }, error = function(e) {
      data.table(
        Disease = disease_label,
        Source = src,
        Anchor = ANCHOR_DRUG,
        Partner_Set = sapply(partner_sets, paste, collapse = " + "),
        Partner_Set_Size = sapply(partner_sets, length),
        Combination = paste(ANCHOR_DRUG, sapply(partner_sets, paste, collapse = " + "), sep = " + "),
        Core_Cosine = NA_real_,
        Combo_Cosine = NA_real_,
        CombinationGain_Cosine = NA_real_,
        Core_Jaccard = NA_real_,
        Combo_Jaccard = NA_real_,
        CombinationGain_Jaccard = NA_real_,
        CoreCoverage = NA_real_,
        ComboCoverage = NA_real_,
        CoverageGain = NA_real_,
        CoreOffTarget = NA_real_,
        ComboOffTarget = NA_real_,
        OffTargetGain = NA_real_,
        MissedCoreRatio = NA_real_,
        Weighted_CosineGain = NA_real_,
        Weighted_JaccardGain = NA_real_,
        Weighted_CoverageGain = NA_real_,
        Weighted_OffTargetPenalty = NA_real_,
        Raw_FinalScore = NA_real_,
        Active_Core_N = NA_integer_,
        Active_Combination_N = NA_integer_,
        Disease_Core_N = NA_integer_,
        Disease_NonCore_N = NA_integer_
      )
    })
    
    source_results[[src]] <- tmp
  }
  
  source_results <- rbindlist(source_results, use.names = TRUE, fill = TRUE)
  
  source_results[
    ,
    Rank_FinalScore := frank(-Raw_FinalScore, ties.method = "average"),
    by = Source
  ]
  
  source_results[
    ,
    Rank_Within_Size := frank(-Raw_FinalScore, ties.method = "average"),
    by = .(Source, Partner_Set_Size)
  ]
  
  ############################################################
  ## 8.3 Best combination and disease-core pathway explanation
  ############################################################
  
  best_row <- main_results[Rank_FinalScore == 1][1]
  best_partner_set <- unlist(strsplit(best_row$Partner_Set, " \\+ "))
  
  disease_core_ids <- main_space$disease_core_ids
  
  best_partner_vec <- make_partner_set_vector(aligned_dt, best_partner_set)
  anchor_vec <- core_vec_all
  combo_vec <- anchor_vec + best_partner_vec
  
  disease_vec <- aligned_dt$Disease
  names(disease_vec) <- aligned_dt$Feature_ID
  
  pathway_explain <- copy(aligned_dt)
  
  pathway_explain[, Disease_Core := Feature_ID %in% disease_core_ids]
  pathway_explain[, Anchor_Value := anchor_vec[Feature_ID]]
  pathway_explain[, Partner_Set_Value := best_partner_vec[Feature_ID]]
  pathway_explain[, Combination_Value := combo_vec[Feature_ID]]
  pathway_explain[, Disease_Value := disease_vec[Feature_ID]]
  pathway_explain[, Added_Core_Pathway := Disease_Core & Anchor_Value <= 0 & Combination_Value > 0]
  pathway_explain[, Enhanced_Core_Pathway := Disease_Core & Anchor_Value > 0 & Partner_Set_Value > 0]
  pathway_explain[, Coverage_Ratio := Combination_Value / pmax(Disease_Value, 1e-12)]
  pathway_explain[, Coverage_Ratio_Capped := pmin(Coverage_Ratio, 2)]
  pathway_explain[, PartnerSpecificContribution := Partner_Set_Value * Disease_Value]
  pathway_explain[, CombinationContribution := Combination_Value * Disease_Value]
  pathway_explain[, Best_Combination := best_row$Combination]
  pathway_explain[, Best_Partner_Set := best_row$Partner_Set]
  pathway_explain[, Best_Partner_Set_Size := best_row$Partner_Set_Size]
  
  added_core_pathways <- pathway_explain[
    Added_Core_Pathway == TRUE
  ][
    order(-PartnerSpecificContribution, -Disease_Value)
  ]
  
  driver_pathways <- pathway_explain[
    Disease_Core == TRUE & Combination_Value > 0
  ][
    order(-PartnerSpecificContribution, -CombinationContribution, -Disease_Value)
  ]
  
  best_summary <- data.table(
    Disease = disease_label,
    Best_Combination = best_row$Combination,
    Best_Partner_Set = best_row$Partner_Set,
    Best_Partner_Set_Size = best_row$Partner_Set_Size,
    Best_FinalScore = best_row$Raw_FinalScore,
    Best_Rank = best_row$Rank_FinalScore,
    Combo_Cosine = best_row$Combo_Cosine,
    Combo_Jaccard = best_row$Combo_Jaccard,
    ComboCoverage = best_row$ComboCoverage,
    ComboOffTarget = best_row$ComboOffTarget,
    CoverageGain = best_row$CoverageGain,
    OffTargetGain = best_row$OffTargetGain,
    Weighted_CosineGain = best_row$Weighted_CosineGain,
    Weighted_JaccardGain = best_row$Weighted_JaccardGain,
    Weighted_CoverageGain = best_row$Weighted_CoverageGain,
    Weighted_OffTargetPenalty = best_row$Weighted_OffTargetPenalty,
    Added_Core_Pathway_N = nrow(added_core_pathways),
    Driver_Core_Pathway_N = nrow(driver_pathways),
    Top_Added_Core_Pathways = paste(head(added_core_pathways$Pathway_Name, 10), collapse = " | "),
    Top_Driver_Core_Pathways = paste(head(driver_pathways$Pathway_Name, 10), collapse = " | ")
  )
  
  ############################################################
  ## 8.4 Per-combination pathway explanation table
  ############################################################
  
  combo_pathway_list <- list()
  
  for (i in seq_along(partner_sets)) {
    pset <- partner_sets[[i]]
    partner_set_name <- paste(pset, collapse = " + ")
    combination_name <- paste(ANCHOR_DRUG, partner_set_name, sep = " + ")
    partner_set_vec <- make_partner_set_vector(aligned_dt, pset)
    combo_vec_i <- core_vec_all + partner_set_vec
    
    tmp <- copy(aligned_dt)
    tmp[, Disease := disease_label]
    tmp[, Disease_Core := Feature_ID %in% disease_core_ids]
    tmp[, Anchor_Value := core_vec_all[Feature_ID]]
    tmp[, Partner_Set_Value := partner_set_vec[Feature_ID]]
    tmp[, Combination_Value := combo_vec_i[Feature_ID]]
    tmp[, Disease_Value := disease_vec[Feature_ID]]
    tmp[, Added_Core_Pathway := Disease_Core & Anchor_Value <= 0 & Combination_Value > 0]
    tmp[, Enhanced_Core_Pathway := Disease_Core & Anchor_Value > 0 & Partner_Set_Value > 0]
    tmp[, Coverage_Ratio := Combination_Value / pmax(Disease_Value, 1e-12)]
    tmp[, Coverage_Ratio_Capped := pmin(Coverage_Ratio, 2)]
    tmp[, PartnerSpecificContribution := Partner_Set_Value * Disease_Value]
    tmp[, CombinationContribution := Combination_Value * Disease_Value]
    tmp[, Partner_Set := partner_set_name]
    tmp[, Partner_Set_Size := length(pset)]
    tmp[, Combination := combination_name]
    
    combo_pathway_list[[i]] <- tmp
  }
  
  combo_pathway_table <- rbindlist(combo_pathway_list, use.names = TRUE, fill = TRUE)
  
  combo_added_summary <- combo_pathway_table[
    Disease_Core == TRUE,
    .(
      Added_Core_Pathway_N = sum(Added_Core_Pathway, na.rm = TRUE),
      Enhanced_Core_Pathway_N = sum(Enhanced_Core_Pathway, na.rm = TRUE),
      Sum_PartnerSpecificContribution = sum(PartnerSpecificContribution, na.rm = TRUE),
      Mean_PartnerSpecificContribution = mean(PartnerSpecificContribution[PartnerSpecificContribution > 0], na.rm = TRUE),
      Sum_CombinationContribution = sum(CombinationContribution, na.rm = TRUE),
      Mean_Coverage_Ratio = mean(Coverage_Ratio, na.rm = TRUE),
      Top_Added_Core_Pathways = paste(
        head(Pathway_Name[Added_Core_Pathway == TRUE][order(-PartnerSpecificContribution[Added_Core_Pathway == TRUE])], 8),
        collapse = " | "
      )
    ),
    by = .(Disease, Partner_Set, Partner_Set_Size, Combination)
  ]
  
  ############################################################
  ## 8.5 Save disease-level tables
  ############################################################
  
  fwrite(
    main_results,
    file.path(table_dir, "01_FullCombination_Ranking_31.tsv"),
    sep = "\t"
  )
  
  fwrite(
    best_summary,
    file.path(table_dir, "02_BestCombination_Summary.tsv"),
    sep = "\t"
  )
  
  fwrite(
    source_results,
    file.path(table_dir, "03_SourceStratified_FullCombination_Scores.tsv"),
    sep = "\t"
  )
  
  fwrite(
    added_core_pathways,
    file.path(table_dir, "04_BestCombination_AddedCorePathways.tsv"),
    sep = "\t"
  )
  
  fwrite(
    driver_pathways,
    file.path(table_dir, "05_BestCombination_TopDriverCorePathways.tsv"),
    sep = "\t"
  )
  
  fwrite(
    combo_added_summary,
    file.path(table_dir, "06_AllCombination_AddedCorePathway_Summary.tsv"),
    sep = "\t"
  )
  
  fwrite(
    combo_pathway_table,
    file.path(table_dir, "07_AllCombination_PathwayContribution_Long.tsv"),
    sep = "\t"
  )
  
  all_disease_results[[disease_label]] <- main_results
  all_source_results[[disease_label]] <- source_results
  all_best_summary[[disease_label]] <- best_summary
  all_added_pathways[[disease_label]] <- added_core_pathways
  all_driver_pathways[[disease_label]] <- driver_pathways
  
  ############################################################
  ## 8.6 Figures for each disease
  ############################################################
  
  size_best <- main_results[
    ,
    .SD[which.max(Raw_FinalScore)],
    by = Partner_Set_Size
  ]
  
  p_landscape <- ggplot(
    main_results,
    aes(
      x = factor(Partner_Set_Size),
      y = Raw_FinalScore,
      color = factor(Partner_Set_Size)
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey72",
      linewidth = 0.4
    ) +
    geom_jitter(
      width = 0.14,
      height = 0,
      size = 2.6,
      alpha = 0.72
    ) +
    geom_point(
      data = size_best,
      aes(
        x = factor(Partner_Set_Size),
        y = Raw_FinalScore
      ),
      size = 4.3,
      shape = 18,
      color = "grey15"
    ) +
    scale_color_manual(values = PAL_SIZE, guide = "none") +
    labs(
      title = disease_label,
      subtitle = "Baicalein-based combination landscape",
      x = "Number of partner compounds added to baicalein",
      y = "Final combination score"
    ) +
    theme_pub(11) +
    theme(plot.margin = margin(14, 24, 16, 16))
  
  save_dual(
    p_landscape,
    "Fig1_CombinationSpace_Landscape",
    7.4,
    5.2,
    fig_dir
  )
  
  top_rank_n <- 15
  
  rank_dt <- copy(main_results)[seq_len(min(top_rank_n, .N))]
  rank_dt[, Combination_Label := wrap_text(Combination, width = 46)]
  rank_dt[, Combination_Label := factor(Combination_Label, levels = rev(Combination_Label))]
  rank_dt[, Size_Factor := factor(Partner_Set_Size)]
  
  fig2_height <- make_dynamic_height(
    n_labels = nrow(rank_dt),
    base = 2.8,
    per_label = 0.38,
    min_height = 6.2,
    max_height = 12.5
  )
  
  p_rank <- ggplot(
    rank_dt,
    aes(
      x = Raw_FinalScore,
      y = Combination_Label,
      color = Size_Factor
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey72",
      linewidth = 0.4
    ) +
    geom_segment(
      aes(
        x = 0,
        xend = Raw_FinalScore,
        y = Combination_Label,
        yend = Combination_Label
      ),
      linewidth = 0.95,
      alpha = 0.72
    ) +
    geom_point(size = 3.8) +
    scale_color_manual(values = PAL_SIZE, name = "Partner\nnumber") +
    labs(
      title = disease_label,
      subtitle = "Top-ranked baicalein-based combinations",
      x = "Final combination score",
      y = NULL
    ) +
    theme_pub(10) +
    theme(
      axis.text.y = element_text(size = 7.8, lineheight = 0.90),
      legend.position = "right",
      plot.margin = margin(14, 26, 16, 16)
    )
  
  save_dual(
    p_rank,
    "Fig2_TopRanked_Combinations",
    9.8,
    fig2_height,
    fig_dir
  )
  
  decomp_dt <- copy(rank_dt)[
    ,
    .(
      Disease,
      Combination,
      Combination_Label,
      Partner_Set_Size,
      Weighted_CosineGain,
      Weighted_JaccardGain,
      Weighted_CoverageGain,
      Weighted_OffTargetPenalty
    )
  ]
  
  decomp_long <- melt(
    decomp_dt,
    id.vars = c("Disease", "Combination", "Combination_Label", "Partner_Set_Size"),
    measure.vars = c(
      "Weighted_CosineGain",
      "Weighted_JaccardGain",
      "Weighted_CoverageGain",
      "Weighted_OffTargetPenalty"
    ),
    variable.name = "Score_Component",
    value.name = "Weighted_Value"
  )
  
  decomp_long[
    ,
    Score_Component := factor(
      Score_Component,
      levels = c(
        "Weighted_CosineGain",
        "Weighted_JaccardGain",
        "Weighted_CoverageGain",
        "Weighted_OffTargetPenalty"
      )
    )
  ]
  
  p_decomp <- ggplot(
    decomp_long,
    aes(
      x = Weighted_Value,
      y = Combination_Label,
      fill = Score_Component
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey72",
      linewidth = 0.4
    ) +
    geom_col(width = 0.68, color = "white", linewidth = 0.15) +
    scale_fill_manual(
      values = PAL_DECOMP,
      labels = c(
        "Cosine gain",
        "Jaccard gain",
        "Coverage gain",
        "Off-target penalty"
      ),
      name = NULL
    ) +
    labs(
      title = disease_label,
      subtitle = "Score decomposition of leading combinations",
      x = "Weighted contribution to final score",
      y = NULL
    ) +
    theme_pub(10) +
    theme(
      axis.text.y = element_text(size = 7.6, lineheight = 0.90),
      legend.position = "bottom",
      plot.margin = margin(14, 28, 18, 16)
    )
  
  save_dual(
    p_decomp,
    "Fig3_ScoreDecomposition_LeadingCombinations",
    10.5,
    fig2_height,
    fig_dir
  )
  
  p_balance <- ggplot(
    main_results,
    aes(
      x = ComboCoverage,
      y = ComboOffTarget,
      size = Raw_FinalScore,
      color = factor(Partner_Set_Size)
    )
  ) +
    geom_point(alpha = 0.76) +
    geom_point(
      data = best_row,
      aes(x = ComboCoverage, y = ComboOffTarget),
      inherit.aes = FALSE,
      shape = 21,
      size = 6.2,
      stroke = 1.1,
      fill = "white",
      color = "grey15"
    ) +
    geom_text(
      data = best_row,
      aes(x = ComboCoverage, y = ComboOffTarget, label = "Best"),
      inherit.aes = FALSE,
      size = 3.2,
      vjust = -1.15,
      color = "grey15"
    ) +
    scale_color_manual(values = PAL_SIZE, name = "Partner\nnumber") +
    scale_size_continuous(range = c(2.5, 7.5), name = "Final\nscore") +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    coord_cartesian(clip = "off") +
    labs(
      title = disease_label,
      subtitle = "Disease-core coverage and non-core expansion",
      x = "Disease-core pathway coverage",
      y = "Non-core pathway activity"
    ) +
    theme_pub(11) +
    theme(plot.margin = margin(14, 24, 16, 16))
  
  save_dual(
    p_balance,
    "Fig4_CoreCoverage_NonCoreExpansion",
    7.8,
    5.8,
    fig_dir
  )
  
  heat_n <- 25
  heat_dt <- copy(driver_pathways)[seq_len(min(heat_n, .N))]
  
  if (nrow(heat_dt) > 0) {
    heat_plot <- heat_dt[
      ,
      .(
        Feature_ID,
        Pathway_Source,
        Pathway_Name,
        Disease_Value,
        Baicalein = Anchor_Value,
        Partner_Set = Partner_Set_Value,
        Combination = Combination_Value
      )
    ]
    
    heat_long <- melt(
      heat_plot,
      id.vars = c("Feature_ID", "Pathway_Source", "Pathway_Name"),
      measure.vars = c("Disease_Value", "Baicalein", "Partner_Set", "Combination"),
      variable.name = "Vector",
      value.name = "PathwayScore"
    )
    
    heat_long[
      ,
      Vector := factor(
        Vector,
        levels = c("Disease_Value", "Baicalein", "Partner_Set", "Combination"),
        labels = c("Disease", "Baicalein", "Partner set", "Combination")
      )
    ]
    
    heat_long[, Pathway_Label := wrap_text(Pathway_Name, width = 42)]
    heat_long[, Pathway_Label := factor(Pathway_Label, levels = rev(unique(Pathway_Label)))]
    
    heat_height <- make_dynamic_height(
      n_labels = length(unique(heat_long$Pathway_Label)),
      base = 2.8,
      per_label = 0.34,
      min_height = 7.4,
      max_height = 14.2
    )
    
    heat_width <- make_dynamic_width(
      labels = heat_long$Pathway_Label,
      base = 8.6,
      per_char = 0.055,
      min_width = 11.0,
      max_width = 14.2
    )
    
    p_driver_heat <- ggplot(
      heat_long,
      aes(
        x = Vector,
        y = Pathway_Label,
        fill = PathwayScore
      )
    ) +
      geom_tile(color = "white", linewidth = 0.22) +
      scale_fill_gradient(
        low = "#F7FBFF",
        high = "#8FB9A8",
        name = "Pathway\nscore"
      ) +
      labs(
        title = disease_label,
        subtitle = "Disease-core pathways driving the top-ranked combination",
        x = "Comparison",
        y = "Pathway"
      ) +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 12,
          hjust = 0,
          color = "grey15",
          lineheight = 0.95,
          margin = margin(b = 3)
        ),
        plot.subtitle = element_text(
          size = 9.5,
          hjust = 0,
          color = "grey30",
          lineheight = 0.95,
          margin = margin(b = 8)
        ),
        axis.title.x = element_text(
          face = "bold",
          size = 9.5,
          color = "grey20",
          margin = margin(t = 7)
        ),
        axis.title.y = element_text(
          face = "bold",
          size = 9.5,
          color = "grey20",
          margin = margin(r = 7)
        ),
        axis.text.x = element_text(
          angle = 0,
          hjust = 0.5,
          size = 9,
          color = "grey20"
        ),
        axis.text.y = element_text(
          size = 7.2,
          color = "grey20",
          lineheight = 0.90
        ),
        axis.ticks = element_blank(),
        legend.title = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 8),
        legend.position = "right",
        plot.margin = margin(18, 56, 22, 20)
      )
    
    save_dual(
      p_driver_heat,
      "Fig5_TopCombination_DriverPathwayHeatmap",
      heat_width,
      heat_height,
      fig_dir
    )
  }
  
  source_plot_dt <- copy(source_results)
  source_plot_dt <- source_plot_dt[
    Source != "All_sources" & Combination %in% rank_dt$Combination
  ]
  
  source_plot_dt[, Combination_Label := wrap_text(Combination, width = 46)]
  source_plot_dt[, Combination_Label := factor(Combination_Label, levels = rev(rank_dt$Combination_Label))]
  source_plot_dt[, Source := factor(Source, levels = c("Reactome_ORA", "GO_BP_ORA", "MSigDB_Hallmark_ORA"))]
  
  p_source <- ggplot(
    source_plot_dt,
    aes(
      x = Raw_FinalScore,
      y = Combination_Label,
      fill = Source
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey72",
      linewidth = 0.4
    ) +
    geom_col(
      position = position_dodge(width = 0.76),
      width = 0.66,
      color = "white",
      linewidth = 0.15
    ) +
    scale_fill_manual(values = PAL_SOURCE, name = NULL) +
    labs(
      title = disease_label,
      subtitle = "Source-stratified support for leading combinations",
      x = "Final combination score",
      y = NULL
    ) +
    theme_pub(9) +
    theme(
      axis.text.y = element_text(size = 6.8, lineheight = 0.88),
      legend.position = "bottom",
      plot.margin = margin(14, 26, 18, 16)
    )
  
  save_dual(
    p_source,
    "Fig6_SourceStratified_LeadingCombinations",
    11.2,
    fig2_height,
    fig_dir
  )
}

############################################################
## 9. Master outputs
############################################################

master_main <- rbindlist(all_disease_results, use.names = TRUE, fill = TRUE)
master_source <- rbindlist(all_source_results, use.names = TRUE, fill = TRUE)
master_best <- rbindlist(all_best_summary, use.names = TRUE, fill = TRUE)
master_added <- rbindlist(all_added_pathways, use.names = TRUE, fill = TRUE)
master_drivers <- rbindlist(all_driver_pathways, use.names = TRUE, fill = TRUE)

master_main[
  ,
  Disease_Rank := frank(-Raw_FinalScore, ties.method = "average"),
  by = Disease
]

master_main[
  ,
  Disease_Size_Rank := frank(-Raw_FinalScore, ties.method = "average"),
  by = .(Disease, Partner_Set_Size)
]

setorder(master_main, Disease, Disease_Rank)

fwrite(
  master_main,
  file.path(OUT_ROOT, "01_MASTER_FullCombination_Ranking_31_perDisease.tsv"),
  sep = "\t"
)

fwrite(
  master_best,
  file.path(OUT_ROOT, "02_MASTER_BestCombination_Summary_perDisease.tsv"),
  sep = "\t"
)

fwrite(
  master_source,
  file.path(OUT_ROOT, "03_MASTER_SourceStratified_FullCombination_Scores.tsv"),
  sep = "\t"
)

fwrite(
  master_added,
  file.path(OUT_ROOT, "04_MASTER_BestCombination_AddedCorePathways.tsv"),
  sep = "\t"
)

fwrite(
  master_drivers,
  file.path(OUT_ROOT, "05_MASTER_BestCombination_DriverCorePathways.tsv"),
  sep = "\t"
)

run_metadata <- data.table(
  Parameter = c(
    "CORE_RATIO",
    "IDEAL_RATIO",
    "PARTIAL_RATIO",
    "REDUNDANT_OVERLAP",
    "DISEASE_CORE_Q",
    "BASE_W_COSINE",
    "BASE_W_JACCARD",
    "BASE_W_COVERAGE",
    "BASE_W_OFFTARGET",
    "N_REPEAT",
    "N_STRESS_REPEAT",
    "N_RANDOM",
    "N_WEIGHT",
    "Anchor_Drug",
    "Partner_Drugs",
    "Partner_Set_Total_N",
    "Disease_N"
  ),
  Value = c(
    CORE_RATIO,
    IDEAL_RATIO,
    PARTIAL_RATIO,
    REDUNDANT_OVERLAP,
    DISEASE_CORE_Q,
    BASE_W_COSINE,
    BASE_W_JACCARD,
    BASE_W_COVERAGE,
    BASE_W_OFFTARGET,
    N_REPEAT,
    N_STRESS_REPEAT,
    N_RANDOM,
    N_WEIGHT,
    ANCHOR_DRUG,
    paste(PARTNER_DRUGS, collapse = "; "),
    length(partner_sets),
    nrow(DISEASE_CONFIG)
  )
)

fwrite(
  run_metadata,
  file.path(OUT_ROOT, "06_RunMetadata_Parameters_Unchanged.tsv"),
  sep = "\t"
)

############################################################
## 10. Master index figures
##     These figures summarize disease-specific results only.
##     They do not imply cross-disease mechanistic comparability.
############################################################

MASTER_FIG_DIR <- file.path(OUT_ROOT, "00_MASTER_Figures")
dir.create(MASTER_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

master_best_plot <- copy(master_best)
master_best_plot[, Disease := factor(Disease, levels = DISEASE_CONFIG$Disease_Label)]
master_best_plot[, Best_Combination_Label := wrap_text(Best_Combination, width = 46)]

p_master_best <- ggplot(
  master_best_plot,
  aes(
    x = Best_FinalScore,
    y = Disease,
    color = factor(Best_Partner_Set_Size)
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Best_FinalScore,
      y = Disease,
      yend = Disease
    ),
    linewidth = 1.0,
    alpha = 0.74
  ) +
  geom_point(size = 4.4) +
  geom_text(
    aes(label = paste0("n=", Best_Partner_Set_Size)),
    hjust = -0.25,
    size = 3.2,
    color = "grey20"
  ) +
  scale_color_manual(values = PAL_SIZE, name = "Partner\nnumber") +
  labs(
    title = "Best baicalein-based combinations",
    subtitle = "Index summary within each disease",
    x = "Final combination score",
    y = NULL
  ) +
  theme_pub(11) +
  theme(
    legend.position = "right",
    plot.margin = margin(14, 30, 16, 16)
  )

save_dual(
  p_master_best,
  "MASTER_BestCombination_PerDisease_Index",
  8.2,
  4.8,
  MASTER_FIG_DIR
)

############################################################
## 11. Final report
############################################################

cat("\n============================================================\n")
cat("Real baicalein-based full-combination analysis completed.\n")
cat("============================================================\n")

cat("\nOutput root:\n")
cat(OUT_ROOT, "\n")

cat("\nMain master tables:\n")
cat(file.path(OUT_ROOT, "01_MASTER_FullCombination_Ranking_31_perDisease.tsv"), "\n")
cat(file.path(OUT_ROOT, "02_MASTER_BestCombination_Summary_perDisease.tsv"), "\n")
cat(file.path(OUT_ROOT, "04_MASTER_BestCombination_AddedCorePathways.tsv"), "\n")
cat(file.path(OUT_ROOT, "05_MASTER_BestCombination_DriverCorePathways.tsv"), "\n")

cat("\nMaster figures:\n")
cat(MASTER_FIG_DIR, "\n")

cat("\nBest-ranked combination per disease:\n")
print(
  master_best[
    ,
    .(
      Disease,
      Best_Combination,
      Best_Partner_Set_Size,
      Best_FinalScore,
      Combo_Cosine,
      Combo_Jaccard,
      ComboCoverage,
      ComboOffTarget,
      CoverageGain,
      OffTargetGain,
      Added_Core_Pathway_N,
      Driver_Core_Pathway_N
    )
  ]
)

cat("============================================================\n")
