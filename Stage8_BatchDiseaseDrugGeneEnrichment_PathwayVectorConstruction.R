############################################################
## Batch pathway enrichment:
## 3 diseases + 6 drugs
## Only reading/output paths changed; enrichment methods unchanged
############################################################

rm(list = ls())

############################################################
## 0. Batch input settings
############################################################

ROOT_DIR <- "/media/desk16/iy15915/中药之开创/富集模型"

INPUT_JOBS <- data.frame(
  Entity_Type = c(
    rep("Disease", 3),
    rep("Drug", 6)
  ),
  Entity_Name = c(
    "类风湿关节炎",
    "糖尿病",
    "骨关节炎",
    "Baicalein",
    "Caffeic Acid",
    "Curcumin",
    "Ferulic Acid",
    "Glycyrrhetinic Acid",
    "Liquiritigenin"
  ),
  Base_Dir = c(
    file.path(ROOT_DIR, "类风湿关节炎"),
    file.path(ROOT_DIR, "糖尿病"),
    file.path(ROOT_DIR, "骨关节炎"),
    file.path(ROOT_DIR, "药物基因", "Baicalein"),
    file.path(ROOT_DIR, "药物基因", "Caffeic Acid"),
    file.path(ROOT_DIR, "药物基因", "Curcumin"),
    file.path(ROOT_DIR, "药物基因", "Ferulic Acid"),
    file.path(ROOT_DIR, "药物基因", "Glycyrrhetinic Acid"),
    file.path(ROOT_DIR, "药物基因", "Liquiritigenin")
  ),
  Input_File = c(
    file.path(ROOT_DIR, "类风湿关节炎", "gene.txt"),
    file.path(ROOT_DIR, "糖尿病", "gene.txt"),
    file.path(ROOT_DIR, "骨关节炎", "gene.txt"),
    file.path(ROOT_DIR, "药物基因", "Baicalein", "gene.txt"),
    file.path(ROOT_DIR, "药物基因", "Caffeic Acid", "gene.txt"),
    file.path(ROOT_DIR, "药物基因", "Curcumin", "gene.txt"),
    file.path(ROOT_DIR, "药物基因", "Ferulic Acid", "gene.txt"),
    file.path(ROOT_DIR, "药物基因", "Glycyrrhetinic Acid", "gene.txt"),
    file.path(ROOT_DIR, "药物基因", "Liquiritigenin", "gene.txt")
  ),
  stringsAsFactors = FALSE
)

RUN_TAG <- format(Sys.time(), "%Y%m%d_%H%M%S")

FDR_SIG <- 0.05
FDR_CAP <- 1e-300
MIN_MAPPED_GENES <- 5
PLOT_TOP_N <- 15
LABEL_WRAP_WIDTH <- 35

############################################################
## 1. Packages
############################################################

pkgs_cran <- c("data.table", "ggplot2", "dplyr", "stringr", "forcats")
pkgs_bioc <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA",
  "msigdbr"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)

############################################################
## 2. Helper functions
############################################################

parse_ratio <- function(x) {
  x <- as.character(x)
  num <- suppressWarnings(as.numeric(sub("/.*", "", x)))
  den <- suppressWarnings(as.numeric(sub(".*/", "", x)))
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}

make_pathway_vector <- function(res_dt, source_name) {
  
  if (is.null(res_dt) || nrow(res_dt) == 0) {
    return(data.table())
  }
  
  x <- copy(res_dt)
  
  required_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  miss <- setdiff(required_cols, colnames(x))
  if (length(miss) > 0) {
    stop("Missing columns in enrichment result: ", paste(miss, collapse = ", "))
  }
  
  x[, Pathway_Source := source_name]
  x[, Pathway_ID := as.character(ID)]
  x[, Pathway_Name := as.character(Description)]
  
  x[, GeneRatio_Value := parse_ratio(GeneRatio)]
  x[, BgRatio_Value := parse_ratio(BgRatio)]
  
  x[, FDR := as.numeric(p.adjust)]
  x[, FDR_Capped := pmax(FDR, FDR_CAP)]
  x[, NegLog10FDR := -log10(FDR_Capped)]
  
  x[, Pvalue_Capped := pmax(as.numeric(pvalue), FDR_CAP)]
  x[, NegLog10P := -log10(Pvalue_Capped)]
  
  x[, Count := as.numeric(Count)]
  
  x[, PathwayScore := GeneRatio_Value * NegLog10FDR]
  x[, Score_NegLog10FDR := NegLog10FDR]
  x[, Score_GeneRatioOnly := GeneRatio_Value]
  x[, Score_CountWeighted := Count * NegLog10FDR]
  
  x[, Is_FDR0p05 := FDR <= FDR_SIG]
  
  x[, .(
    Pathway_Source,
    Pathway_ID,
    Pathway_Name,
    GeneRatio,
    BgRatio,
    GeneRatio_Value,
    BgRatio_Value,
    Count,
    pvalue,
    p.adjust,
    qvalue,
    FDR,
    NegLog10FDR,
    NegLog10P,
    PathwayScore,
    Score_NegLog10FDR,
    Score_GeneRatioOnly,
    Score_CountWeighted,
    Is_FDR0p05,
    GeneID = geneID
  )]
}

############################################################
## 3. Batch runner
############################################################

run_one_entity <- function(entity_type, entity_name, base_dir, input_file) {
  
  if (!dir.exists(base_dir)) {
    warning("Base directory not found: ", base_dir)
    return(NULL)
  }
  
  if (!file.exists(input_file)) {
    warning("Input gene file not found: ", input_file)
    return(NULL)
  }
  
  safe_name <- gsub("[/\\\\:*?\"<>| ]+", "_", entity_name)
  
  OUT_DIR <- file.path(
    base_dir,
    paste0(entity_type, "_PathwayVector_", safe_name, "_", RUN_TAG)
  )
  
  TABLE_DIR <- file.path(OUT_DIR, "01_Tables")
  VECTOR_DIR <- file.path(OUT_DIR, "02_ModelInput_PathwayVectors")
  FIG_DIR <- file.path(OUT_DIR, "03_QC_Figures")
  RDS_DIR <- file.path(OUT_DIR, "04_RDS")
  
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(VECTOR_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
  
  save_enrich_result <- function(enrich_obj, prefix) {
    
    saveRDS(enrich_obj, file.path(RDS_DIR, paste0(prefix, ".rds")))
    
    res <- as.data.table(enrich_obj@result)
    
    if (nrow(res) == 0) {
      fwrite(data.table(), file.path(TABLE_DIR, paste0(prefix, "_FULL.tsv")), sep = "\t")
      fwrite(data.table(), file.path(TABLE_DIR, paste0(prefix, "_FDR0p05.tsv")), sep = "\t")
      return(res)
    }
    
    res[, Source_File_Prefix := prefix]
    
    fwrite(res, file.path(TABLE_DIR, paste0(prefix, "_FULL.tsv")), sep = "\t")
    fwrite(res[p.adjust <= FDR_SIG], file.path(TABLE_DIR, paste0(prefix, "_FDR0p05.tsv")), sep = "\t")
    
    return(res)
  }
  
  plot_qc_bubble <- function(vector_dt, prefix, title_text, top_n = 15) {
    
    if (is.null(vector_dt) || nrow(vector_dt) == 0) {
      message("No result for plotting: ", prefix)
      return(NULL)
    }
    
    plot_dt <- copy(vector_dt)
    plot_dt <- plot_dt[is.finite(PathwayScore)]
    
    if (nrow(plot_dt) == 0) {
      message("No finite pathway scores for plotting: ", prefix)
      return(NULL)
    }
    
    setorder(plot_dt, p.adjust, -PathwayScore)
    plot_dt <- plot_dt[seq_len(min(top_n, .N))]
    
    plot_dt[, Pathway_Label := stringr::str_wrap(Pathway_Name, width = LABEL_WRAP_WIDTH)]
    plot_dt[, Pathway_Label := factor(Pathway_Label, levels = rev(Pathway_Label))]
    
    fig_height <- max(5.5, 0.36 * nrow(plot_dt) + 2.0)
    
    p <- ggplot(
      plot_dt,
      aes(
        x = GeneRatio_Value,
        y = Pathway_Label,
        size = Count,
        color = NegLog10FDR
      )
    ) +
      geom_point(alpha = 0.85) +
      scale_size_continuous(range = c(2.5, 8), name = "Count") +
      scale_color_gradient(name = "-log10(FDR)", low = "#6BAED6", high = "#DE2D26") +
      labs(
        title = title_text,
        x = "Gene ratio",
        y = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0),
        axis.text.y = element_text(size = 8.5, color = "black", lineheight = 0.92),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(face = "bold", size = 11),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 20, 10, 15)
      )
    
    ggsave(
      file.path(FIG_DIR, paste0(prefix, "_QC_Bubble.png")),
      p,
      width = 10.5,
      height = fig_height,
      dpi = 600,
      bg = "white",
      limitsize = FALSE
    )
    
    ggsave(
      file.path(FIG_DIR, paste0(prefix, "_QC_Bubble.pdf")),
      p,
      width = 10.5,
      height = fig_height,
      bg = "white",
      limitsize = FALSE
    )
    
    return(p)
  }
  
  make_wide_vector <- function(vector_dt, value_col, out_name) {
    
    if (nrow(vector_dt) == 0) {
      fwrite(data.table(), file.path(VECTOR_DIR, out_name), sep = "\t")
      return(invisible(NULL))
    }
    
    x <- copy(vector_dt)
    
    x[, Feature_ID := paste(Pathway_Source, Pathway_ID, sep = "__")]
    x[, Feature_Name := paste(Pathway_Source, Pathway_Name, sep = "__")]
    
    wide <- x[, .(
      Feature_ID,
      Feature_Name,
      Pathway_Source,
      Pathway_ID,
      Pathway_Name,
      Value = get(value_col)
    )]
    
    setnames(wide, "Value", entity_name)
    
    fwrite(wide, file.path(VECTOR_DIR, out_name), sep = "\t")
    
    invisible(wide)
  }
  
  ############################################################
  ## Read and clean gene list
  ############################################################
  
  genes_raw <- fread(
    input_file,
    header = FALSE,
    data.table = FALSE,
    stringsAsFactors = FALSE
  )[[1]]
  
  genes_symbol <- genes_raw |>
    as.character() |>
    stringr::str_trim() |>
    toupper() |>
    unique()
  
  genes_symbol <- genes_symbol[
    !is.na(genes_symbol) &
      genes_symbol != "" &
      genes_symbol != "NA"
  ]
  
  fwrite(
    data.table(Gene_Symbol = genes_symbol),
    file.path(TABLE_DIR, "00_Input_Cleaned_GeneSymbols.tsv"),
    sep = "\t"
  )
  
  cat("\n============================================================\n")
  cat("Running:", entity_type, "-", entity_name, "\n")
  cat("Input file:", input_file, "\n")
  cat("Input genes:", length(genes_raw), "\n")
  cat("Unique cleaned gene symbols:", length(genes_symbol), "\n")
  
  ############################################################
  ## SYMBOL to ENTREZID mapping
  ############################################################
  
  gene_map <- bitr(
    genes_symbol,
    fromType = "SYMBOL",
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Hs.eg.db
  )
  
  gene_map <- unique(as.data.table(gene_map))
  mapped_entrez <- unique(gene_map$ENTREZID)
  unmapped_genes <- setdiff(genes_symbol, gene_map$SYMBOL)
  
  fwrite(
    gene_map,
    file.path(TABLE_DIR, "00_Gene_SYMBOL_to_ENTREZID_Mapping.tsv"),
    sep = "\t"
  )
  
  fwrite(
    data.table(Unmapped_Gene = unmapped_genes),
    file.path(TABLE_DIR, "00_Unmapped_Genes.tsv"),
    sep = "\t"
  )
  
  cat("Mapped ENTREZ genes:", length(mapped_entrez), "\n")
  cat("Unmapped genes:", length(unmapped_genes), "\n")
  
  if (length(mapped_entrez) < MIN_MAPPED_GENES) {
    warning("Too few mapped genes for enrichment analysis: ", entity_name)
    return(NULL)
  }
  
  ############################################################
  ## Reactome ORA
  ############################################################
  
  reactome_ora <- enrichPathway(
    gene = mapped_entrez,
    organism = "human",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  reactome_res <- save_enrich_result(reactome_ora, "01_Reactome_ORA")
  reactome_vec <- make_pathway_vector(reactome_res, "Reactome_ORA")
  
  fwrite(
    reactome_vec,
    file.path(VECTOR_DIR, "01_Reactome_ORA_PathwayVector.tsv"),
    sep = "\t"
  )
  
  plot_qc_bubble(
    reactome_vec,
    "01_Reactome_ORA",
    paste0(entity_name, " | Reactome ORA pathway vector QC"),
    top_n = PLOT_TOP_N
  )
  
  ############################################################
  ## GO Biological Process ORA
  ############################################################
  
  go_bp_ora <- enrichGO(
    gene = mapped_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  go_bp_res <- save_enrich_result(go_bp_ora, "02_GO_BP_ORA")
  go_bp_vec <- make_pathway_vector(go_bp_res, "GO_BP_ORA")
  
  fwrite(
    go_bp_vec,
    file.path(VECTOR_DIR, "02_GO_BP_ORA_PathwayVector.tsv"),
    sep = "\t"
  )
  
  plot_qc_bubble(
    go_bp_vec,
    "02_GO_BP_ORA",
    paste0(entity_name, " | GO Biological Process ORA pathway vector QC"),
    top_n = PLOT_TOP_N
  )
  
  ############################################################
  ## MSigDB Hallmark ORA
  ############################################################
  
  hallmark_df <- msigdbr(
    species = "Homo sapiens",
    category = "H"
  )
  
  hallmark_term2gene <- hallmark_df |>
    dplyr::select(gs_name, entrez_gene) |>
    dplyr::filter(!is.na(entrez_gene)) |>
    dplyr::distinct()
  
  hallmark_term2name <- hallmark_df |>
    dplyr::select(gs_name, gs_description) |>
    dplyr::distinct()
  
  hallmark_ora <- enricher(
    gene = mapped_entrez,
    TERM2GENE = hallmark_term2gene,
    TERM2NAME = hallmark_term2name,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH"
  )
  
  hallmark_res <- save_enrich_result(hallmark_ora, "03_MSigDB_Hallmark_ORA")
  hallmark_vec <- make_pathway_vector(hallmark_res, "MSigDB_Hallmark_ORA")
  
  fwrite(
    hallmark_vec,
    file.path(VECTOR_DIR, "03_MSigDB_Hallmark_ORA_PathwayVector.tsv"),
    sep = "\t"
  )
  
  plot_qc_bubble(
    hallmark_vec,
    "03_MSigDB_Hallmark_ORA",
    paste0(entity_name, " | MSigDB Hallmark ORA pathway vector QC"),
    top_n = PLOT_TOP_N
  )
  
  ############################################################
  ## Unified model-input pathway vector
  ############################################################
  
  entity_pathway_vector <- rbindlist(
    list(
      reactome_vec,
      go_bp_vec,
      hallmark_vec
    ),
    use.names = TRUE,
    fill = TRUE
  )
  
  entity_pathway_vector[, Entity_Type := entity_type]
  entity_pathway_vector[, Entity_Name := entity_name]
  
  setorder(entity_pathway_vector, Pathway_Source, p.adjust, -PathwayScore)
  
  fwrite(
    entity_pathway_vector,
    file.path(VECTOR_DIR, "04_Entity_PathwayVector_ForModelInput_LONG.tsv"),
    sep = "\t"
  )
  
  make_wide_vector(
    entity_pathway_vector,
    value_col = "PathwayScore",
    out_name = "05_Entity_PathwayVector_PathwayScore_WIDE.tsv"
  )
  
  make_wide_vector(
    entity_pathway_vector,
    value_col = "Score_NegLog10FDR",
    out_name = "06_Entity_PathwayVector_NegLog10FDR_WIDE.tsv"
  )
  
  make_wide_vector(
    entity_pathway_vector,
    value_col = "Score_GeneRatioOnly",
    out_name = "07_Entity_PathwayVector_GeneRatio_WIDE.tsv"
  )
  
  ############################################################
  ## Source-level QC summary
  ############################################################
  
  source_qc <- entity_pathway_vector[
    ,
    .(
      Total_Term_N = .N,
      FDR0p05_Term_N = sum(Is_FDR0p05, na.rm = TRUE),
      Median_PathwayScore = median(PathwayScore, na.rm = TRUE),
      Max_PathwayScore = max(PathwayScore, na.rm = TRUE),
      Median_GeneRatio = median(GeneRatio_Value, na.rm = TRUE),
      Max_GeneRatio = max(GeneRatio_Value, na.rm = TRUE)
    ),
    by = Pathway_Source
  ]
  
  fwrite(
    source_qc,
    file.path(TABLE_DIR, "04_SourceLevel_QC_Summary.tsv"),
    sep = "\t"
  )
  
  top_terms_for_qc <- entity_pathway_vector[
    order(Pathway_Source, p.adjust, -PathwayScore)
  ][
    ,
    head(.SD, 30),
    by = Pathway_Source
  ]
  
  fwrite(
    top_terms_for_qc,
    file.path(TABLE_DIR, "05_Top30_Terms_PerSource_ForQC.tsv"),
    sep = "\t"
  )
  
  ############################################################
  ## Metadata
  ############################################################
  
  metadata <- data.table(
    Item = c(
      "Entity_Type",
      "Entity_Name",
      "Input_File",
      "Base_Dir",
      "Output_Dir",
      "Run_Tag",
      "Input_Gene_N",
      "Unique_Cleaned_Gene_N",
      "Mapped_ENTREZ_N",
      "Unmapped_Gene_N",
      "Reactome_Term_N",
      "Reactome_FDR0p05_N",
      "GO_BP_Term_N",
      "GO_BP_FDR0p05_N",
      "Hallmark_Term_N",
      "Hallmark_FDR0p05_N",
      "Main_Model_Input_Long",
      "Main_Model_Input_Wide"
    ),
    Value = c(
      entity_type,
      entity_name,
      input_file,
      base_dir,
      OUT_DIR,
      RUN_TAG,
      length(genes_raw),
      length(genes_symbol),
      length(mapped_entrez),
      length(unmapped_genes),
      nrow(reactome_vec),
      sum(reactome_vec$Is_FDR0p05, na.rm = TRUE),
      nrow(go_bp_vec),
      sum(go_bp_vec$Is_FDR0p05, na.rm = TRUE),
      nrow(hallmark_vec),
      sum(hallmark_vec$Is_FDR0p05, na.rm = TRUE),
      file.path(VECTOR_DIR, "04_Entity_PathwayVector_ForModelInput_LONG.tsv"),
      file.path(VECTOR_DIR, "05_Entity_PathwayVector_PathwayScore_WIDE.tsv")
    )
  )
  
  fwrite(
    metadata,
    file.path(OUT_DIR, "00_Run_Metadata.tsv"),
    sep = "\t"
  )
  
  cat("Completed:", entity_name, "\n")
  cat("Output folder:", OUT_DIR, "\n")
  cat("============================================================\n")
  
  return(data.table(
    Entity_Type = entity_type,
    Entity_Name = entity_name,
    Input_File = input_file,
    Output_Dir = OUT_DIR,
    Input_Gene_N = length(genes_raw),
    Cleaned_Gene_N = length(genes_symbol),
    Mapped_ENTREZ_N = length(mapped_entrez),
    Unmapped_Gene_N = length(unmapped_genes),
    Reactome_Term_N = nrow(reactome_vec),
    GO_BP_Term_N = nrow(go_bp_vec),
    Hallmark_Term_N = nrow(hallmark_vec)
  ))
}

############################################################
## 4. Run all jobs
############################################################

all_summary <- list()

for (i in seq_len(nrow(INPUT_JOBS))) {
  
  one <- INPUT_JOBS[i, ]
  
  res <- tryCatch(
    run_one_entity(
      entity_type = one$Entity_Type,
      entity_name = one$Entity_Name,
      base_dir = one$Base_Dir,
      input_file = one$Input_File
    ),
    error = function(e) {
      warning("Failed: ", one$Entity_Name, " | ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (!is.null(res)) {
    all_summary[[length(all_summary) + 1]] <- res
  }
}

batch_summary <- rbindlist(all_summary, use.names = TRUE, fill = TRUE)

BATCH_OUT <- file.path(ROOT_DIR, paste0("Batch_PathwayVector_Summary_", RUN_TAG))
dir.create(BATCH_OUT, recursive = TRUE, showWarnings = FALSE)

fwrite(
  INPUT_JOBS,
  file.path(BATCH_OUT, "00_Batch_Input_Jobs.tsv"),
  sep = "\t"
)

fwrite(
  batch_summary,
  file.path(BATCH_OUT, "01_Batch_Run_Summary.tsv"),
  sep = "\t"
)

cat("\n============================================================\n")
cat("Batch pathway-vector construction completed.\n")
cat("Batch summary folder:\n")
cat(BATCH_OUT, "\n")
cat("============================================================\n")