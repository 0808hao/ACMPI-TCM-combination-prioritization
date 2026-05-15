############################################################
## ANPS Step 1
## Multi-disease genes -> STRING PPI mapping
## RA / OA / DM disease modules and node-to-disease distance
############################################################

rm(list = ls())

############################################################
## 0. Paths and parameters
############################################################

DISEASE_ROOT <- "/media/desk16/iy15915/中药之开创/疾病验证"

## STRING background files: DO NOT CHANGE
BG_DIR <- "/media/desk16/iy15915/中药之开创/脑卒中疾病基因/背景基因"

STRING_ALIAS_FILE <- file.path(BG_DIR, "9606.protein.aliases.v12.0.txt")
STRING_INFO_FILE  <- file.path(BG_DIR, "9606.protein.info.v12.0.txt")
STRING_LINK_FILE  <- file.path(BG_DIR, "9606.protein.links.v12.0.txt")

STRING_SCORE_MIN <- 700

disease_jobs <- data.table::data.table(
  Disease_ID = c("RA", "OA", "DM"),
  Disease_Label = c(
    "Rheumatoid_Arthritis",
    "Osteoarthritis",
    "Diabetes_Mellitus"
  ),
  Gene_File = c(
    file.path(DISEASE_ROOT, "类风湿性关节炎", "RA_vs_Control_Candidate_Genes.txt"),
    file.path(DISEASE_ROOT, "类风湿性关节炎", "OA_vs_Control_Candidate_Genes.txt"),
    file.path(DISEASE_ROOT, "糖尿病", "DM_Candidate_Genes_Final.txt")
  ),
  Output_Parent = c(
    file.path(DISEASE_ROOT, "类风湿性关节炎"),
    file.path(DISEASE_ROOT, "类风湿性关节炎"),
    file.path(DISEASE_ROOT, "糖尿病")
  )
)

RUN_TAG <- format(Sys.time(), "%Y%m%d_%H%M%S")

cat("\n============================================================\n")
cat("ANPS Step 1: Multi-disease STRING PPI module\n")
cat("STRING background directory:\n")
cat(BG_DIR, "\n")
cat("============================================================\n")

############################################################
## 1. Packages
############################################################

pkgs <- c("data.table", "igraph")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(data.table)
library(igraph)

############################################################
## 2. Input checks
############################################################

string_files <- c(
  STRING_ALIAS_FILE,
  STRING_INFO_FILE,
  STRING_LINK_FILE
)

missing_string_files <- string_files[!file.exists(string_files)]

if (length(missing_string_files) > 0) {
  stop(
    "Missing STRING background files:\n",
    paste(missing_string_files, collapse = "\n"),
    call. = FALSE
  )
}

missing_gene_files <- disease_jobs$Gene_File[!file.exists(disease_jobs$Gene_File)]

if (length(missing_gene_files) > 0) {
  stop(
    "Missing disease gene files:\n",
    paste(missing_gene_files, collapse = "\n"),
    call. = FALSE
  )
}

############################################################
## 3. Read STRING protein info
############################################################

cat("\nReading STRING protein info...\n")

string_info <- data.table::fread(
  STRING_INFO_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  data.table = TRUE,
  showProgress = TRUE
)

setnames(string_info, colnames(string_info)[1], "string_id")

required_info_cols <- c("string_id", "preferred_name")
missing_info_cols <- setdiff(required_info_cols, colnames(string_info))

if (length(missing_info_cols) > 0) {
  stop(
    "STRING info file missing columns:\n",
    paste(missing_info_cols, collapse = "\n"),
    call. = FALSE
  )
}

string_info <- unique(
  string_info[, .(
    string_id = as.character(string_id),
    preferred_name = as.character(preferred_name),
    protein_size = if ("protein_size" %in% colnames(string_info)) protein_size else NA,
    annotation = if ("annotation" %in% colnames(string_info)) annotation else NA
  )]
)

############################################################
## 4. Read STRING aliases once
############################################################

cat("\nReading STRING aliases...\n")

string_alias <- data.table::fread(
  STRING_ALIAS_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  data.table = TRUE,
  showProgress = TRUE
)

setnames(string_alias, colnames(string_alias)[1], "string_id")

required_alias_cols <- c("string_id", "alias", "source")
missing_alias_cols <- setdiff(required_alias_cols, colnames(string_alias))

if (length(missing_alias_cols) > 0) {
  stop(
    "STRING alias file missing columns:\n",
    paste(missing_alias_cols, collapse = "\n"),
    call. = FALSE
  )
}

string_alias <- string_alias[, .(
  string_id = as.character(string_id),
  alias = as.character(alias),
  source = as.character(source)
)]

############################################################
## 5. Read STRING PPI links and build background graph once
############################################################

cat("\nReading STRING PPI links...\n")
cat("This may take several minutes for the full STRING file.\n")

ppi_links <- data.table::fread(
  STRING_LINK_FILE,
  sep = " ",
  header = TRUE,
  quote = "",
  data.table = TRUE,
  showProgress = TRUE
)

required_link_cols <- c("protein1", "protein2", "combined_score")
missing_link_cols <- setdiff(required_link_cols, colnames(ppi_links))

if (length(missing_link_cols) > 0) {
  stop(
    "STRING link file missing columns:\n",
    paste(missing_link_cols, collapse = "\n"),
    call. = FALSE
  )
}

ppi_links[, protein1 := as.character(protein1)]
ppi_links[, protein2 := as.character(protein2)]
ppi_links[, combined_score := as.numeric(combined_score)]

ppi_filtered <- ppi_links[
  !is.na(combined_score) &
    combined_score >= STRING_SCORE_MIN &
    protein1 != protein2,
  .(protein1, protein2, combined_score)
]

rm(ppi_links)
gc()

cat("\nSTRING PPI background after filtering:\n")
cat("Edges with combined_score >=", STRING_SCORE_MIN, ":", nrow(ppi_filtered), "\n")

cat("\nBuilding full STRING PPI graph...\n")

g_ppi <- igraph::graph_from_data_frame(
  d = ppi_filtered[, .(
    from = protein1,
    to = protein2,
    combined_score = combined_score
  )],
  directed = FALSE
)

g_ppi <- igraph::simplify(
  g_ppi,
  remove.multiple = TRUE,
  remove.loops = TRUE,
  edge.attr.comb = list(combined_score = "max")
)

rm(ppi_filtered)
gc()

ppi_nodes <- V(g_ppi)$name

cat("Full PPI nodes:", igraph::vcount(g_ppi), "\n")
cat("Full PPI edges:", igraph::ecount(g_ppi), "\n")

############################################################
## 6. Helper function: one disease
############################################################

run_one_disease <- function(
    disease_id,
    disease_label,
    disease_gene_file,
    output_parent
) {
  
  cat("\n============================================================\n")
  cat("Processing disease:", disease_id, "-", disease_label, "\n")
  cat("Gene file:\n")
  cat(disease_gene_file, "\n")
  cat("============================================================\n")
  
  disease_genes_raw <- readLines(disease_gene_file, warn = FALSE)
  disease_genes_raw <- trimws(disease_genes_raw)
  disease_genes_raw <- disease_genes_raw[disease_genes_raw != ""]
  disease_genes <- unique(disease_genes_raw)
  
  if (length(disease_genes) == 0) {
    stop("No disease genes found in: ", disease_gene_file, call. = FALSE)
  }
  
  cat("\nDisease genes loaded:\n")
  cat("Raw non-empty lines:", length(disease_genes_raw), "\n")
  cat("Unique genes:", length(disease_genes), "\n")
  
  ##########################################################
  ## Mapping disease genes to STRING IDs
  ##########################################################
  
  disease_dt <- data.table(
    Input_Gene = disease_genes,
    Input_Gene_Upper = toupper(disease_genes)
  )
  
  info_direct <- string_info[
    preferred_name %in% disease_genes,
    .(
      Input_Gene = preferred_name,
      string_id,
      preferred_name,
      Mapping_Method = "preferred_name",
      Mapping_Source = "protein.info"
    )
  ]
  
  info_direct <- info_direct[
    order(Input_Gene, string_id)
  ][
    ,
    .SD[1],
    by = Input_Gene
  ]
  
  unmapped_after_direct <- setdiff(disease_genes, info_direct$Input_Gene)
  
  cat("Direct preferred_name mapped genes:", nrow(info_direct), "\n")
  cat("Genes requiring alias mapping:", length(unmapped_after_direct), "\n")
  
  alias_selected <- data.table()
  
  if (length(unmapped_after_direct) > 0) {
    
    alias_hits <- string_alias[
      alias %in% unmapped_after_direct
    ]
    
    if (nrow(alias_hits) > 0) {
      
      alias_hits <- merge(
        alias_hits,
        string_info[, .(string_id, preferred_name)],
        by = "string_id",
        all.x = TRUE
      )
      
      alias_hits[, source_lower := tolower(source)]
      
      alias_hits[, source_priority := fifelse(
        grepl("hgnc", source_lower) & grepl("symbol", source_lower), 1L,
        fifelse(
          grepl("hgnc", source_lower), 2L,
          fifelse(
            grepl("uniprot", source_lower) & grepl("gn", source_lower), 3L,
            fifelse(
              grepl("geneid", source_lower) | grepl("entrez", source_lower), 4L,
              9L
            )
          )
        )
      )]
      
      alias_hits[, preferred_match := as.integer(alias == preferred_name)]
      
      setorder(alias_hits, alias, source_priority, -preferred_match, string_id)
      
      alias_selected <- alias_hits[
        ,
        .SD[1],
        by = alias
      ][
        ,
        .(
          Input_Gene = alias,
          string_id,
          preferred_name,
          Mapping_Method = "alias",
          Mapping_Source = source
        )
      ]
    }
  }
  
  mapping_combined <- rbindlist(
    list(info_direct, alias_selected),
    use.names = TRUE,
    fill = TRUE
  )
  
  if (nrow(mapping_combined) > 0) {
    mapping_combined[, method_priority := fifelse(Mapping_Method == "preferred_name", 1L, 2L)]
    setorder(mapping_combined, Input_Gene, method_priority, string_id)
    
    mapping_selected <- mapping_combined[
      ,
      .SD[1],
      by = Input_Gene
    ][
      ,
      method_priority := NULL
    ]
  } else {
    mapping_selected <- data.table(
      Input_Gene = character(),
      string_id = character(),
      preferred_name = character(),
      Mapping_Method = character(),
      Mapping_Source = character()
    )
  }
  
  mapping_table <- merge(
    disease_dt[, .(Input_Gene)],
    mapping_selected,
    by = "Input_Gene",
    all.x = TRUE,
    sort = FALSE
  )
  
  mapping_table[, Mapped_to_STRING := !is.na(string_id)]
  mapping_table[, Present_in_STRING700_PPI := string_id %in% ppi_nodes]
  
  cat("\nFinal mapping summary:\n")
  cat("Input disease genes:", nrow(mapping_table), "\n")
  cat("Mapped to STRING:", sum(mapping_table$Mapped_to_STRING), "\n")
  cat("Present in STRING700 PPI:", sum(mapping_table$Present_in_STRING700_PPI, na.rm = TRUE), "\n")
  cat("Unmapped:", sum(!mapping_table$Mapped_to_STRING), "\n")
  
  ##########################################################
  ## Disease nodes in PPI
  ##########################################################
  
  mapped_genes <- mapping_table[
    Mapped_to_STRING == TRUE &
      string_id %in% ppi_nodes
  ]
  
  disease_string_ids <- unique(mapped_genes$string_id)
  
  if (length(disease_string_ids) == 0) {
    stop(
      "No mapped disease genes are present in the filtered STRING PPI network for: ",
      disease_id,
      call. = FALSE
    )
  }
  
  ##########################################################
  ## Disease induced subgraph
  ##########################################################
  
  g_disease <- igraph::induced_subgraph(
    graph = g_ppi,
    vids = disease_string_ids
  )
  
  disease_components <- igraph::components(g_disease)
  
  largest_component_size <- if (igraph::vcount(g_disease) > 0) {
    max(disease_components$csize)
  } else {
    0
  }
  
  largest_component_fraction <- if (igraph::vcount(g_disease) > 0) {
    largest_component_size / igraph::vcount(g_disease)
  } else {
    NA_real_
  }
  
  cat("\nDisease PPI module:\n")
  cat("Disease PPI nodes:", igraph::vcount(g_disease), "\n")
  cat("Disease PPI edges:", igraph::ecount(g_disease), "\n")
  cat("Disease PPI components:", disease_components$no, "\n")
  cat("Largest component size:", largest_component_size, "\n")
  cat("Largest component fraction:", round(largest_component_fraction, 4), "\n")
  
  ##########################################################
  ## Node-to-disease shortest-path distance
  ##########################################################
  
  cat("\nPrecomputing node-to-disease shortest-path distance...\n")
  
  SUPER_NODE <- paste0("__", disease_id, "_DISEASE_SUPER_NODE__")
  
  g_tmp <- igraph::add_vertices(
    g_ppi,
    nv = 1,
    name = SUPER_NODE
  )
  
  super_edges <- as.vector(rbind(
    rep(SUPER_NODE, length(disease_string_ids)),
    disease_string_ids
  ))
  
  g_tmp <- igraph::add_edges(
    g_tmp,
    edges = super_edges
  )
  
  dist_vec <- igraph::distances(
    graph = g_tmp,
    v = SUPER_NODE,
    to = ppi_nodes,
    mode = "all",
    weights = NA
  )
  
  dist_vec <- as.numeric(dist_vec[1, ])
  distance_to_disease <- dist_vec - 1
  distance_to_disease[is.infinite(distance_to_disease)] <- NA_real_
  
  rm(g_tmp, dist_vec)
  gc()
  
  ##########################################################
  ## Output tables
  ##########################################################
  
  mapping_out <- copy(mapping_table)
  
  setcolorder(
    mapping_out,
    c(
      "Input_Gene",
      "Mapped_to_STRING",
      "Present_in_STRING700_PPI",
      "string_id",
      "preferred_name",
      "Mapping_Method",
      "Mapping_Source"
    )
  )
  
  degree_vec <- igraph::degree(g_ppi, v = ppi_nodes)
  
  node_distance <- data.table(
    string_id = ppi_nodes,
    distance_to_disease = distance_to_disease,
    degree_STRING700 = as.numeric(degree_vec),
    is_disease_node = ppi_nodes %in% disease_string_ids
  )
  
  setnames(
    node_distance,
    "distance_to_disease",
    paste0("distance_to_", disease_id, "_disease")
  )
  
  distance_col <- paste0("distance_to_", disease_id, "_disease")
  
  node_distance <- merge(
    node_distance,
    string_info[, .(string_id, preferred_name)],
    by = "string_id",
    all.x = TRUE
  )
  
  gene_collapse <- mapping_out[
    Mapped_to_STRING == TRUE,
    .(
      Input_Disease_Genes = paste(unique(Input_Gene), collapse = ";")
    ),
    by = string_id
  ]
  
  node_distance <- merge(
    node_distance,
    gene_collapse,
    by = "string_id",
    all.x = TRUE
  )
  
  setcolorder(
    node_distance,
    c(
      "string_id",
      "preferred_name",
      distance_col,
      "degree_STRING700",
      "is_disease_node",
      "Input_Disease_Genes"
    )
  )
  
  setorderv(node_distance, c(distance_col, "degree_STRING700", "preferred_name"), c(1, -1, 1))
  
  disease_edges <- igraph::as_data_frame(
    g_disease,
    what = "edges"
  )
  
  if (nrow(disease_edges) > 0) {
    disease_edges <- as.data.table(disease_edges)
    
    disease_edges <- merge(
      disease_edges,
      string_info[, .(from = string_id, from_gene = preferred_name)],
      by = "from",
      all.x = TRUE
    )
    
    disease_edges <- merge(
      disease_edges,
      string_info[, .(to = string_id, to_gene = preferred_name)],
      by = "to",
      all.x = TRUE
    )
    
    setcolorder(
      disease_edges,
      c("from", "from_gene", "to", "to_gene", "combined_score")
    )
    
    setorder(disease_edges, -combined_score, from_gene, to_gene)
  } else {
    disease_edges <- data.table(
      from = character(),
      from_gene = character(),
      to = character(),
      to_gene = character(),
      combined_score = numeric()
    )
  }
  
  ##########################################################
  ## QC summary
  ##########################################################
  
  mapped_n <- sum(mapping_out$Mapped_to_STRING)
  ppi_present_n <- sum(mapping_out$Present_in_STRING700_PPI, na.rm = TRUE)
  
  distance_stats <- node_distance[
    !is.na(get(distance_col)),
    .(
      min_distance = min(get(distance_col)),
      median_distance = median(get(distance_col)),
      mean_distance = mean(get(distance_col)),
      max_distance = max(get(distance_col))
    )
  ]
  
  disease_degree <- node_distance[
    is_disease_node == TRUE,
    degree_STRING700
  ]
  
  background_degree <- node_distance[
    ,
    degree_STRING700
  ]
  
  qc_summary <- data.table(
    Metric = c(
      "Disease ID",
      "Disease label",
      "Input disease genes",
      "Mapped genes to STRING",
      "Unmapped genes",
      "Mapping rate to STRING",
      "Mapped genes present in STRING700 PPI",
      "PPI-present mapping rate",
      "Unique disease STRING nodes",
      "STRING score threshold",
      "Full PPI nodes",
      "Full PPI edges",
      "Disease PPI nodes",
      "Disease PPI edges",
      "Disease PPI components",
      "Largest disease component size",
      "Largest disease component fraction",
      "Median degree of disease nodes",
      "Median degree of background nodes",
      "Mean distance among reachable PPI nodes to disease module",
      "Median distance among reachable PPI nodes to disease module"
    ),
    Value = as.character(c(
      disease_id,
      disease_label,
      length(disease_genes),
      mapped_n,
      length(disease_genes) - mapped_n,
      round(mapped_n / length(disease_genes), 4),
      ppi_present_n,
      round(ppi_present_n / length(disease_genes), 4),
      length(disease_string_ids),
      STRING_SCORE_MIN,
      igraph::vcount(g_ppi),
      igraph::ecount(g_ppi),
      igraph::vcount(g_disease),
      igraph::ecount(g_disease),
      disease_components$no,
      largest_component_size,
      round(largest_component_fraction, 4),
      round(median(disease_degree, na.rm = TRUE), 4),
      round(median(background_degree, na.rm = TRUE), 4),
      round(distance_stats$mean_distance, 4),
      round(distance_stats$median_distance, 4)
    ))
  )
  
  ##########################################################
  ## Save
  ##########################################################
  
  out_dir <- file.path(
    output_parent,
    paste0("ANPS_", disease_id, "_STRING", STRING_SCORE_MIN, "_", RUN_TAG)
  )
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_mapping <- file.path(out_dir, paste0("01_", disease_id, "_GeneMapping_STRING.tsv"))
  out_distance <- file.path(out_dir, paste0("02_", disease_id, "_STRING700_NodeDistanceToDisease.tsv"))
  out_edges <- file.path(out_dir, paste0("03_", disease_id, "_DiseaseModuleEdges_STRING700.tsv"))
  out_qc <- file.path(out_dir, paste0("04_", disease_id, "_PPI_QC.tsv"))
  
  data.table::fwrite(mapping_out, out_mapping, sep = "\t")
  data.table::fwrite(node_distance, out_distance, sep = "\t")
  data.table::fwrite(disease_edges, out_edges, sep = "\t")
  data.table::fwrite(qc_summary, out_qc, sep = "\t")
  
  cat("\nEssential outputs saved:\n")
  cat(out_mapping, "\n")
  cat(out_distance, "\n")
  cat(out_edges, "\n")
  cat(out_qc, "\n")
  
  cat("\nTop 20 closest non-disease PPI nodes:\n")
  print(
    node_distance[
      is_disease_node == FALSE &
        !is.na(get(distance_col))
    ][
      order(get(distance_col), -degree_STRING700)
    ][
      1:min(.N, 20),
      .(
        preferred_name,
        string_id,
        distance = get(distance_col),
        degree_STRING700
      )
    ]
  )
  
  return(data.table(
    Disease_ID = disease_id,
    Disease_Label = disease_label,
    Output_Dir = out_dir,
    Input_Genes = length(disease_genes),
    Mapped_to_STRING = mapped_n,
    Present_in_STRING700_PPI = ppi_present_n,
    Disease_PPI_Nodes = igraph::vcount(g_disease),
    Disease_PPI_Edges = igraph::ecount(g_disease),
    Largest_Component_Fraction = largest_component_fraction
  ))
}

############################################################
## 7. Run all diseases
############################################################

all_summary <- rbindlist(
  lapply(seq_len(nrow(disease_jobs)), function(i) {
    run_one_disease(
      disease_id = disease_jobs$Disease_ID[i],
      disease_label = disease_jobs$Disease_Label[i],
      disease_gene_file = disease_jobs$Gene_File[i],
      output_parent = disease_jobs$Output_Parent[i]
    )
  }),
  use.names = TRUE,
  fill = TRUE
)

summary_out <- file.path(
  DISEASE_ROOT,
  paste0("ANPS_AllDisease_STRING", STRING_SCORE_MIN, "_RunSummary_", RUN_TAG, ".tsv")
)

data.table::fwrite(all_summary, summary_out, sep = "\t")

cat("\n============================================================\n")
cat("All disease PPI mapping completed.\n")
cat("Run summary saved:\n")
cat(summary_out, "\n")
cat("============================================================\n")

print(all_summary)