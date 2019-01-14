#### This script identifies and annotates genes of interest from clustered data.
#### It requires the working environment of cart19_intsite_analysis.Rmd.
output_files <- c("cart19_goi_impact.csv", "cart19_gene_stats.csv")

if(
  !all(sapply(output_files, function(x) file.exists(file.path(outputDir, x))))
){

  # Load genes of interest -------------------------------------------------------
  top_one_pc_sonicAbund <- as.character(read.table(
    file.path(outputDir, "top_one_percent_of_gene_ids_by_sonicAbund.txt"),
    header = FALSE, sep = "\t")[1,1])
  
  top_ten_pc_sonicAbund <- as.character(read.table(
    file.path(outputDir, "top_ten_percent_of_gene_ids_by_sonicAbund.txt"),
    header = FALSE, sep = "\t")[1,1])
  
  ## Identify genes of interest from clusters and clonal expansions --------------
  cluster_gene_list <- unique(unlist(sapply(
    red_clusters$genes_in_cluster, 
    function(x) unlist(strsplit(x, ", "))), use.names = FALSE))
  
  top_clonal_genes <- unique(unlist(strsplit(top_one_pc_sonicAbund, ", ")))
  top_ten_pc_clonal_genes <- unique(unlist(strsplit(top_ten_pc_sonicAbund, ", ")))
  consolidated_genes <- union(cluster_gene_list, top_clonal_genes)
  
  consol_gene_ranges <- unlist(GRangesList(lapply(consolidated_genes, function(x){
    gene <- GenomicRanges::reduce(refGenes[which(refGenes$name2 == x)])
    gene$name <- x
    gene
  })))
  
  tdn_gr <- format_sites(tdn_sites)
  tp_gr <- format_sites(timepoint_sites)
  
  gene_tdn_sites <- GRangesList(lapply(consol_gene_ranges, function(x){
    tdn_gr[subjectHits(findOverlaps(
      x, tdn_gr, maxgap = 5000, ignore.strand = TRUE))]
  }))
  
  gene_tp_sites <- GRangesList(lapply(consol_gene_ranges, function(x){
    tp_gr[subjectHits(findOverlaps(
      x, tp_gr, maxgap = 5000, ignore.strand = TRUE))]
  }))
  
  goi_stats <- data.frame(
    "gene_name" = consol_gene_ranges$name,
    "gene_ort" = as.character(strand(consol_gene_ranges)),
    "TDN_num_patients" = sapply(
      gene_tdn_sites, function(x) length(unique(x$patient))
    ),
    "TP_num_patients" = sapply(
      gene_tp_sites, function(x) length(unique(x$patient))
    ),
    "TDN_num_sites" = sapply(gene_tdn_sites, length),
    "TP_num_sites" = sapply(gene_tp_sites, length),
    "TDN_sum_abund" = sapply(gene_tdn_sites, function(x) sum(x$estAbund)),
    "TP_sum_abund" = sapply(gene_tp_sites, function(x) sum(x$estAbund))) %>%
    dplyr::mutate("TDN_same_ort" = mapply(function(x, ort){
      length(which(as.character(strand(x)) == ort))}, 
      gene_tdn_sites, gene_ort)) %>%
    dplyr::mutate("TDN_oppo_ort" = mapply(function(x, ort){
      length(which(as.character(strand(x)) != ort))}, 
      gene_tdn_sites, gene_ort)) %>%
    dplyr::mutate("TP_same_ort" = mapply(function(x, ort){
      length(which(as.character(strand(x)) == ort))}, 
      gene_tp_sites, gene_ort)) %>%
    dplyr::mutate("TP_oppo_ort" = mapply(function(x, ort){
      length(which(as.character(strand(x)) != ort))}, 
      gene_tp_sites, gene_ort))
  
  goi_stats$ort_fisher_test <- sapply(seq_len(nrow(goi_stats)), function(i){
    m <- matrix(as.numeric(goi_stats[i,9:12]), nrow = 2, ncol = 2)
    f <- fisher.test(m)
    f$p.value
  })
  
  goi_stats$abund_gini <- sapply(gene_tp_sites, function(x){
    df <- GenomicRanges::as.data.frame(x, row.names = NULL) %>%
      dplyr::select(patient, posid, estAbund) %>%
      dplyr::group_by(patient, posid) %>%
      dplyr::summarise(abund = sum(estAbund)) %>%
      dplyr::ungroup(.) %>% 
      as.data.frame(.)
    pop_calcs(df$abund, calc = "gini")
  })
  
  # Annotate which genes are overlapping -----------------------------------------
  goi_stats$Overlapping_Genes <- sapply(goi_stats$gene_name, function(gene){
    gene_range <- refGenes[which(refGenes$name2 == gene)]
    ovlp_hits <- refGenes[subjectHits(findOverlaps(
      gene_range, refGenes, ignore.strand = TRUE))]$name2
    paste(unique(ovlp_hits[which(ovlp_hits != gene)]), collapse = ", ")
  })
  
  # Annotate with which list the gene belongs ------------------------------------
  goi_stats <- goi_stats %>%
    dplyr::mutate(
      "On_Onco_List" = gene_name %in% oncoGenes,
      "Top_1pc_Abund" = gene_name %in% top_clonal_genes,
      "Top_10pc_Abund" = gene_name %in% top_ten_pc_clonal_genes,
      "Within_Cluster" = gene_name %in% cluster_gene_list
    )
  
  # Identify which clusters the genes belong to ----------------------------------
  all_clusters <- cart19_clusters
  all_clusters$genes_in_cluster <- scan_genes(all_clusters, refGenes)
  
  goi_stats$Cluster_target.min <- sapply(goi_stats$gene_name, function(gene){
    clusters <- all_clusters[grep(gene, all_clusters$genes_in_cluster)]
    if(length(clusters) == 0){
      target <- NA
    }else{
      target <- min(clusters$target.min)
    }
    target
  })
  
  goi_stats <- dplyr::filter(goi_stats, TP_num_sites > 0)
  
  goi_stats <- goi_stats[order(goi_stats$abund_gini, decreasing = TRUE),]
  
  goi_stats$gene_name <- paste0("'", goi_stats$gene_name)
  
  # Nearest gene approach for binning the data -----------------------------------
  refGenesOrt <- as.data.frame(refGenes, row.names = NULL) %>%
    dplyr::select(name2, strand) %>%
    dplyr::group_by(name2) %>%
    dplyr::summarise(
      gene_ort = names(sort(table(strand), decreasing = TRUE)[1])
    ) %>%
    dplyr::ungroup() %>% 
    dplyr::rename(gene = name2) %>%
    as.data.frame()
  
  gene_stats <- as.data.frame(cond_uniq_sites, row.names = NULL) %>%
    dplyr::select(
      specimen, patient, celltype, timepoint, estAbund, relAbund, posid, 
      gene_id_wo_annot, in_gene, nearest_geneDist, start, strand, seqnames) %>%
    dplyr::mutate(
      type = ifelse(timepoint == "d0", "TDN", "TP"),
      gene_name = ifelse(is.na(gene_id_wo_annot), posid, gene_id_wo_annot),
      nearest_geneDist = ifelse(is.na(nearest_geneDist), 0, nearest_geneDist)) %>%
    dplyr::rename(ort = strand) %>%
    dplyr::group_by(type, patient, gene_name, posid, ort) %>%
    dplyr::summarise(
      long_count = n_distinct(timepoint),
      peak_abund = max(estAbund),
      peak_relAbund = max(relAbund),
      max_geneDist = max(ifelse(in_gene == FALSE, abs(nearest_geneDist), 0)),
      loci = unique(start)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(., refGenesOrt, by = c("gene_name" = "gene")) %>%
    dplyr::mutate(gene_ort = ifelse(is.na(gene_ort), "*", gene_ort)) %>%
    dplyr::group_by(type, gene_name, gene_ort) %>%
    dplyr::summarise(
      num_patients = n_distinct(patient),
      num_sites = n_distinct(posid),
      long_count = max(long_count),
      abund_gini = pop_calcs(peak_abund, calc = "gini"),
      peak_abund = max(peak_abund),
      peak_relAbund = max(peak_relAbund),
      max_geneDist = max(max_geneDist),
      min_loci = min(loci),
      max_loci = max(loci),
      same_ort = length(which(ort == gene_ort)),
      oppo_ort = length(which(ort != gene_ort))) %>%
    dplyr::ungroup() %>% 
    as.data.frame() %>%
    melt(
      id.vars = c("type", "gene_name", "gene_ort"), 
      measure.vars = c(
        "num_patients", "num_sites", "long_count", 
        "peak_abund", "peak_relAbund", "abund_gini", "max_geneDist",
        "min_loci", "max_loci", "same_ort", "oppo_ort")) %>%
    dcast(gene_name + gene_ort ~ type + variable, fill = 0) %>%
    dplyr::group_by(gene_name, gene_ort) %>%
    dplyr::mutate(
      long_count = TP_long_count,
      abund_gini = TP_abund_gini,
      max_geneDist = max(TDN_max_geneDist, TP_max_geneDist),
      genomic_range = max(c(TDN_max_loci, TP_max_loci)) - 
        min(c(TDN_min_loci, TP_min_loci)),
      ort_fisher_test = fisher.test(
        matrix(
          c(TDN_same_ort, TDN_oppo_ort, TP_same_ort, TP_oppo_ort), 
          nrow = 2, ncol = 2))$p.value) %>%
    dplyr::select(
      gene_name, gene_ort, TDN_num_patients, TP_num_patients, TDN_num_sites, 
      TP_num_sites, TDN_peak_abund, TP_peak_abund, TDN_peak_relAbund, 
      TP_peak_relAbund, abund_gini, long_count, max_geneDist, genomic_range, 
      ort_fisher_test) %>%
    dplyr::ungroup() %>% 
    as.data.frame()
  
  ## Annotate with which list the gene belongs ===================================
  gene_stats <- gene_stats %>%
    dplyr::mutate(
      "On_Onco_List" = gene_name %in% oncoGenes,
      "Top_1pc_Abund" = gene_name %in% top_clonal_genes,
      "Top_10pc_Abund" = gene_name %in% top_ten_pc_clonal_genes,
      "Within_Cluster" = gene_name %in% cluster_gene_list
    )
  
  ## Identify which clusters the genes belong to =================================
  gene_stats$Cluster_target.min <- sapply(
    gene_stats$gene_name, function(gene){
      clusters <- all_clusters[grep(gene, all_clusters$genes_in_cluster)]
      if(length(clusters) == 0){
        target <- NA
      }else{
        target <- min(clusters$target.min)
      }
      target
    })
  
  write.csv(
    goi_stats,
    file = file.path(outputDir, "cart19_goi_impact.csv"),
    quote = TRUE,
    row.names = FALSE
  )
  
  write.csv(
    gene_stats,
    file = file.path(outputDir, "cart19_gene_stats.csv"),
    quote = TRUE,
    row.names = FALSE
  )

}