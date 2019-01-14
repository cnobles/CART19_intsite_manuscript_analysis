#### This R-script is for running scan statistics analysises across the CART19
#### data set. It depends on the setup environment of 
#### cart19_intsite_analysis.Rmd.
output_files <- c(
  "cart19_tdn_enriched_clusters.rds", "cart19_timepoint_enriched_clusters.rds",
  "cart19_clusters.rds", "cart19_cluster_gene_info.csv"
)

if(
  !all(sapply(output_files, function(x) file.exists(file.path(outputDir, x))))
){
  
  # Transduction vs All Timepoints -----------------------------------------------
  scan_sites <- gintools:::scan_format(
    tdn_sites, timepoint_sites, grouping = "patient")
  
  scanned_clus <- gRxCluster(
    object = scan_sites$chr, 
    starts = scan_sites$pos, 
    group = scan_sites$grp,
    kvals = c(10L:50L),
    nperm = 100L,
    cutpt.tail.expr = critVal.target(k, n, target = 7.5, posdiff = x),
    cutpt.filter.expr = apply(x, 2, quantile, probs = 0.35, na.rm = TRUE))
  
  scan_summary <- gRxSummary(scanned_clus)
  
  scanned_clus <- annotate_scan_clusters(
    scanned_clus, tdn_sites, timepoint_sites, refGenes)
  
  scanned_clus <- scanned_clus[
    scanned_clus$n_patients_tdn > 1 | scanned_clus$n_patients_tp > 1]
  
  enriched_tdn_clus <- scanned_clus[
    (scanned_clus$n_sites_tdn/length(tdn_sites)) > 
      (scanned_clus$n_sites_tp/length(timepoint_sites))]
  
  enriched_timepoint_clus <- scanned_clus[
    (scanned_clus$n_sites_tdn/length(tdn_sites)) < 
      (scanned_clus$n_sites_tp/length(timepoint_sites))]
  
  saveRDS(
    enriched_tdn_clus, 
    file = file.path(outputDir, "cart19_tdn_enriched_clusters.rds"))
  
  write.table(
    paste(
      unlist(strsplit(enriched_tdn_clus$genes_in_cluster, ", ")), 
      collapse = ", "),
    file = file.path(outputDir, "genes_in_tdn_clusters.txt"), 
    quote = FALSE, 
    row.names = FALSE, 
    col.names = FALSE
  )
  
  saveRDS(
    enriched_timepoint_clus,
    file = file.path(outputDir, "cart19_timepoint_enriched_clusters.rds")
  )
  
  write.table(
    paste(
      unlist(strsplit(enriched_timepoint_clus$genes_in_cluster, ", ")), 
      collapse = ", "),
    file = file.path(outputDir, "genes_in_timepoint_clusters.txt"), 
    quote = FALSE, 
    row.names = FALSE, 
    col.names = FALSE
  )
  
  # Scans for clusters with orientation biases -----------------------------------
  pos_strand_tdn_sites <- tdn_sites[strand(tdn_sites) == "+"]
  pos_strand_timepoint_sites <- timepoint_sites[strand(timepoint_sites) == "+"]
  
  neg_strand_tdn_sites <- tdn_sites[strand(tdn_sites) == "-"]
  neg_strand_timepoint_sites <- timepoint_sites[strand(timepoint_sites) == "-"]
  
  pos_strand_scan_sites <- gintools:::scan_format(
    pos_strand_tdn_sites, pos_strand_timepoint_sites, grouping = "patient")
  neg_strand_scan_sites <- gintools:::scan_format(
    neg_strand_tdn_sites, neg_strand_timepoint_sites, grouping = "patient")
  
  pos_scanned_clus <- gRxCluster(
    object = pos_strand_scan_sites$chr, 
    starts = pos_strand_scan_sites$pos, 
    group = pos_strand_scan_sites$grp,
    kvals = c(10L:50L),
    nperm = 100L,
    cutpt.tail.expr = critVal.target(k, n, target = 10, posdiff = x),
    cutpt.filter.expr = apply(x, 2, quantile, probs = 0.25, na.rm = TRUE))
  
  pos_scan_summary <- gRxSummary(pos_scanned_clus)
  
  neg_scanned_clus <- gRxCluster(
    object = neg_strand_scan_sites$chr, 
    starts = neg_strand_scan_sites$pos, 
    group = neg_strand_scan_sites$grp,
    kvals = c(10L:50L),
    nperm = 100L,
    cutpt.tail.expr = critVal.target(k, n, target = 10, posdiff = x),
    cutpt.filter.expr = apply(x, 2, quantile, probs = 0.25, na.rm = TRUE))
  
  neg_scan_summary <- gRxSummary(neg_scanned_clus)
  
  ## Clusters are compared to all sites, not just the same strand
  pos_scanned_clus <- annotate_scan_clusters(
    pos_scanned_clus, tdn_sites[strand(tdn_sites) == "+"], 
    timepoint_sites[strand(timepoint_sites) == "+"], refGenes)
  
  neg_scanned_clus <- annotate_scan_clusters(
    neg_scanned_clus, tdn_sites[strand(tdn_sites) == "-"], 
    timepoint_sites[strand(timepoint_sites) == "-"], refGenes)
  
  pos_scanned_clus <- pos_scanned_clus[
    pos_scanned_clus$n_patients_tdn > 1 | pos_scanned_clus$n_patients_tp > 1]
  
  neg_scanned_clus <- neg_scanned_clus[
    neg_scanned_clus$n_patients_tdn > 1 | neg_scanned_clus$n_patients_tp > 1]
  
  pos_strand_timepoint_enriched_clus <- pos_scanned_clus[
    (pos_scanned_clus$n_sites_tp/length(pos_strand_timepoint_sites)) > 
      (pos_scanned_clus$n_sites_tdn/length(pos_strand_tdn_sites))]
  
  neg_strand_timepoint_enriched_clus <- neg_scanned_clus[
    (neg_scanned_clus$n_sites_tp/length(neg_strand_timepoint_sites)) > 
      (neg_scanned_clus$n_sites_tdn/length(neg_strand_tdn_sites))]
  
  # Scans for clusters in higher abundance ---------------------------------------
  # over lower abundance in other timepoints
  cutoff <- 2.5 # percent
  
  tp_high_abund_posids <- get_top_sites(timepoint_sites, cutoff, "estAbund")
  timepoint_sites$abund_status <- ifelse(
    timepoint_sites$posid %in% tp_high_abund_posids,
    "High Abundance", "Low Abundance")
  
  high_sites <- timepoint_sites[
    timepoint_sites$abund_status == "High Abundance"]
  low_sites <- timepoint_sites[
    timepoint_sites$abund_status == "Low Abundance"]
  
  high_low_sites <- gintools:::scan_format(
    low_sites, high_sites, grouping = "patient")
  
  high_low_clus <- gRxCluster(
    object = high_low_sites$chr, 
    starts = high_low_sites$pos, 
    group = high_low_sites$grp,
    kvals = c(10L:35L),
    nperm = 100L,
    cutpt.tail.expr = critVal.target(k, n, target = 2, posdiff = x),
    cutpt.filter.expr = apply(x, 2, quantile, probs = 0.15, na.rm = TRUE))
  
  high_low_summary <- gRxSummary(high_low_clus)
  
  high_low_clus <- high_low_clus[width(high_low_clus) > 1]
  high_low_clus$n_sites_high <- scan_sites_count(high_low_clus, high_sites)
  high_low_clus$n_sites_low <- scan_sites_count(high_low_clus, low_sites)
  high_low_clus$n_patients_high <- scan_patients(high_low_clus, high_sites)
  high_low_clus$n_patients_low <- scan_patients(high_low_clus, low_sites)
  high_low_clus$genes_in_cluster <- scan_genes(high_low_clus, refGenes)
  high_low_clus$in_gene_ort_high <- scan_orientation(high_low_clus, high_sites)
  high_low_clus$in_gene_ort_low <- scan_orientation(high_low_clus, low_sites)
  high_low_clus$ort_fisher_test <- scan_fisher_test_ort(
    high_low_clus$in_gene_ort_high, high_low_clus$in_gene_ort_low)
  enriched_high_clus <- high_low_clus[
    (high_low_clus$n_sites_high/length(high_sites)) > 
      (high_low_clus$n_sites_low/length(low_sites))]
  enriched_low_clus <- high_low_clus[
    (high_low_clus$n_sites_high/length(high_sites)) < 
      (high_low_clus$n_sites_low/length(low_sites))]
  
  # Combine all clusters into the same frame 
  # and reduce to find all unique clusters
  clusters <- list(
    "timepoint_clusters" = enriched_timepoint_clus,
    "positive_strand_clusters" = pos_strand_timepoint_enriched_clus,
    "negative_strand_clusters" = neg_strand_timepoint_enriched_clus,
    "high_abund_clusters" = enriched_high_clus)
  
  cart19_clusters <- unlist(GRangesList(lapply(
    1:4,
    function(i){
      clus_name <- c("Timepoint", "Pos-Strand", "Neg-Strand", "Abundance")[i]
      cluster_group <- granges(clusters[[i]])
      cluster_group$clus.origin <- clus_name
      cluster_group$target.min <- clusters[[i]]$target.min
      cluster_group
    }
  )))
  
  red_clusters <- GenomicRanges::reduce(cart19_clusters, with.revmap = TRUE)
  
  red_clusters$cluster_origin <- sapply(red_clusters$revmap, function(x){
    paste(unique(cart19_clusters[x]$clus.origin), collapse = ", ")})
  
  red_clusters$revmap <- NULL
  
  red_clusters <- annotate_scan_clusters(
    red_clusters, tdn_sites, timepoint_sites, refGenes)
  
  red_df <- dplyr::select(as.data.frame(red_clusters), -strand)
  
  names(red_df) <- c(
    "Chr", "Start", "End", "Width", "Cluster Origin", "TDN: Num. Sites", 
    "TP: Num. Sites", "TDN: Abundance Sum", "TP: Abundance Sum", 
    "TDN: Num. Patients", "TP: Num. Patients", "Genes In Cluster", 
    "TDN: Within Gene Ort.", "TP: Within Gene Ort.", 
    "Orientation Fisher Test P-value")
  
  saveRDS(
    c(clusters, 
      list("all_clusters" = cart19_clusters, "red_clusters" = red_clusters)), 
    file = file.path(outputDir, "cart19_clusters.rds"))
  
  write.csv(
    red_df,
    file = file.path(outputDir, "cart19_clusters.csv"),
    quote = TRUE,
    row.names = FALSE
  )
  
  # Compile Fisher Exact tests ---------------------------------------------------
  # each gene found within clusters to determine orientation bias
  cluster_gene_list <- unique(unlist(sapply(
    red_clusters$genes_in_cluster, function(x)
      {unlist(strsplit(x, ", "))}
    ), use.names = FALSE))
  
  gene_ranges <- unlist(GRangesList(lapply(cluster_gene_list, function(x){
    gene <- GenomicRanges::reduce(refGenes[which(refGenes$name2 == x)])
    gene$name <- x
    gene
  })))
  
  mod_gene_list <- gene_ranges$name
  
  tdn_gr <- format_sites(tdn_sites)
  tp_gr <- format_sites(timepoint_sites)
  
  gene_tdn_sites <- GRangesList(lapply(gene_ranges, function(x){
    tdn_gr[subjectHits(findOverlaps(x, tdn_gr, ignore.strand = TRUE))]
  }))
  
  gene_tp_sites <- GRangesList(lapply(gene_ranges, function(x){
    tp_gr[subjectHits(findOverlaps(x, tp_gr, ignore.strand = TRUE))]
  }))
  
  cluster_genes <- data.frame(
    "gene_name" = mod_gene_list,
    "gene_ort" = as.character(strand(gene_ranges)),
    "num_tdn_patients" = sapply(
      gene_tdn_sites, function(x) length(unique(x$patient))
    ),
    "num_tp_patients" = sapply(
      gene_tp_sites, function(x) length(unique(x$patient))
    ),
    "num_tdn_sites" = sapply(gene_tdn_sites, length),
    "num_tp_sites" = sapply(gene_tp_sites, length),
    "tdn_sum_abund" = sapply(gene_tdn_sites, function(x) sum(x$estAbund)),
    "tp_sum_abund" = sapply(gene_tp_sites, function(x) sum(x$estAbund))) %>%
    dplyr::mutate(tdn_same = mapply(function(x, ort){
      length(which(as.character(strand(x)) == ort))}, 
      gene_tdn_sites, gene_ort)) %>%
    dplyr::mutate(tdn_oppo = mapply(function(x, ort){
      length(which(as.character(strand(x)) != ort))}, 
      gene_tdn_sites, gene_ort)) %>%
    dplyr::mutate(tp_same = mapply(function(x, ort){
      length(which(as.character(strand(x)) == ort))}, 
      gene_tp_sites, gene_ort)) %>%
    dplyr::mutate(tp_oppo = mapply(function(x, ort){
      length(which(as.character(strand(x)) != ort))}, 
      gene_tp_sites, gene_ort))
  
  cluster_genes$fisher_p_value <- sapply(
    seq_len(nrow(cluster_genes)), function(i){
      m <- matrix(as.numeric(cluster_genes[i,9:12]), nrow = 2, ncol = 2)
      f <- fisher.test(m)
      f$p.value
    }
  )
  
  cluster_genes$tp_abund_gini <- sapply(gene_tp_sites, function(x){
    df <- GenomicRanges::as.data.frame(x, row.names = NULL) %>%
      dplyr::select(patient, posid, estAbund) %>%
      dplyr::group_by(patient, posid) %>%
      dplyr::summarise(abund = sum(estAbund)) %>%
      ungroup(.) %>% 
      as.data.frame(.)
    pop_calcs(df$abund, calc = "gini")
  })
  
  write.csv(
    cluster_genes,
    file = file.path(outputDir, "cart19_cluster_gene_info.csv"),
    quote = TRUE,
    row.names = FALSE
  )

}