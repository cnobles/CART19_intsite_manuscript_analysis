#### This script is meant to tabulate the intergenic impact of integration 
#### sites. This only considers integrations outside genes, but every region 
#### with an integation across the entire dataset.

output_files <- c("cart19_intergene_impact.rds", "cart19_intergene_impact.csv")

if(
  !all(sapply(output_files, function(x) file.exists(file.path(outputDir, x))))
){

  # Compile Data for each gene ---------------------------------------------------
  mod_refGenes <- unlist(GenomicRanges::reduce(split(
    refGenes, refGenes$name2), min.gapwidth = 100000L))
  mod_refGenes$name <- names(mod_refGenes)
  strand(mod_refGenes) <- "*"
  mod_gapRegions <- gaps(mod_refGenes)
  mod_gapRegions <- mod_gapRegions[
    seqnames(mod_gapRegions) %in% paste0("chr", c(1:22, "X", "Y", "M"))]
  mod_gapRegions <- mod_gapRegions[strand(mod_gapRegions) == "*"]
  mod_gapRegions_names <- getNearestFeature(
    mod_gapRegions, mod_refGenes, colnam = "nearestTSS", side = "5p", 
    feature.colnam = "name")
  mod_gapRegions_names <- CharacterList(
    strsplit(mod_gapRegions_names$X5pnearestTSS, ","))
  mod_gapRegions$name <- unlist(mod_gapRegions_names)[
    start(mod_gapRegions_names@partitioning)]
  mod_gapRegions <- unique_granges(mod_gapRegions)
  
  ## Annotate Cluster FDR across gene set ========================================
  mod_gapXclusters <- suppressWarnings(
    findOverlaps(mod_gapRegions, cart19_clusters))
  mod_gapRegions$Within_Cluster <- FALSE
  mod_gapRegions$Cluster_target.min <- as.numeric(NA)
  mod_gapRegions[queryHits(mod_gapXclusters)]$Within_Cluster <- TRUE
  mod_gapRegions[queryHits(mod_gapXclusters)]$Cluster_target.min <- 
    cart19_clusters$target.min[subjectHits(mod_gapXclusters)]
  
  ## List of all genes ===========================================================
  region_list <- get_loci_id(mod_gapRegions)
  region_onco_hits <- findOverlaps(
    mod_gapRegions, mod_refGenes[mod_refGenes$name %in% oncoGenes], 
    maxgap = 1L, ignore.strand = TRUE)
  region_near_onco <- region_list[unique(queryHits(region_onco_hits))]
  region_1pct_hits <- findOverlaps(
    mod_gapRegions, 
    timepoint_sites[
      timepoint_sites$estAbund >= quantile(
        timepoint_sites$estAbund, probs = 0.99)],
    ignore.strand = TRUE)
  region_1pct <- region_list[unique(queryHits(region_1pct_hits))]
  region_10pct_hits <- findOverlaps(
    mod_gapRegions, 
    timepoint_sites[
      timepoint_sites$estAbund >= quantile(
        timepoint_sites$estAbund, probs = 0.90)],
    ignore.strand = TRUE)
  region_10pct <- region_list[unique(queryHits(region_10pct_hits))]
  region_names <- mod_gapRegions$name
  
  ## Format TP and TDN sites for analysis ========================================
  tdn_gr <- format_sites(tdn_sites)
  tp_gr <- format_sites(timepoint_sites)
  
  # Determine which integrations were within a transcriptional unit --------------
  tdn_in_region <- findOverlaps(
    mod_gapRegions, tdn_gr, maxgap = 0L, ignore.strand = TRUE)
  
  tp_in_region <- findOverlaps(
    mod_gapRegions, tp_gr, maxgap = 0L, ignore.strand = TRUE)
  
  # Assemble stats for each transcriptional unit ---------------------------------
  stats_tdn <- as.data.frame(
      tdn_gr[subjectHits(tdn_in_region)], row.names = NULL) %>%
    dplyr::mutate(
      loci = get_loci_id(mod_gapRegions[queryHits(tdn_in_region)]),
      region_name = as.character(mod_gapRegions$name[queryHits(tdn_in_region)]),
      strand = as.character(strand)) %>%
    dplyr::select(
      loci, region_name, posid, strand, estAbund, relAbund, patient) %>%
    dplyr::group_by(loci, region_name) %>%
    dplyr::summarise(
      "TDN_num_patients" = n_distinct(patient),
      "TDN_num_sites" = n_distinct(posid),
      "TDN_sum_abund" = sum(estAbund),
      "TDN_peak_abund" = max(estAbund),
      "TDN_peak_relAbund" = max(relAbund),
      "tdn_pos" = sum(as.integer(strand == "+")),
      "tdn_neg" = sum(as.integer(strand == "-")))
  
  stats_tp <- as.data.frame(
      tp_gr[subjectHits(tp_in_region)], row.names = NULL) %>%
    dplyr::mutate(
      loci = get_loci_id(mod_gapRegions[queryHits(tp_in_region)]),
      region_name = as.character(mod_gapRegions$name[queryHits(tp_in_region)]),
      strand = as.character(strand),
      timepoint = convert_time(as.character(timepoint))) %>%
    dplyr::select(
      loci, region_name, posid, strand, 
      estAbund, relAbund, patient, timepoint) %>%
    dplyr::group_by(loci, region_name, patient, posid, strand) %>%
    dplyr::summarise(
      "long_count" = n_distinct(timepoint),
      "first_time" = min(timepoint),
      "last_time" = max(timepoint),
      "sum_abund" = sum(estAbund),
      "peak_abund" = max(estAbund),
      "peak_relAbund" = max(relAbund)) %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(loci, region_name) %>%
    dplyr::summarise(
      "TP_num_patients" = n_distinct(patient),
      "long_count" = max(long_count),
      "max_time" = max(last_time),
      "max_span" = max(last_time - first_time),
      "TP_num_sites" = n_distinct(posid),
      "TP_sum_abund" = sum(sum_abund),
      "TP_peak_abund" = max(peak_abund),
      "TP_peak_relAbund" = max(peak_relAbund),
      "abund_gini" = pop_calcs(peak_abund, calc = "gini"),
      "tp_pos" = sum(as.integer(strand == "+")),
      "tp_neg" = sum(as.integer(strand == "-")))
  
  region_stats <- data.frame(
      "loci" = get_loci_id(mod_gapRegions),
      "region_name" = as.character(mod_gapRegions$name),
      "region_width" = (width(mod_gapRegions) + 10000) / 1000,
      "Within_Cluster" = mod_gapRegions$Within_Cluster,
      "Cluster_target.min" = mod_gapRegions$Cluster_target.min) %>%
    dplyr::full_join(., stats_tdn, by = c("loci", "region_name")) %>%
    dplyr::full_join(., stats_tp, by = c("loci", "region_name")) %>%
    dplyr::mutate(
      "Near_Onco" = loci %in% region_near_onco,
      "Top_1pc_Abund" = loci %in% region_1pct,
      "Top_10pc_Abund" = loci %in% region_10pct
    )
  
  region_stats[is.na(region_stats)] <- 0
  total_tdn_sites <- summaries$timepoint$num_UniqueSites[
    summaries$timepoint$timepoint == "d0"]
  total_pat_sites <- sum(summaries$timepoint$num_UniqueSites[
    summaries$timepoint$timepoint != "d0"])
  
  region_stats <- dplyr::filter(
      region_stats, TDN_num_patients > 0 | TP_num_patients > 0) %>%
    dplyr::group_by(loci, region_name) %>%
    dplyr::mutate(
      "ort_fisher_test" = fisher.test(matrix(
          c(tdn_pos, tdn_neg, tp_pos, tp_neg), 
        nrow = 2, ncol = 2))$p.value
    ) %>% 
    dplyr::ungroup() %>%
    dplyr::select(
      loci, region_name, region_width, TDN_num_patients, TP_num_patients, 
      TDN_num_sites, TP_num_sites, TDN_peak_abund, TP_peak_abund, 
      TDN_peak_relAbund, TP_peak_relAbund, TDN_sum_abund, TP_sum_abund, 
      long_count, max_time, max_span, abund_gini, ort_fisher_test, Near_Onco,
      Top_1pc_Abund, Top_10pc_Abund, Within_Cluster, Cluster_target.min) %>%
    dplyr::mutate(
      TDN_freq = TDN_num_sites / (region_width * total_tdn_sites),
      TP_freq = TP_num_sites/(region_width * total_pat_sites),
      freq_diff = (TP_freq - TDN_freq),
      pct_chg = 100 * freq_diff / TDN_freq) %>%
    as.data.frame()
  
  saveRDS(
    region_stats,
    file = file.path(outputDir, "cart19_intergene_impact.rds"))
  
  write.csv(
    region_stats,
    file = file.path(outputDir, "cart19_intergene_impact.csv"),
    quote = TRUE,
    row.names = FALSE
  )

}