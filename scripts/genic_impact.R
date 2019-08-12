#### This only considers integrations within genes, but every gene with an 
#### integation across the entire dataset.

output_files <- c(
  "cart19_gene_impact.rds", "cart19_gene_impact.csv", 
  "cart19_cll_gene_impact.csv", "cart19_cll_gene_impact.rds",
  "cart19_all_gene_impact.csv", "cart19_all_gene_impact.rds",
  "cart19_cr_gene_impact.csv", "cart19_cr_gene_impact.rds"
)

if(
  !all(sapply(output_files, function(x) file.exists(file.path(outputDir, x))))
){

  # Compile Data for each gene -------------------------------------------------
  mod_refGenes <- unlist(GenomicRanges::reduce(split(
    refGenes, refGenes$name2), min.gapwidth = 100000L))
  
  mod_refGenes$name <- names(mod_refGenes)
  
  ## Annotate Cluster FDR across gene set ======================================
  if(trial == "CART19"){
    mod_refXclusters <- suppressWarnings(
      findOverlaps(mod_refGenes, cart19_clusters))
    mod_refGenes$Within_Cluster <- FALSE
    mod_refGenes$Cluster_target.min <- as.numeric(NA)
    mod_refGenes[queryHits(mod_refXclusters)]$Within_Cluster <- TRUE
    mod_refGenes[queryHits(mod_refXclusters)]$Cluster_target.min <- 
      cart19_clusters$target.min[subjectHits(mod_refXclusters)]
  }
    
  ## List of all genes =========================================================
  gene_list <- mod_refGenes$name
  
  ## Format TP and TDN sites for analysis ======================================
  tdn_gr <- format_sites(tdn_sites)
  tp_gr <- format_sites(timepoint_sites)
  
  # Determine which integrations were within a transcriptional unit ------------
  tdn_in_gene <- findOverlaps(
    mod_refGenes, tdn_gr, maxgap = 5000, ignore.strand = TRUE)
  
  tp_in_gene <- findOverlaps(
    mod_refGenes, tp_gr, maxgap = 5000, ignore.strand = TRUE)
  
  # Assemble stats for each transcriptional unit -------------------------------
  ## From all patients
  stats_tdn <- as.data.frame(
      tdn_gr[subjectHits(tdn_in_gene)], row.names = NULL) %>%
    dplyr::mutate(
      loci = get_loci_id(mod_refGenes[queryHits(tdn_in_gene)]),
      gene_name = as.character(mod_refGenes$name[queryHits(tdn_in_gene)]),
      gene_ort = as.character(strand(mod_refGenes[queryHits(tdn_in_gene)])),
      strand = as.character(strand)) %>%
    dplyr::select(
      loci, gene_name, gene_ort, posid, strand, estAbund, relAbund, patient) %>%
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::summarise(
      "TDN_num_patients" = n_distinct(patient),
      "TDN_num_pats_CR" = n_distinct(patient[patient %in% CR_pats]),
      "TDN_num_pats_NR" = n_distinct(patient[patient %in% NR_pats]),
      "TDN_num_sites" = n_distinct(posid),
      "TDN_num_sites_CR" = n_distinct(posid[patient %in% CR_pats]),
      "TDN_num_sites_NR" = n_distinct(posid[patient %in% NR_pats]),
      "TDN_sum_abund" = sum(estAbund),
      "TDN_sum_abund_CR" = sum(estAbund[patient %in% CR_pats]),
      "TDN_sum_abund_NR" = sum(estAbund[patient %in% NR_pats]),
      "TDN_peak_abund" = max(estAbund),
      "TDN_peak_abund_CR" = max(estAbund[patient %in% CR_pats]),
      "TDN_peak_abund_CR" = ifelse(TDN_peak_abund_CR < 0, 0, TDN_peak_abund_CR),
      "TDN_peak_abund_NR" = max(estAbund[patient %in% NR_pats]),
      "TDN_peak_abund_NR" = ifelse(TDN_peak_abund_NR < 0, 0, TDN_peak_abund_NR),
      "TDN_peak_relAbund" = max(relAbund),
      "tdn_same" = sum(as.integer(strand == gene_ort)),
      "tdn_same_CR" = sum(as.integer(strand == gene_ort)[patient %in% CR_pats]),
      "tdn_same_NR" = sum(as.integer(strand == gene_ort)[patient %in% NR_pats]),
      "tdn_oppo" = sum(as.integer(strand != gene_ort)),
      "tdn_oppo_CR" = sum(as.integer(strand != gene_ort)[patient %in% CR_pats]),
      "tdn_oppo_NR" = sum(as.integer(strand != gene_ort)[patient %in% NR_pats]))
  
  stats_tp <- as.data.frame(
      tp_gr[subjectHits(tp_in_gene)], row.names = NULL) %>%
    dplyr::mutate(
      loci = get_loci_id(mod_refGenes[queryHits(tp_in_gene)]),
      gene_name = as.character(mod_refGenes$name[queryHits(tp_in_gene)]),
      gene_ort = as.character(strand(mod_refGenes[queryHits(tp_in_gene)])),
      strand = as.character(strand),
      timepoint = convert_time(as.character(timepoint))) %>%
    dplyr::select(
      loci, gene_name, gene_ort, posid, strand, 
      estAbund, relAbund, patient, timepoint) %>%
    dplyr::group_by(loci, gene_name, gene_ort, patient, posid, strand) %>%
    dplyr::summarise(
      "long_count" = n_distinct(timepoint),
      "first_time" = min(timepoint),
      "last_time" = max(timepoint),
      "sum_abund" = sum(estAbund),
      "peak_abund" = max(estAbund),
      "peak_relAbund" = max(relAbund)) %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::summarise(
      "TP_num_patients" = n_distinct(patient),
      "TP_num_pats_CR" = n_distinct(patient[patient %in% CR_pats]),
      "TP_num_pats_NR" = n_distinct(patient[patient %in% NR_pats]),
      "long_count" = max(long_count),
      "max_time" = max(last_time),
      "max_span" = max(last_time - first_time),
      "TP_num_sites" = n_distinct(posid),
      "TP_num_sites_CR" = n_distinct(posid[patient %in% CR_pats]),
      "TP_num_sites_NR" = n_distinct(posid[patient %in% NR_pats]),
      "TP_sum_abund" = sum(sum_abund),
      "TP_sum_abund_CR" = sum(sum_abund[patient %in% CR_pats]),
      "TP_sum_abund_NR" = sum(sum_abund[patient %in% NR_pats]),
      "TP_peak_abund" = max(peak_abund),
      "TP_peak_abund_CR" = max(peak_abund[patient %in% CR_pats]),
      "TP_peak_abund_CR" = ifelse(TP_peak_abund_CR < 0, 0, TP_peak_abund_CR),
      "TP_peak_abund_NR" = max(peak_abund[patient %in% NR_pats]),
      "TP_peak_abund_NR" = ifelse(TP_peak_abund_NR < 0, 0, TP_peak_abund_NR),
      "TP_peak_relAbund" = max(peak_relAbund),
      "abund_gini" = pop_calcs(peak_abund, calc = "gini"),
      "tp_same" = sum(as.integer(strand == gene_ort)),
      "tp_same_CR" = sum(as.integer(strand == gene_ort)[patient %in% CR_pats]),
      "tp_same_NR" = sum(as.integer(strand == gene_ort)[patient %in% NR_pats]),
      "tp_oppo" = sum(as.integer(strand != gene_ort)),
      "tp_oppo_CR" = sum(as.integer(strand != gene_ort)[patient %in% CR_pats]),
      "tp_oppo_NR" = sum(as.integer(strand != gene_ort)[patient %in% NR_pats]))
  
  gene_stats <- data.frame(
      "loci" = get_loci_id(mod_refGenes),
      "gene_name" = as.character(mod_refGenes$name),
      "gene_ort" = as.character(strand(mod_refGenes)),
      "gene_width" = (width(mod_refGenes) + 10000) / 1000,
      "Within_Cluster" = mod_refGenes$Within_Cluster,
      "Cluster_target.min" = mod_refGenes$Cluster_target.min) %>%
    dplyr::full_join(., stats_tdn, by = c("loci", "gene_name", "gene_ort")) %>%
    dplyr::full_join(., stats_tp, by = c("loci", "gene_name", "gene_ort")) %>%
    dplyr::mutate(
      "On_Onco_List" = gene_name %in% oncoGenes,
      "Top_1pc_Abund" = gene_name %in% top_one_pc_sonicAbund,
      "Top_10pc_Abund" = gene_name %in% top_ten_pc_sonicAbund)
  
  gene_stats[is.na(gene_stats)] <- 0
  total_tdn_sites <- length(unique(paste(tdn_gr$patient, tdn_gr$posid)))
  total_tdn_sites_CR <- length(unique(paste(
    tdn_gr$patient, tdn_gr$posid)[tdn_gr$patient %in% CR_pats]))
  total_tdn_sites_NR <- length(unique(paste(
    tdn_gr$patient, tdn_gr$posid)[tdn_gr$patient %in% NR_pats]))
  total_pat_sites <- length(unique(paste(tp_gr$patient, tp_gr$posid)))
  total_pat_sites_CR <- length(unique(paste(
    tp_gr$patient, tp_gr$posid)[tp_gr$patient %in% CR_pats]))
  total_pat_sites_NR <- length(unique(paste(
    tp_gr$patient, tp_gr$posid)[tp_gr$patient %in% NR_pats]))
  
  gene_stats <- dplyr::filter(
      gene_stats, TDN_num_patients > 0 | TP_num_patients > 0
    ) %>%
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::mutate(
      "ort_fisher_test" = fisher.test(matrix(
          c(tdn_same, tdn_oppo, tp_same, tp_oppo), 
        nrow = 2, ncol = 2))$p.value,
      "ort_fisher_test_CR" = fisher.test(matrix(
        c(tdn_same_CR, tdn_oppo_CR, tp_same_CR, tp_oppo_CR), 
        nrow = 2, ncol = 2))$p.value,
      "ort_fisher_test_NR" = fisher.test(matrix(
        c(tdn_same_NR, tdn_oppo_NR, tp_same_NR, tp_oppo_NR), 
        nrow = 2, ncol = 2))$p.value) %>% 
    dplyr::ungroup() %>%
    dplyr::select(
      loci, gene_name, gene_ort, gene_width, TDN_num_patients, TDN_num_pats_CR, 
      TDN_num_pats_NR, TP_num_patients, TP_num_pats_CR, TP_num_pats_NR,
      TDN_num_sites, TDN_num_sites_CR, TDN_num_sites_NR, TP_num_sites, 
      TP_num_sites_CR, TP_num_sites_NR, TDN_peak_abund, TDN_peak_abund_CR, 
      TDN_peak_abund_NR, TP_peak_abund, TP_peak_abund_CR, TP_peak_abund_NR,
      TDN_peak_relAbund, TP_peak_relAbund, TDN_sum_abund, TDN_sum_abund_CR, 
      TDN_sum_abund_NR, TP_sum_abund, TP_sum_abund_CR, TP_sum_abund_NR, 
      long_count, max_time, max_span, abund_gini, ort_fisher_test, 
      ort_fisher_test_CR, ort_fisher_test_NR, On_Onco_List, Top_1pc_Abund, 
      Top_10pc_Abund, Within_Cluster, Cluster_target.min) %>%
    dplyr::mutate(
      TDN_freq = TDN_num_sites / (gene_width * total_tdn_sites),
      TDN_freq_CR = TDN_num_sites_CR / (gene_width * total_tdn_sites_CR),
      TDN_freq_NR = TDN_num_sites_NR / (gene_width * total_tdn_sites_NR),
      TP_freq = TP_num_sites/(gene_width * total_pat_sites),
      TP_freq_CR = TP_num_sites_CR/(gene_width * total_pat_sites_CR),
      TP_freq_NR = TP_num_sites_NR/(gene_width * total_pat_sites_NR),
      freq_diff = (TP_freq - TDN_freq),
      freq_diff_CR = (TP_freq_CR - TDN_freq_CR),
      freq_diff_NR = (TP_freq_NR - TDN_freq_NR),
      pct_chg = 100 * freq_diff / TDN_freq,
      pct_chg_CR = 100 * freq_diff_CR / TDN_freq_CR,
      pct_chg_NR = 100 * freq_diff_NR / TDN_freq_NR) %>%
    as.data.frame()
  
  saveRDS(
    gene_stats,
    file = file.path(outputDir, "cart19_gene_impact.rds"))
  
  write.csv(
    gene_stats,
    file = file.path(outputDir, "cart19_gene_impact.csv"),
    quote = TRUE,
    row.names = FALSE
  )

  # Assemble gene specific stats for each disease type independently -----------
  ## CLL ----
  CLL_stats_tdn <- as.data.frame(
      tdn_gr[subjectHits(tdn_in_gene)], row.names = NULL
    ) %>%
    dplyr::mutate(
      loci = get_loci_id(mod_refGenes[queryHits(tdn_in_gene)]),
      gene_name = as.character(mod_refGenes$name[queryHits(tdn_in_gene)]),
      gene_ort = as.character(strand(mod_refGenes[queryHits(tdn_in_gene)])),
      strand = as.character(strand)) %>%
    dplyr::select(
      loci, gene_name, gene_ort, posid, strand, estAbund, relAbund, patient
    ) %>%
    dplyr::filter(patient %in% CLL_pats) %>%                                    # Disease specificity
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::summarise(
      "TDN_num_patients" = n_distinct(patient),
      "TDN_num_pats_CR" = n_distinct(patient[patient %in% CR_pats]),
      "TDN_num_pats_NR" = n_distinct(patient[patient %in% NR_pats]),
      "TDN_num_sites" = n_distinct(posid),
      "TDN_num_sites_CR" = n_distinct(posid[patient %in% CR_pats]),
      "TDN_num_sites_NR" = n_distinct(posid[patient %in% NR_pats]),
      "TDN_sum_abund" = sum(estAbund),
      "TDN_sum_abund_CR" = sum(estAbund[patient %in% CR_pats]),
      "TDN_sum_abund_NR" = sum(estAbund[patient %in% NR_pats]),
      "TDN_peak_abund" = max(estAbund),
      "TDN_peak_abund_CR" = max(estAbund[patient %in% CR_pats]),
      "TDN_peak_abund_CR" = ifelse(TDN_peak_abund_CR < 0, 0, TDN_peak_abund_CR),
      "TDN_peak_abund_NR" = max(estAbund[patient %in% NR_pats]),
      "TDN_peak_abund_NR" = ifelse(TDN_peak_abund_NR < 0, 0, TDN_peak_abund_NR),
      "TDN_peak_relAbund" = max(relAbund),
      "tdn_same" = sum(as.integer(strand == gene_ort)),
      "tdn_same_CR" = sum(as.integer(strand == gene_ort)[patient %in% CR_pats]),
      "tdn_same_NR" = sum(as.integer(strand == gene_ort)[patient %in% NR_pats]),
      "tdn_oppo" = sum(as.integer(strand != gene_ort)),
      "tdn_oppo_CR" = sum(as.integer(strand != gene_ort)[patient %in% CR_pats]),
      "tdn_oppo_NR" = sum(as.integer(strand != gene_ort)[patient %in% NR_pats])
    )
  
  CLL_stats_tp <- as.data.frame(
      tp_gr[subjectHits(tp_in_gene)], row.names = NULL
    ) %>%
    dplyr::mutate(
      loci = get_loci_id(mod_refGenes[queryHits(tp_in_gene)]),
      gene_name = as.character(mod_refGenes$name[queryHits(tp_in_gene)]),
      gene_ort = as.character(strand(mod_refGenes[queryHits(tp_in_gene)])),
      strand = as.character(strand),
      timepoint = convert_time(as.character(timepoint))) %>%
    dplyr::select(
      loci, gene_name, gene_ort, posid, strand, 
      estAbund, relAbund, patient, timepoint
    ) %>%
    dplyr::filter(patient %in% CLL_pats) %>%                                    # Disease specificity
    dplyr::group_by(loci, gene_name, gene_ort, patient, posid, strand) %>%
    dplyr::summarise(
      "long_count" = n_distinct(timepoint),
      "first_time" = min(timepoint),
      "last_time" = max(timepoint),
      "sum_abund" = sum(estAbund),
      "peak_abund" = max(estAbund),
      "peak_relAbund" = max(relAbund)
    ) %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::summarise(
      "TP_num_patients" = n_distinct(patient),
      "TP_num_pats_CR" = n_distinct(patient[patient %in% CR_pats]),
      "TP_num_pats_NR" = n_distinct(patient[patient %in% NR_pats]),
      "long_count" = max(long_count),
      "max_time" = max(last_time),
      "max_span" = max(last_time - first_time),
      "TP_num_sites" = n_distinct(posid),
      "TP_num_sites_CR" = n_distinct(posid[patient %in% CR_pats]),
      "TP_num_sites_NR" = n_distinct(posid[patient %in% NR_pats]),
      "TP_sum_abund" = sum(sum_abund),
      "TP_sum_abund_CR" = sum(sum_abund[patient %in% CR_pats]),
      "TP_sum_abund_NR" = sum(sum_abund[patient %in% NR_pats]),
      "TP_peak_abund" = max(peak_abund),
      "TP_peak_abund_CR" = max(peak_abund[patient %in% CR_pats]),
      "TP_peak_abund_CR" = ifelse(TP_peak_abund_CR < 0, 0, TP_peak_abund_CR),
      "TP_peak_abund_NR" = max(peak_abund[patient %in% NR_pats]),
      "TP_peak_abund_NR" = ifelse(TP_peak_abund_NR < 0, 0, TP_peak_abund_NR),
      "TP_peak_relAbund" = max(peak_relAbund),
      "abund_gini" = pop_calcs(peak_abund, calc = "gini"),
      "tp_same" = sum(as.integer(strand == gene_ort)),
      "tp_same_CR" = sum(as.integer(strand == gene_ort)[patient %in% CR_pats]),
      "tp_same_NR" = sum(as.integer(strand == gene_ort)[patient %in% NR_pats]),
      "tp_oppo" = sum(as.integer(strand != gene_ort)),
      "tp_oppo_CR" = sum(as.integer(strand != gene_ort)[patient %in% CR_pats]),
      "tp_oppo_NR" = sum(as.integer(strand != gene_ort)[patient %in% NR_pats])
    )
  
  CLL_gene_stats <- data.frame(
    "loci" = get_loci_id(mod_refGenes),
    "gene_name" = as.character(mod_refGenes$name),
    "gene_ort" = as.character(strand(mod_refGenes)),
    "gene_width" = (width(mod_refGenes) + 10000) / 1000,
    "Within_Cluster" = mod_refGenes$Within_Cluster,
    "Cluster_target.min" = mod_refGenes$Cluster_target.min) %>%
    dplyr::full_join(CLL_stats_tdn, by = c("loci", "gene_name", "gene_ort")) %>%
    dplyr::full_join(CLL_stats_tp, by = c("loci", "gene_name", "gene_ort")) %>%
    dplyr::mutate(
      "On_Onco_List" = gene_name %in% oncoGenes,
      "Top_1pc_Abund" = gene_name %in% top_one_pc_sonicAbund,
      "Top_10pc_Abund" = gene_name %in% top_ten_pc_sonicAbund
    )
  
  CLL_gene_stats[is.na(CLL_gene_stats)] <- 0
  
  CLL_tdn_sites <- length(unique(
    paste(tdn_gr$patient, tdn_gr$posid)[tdn_gr$patient %in% CLL_pats]           # Disease specific
  ))
  
  CLL_tdn_sites_CR <- length(unique(paste(
    tdn_gr$patient, tdn_gr$posid)[
      tdn_gr$patient %in% intersect(CR_pats, CLL_pats)                          # Disease specific
    ]
  ))
  
  CLL_tdn_sites_NR <- length(unique(paste(
    tdn_gr$patient, tdn_gr$posid)[
      tdn_gr$patient %in% intersect(NR_pats, CLL_pats)                          # Disease specific
    ]
  ))
  
  CLL_pat_sites <- length(unique(paste(tp_gr$patient, tp_gr$posid)[
    tp_gr$patient %in% CLL_pats                                                 # Disease specific
  ]))
  
  CLL_pat_sites_CR <- length(unique(paste(
    tp_gr$patient, tp_gr$posid)[
      tp_gr$patient %in% intersect(CR_pats, CLL_pats)                           # Disease specific
    ]
  ))
  
  CLL_pat_sites_NR <- length(unique(paste(
    tp_gr$patient, tp_gr$posid)[
      tp_gr$patient %in% intersect(NR_pats, CLL_pats)                           # Disease specific
    ]
  ))
  
  
  CLL_gene_stats <- dplyr::filter(
      CLL_gene_stats, TDN_num_patients > 0 | TP_num_patients > 0
    ) %>%
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::mutate(
      "ort_fisher_test" = fisher.test(matrix(
        c(tdn_same, tdn_oppo, tp_same, tp_oppo), 
        nrow = 2, ncol = 2))$p.value,
      "ort_fisher_test_CR" = fisher.test(matrix(
        c(tdn_same_CR, tdn_oppo_CR, tp_same_CR, tp_oppo_CR), 
        nrow = 2, ncol = 2))$p.value,
      "ort_fisher_test_NR" = fisher.test(matrix(
        c(tdn_same_NR, tdn_oppo_NR, tp_same_NR, tp_oppo_NR), 
        nrow = 2, ncol = 2))$p.value) %>% 
    dplyr::ungroup() %>%
    dplyr::select(
      loci, gene_name, gene_ort, gene_width, TDN_num_patients, TDN_num_pats_CR, 
      TDN_num_pats_NR, TP_num_patients, TP_num_pats_CR, TP_num_pats_NR,
      TDN_num_sites, TDN_num_sites_CR, TDN_num_sites_NR, TP_num_sites, 
      TP_num_sites_CR, TP_num_sites_NR, TDN_peak_abund, TDN_peak_abund_CR, 
      TDN_peak_abund_NR, TP_peak_abund, TP_peak_abund_CR, TP_peak_abund_NR,
      TDN_peak_relAbund, TP_peak_relAbund, TDN_sum_abund, TDN_sum_abund_CR, 
      TDN_sum_abund_NR, TP_sum_abund, TP_sum_abund_CR, TP_sum_abund_NR, 
      long_count, max_time, max_span, abund_gini, ort_fisher_test, 
      ort_fisher_test_CR, ort_fisher_test_NR, On_Onco_List, Top_1pc_Abund, 
      Top_10pc_Abund, Within_Cluster, Cluster_target.min) %>%
    dplyr::mutate(
      TDN_freq = TDN_num_sites / (gene_width * CLL_tdn_sites),
      TDN_freq_CR = TDN_num_sites_CR / (gene_width * CLL_tdn_sites_CR),
      TDN_freq_NR = TDN_num_sites_NR / (gene_width * CLL_tdn_sites_NR),
      TP_freq = TP_num_sites/(gene_width * CLL_pat_sites),
      TP_freq_CR = TP_num_sites_CR/(gene_width * CLL_pat_sites_CR),
      TP_freq_NR = TP_num_sites_NR/(gene_width * CLL_pat_sites_NR),
      freq_diff = (TP_freq - TDN_freq),
      freq_diff_CR = (TP_freq_CR - TDN_freq_CR),
      freq_diff_NR = (TP_freq_NR - TDN_freq_NR),
      pct_chg = 100 * freq_diff / TDN_freq,
      pct_chg_CR = 100 * freq_diff_CR / TDN_freq_CR,
      pct_chg_NR = 100 * freq_diff_NR / TDN_freq_NR) %>%
    as.data.frame()
  
  saveRDS(
    CLL_gene_stats,
    file = file.path(outputDir, "cart19_cll_gene_impact.rds"))
  
  write.csv(
    CLL_gene_stats,
    file = file.path(outputDir, "cart19_cll_gene_impact.csv"),
    quote = TRUE,
    row.names = FALSE
  )
  
  ## Adult ALL ----
  ALL_stats_tdn <- as.data.frame(
      tdn_gr[subjectHits(tdn_in_gene)], row.names = NULL
    ) %>%
    dplyr::mutate(
      loci = get_loci_id(mod_refGenes[queryHits(tdn_in_gene)]),
      gene_name = as.character(mod_refGenes$name[queryHits(tdn_in_gene)]),
      gene_ort = as.character(strand(mod_refGenes[queryHits(tdn_in_gene)])),
      strand = as.character(strand)) %>%
    dplyr::select(
      loci, gene_name, gene_ort, posid, strand, estAbund, relAbund, patient
    ) %>%
    dplyr::filter(patient %in% ALL_pats) %>%                                    # Disease specificity
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::summarise(
      "TDN_num_patients" = n_distinct(patient),
      "TDN_num_pats_CR" = n_distinct(patient[patient %in% CR_pats]),
      "TDN_num_pats_NR" = n_distinct(patient[patient %in% NR_pats]),
      "TDN_num_sites" = n_distinct(posid),
      "TDN_num_sites_CR" = n_distinct(posid[patient %in% CR_pats]),
      "TDN_num_sites_NR" = n_distinct(posid[patient %in% NR_pats]),
      "TDN_sum_abund" = sum(estAbund),
      "TDN_sum_abund_CR" = sum(estAbund[patient %in% CR_pats]),
      "TDN_sum_abund_NR" = sum(estAbund[patient %in% NR_pats]),
      "TDN_peak_abund" = max(estAbund),
      "TDN_peak_abund_CR" = max(estAbund[patient %in% CR_pats]),
      "TDN_peak_abund_CR" = ifelse(TDN_peak_abund_CR < 0, 0, TDN_peak_abund_CR),
      "TDN_peak_abund_NR" = max(estAbund[patient %in% NR_pats]),
      "TDN_peak_abund_NR" = ifelse(TDN_peak_abund_NR < 0, 0, TDN_peak_abund_NR),
      "TDN_peak_relAbund" = max(relAbund),
      "tdn_same" = sum(as.integer(strand == gene_ort)),
      "tdn_same_CR" = sum(as.integer(strand == gene_ort)[patient %in% CR_pats]),
      "tdn_same_NR" = sum(as.integer(strand == gene_ort)[patient %in% NR_pats]),
      "tdn_oppo" = sum(as.integer(strand != gene_ort)),
      "tdn_oppo_CR" = sum(as.integer(strand != gene_ort)[patient %in% CR_pats]),
      "tdn_oppo_NR" = sum(as.integer(strand != gene_ort)[patient %in% NR_pats])
    )
  
  ALL_stats_tp <- as.data.frame(
      tp_gr[subjectHits(tp_in_gene)], row.names = NULL
    ) %>%
    dplyr::mutate(
      loci = get_loci_id(mod_refGenes[queryHits(tp_in_gene)]),
      gene_name = as.character(mod_refGenes$name[queryHits(tp_in_gene)]),
      gene_ort = as.character(strand(mod_refGenes[queryHits(tp_in_gene)])),
      strand = as.character(strand),
      timepoint = convert_time(as.character(timepoint))) %>%
    dplyr::select(
      loci, gene_name, gene_ort, posid, strand, 
      estAbund, relAbund, patient, timepoint
    ) %>%
    dplyr::filter(patient %in% ALL_pats) %>%                                    # Disease specificity
    dplyr::group_by(loci, gene_name, gene_ort, patient, posid, strand) %>%
    dplyr::summarise(
      "long_count" = n_distinct(timepoint),
      "first_time" = min(timepoint),
      "last_time" = max(timepoint),
      "sum_abund" = sum(estAbund),
      "peak_abund" = max(estAbund),
      "peak_relAbund" = max(relAbund)
    ) %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::summarise(
      "TP_num_patients" = n_distinct(patient),
      "TP_num_pats_CR" = n_distinct(patient[patient %in% CR_pats]),
      "TP_num_pats_NR" = n_distinct(patient[patient %in% NR_pats]),
      "long_count" = max(long_count),
      "max_time" = max(last_time),
      "max_span" = max(last_time - first_time),
      "TP_num_sites" = n_distinct(posid),
      "TP_num_sites_CR" = n_distinct(posid[patient %in% CR_pats]),
      "TP_num_sites_NR" = n_distinct(posid[patient %in% NR_pats]),
      "TP_sum_abund" = sum(sum_abund),
      "TP_sum_abund_CR" = sum(sum_abund[patient %in% CR_pats]),
      "TP_sum_abund_NR" = sum(sum_abund[patient %in% NR_pats]),
      "TP_peak_abund" = max(peak_abund),
      "TP_peak_abund_CR" = max(peak_abund[patient %in% CR_pats]),
      "TP_peak_abund_CR" = ifelse(TP_peak_abund_CR < 0, 0, TP_peak_abund_CR),
      "TP_peak_abund_NR" = max(peak_abund[patient %in% NR_pats]),
      "TP_peak_abund_NR" = ifelse(TP_peak_abund_NR < 0, 0, TP_peak_abund_NR),
      "TP_peak_relAbund" = max(peak_relAbund),
      "abund_gini" = pop_calcs(peak_abund, calc = "gini"),
      "tp_same" = sum(as.integer(strand == gene_ort)),
      "tp_same_CR" = sum(as.integer(strand == gene_ort)[patient %in% CR_pats]),
      "tp_same_NR" = sum(as.integer(strand == gene_ort)[patient %in% NR_pats]),
      "tp_oppo" = sum(as.integer(strand != gene_ort)),
      "tp_oppo_CR" = sum(as.integer(strand != gene_ort)[patient %in% CR_pats]),
      "tp_oppo_NR" = sum(as.integer(strand != gene_ort)[patient %in% NR_pats])
    )
  
  ALL_gene_stats <- data.frame(
    "loci" = get_loci_id(mod_refGenes),
    "gene_name" = as.character(mod_refGenes$name),
    "gene_ort" = as.character(strand(mod_refGenes)),
    "gene_width" = (width(mod_refGenes) + 10000) / 1000,
    "Within_Cluster" = mod_refGenes$Within_Cluster,
    "Cluster_target.min" = mod_refGenes$Cluster_target.min) %>%
    dplyr::full_join(ALL_stats_tdn, by = c("loci", "gene_name", "gene_ort")) %>%
    dplyr::full_join(ALL_stats_tp, by = c("loci", "gene_name", "gene_ort")) %>%
    dplyr::mutate(
      "On_Onco_List" = gene_name %in% oncoGenes,
      "Top_1pc_Abund" = gene_name %in% top_one_pc_sonicAbund,
      "Top_10pc_Abund" = gene_name %in% top_ten_pc_sonicAbund
  )
  
  ALL_gene_stats[is.na(ALL_gene_stats)] <- 0
  
  ALL_tdn_sites <- length(unique(
    paste(tdn_gr$patient, tdn_gr$posid)[tdn_gr$patient %in% ALL_pats]           # Disease specific
  ))
  
  ALL_tdn_sites_CR <- length(unique(paste(
    tdn_gr$patient, tdn_gr$posid)[
      tdn_gr$patient %in% intersect(CR_pats, ALL_pats)                          # Disease specific
    ]
  ))
  
  ALL_tdn_sites_NR <- length(unique(paste(
    tdn_gr$patient, tdn_gr$posid)[
      tdn_gr$patient %in% intersect(NR_pats, ALL_pats)                          # Disease specific
    ]
  ))
  
  ALL_pat_sites <- length(unique(paste(tp_gr$patient, tp_gr$posid)[
    tp_gr$patient %in% ALL_pats                                                 # Disease specific
  ]))
  
  ALL_pat_sites_CR <- length(unique(paste(
    tp_gr$patient, tp_gr$posid)[
      tp_gr$patient %in% intersect(CR_pats, ALL_pats)                           # Disease specific
    ]
  ))
  
  ALL_pat_sites_NR <- length(unique(paste(
    tp_gr$patient, tp_gr$posid)[
      tp_gr$patient %in% intersect(NR_pats, ALL_pats)                           # Disease specific
    ]
  ))
  
  ALL_gene_stats <- dplyr::filter(
      ALL_gene_stats, TDN_num_patients > 0 | TP_num_patients > 0
    ) %>%
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::mutate(
      "ort_fisher_test" = fisher.test(matrix(
        c(tdn_same, tdn_oppo, tp_same, tp_oppo), 
        nrow = 2, ncol = 2))$p.value,
      "ort_fisher_test_CR" = fisher.test(matrix(
        c(tdn_same_CR, tdn_oppo_CR, tp_same_CR, tp_oppo_CR), 
        nrow = 2, ncol = 2))$p.value,
      "ort_fisher_test_NR" = fisher.test(matrix(
        c(tdn_same_NR, tdn_oppo_NR, tp_same_NR, tp_oppo_NR), 
        nrow = 2, ncol = 2))$p.value) %>% 
    dplyr::ungroup() %>%
    dplyr::select(
      loci, gene_name, gene_ort, gene_width, TDN_num_patients, TDN_num_pats_CR, 
      TDN_num_pats_NR, TP_num_patients, TP_num_pats_CR, TP_num_pats_NR,
      TDN_num_sites, TDN_num_sites_CR, TDN_num_sites_NR, TP_num_sites, 
      TP_num_sites_CR, TP_num_sites_NR, TDN_peak_abund, TDN_peak_abund_CR, 
      TDN_peak_abund_NR, TP_peak_abund, TP_peak_abund_CR, TP_peak_abund_NR,
      TDN_peak_relAbund, TP_peak_relAbund, TDN_sum_abund, TDN_sum_abund_CR, 
      TDN_sum_abund_NR, TP_sum_abund, TP_sum_abund_CR, TP_sum_abund_NR, 
      long_count, max_time, max_span, abund_gini, ort_fisher_test, 
      ort_fisher_test_CR, ort_fisher_test_NR, On_Onco_List, Top_1pc_Abund, 
      Top_10pc_Abund, Within_Cluster, Cluster_target.min) %>%
    dplyr::mutate(
      TDN_freq = TDN_num_sites / (gene_width * ALL_tdn_sites),
      TDN_freq_CR = TDN_num_sites_CR / (gene_width * ALL_tdn_sites_CR),
      TDN_freq_NR = TDN_num_sites_NR / (gene_width * ALL_tdn_sites_NR),
      TP_freq = TP_num_sites/(gene_width * ALL_pat_sites),
      TP_freq_CR = TP_num_sites_CR/(gene_width * ALL_pat_sites_CR),
      TP_freq_NR = TP_num_sites_NR/(gene_width * ALL_pat_sites_NR),
      freq_diff = (TP_freq - TDN_freq),
      freq_diff_CR = (TP_freq_CR - TDN_freq_CR),
      freq_diff_NR = (TP_freq_NR - TDN_freq_NR),
      pct_chg = 100 * freq_diff / TDN_freq,
      pct_chg_CR = 100 * freq_diff_CR / TDN_freq_CR,
      pct_chg_NR = 100 * freq_diff_NR / TDN_freq_NR) %>%
    as.data.frame()
  
  saveRDS(
    ALL_gene_stats,
    file = file.path(outputDir, "cart19_all_gene_impact.rds"))
  
  write.csv(
    ALL_gene_stats,
    file = file.path(outputDir, "cart19_all_gene_impact.csv"),
    quote = TRUE,
    row.names = FALSE
  )
  
  ## Responders only ----
  CR_stats_tdn <- as.data.frame(
      tdn_gr[subjectHits(tdn_in_gene)], row.names = NULL
    ) %>%
    dplyr::mutate(
      loci = get_loci_id(mod_refGenes[queryHits(tdn_in_gene)]),
      gene_name = as.character(mod_refGenes$name[queryHits(tdn_in_gene)]),
      gene_ort = as.character(strand(mod_refGenes[queryHits(tdn_in_gene)])),
      strand = as.character(strand)) %>%
    dplyr::select(
      loci, gene_name, gene_ort, posid, strand, estAbund, relAbund, patient
    ) %>%
    dplyr::filter(patient %in% Std_CR_pats) %>%                                 # patient specificity
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::summarise(
      "TDN_num_patients" = n_distinct(patient),
      "TDN_num_sites" = n_distinct(posid),
      "TDN_sum_abund" = sum(estAbund),
      "TDN_peak_abund" = max(estAbund),
      "TDN_peak_relAbund" = max(relAbund),
      "tdn_same" = sum(as.integer(strand == gene_ort)),
      "tdn_oppo" = sum(as.integer(strand != gene_ort))
    )
  
  CR_stats_tp <- as.data.frame(
      tp_gr[subjectHits(tp_in_gene)], row.names = NULL
    ) %>%
    dplyr::mutate(
      loci = get_loci_id(mod_refGenes[queryHits(tp_in_gene)]),
      gene_name = as.character(mod_refGenes$name[queryHits(tp_in_gene)]),
      gene_ort = as.character(strand(mod_refGenes[queryHits(tp_in_gene)])),
      strand = as.character(strand),
      timepoint = convert_time(as.character(timepoint))) %>%
    dplyr::select(
      loci, gene_name, gene_ort, posid, strand, 
      estAbund, relAbund, patient, timepoint
    ) %>%
    dplyr::filter(patient %in% Std_CR_pats) %>%                                 # patient specificity
    dplyr::group_by(loci, gene_name, gene_ort, patient, posid, strand) %>%
    dplyr::summarise(
      "long_count" = n_distinct(timepoint),
      "first_time" = min(timepoint),
      "last_time" = max(timepoint),
      "sum_abund" = sum(estAbund),
      "peak_abund" = max(estAbund),
      "peak_relAbund" = max(relAbund)
    ) %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(loci, gene_name, gene_ort) %>%
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
      "tp_same" = sum(as.integer(strand == gene_ort)),
      "tp_oppo" = sum(as.integer(strand != gene_ort))
    )
  
  CR_gene_stats <- data.frame(
    "loci" = get_loci_id(mod_refGenes),
    "gene_name" = as.character(mod_refGenes$name),
    "gene_ort" = as.character(strand(mod_refGenes)),
    "gene_width" = (width(mod_refGenes) + 10000) / 1000,
    "Within_Cluster" = mod_refGenes$Within_Cluster,
    "Cluster_target.min" = mod_refGenes$Cluster_target.min) %>%
    dplyr::full_join(CLL_stats_tdn, by = c("loci", "gene_name", "gene_ort")) %>%
    dplyr::full_join(CLL_stats_tp, by = c("loci", "gene_name", "gene_ort")) %>%
    dplyr::mutate(
      "On_Onco_List" = gene_name %in% oncoGenes,
      "Top_1pc_Abund" = gene_name %in% top_one_pc_sonicAbund,
      "Top_10pc_Abund" = gene_name %in% top_ten_pc_sonicAbund
    )
  
  CR_gene_stats[is.na(CR_gene_stats)] <- 0
  
  CR_tdn_sites <- length(unique(
    paste(tdn_gr$patient, tdn_gr$posid)[tdn_gr$patient %in% Std_CR_pats]        # patient specific
  ))
  
  CR_pat_sites <- length(unique(paste(tp_gr$patient, tp_gr$posid)[
    tp_gr$patient %in% Std_CR_pats                                              # patient specific
  ]))
  

  CR_gene_stats <- dplyr::filter(
      CR_gene_stats, TDN_num_patients > 0 | TP_num_patients > 0
    ) %>%
    dplyr::group_by(loci, gene_name, gene_ort) %>%
    dplyr::mutate(
      "ort_fisher_test" = fisher.test(matrix(
        c(tdn_same, tdn_oppo, tp_same, tp_oppo), 
        nrow = 2, ncol = 2))$p.value
    ) %>% 
    dplyr::ungroup() %>%
    dplyr::select(
      loci, gene_name, gene_ort, gene_width, TDN_num_patients, 
      TP_num_patients, TDN_num_sites, TP_num_sites, TDN_peak_abund, 
      TP_peak_abund, TDN_peak_relAbund, TP_peak_relAbund, TDN_sum_abund, 
      TP_sum_abund, long_count, max_time, max_span, abund_gini, ort_fisher_test, 
      On_Onco_List, Top_1pc_Abund, Top_10pc_Abund, Within_Cluster, 
      Cluster_target.min
    ) %>%
    dplyr::mutate(
      TDN_freq = TDN_num_sites / (gene_width * CR_tdn_sites),
      TP_freq = TP_num_sites/(gene_width * CR_pat_sites),
      freq_diff = (TP_freq - TDN_freq),
      pct_chg = 100 * freq_diff / TDN_freq
    ) %>%
    as.data.frame()
  
  saveRDS(
    CR_gene_stats,
    file = file.path(outputDir, "cart19_cr_gene_impact.rds"))
  
  write.csv(
    CR_gene_stats,
    file = file.path(outputDir, "cart19_cr_gene_impact.csv"),
    quote = TRUE,
    row.names = FALSE
  )

}else{
  
  gene_stats <- readRDS(file.path(outputDir, "cart19_gene_impact.rds"))
  CLL_gene_stats <- readRDS(file.path(outputDir, "cart19_cll_gene_impact.rds"))
  ALL_gene_stats <- readRDS(file.path(outputDir, "cart19_all_gene_impact.rds"))
  CR_gene_stats <- readRDS(file.path(outputDir, "cart19_cr_gene_impact.rds"))
  
}