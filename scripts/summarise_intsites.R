#### This script is for summarizing the data found in the CART19 unique 
#### integration sites. This script is dependent on the 
#### cart19_intsite_analysis.Rmd setup environment
output_files <- c(
  "cart19_specimen_summary.csv", "cart19_annotated_uniq_sites.rds", 
  "cart19_timepoint_summary.csv", "cart19_patient_summary.csv",
  "cart19_celltype_summary.csv", "cart19_summaries.rds",
  "top_one_percent_of_gene_ids_by_sonicAbund.txt",
  "top_ten_percent_of_gene_ids_by_relAbund.txt"
)

if(
  !all(sapply(output_files, function(x) file.exists(file.path(outputDir, x))))
){
  
  packs <- c(
    "hiAnnotator", "intSiteRetriever", "GCcontent", "BSgenome")
  null <- suppressMessages(sapply(packs, require, character.only = TRUE)) 
  
  # Load genomic and epigenetic features -----------------------------------------
  genome_sequence <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  genome_sequence@user_seqnames <- genome_sequence@user_seqnames[
    genome_sequence@user_seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))]
  genome_sequence@seqinfo <- genome_sequence@seqinfo[
    paste0("chr", c(1:22, "X", "Y", "M"))]
  
  #CpG_islands <- getCpG_islands(genomicFreeze)
  CpG_data <- cpg <- getUCSCtable(
      "cpgIslandExt", "CpG Islands", freeze = "hg38"
    ) %>%
    dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))
  
  CpG_islands <- GenomicRanges::GRanges(
    seqnames = CpG_data$chrom,
    ranges = IRanges::IRanges(
      start = CpG_data$chromStart, end = CpG_data$chromEnd
    ),
    strand = "*",
    seqinfo = GenomeInfoDb::seqinfo(genome_sequence)
  )
  
  mcols(CpG_islands) <- CpG_data
  
  #DNaseI <- suppressWarnings(getDNaseI(genomicFreeze))
  DNaseI_data <- getUCSCtable(
      "wgEncodeRegDnaseClustered", "DNase Clusters", freeze = "hg38"
    ) %>%
    dplyr::filter(chrom %in% paste0("chr", c(1:22, "X", "Y", "M")))
  
  DNaseI <- GenomicRanges::GRanges(
    seqnames = DNaseI_data$chrom,
    ranges = IRanges::IRanges(
      start = DNaseI_data$chromStart, end = DNaseI_data$chromEnd
    ),
    strand = "*",
    seqinfo = GenomeInfoDb::seqinfo(genome_sequence)
  )
  
  mcols(DNaseI) <- DNaseI_data
  
  ## windows
  window_size_refSeq <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
  window_size_CpG_counts <- c("1k"=1e3, "10k"=1e4)
  window_size_CpG_density <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
  window_size_GC <- c("100"=100, "1k"=1000, "10k"=1e4, "100k"=1e5, "1M"=1e6)
  window_size_DNaseI <- c("1k"=1e3, "10k"=1e4, "100k"=1e5, "1M"=1e6)
  window_size_epi <- c("10k"=1e4)
  
  ## Load epigenetic features
  epi_names <- str_split(epigenetic_features_files, "/")
  epi_names <- sapply(epi_names, "[[", 6)
  epi_names <- str_remove(epi_names, ".RData") %>% str_replace_all("-", "_")
  
  epi_env <- new.env()
  epi_env$epi_features <- list()
  
  null <- lapply(epigenetic_features_files, function(path){
    load(path)
    epi_env$epi_features <- c(epi_env$epi_features, list(epigenData))
  })
  
  epi_features <- epi_env$epi_features
  names(epi_features) <- epi_names
  
  # Supporting functions ---------------------------------------------------------
  jackIID <- function(ids, jrep = NULL, nrep = 10L){
    if(is.null(jrep)) jrep <- sample(rep(1:nrep,length=length(ids)))
    est0 <- vegan::estimateR(table(ids))
    jackrep <- jrep
    urepl <- unique(jrep)
    jackmat <- sapply(
      urepl, function(x) vegan::estimateR(table(ids[jackrep!=x])))
    pseudo <- length(urepl)*est0 - (length(urepl)-1)*jackmat
    rowMeans(pseudo)
  }
  
  #' jackknife biased or unbiased Chao estimator
  #'
  #' @param replicatedSites df with column posid
  #' @return number population size estimate
  calculateChao <- function(ids, biased=TRUE){
    if ( !biased ) { #regular Chao
      cluster.tab <- table(ids)
      return(round(estimateR(cluster.tab)["S.chao1"]))
    }
    round(jackIID(ids)["S.chao1"])
  }
  
  # Generate summary stats for each sample ---------------------------------------
  ## Append genomic features
  cond_uniq_sites <- cond_uniq_sites[
    seqnames(cond_uniq_sites) %in% paste0("chr", c(1:22, "X", "Y", "M"))]
  seqlevels(cond_uniq_sites) <- paste0("chr", c(1:22, "X", "Y", "M"))
  seqinfo(cond_uniq_sites) <- seqinfo(genome_sequence)
  
  cond_uniq_list <- split(
    cond_uniq_sites, 
    ceiling(seq_along(cond_uniq_sites) / 
              (length(cond_uniq_sites)/length(unique(cond_uniq_sites$specimen))))
  )
  
  buster <- parallel::makeCluster(numCores) 
  # Memory requirement ~2.5 GB per core
  
  parallel::clusterExport(
    buster, 
    varlist = c(
      "from_counts_to_density", "refGenes", "CpG_islands", "DNaseI", 
      "window_size_refSeq", "window_size_GC", "window_size_CpG_counts", 
      "window_size_CpG_density", "window_size_DNaseI", "genome_sequence"
    ), 
    envir = environment()
  )
  
  cond_uniq_annot_gen <- unname(unlist(GRangesList(parallel::parLapply(
    buster, 
    cond_uniq_list,
    function(gr){
  
      library(magrittr)
      library(GenomicRanges)
  
      gr %>%
        hiAnnotator::getFeatureCounts(
          refGenes, "refSeq_counts", width = window_size_refSeq) %>%
        GCcontent::getGCpercentage(
          "GC", window_size_GC, genome_sequence) %>%
        hiAnnotator::getFeatureCounts(
          CpG_islands, "CpG_counts", width = window_size_CpG_counts) %>%
        hiAnnotator::getFeatureCounts(
          CpG_islands, "CpG_density", width = window_size_CpG_density) %>%
        from_counts_to_density(
          "CpG_density", window_size_CpG_density) %>%
        hiAnnotator::getFeatureCounts(
          DNaseI, "DNaseI_count", width = window_size_DNaseI)
  
    }
  ))))
  
  parallel::stopCluster(buster)
  
  
  ## Append epigenetic features
  
  buster <- parallel::makeCluster(numCores) 
  # Memory requirement ~1.5 GB per core, max ~ 6GB
   
  parallel::clusterExport(
    buster, 
    varlist = c("cond_uniq_sites", "window_size_epi"), 
    envir = environment()
  )
  
  epi_annots <- dplyr::bind_cols(parallel::clusterMap(
    buster,
    function(epi, name){
      
      library(GenomicRanges)
      
      annot_gr <- hiAnnotator::getFeatureCounts(
        cond_uniq_sites, epi, name, width = window_size_epi
      )
      
      feat_cols <- names(mcols(annot_gr))[
        !names(mcols(annot_gr)) %in% names(mcols(cond_uniq_sites))
      ]
      
      as.data.frame(mcols(annot_gr)[, feat_cols, drop = FALSE])
      
    },
    epi = epi_features,
    name = names(epi_features)
  ))
  
  parallel::stopCluster(buster)
  
  
  ## Combine and summarise features
  cond_uniq_annot <- cond_uniq_annot_gen
  mcols(cond_uniq_annot) <- bind_cols(
    as.data.frame(mcols(cond_uniq_annot)), 
    epi_annots
  )
  
  if( !all(cond_uniq_annot$posid == cond_uniq_sites$posid) ){
    stop("Indexing error occured during parallel processing.\n")
  }
  
  gen_epi_stats <- as.data.frame(cond_uniq_annot, row.names = NULL) %>%
    dplyr::select(
      -seqnames, -start, -end, -width, -refgenome, -relRank, -nearest_geneDist, 
      -nearest_gene, -nearest_geneOrt, -patient, -celltype, -timepoint, 
      -gene_id_wo_annot, -gene_id, -strand, -in_gene, -in_geneOrt)
  
  cond_uniq_df <- as.data.frame(cond_uniq_sites, row.names = NULL) %>% 
    dplyr::mutate(
      within_gene = in_gene != FALSE,
      same_ort = same_ort(strand, in_geneOrt),
      same_ort = ifelse(is.na(in_geneOrt), NA, same_ort))
  
  stats <- cond_uniq_df %>%
    dplyr::group_by(specimen) %>%
    dplyr::summarise(
      "numUniqSites" = n(), 
      "ShannonIndex" = pop_calcs(estAbund, calc = "shannon"),
      "GiniIndex" = pop_calcs(estAbund, calc = "gini"),
      "Chao1" = calculateChao(as.character(Rle(posid, estAbund))),
      "UC50" = pop_calcs(estAbund, calc = "uc50"),
      "pctTxnUnit" = 100 * sum(within_gene, na.rm = TRUE)/n(),
      "pctSameOrt" = 100 * sum(same_ort, na.rm = TRUE) / 
        sum(within_gene, na.rm = TRUE),
      "pctNearTxnUn" = 100 * sum(
        in_gene == FALSE & abs(nearest_geneDist) <= 5000, na.rm = TRUE) / 
        sum(in_gene == FALSE, na.rm = TRUE),
      "pctInOnco" = 100 * sum(in_gene %in% oncoGenes, na.rm = TRUE) / n())
  
  gen_epi_summary <- gen_epi_stats %>%
    dplyr::select(-posid) %>%
    dplyr::group_by(specimen) %>%
    dplyr::summarise_all(mean, na.rm = TRUE)
  
  ## Specimen Summary ------------------------------------------------------------
  
  specimen_data <- patient_data %>%
    dplyr::full_join(specimen_data, by = "patient") %>%
    dplyr::rename(specimen = specimenaccnum) %>%  
    dplyr::inner_join(stats, by = "specimen") %>%
    dplyr::left_join(gen_epi_summary, by = "specimen")
  
  specimen_data$clustersRepresented <- sapply(
    specimen_data$specimen, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$specimen == x]
      hits <- findOverlaps(sites, red_clusters)
      length(unique(subjectHits(hits)))
    })

  specimen_data$numSitesInClusters <- sapply(
    specimen_data$specimen, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$specimen == x]
      hits <- findOverlaps(sites, red_clusters)
      length(unique(queryHits(hits)))
    })

  specimen_data$abundInClusters <- sapply(
    specimen_data$specimen, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$specimen == x]
      hits <- findOverlaps(sites, red_clusters)
      sum(sites[queryHits(hits)]$estAbund)
    })

  write.csv(
    specimen_data,
    file = file.path(outputDir, "cart19_specimen_summary.csv"),
    quote = TRUE,
    row.names = FALSE)
  
  saveRDS(
    cond_uniq_annot,
    file = file.path(outputDir, "cart19_annotated_uniq_sites.rds"))
  

  # Timepoint Summary ----------------------------------------------------------
  cond_uniq_sites$timepoint <- factor(
    cond_uniq_sites$timepoint, levels = timepointLevels)
  
  gen_epi_timepoint_summary <- gen_epi_stats %>%
    dplyr::select(-posid) %>%
    dplyr::left_join(
      dplyr::select(specimen_data, specimen, timepoint), 
      by = "specimen"
    ) %>%
    dplyr::select(-specimen) %>%
    dplyr::group_by(timepoint) %>%
    dplyr::summarise_all(mean, na.rm = TRUE)
  
  timepoint_summary <- data.frame(
    "timepoint" = names(split(cond_uniq_sites, cond_uniq_sites$timepoint)),
    "num_UniqueSites" = sapply(
      split(cond_uniq_sites, cond_uniq_sites$timepoint), 
      function(x){length(unique(x$posid))}),
    "num_UniqueSites_CR" = sapply(
      split(cond_uniq_sites, cond_uniq_sites$timepoint), 
      function(x){length(unique(x$posid[x$patient %in% CR_pats]))}),
    "num_UniqueSites_NR" = sapply(
      split(cond_uniq_sites, cond_uniq_sites$timepoint), 
      function(x){length(unique(x$posid[x$patient %in% NR_pats]))}),
    "num_Patients" = sapply(
      split(cond_uniq_sites, cond_uniq_sites$timepoint), 
      function(x) length(unique(x$patient))),
    "num_Patients_CR" = sapply(
      split(cond_uniq_sites, cond_uniq_sites$timepoint), 
      function(x) length(unique(x$patient[x$patients %in% CR_pats]))),
    "num_Patients_NR" = sapply(
      split(cond_uniq_sites, cond_uniq_sites$timepoint), 
      function(x) length(unique(x$patient[x$patients %in% NR_pats])))) %>%
    left_join(gen_epi_timepoint_summary, by = "timepoint")
  
  timepoint_summary <- filter(timepoint_summary, num_Patients > 0)
  
  timepoint_summary$clustersRepresented <- sapply(
    timepoint_summary$timepoint, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$timepoint == x]
      hits <- findOverlaps(sites, red_clusters)
      length(unique(subjectHits(hits)))
    })
  
  timepoint_summary$numSitesInClusters <- sapply(
    timepoint_summary$timepoint, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$timepoint == x]
      hits <- findOverlaps(sites, red_clusters)
      length(unique(queryHits(hits)))
    })
  
  timepoint_summary$abundInClusters <- sapply(
    timepoint_summary$timepoint, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$timepoint == x]
      hits <- findOverlaps(sites, red_clusters)
      sum(sites[queryHits(hits)]$estAbund)
    })
  
  write.csv(
    timepoint_summary,
    file = file.path(outputDir, "cart19_timepoint_summary.csv"),
    quote = TRUE,
    row.names = FALSE)
  
  
  # Patient Summary --------------------------------------------------------------
  patient_summary <- group_by(specimen_data, patient) %>% 
    dplyr::arrange(timepoint) %>%
    dplyr::summarise(
      numTimepoints = n_distinct(timepoint), 
      numCellTypes = n_distinct(celltype), 
      numSpecimens = n_distinct(specimen), 
      tdnSample = any(timepoint == "d0"), 
      d28Sample = any(timepoint == "d28"), 
      lastTimepoint = last(timepoint)) %>% 
    as.data.frame()
  
  patient_summary$numUniqSites <- sapply(
    patient_summary$patient, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$patient == x]
      length(unique(sites$posid))
    })
  
  patient_summary$clustersRepresented <- sapply(
    patient_summary$patient, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$patient == x]
      hits <- findOverlaps(sites, red_clusters)
      length(unique(subjectHits(hits)))
    })
  
  patient_summary$numSitesInClusters <- sapply(
    patient_summary$patient, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$patient == x]
      hits <- findOverlaps(sites, red_clusters)
      length(unique(queryHits(hits)))
    })
  
  patient_summary$abundInClusters <- sapply(
    patient_summary$patient, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$patient == x]
      hits <- findOverlaps(sites, red_clusters)
      sum(sites[queryHits(hits)]$estAbund)
    })
  
  gen_epi_patient_summary <- gen_epi_stats %>%
    dplyr::select(-posid) %>%
    dplyr::left_join(
      dplyr::select(specimen_data, specimen, patient), 
      by = "specimen"
    ) %>% 
    dplyr::select(-specimen) %>%
    dplyr::group_by(patient) %>%
    dplyr::summarise_all(mean, na.rm = TRUE)
  
  patient_summary <- dplyr::left_join(
    patient_summary, gen_epi_patient_summary, by = "patient")
  
  write.csv(
    patient_summary,
    file = file.path(outputDir, "cart19_patient_summary.csv"),
    quote = TRUE,
    row.names = FALSE
  )

  # Celltype Summary -------------------------------------------------------------
  cond_uniq_sites$celltype <- factor(
    cond_uniq_sites$celltype, levels = celltypeLevels)
  
  celltype_summary <- data.frame(
    "celltype" = names(split(cond_uniq_sites, cond_uniq_sites$celltype)),
    "num_UniqueSites" = sapply(
      split(cond_uniq_sites, cond_uniq_sites$celltype), 
      function(x) length(unique(x$posid)) ),
    "num_Patients" = sapply(
      split(cond_uniq_sites, cond_uniq_sites$celltype), 
      function(x) length(unique(x$patient)) )
  )
  
  celltype_summary <- filter(celltype_summary, num_Patients > 0)
  
  celltype_summary$clustersRepresented <- sapply(
    celltype_summary$celltype, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$celltype == x]
      hits <- findOverlaps(sites, red_clusters)
      length(unique(subjectHits(hits)))
    })
  
  celltype_summary$numSitesInClusters <- sapply(
    celltype_summary$celltype, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$celltype == x]
      hits <- findOverlaps(sites, red_clusters)
      length(unique(queryHits(hits)))
    })
  
  celltype_summary$abundInClusters <- sapply(
    celltype_summary$celltype, function(x){
      sites <- cond_uniq_sites[cond_uniq_sites$celltype == x]
      hits <- findOverlaps(sites, red_clusters)
      sum(sites[queryHits(hits)]$estAbund)
    })
    
  gen_epi_celltype_summary <- gen_epi_stats %>%
    dplyr::select(-posid) %>%
    dplyr::left_join(
      dplyr::select(specimen_data, specimen, celltype), 
      by = "specimen"
    ) %>% 
    dplyr::select(-specimen) %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise_all(mean, na.rm = TRUE)
  
  celltype_summary <- dplyr::left_join(
    celltype_summary, gen_epi_celltype_summary, by = "celltype")
  
  write.csv(
    celltype_summary,
    file = file.path(outputDir, "cart19_celltype_summary.csv"),
    quote = TRUE,
    row.names = FALSE
  )
  
  # Save all summaries as a single object ----------------------------------------
  summaries <- list(
    "timepoint" = timepoint_summary,
    "specimen" = specimen_data,
    "patient" = patient_summary,
    "celltype" = celltype_summary)
  
  saveRDS(summaries, file = file.path(outputDir, "cart19_summaries.rds"))
  
  # Top ranked genes from each patient (by estAbund and relAbund) ----------------
  sites_by_patient <- split(
    cond_uniq_sites[cond_uniq_sites$timepoint != "d0"], 
    cond_uniq_sites[cond_uniq_sites$timepoint != "d0"]$patient)
  
  top_one_pc_sonicAbund <- paste(
    unlist(lapply(
      sites_by_patient,
      get_top_genes,
      genes = refGenes,
      percent = 1,
      rank_by = "estAbund")),
    collapse = ", ")
  
  write.table(
    top_one_pc_sonicAbund, 
    file = file.path(
      outputDir, "top_one_percent_of_gene_ids_by_sonicAbund.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE)
  
  top_ten_pc_sonicAbund <- paste(
    unlist(lapply(
      sites_by_patient,
      get_top_genes,
      genes = refGenes,
      percent = 10,
      rank_by = "estAbund")),
    collapse = ", ")
  
  write.table(
    top_ten_pc_sonicAbund, 
    file = file.path(
      outputDir, "top_ten_percent_of_gene_ids_by_sonicAbund.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE)
  
  top_one_pc_relAbund <- paste(
    unlist(lapply(
      sites_by_patient,
      get_top_genes,
      genes = refGenes,
      percent = 1,
      rank_by = "relAbund")),
    collapse = ", ")
  
  write.table(
    top_one_pc_relAbund, 
    file = file.path(
      outputDir, "top_one_percent_of_gene_ids_by_relAbund.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE)
  
  top_ten_pc_relAbund <- paste(
    unlist(lapply(
      sites_by_patient,
      get_top_genes,
      genes = refGenes,
      percent = 10,
      rank_by = "relAbund")),
    collapse = ", ")
  
  write.table(
    top_ten_pc_relAbund, 
    file = file.path(
      outputDir, "top_ten_percent_of_gene_ids_by_relAbund.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE)

}else{
  
  cond_uniq_annot <- readRDS(
    file.path(outputDir, "cart19_annotated_uniq_sites.rds")
  )
  
  summaries <- readRDS(file.path(outputDir, "cart19_summaries.rds"))
  
}