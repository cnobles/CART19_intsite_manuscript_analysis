#### This R-script is for processing and annotating the unique integration sites
#### This script needs to be run in the proper environment and depends on
#### cart19_intsite_analysis.Rmd.
output_files <- c(
  "condensed_intsites.rds", "cart19_condensed_intsites.csv",
  "cart19_tdn_sites.rds", "cart19_timepoint_sites.rds"
)

if(
  !all(sapply(output_files, function(x) file.exists(file.path(outputDir, x))))
){
  # Load data files ------------------------------------------------------------
  data_files <- list.files(outputDir)
  file_uniq_sites <- grep("filtered_unique_intsites", data_files, value = TRUE)
  
  uniq_sites <- readRDS(file.path(outputDir, file_uniq_sites))
  uniq_sites$celltype <- specimen_data$celltype[
    match(uniq_sites$specimen, specimen_data$specimenaccnum)]
  
  # Condense the integration site ranges to unique locations and abundances ----
  uniq_list <- split(uniq_sites, uniq_sites$gtsp)
  
  cond_uniq_sites <- unlist(GRangesList(lapply(
    uniq_list,
    function(x){
      sites <- condense_intsites(
        sites.to.condense = x,
        return.abundance = TRUE,
        method = "fragLen",
        replicates = "samplename"
      )
      sites$samplename <- NULL
      sites$sampleid <- NULL
      sites$count <- NULL
      sites$miseqid <- NULL
      sites$gender <- NULL
      sites$siteid <- NULL
      sites$clusID <- NULL
      sites$called.pos <- NULL
      sites$adj.pos <- NULL
      sites$gtsp <- NULL
      sites
    }
  )))
  
  # Annotate Sites with Gene names for within and nearest genes ----------------
  cond_uniq_sites <- getSitesInFeature(
    cond_uniq_sites, 
    refGenes, 
    colnam = "in_gene", 
    feature.colnam = "name2"
  )
  cond_uniq_sites <- getNearestFeature(
    cond_uniq_sites,
    refGenes,
    colnam = "nearest_gene",
    feature.colnam = "name2"
  )
  
  ## Add gene marks ============================================================
  ## ("*" for in_gene, "~" for onco gene, and "!" for badActors)
  cond_uniq_sites$gene_id_wo_annot <- ifelse(
    cond_uniq_sites$in_gene == "FALSE",
    cond_uniq_sites$nearest_gene,
    cond_uniq_sites$in_gene
  )
  
  cond_uniq_sites$gene_id_wo_annot <- sapply(
    strsplit(cond_uniq_sites$gene_id_wo_annot, ","), "[[", 1)
  
  cond_uniq_sites$gene_id <- paste0(cond_uniq_sites$gene_id_wo_annot, " ")
  
  cond_uniq_sites$gene_id <- ifelse(
    cond_uniq_sites$in_gene == "FALSE",
    cond_uniq_sites$gene_id,
    paste0(cond_uniq_sites$gene_id, "*")
  )
  
  cond_uniq_sites$gene_id <- ifelse(
    cond_uniq_sites$gene_id_wo_annot %in% oncoGenes,
    paste0(cond_uniq_sites$gene_id, "~"),
    cond_uniq_sites$gene_id
  )
  
  cond_uniq_sites$gene_id <- ifelse(
    cond_uniq_sites$gene_id_wo_annot %in% badActors,
    paste0(cond_uniq_sites$gene_id, "!"),
    cond_uniq_sites$gene_id
  )
  
  cond_uniq_sites$gene_id <- ifelse(
    stringr::str_detect(cond_uniq_sites$gene_id, "[\\s]$"),
    cond_uniq_sites$gene_id_wo_annot,
    cond_uniq_sites$gene_id
  )
  
  ## Remove artifact sites =====================================================
  ## PCCA and NRP1 annotated expansions were found in Bcell malignacies, not 
  ## Tcells. These two expansions will be removed from the analysis for now.
  excluded_sites <- c("PCCA" = "chr13-100350286", "NRP1" = "chr10+33114980")
  cond_uniq_sites <- cond_uniq_sites[!cond_uniq_sites$posid %in% excluded_sites]
  
  excluded_specimens <- specimen_data %>%
    dplyr::filter(patient %in% c("p04409-01", "p04409-02")) %>%
    dplyr::arrange(timepoint, celltype) %>%
    dplyr::mutate(fct_timepoint = as.integer(timepoint)) %>%
    dplyr::filter(fct_timepoint > grep("m12", timepointLevels)) %$%
    specimenaccnum
  
  cond_uniq_sites <- cond_uniq_sites[
    !cond_uniq_sites$specimen %in% excluded_specimens]
  
  # Write condensed sites to file ----------------------------------------------
  saveRDS(
    cond_uniq_sites,
    file = file.path(outputDir, "condensed_intsites.rds")
  )
  
  write.csv(
    as.data.frame(cond_uniq_sites, row.names = NULL),
    file = file.path(outputDir, "cart19_condensed_intsites.csv"),
    quote = TRUE,
    row.names = FALSE
  )
  
  # Split unique sites to pre and post transduction sets
  tdn_sites <- cond_uniq_sites[cond_uniq_sites$timepoint == "d0"]
  tdn_sites <- tdn_sites[tdn_sites$patient %in% std_clin_patients]
  timepoint_sites <- cond_uniq_sites[cond_uniq_sites$timepoint != "d0"]
  timepoint_sites <- timepoint_sites[
    timepoint_sites$patient %in% std_clin_patients]
  
  saveRDS(
    tdn_sites, file = file.path(outputDir, "cart19_tdn_sites.rds"))
  
  saveRDS(
    timepoint_sites, file = file.path(outputDir, "cart19_timepoint_sites.rds"))
}

