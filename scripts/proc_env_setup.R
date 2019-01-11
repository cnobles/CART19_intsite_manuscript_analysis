# proc_env_setup.R
# This script loads required packages and reference data into a processing
# environment. It is sourced by initiate_analysis.R.

options(stringsAsFactors = FALSE, scipen = 99, width = 999)

packs <- c(
  "BiocGenerics", "Biostrings", "GenomicRanges", "geneRxCluster", "hiAnnotator",
  "igraph", "Matrix", "parallel", "reldist", "reshape2", "RMySQL", "data.table",
  "sonicLength", "vegan", "lubridate", "gintools", "spraphal", "ggrepel", 
  "scales", "grid", "gridExtra", "RColorBrewer", "knitr", "magrittr", "pander",
  "KEGGREST", "GO.db", "UpSetR", "BiasedUrn", "foreach", "tidyverse")

packsLoaded <- suppressMessages(
  sapply(packs, require, character.only = TRUE)
)

if( !all(packsLoaded) ){
  
  pandoc.table(data.frame(
    "R-Packages" = names(packsLoaded), 
    "Loaded" = packsLoaded, 
    row.names = NULL
  ))
  
  stop("Check dependancies.")
  
}

## Processing information ----
genomicFreeze <- "hg38"

analysisDate <- Sys.Date()

trial <- "CART19"
std_clin_trials <- c("959", "04409", "03712")

set.seed(1234)

## Reference and supporting files ----
source(file.path(scriptDir, "supporting_functions.R"))
source(file.path(scriptDir, "supporting_goa_functions.R"))
source(file.path(scriptDir, "GenomicHeatmapMaker/utils.R"))
source(file.path(scriptDir, "EpigeneticHeatmapMaker/utils.R"))

key_gois <- c("TET2")

if( file.exists(file.path(outputDir, "cart19_intsite_sample_list.csv")) ){
  
  intsiteSpecimenMetadata <- read.csv(
    file.path(workingDir, "data/cart19_intsite_sample_list.csv")) %>%
    mutate(clin_trial = str_extract(Patient_ID, "[0-9]+"))
  
}else{
  
  stop("Cannot find 'data/cart19_intsite_sample_list.csv' file.")
  
}

refGenes <- readRDS(file.path(utilsDir, "hg38.refSeq.rds"))

refGenes <- refGenes[
  seqnames(refGenes) %in% paste0("chr", c(1:22, "X", "Y", "M"))
  ]

oncoGenesData <- read.delim(
  file.path(utilsDir, "allOnco.human.v3.tsv"),
  header = TRUE, 
  sep = "\t",
  stringsAsFactors = FALSE
)

oncoGenes <- unique(oncoGenesData[,"symbol"])
nonOncoGenes <- unique(refGenes$name2[!refGenes$name2 %in% oncoGenes])

badActors <- read.delim(
  file.path(utilsDir, "humanLymph.v1.list"),
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)[,1]

genomic_features_path <- file.path(scriptDir, "GenomicHeatmapMaker")
epigenetic_features_path <- "/data/internal/epigeneticData/hg38" # Not archivable

epigenetic_features_files <- list.files(epigenetic_features_path) %>%
  grep("HeLa", ., invert = TRUE, value = TRUE) %>%
  grep("CD133", ., invert = TRUE, value = TRUE) %>%
  grep("HEK293T", ., invert = TRUE, value = TRUE) %>%
  grep("Jurkat", ., invert = TRUE, value = TRUE) %>%
  file.path(epigenetic_features_path, .)

# Import HGNC reference data for annotation and consistency ----
hgnc_complete <- fread(
  paste0("zcat ", file.path(utilsDir, "hgnc_complete_set.180207.txt.gz")),
  sep = "\t", header = TRUE, 
  select = c(
    "HGNC ID", "Approved Symbol", "Approved Name", "Locus Group", "Locus Type", 
    "Synonyms", "Previous Symbols", "Entrez Gene ID", "Ensembl Gene ID", 
    "RefSeq (supplied by NCBI)", "UniProt ID (supplied by UniProt)"),
  data.table = FALSE
)

names(hgnc_complete) <- c(
  "hgnc_id", "symbol", "name", "locus_group", "locus_type", "alias_symbol", 
  "prev_symbol", "entrez_id", "ensembl_gene_id", "refseq_accession", 
  "uniprot_ids"
)

hgnc_complete <- dplyr::filter(hgnc_complete, !grepl("withdrawn", symbol)) %>%
  dplyr::mutate(
    kegg_id = paste0("hsa:", entrez_id),
    entrez_id = paste0(entrez_id, ":EZID"),
    alias_symbol = gsub(", ", "|", alias_symbol),
    prev_symbol = gsub(", ", "|", prev_symbol),
    extended_alias = paste0(
      alias_symbol, "|", prev_symbol, "|", ensembl_gene_id, "|", 
      refseq_accession, "|", uniprot_ids)
  )

## Develop standard factor scales for celltypes and timepoints ----
celltypeLevels <- c(
  "PB", "PBMC", "PBL", "Whole Blood", "Tcells", "Tcells:CAR+", 
  "Tcells:CAR+CD4+", "Tcells:CAR+CD8+", "Tcells:CAR+CD8-", "Tcells:CAR-CD4+", 
  "Tcells:CD4+SP", "Tcells:CAR-CD8+", "Tcells:CAR-CD8-", "Tcells:CD4+", 
  "Tcells:CD8+", "Tcells:CD8+Naive", "Tcells:CD8+Tscm", "Tcells:CD8+Tcm", 
  "Tcells:CD8+Tem", "Tcells:CD8+Te", "Tcells:CD8+Tm", "Tcells:CD4+CD8+", 
  "Tcells:CD4+CD8+DP", "Tcells:CD4+CD8+DN", "Bone Marrow", "BM", "BMMC", 
  "BM:CAR+", "CD3-"
)

timepointLevels <- c(
  "d-10", "d-1", "d0", "d1", "d5", "d7", "d9", "d10", "d11", "d13", "d14", 
  "d15", "d17", "d21", "d23", "d25", "d28", "d30", "d35", "d36", "d42", "d49",
  "d50", "m2", "d63", "d75", "d90", "m3", "d92", "d120", "d121", "m4", "d133",
  "d147", "m5", "d169", "m6", "d204", "m9", "m12", "d442", "m15", "m18", 
  "y1.5", "m20", "m21", "d720", "m24", "y2", "d801",  "y2.5", "m32", "y3", 
  "y4", "d1584", "y4.5", "m60", "y5", "y5.5", "y6", "y6.5", "y7", "y8"
) 

# Load patient data ----
patient_data <- read.csv(
    file.path(outputDir, "cart19_intsite_sample_list.csv")
  ) %>%
  dplyr::distinct(Patient_ID, Disease, Response_to_Treatment) %>%
  dplyr::rename(
    "patient" = Patient_ID, 
    "disease" = Disease, 
    "response" = Response_to_Treatment) %>%
  dplyr::mutate(
    clin_trial = str_extract(patient, "[0-9]+"),
    patient = gsub("_", "-", patient),
    disease = gsub("Adult ALL", "aALL", disease),
    disease = gsub("Pediatric ALL", "pALL", disease),
    disease = factor(disease, levels = c("pALL", "aALL", "CLL")),
    response = gsub(" but relapse", "\nw relapse", response),
    response = gsub(" - with transformed dz", "\nw TnDz", response),
    response = factor(
      response, 
      levels = c(
        "None", "Partial", "Partial\nw TnDz", 
        "Complete", "Complete\nw relapse")),
    simple_response = factor(
      str_extract(as.character(response), "[\\w]+"), 
      levels = c("None", "Partial", "Complete")),
    general_response = factor(ifelse(
      response == "None", "Non-responder", "Responder"), 
      levels = c("Non-responder", "Responder")),
    determinant_response = factor(ifelse(
      response %in% c("Complete", "Partial\nw TnDz", "Complete\nw relapse"), 
      "CR_PRtd", "PR_NR"), levels = c("CR_PRtd", "PR_NR"))) %>%
  dplyr::arrange(disease, response)

patient_data[patient_data$patient == "p03712-12",]$response <- "None"
patient_data[patient_data$patient == "p03712-12",]$simple_response <- "None"
patient_data[patient_data$patient == "p03712-12",]$general_response <- "Non-responder"

std_clin_patients <- patient_data %>%
  dplyr::filter(clin_trial %in% std_clin_trials) %$%
  patient

CR_pats <- patient_data$patient[patient_data$determinant_response == "CR_PRtd"]
NR_pats <- patient_data$patient[patient_data$determinant_response == "PR_NR"]

## Load integration site analysis data ----
data_files <- list.files(outputDir)

### Load specimen_data
file <- grep("specimen_data", data_files, value = TRUE)
specimen_data <- readRDS(file.path(outputDir, file)) %>%
  dplyr::mutate(
    celltype = ifelse(
      as.character(celltype) == "PBMC", "PBL", as.character(celltype)),
    celltype = factor(celltype, levels = celltypeLevels))

if( !all(patient_data$patient %in% specimen_data$patient) ){
  warning("Not all specimens in database, using custom data.")
  specimen_data <- dplyr::select(
      intsiteSpecimenMetadata, 
      Trial, Patient_ID, Abv_Cell_Type, 
      Sorting_Parameters, Timepoint, GTSP
    ) %>%
    dplyr::rename(
      "trial" = Trial, "patient" = Patient_ID, 
      "timepoint" = Timepoint, "specimenaccnum" = GTSP) %>%       
    dplyr::mutate(
      celltype = ifelse(
        grepl("^s", Abv_Cell_Type),
        paste0(gsub("^s", "", Abv_Cell_Type), ":", Sorting_Parameters),
        Abv_Cell_Type),
      celltype = ifelse(
        as.character(celltype) == "PBMC", "PBL", as.character(celltype)),
      celltype = ifelse(
        as.character(celltype) == "Whole Blood", 
        "PBL", as.character(celltype)),
      celltype = factor(celltype, levels = celltypeLevels),
      timepoint = factor(timepoint, levels = timepointLevels)) %>%
    dplyr::select(trial, patient, celltype, timepoint, specimenaccnum) %>% 
    dplyr::arrange(specimenaccnum)
}

