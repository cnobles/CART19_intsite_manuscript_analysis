#' Initial processing for CART19 Study
#' This script should only be used for pulling data from a SQL database, as
#' integration sites are housed in the Bushman lab, and outputing the 
#' unprocessed integration sites, unique and multihit integration sites, and 
#' contamination filtered integration sites. Additionally, the specimen info
#' will be generated and the vcn data will be assembled into a single object.
#' 
#' All data will be saved in the data directory of the cloned repository.
#' Additional analysis should rely on that data. Trying to reprocess from the 
#' unprocessed data may take some time as this script is paralleled with 
#' r-parallel. If reprocessing, make sure to change the number of cores 
#' requested.
#' 

## Initial check to see if any work needs to be done.
outputFiles <- list.files(outputDir)

if( all(c(
  any(grepl("specimen_data.rds", outputFiles)),
  any(grepl("unprocessed_intsites.rds", outputFiles)),
  any(grepl("[0-9]_unique_intsites.rds", outputFiles, perl = TRUE)),
  any(grepl("multihits.rds", outputFiles)),
  any(grepl("vcn_data.rds", outputFiles)),
  any(grepl("filtered_unique_intsites.rds", outputFiles))
)) ){
  
  cat("[", paste(Sys.time()), "] All initial processing files identified, ending initial processing.\n")
  
}else{

  cat("[", paste(Sys.time()), "]\n",
      "Initial Processing for CART19 Specimens\n",
      "Author: Christopher Nobles, Ph.D.\n",
      "Process Date: ", as.character(Sys.Date()),
      "\n\n***************************************\n\n"
  )
  
  # Processing setup ----
  ## R-packages required for analysis ----
  options(stringsAsFactors = FALSE, scipen = 99)
  suppressMessages(library("pander"))
  
  cat("[", paste(Sys.time()), "] Loading R-packages.\n")
  
  packs <- packs <- c(
    "BiocGenerics", "Biostrings", "GenomicRanges", "gintools", "hiAnnotator", 
    "igraph", "Matrix", "parallel", "reldist", "reshape2", "RMySQL", 
    "sonicLength", "vegan", "lubridate", "forcats", "magrittr", "dplyr"
  )
  
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
    
  }else{
    
    cat("[", paste(Sys.time()), "] All R-packages loaded successfully.\n")
    
  }
  
  ## Processing information ----
  genomicFreeze <- "hg38"
  
  connections <- list(
    specimen_db_group = "specimen_management",
    intsites_db_group = "intsites_cart_20180808",
    intsites_db_name = "intsites_cart_20180808"
  )
  
  analysisDate <- as.character(Sys.Date())
  
  trial <- "CART19"
  include_TCR_IGH_Data <- FALSE
  
  set.seed(1234)
  
  paramTable <- data.frame(
    "Parameters" = c(
      "Genome:", "Specimen_DB Group:", "IntSites_DB Group:", 
      "IntSites_DB Name:", "Maximum Cores:", "Analysis Date:", 
      "Clinical Trial:", "Include TCR/IGH:"
    ),
    "Values" = c(
      genomicFreeze, unlist(connections), numCores, analysisDate, 
      trial, include_TCR_IGH_Data
    ),
    row.names = NULL
  )
  
  pandoc.table(
    paramTable,
    style = "simple",
    emphasize.rownames = FALSE,
    justify = "left"
  )
  
  rm(paramTable)
  
  pander(
    data.frame(
      "Directories" = c("Working:", "Utilities:", "Output:", "VCN Data:"),
      "Paths" = c(workingDir, utilsDir, outputDir, vcnDir),
      row.names = NULL
    ),
    style = "simple",
    emphasize.rownames = FALSE,
    justify = "left"
  )
  
  ## Reference and supporting files ----
  intsiteSpecimenMetadata <- read.csv(
    file.path(outputDir, "cart19_intsite_sample_list.csv")
  )
  
  refGenes <- readRDS(
    file.path(utilsDir, "hg38.refSeq.rds")
  )
  
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
    "d147", "m5", "d169", "m6", "d204", "m9", "m12", "y1", "d442", "m15", "m18", 
    "y1.5", "m20", "m21", "d720", "m24", "y2", "d801",  "y2.5", "m32", "y3", 
    "y4", "d1584", "y4.5", "m60", "y5", "y5.5", "y6", "y6.5", "y7", "y8"
  ) 
  
  # Process data from database and acquired files ----
  ## Acquire all patient data available in the specimen management database ----
  if( any(grepl("specimen_data.rds", list.files(outputDir))) ){
    
    cat("[", paste(Sys.time()), "] Specimen data found, loading...\n")
    specimenData <- readRDS(file.path(
      outputDir, 
      grep("specimen_data.rds", list.files(outputDir), value = TRUE)
    ))
    
  }else{
  
    cat("[", paste(Sys.time()), "] Connecting to Specimen Database...\n")
    junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
    dbConn <- dbConnect(MySQL(), group = connections$specimen_db_group)
    stopifnot(dbGetQuery(dbConn, "SELECT 1") == 1)
    
    metaCols <- c(
      "Trial", "Patient", "CellType", "Timepoint", "VCN", "SpecimenAccNum"
    )
    
    querySelection <- paste(
      "SELECT", 
      paste(metaCols, collapse = ", "),
      "FROM specimen_management.gtsp",
      sep = " "
    )
    
    queryCondition <- paste0("WHERE Trial = '", trial, "'")
    query <- paste(querySelection, queryCondition, sep = " ")
    specimenData <- dbGetQuery(dbConn, query)
    names(specimenData) <- tolower(names(specimenData))
    
    cat(
      "[", paste(Sys.time()), "] Specimen information includes data for ", 
      nrow(specimenData), 
      " samples.\n"
    )
    
    # Special condition, same patient with two identifiers.
    specimenData$patient <- gsub(
      "p36908-115", "p04409-10", specimenData$patient
    )
    
    specimenData$timepoint <- factor(
      specimenData$timepoint, levels = timepointLevels
    )
    
    specimenData$celltype <- factor(
      specimenData$celltype, levels = celltypeLevels
    )
    
    null <- dbDisconnect(dbConn)
    rm(query, queryCondition, querySelection, metaCols, dbConn, junk)
    
    if( exists("specimenData") ){
      
      if( all(intsiteSpecimenMetadata$GTSP %in% specimenData$specimenaccnum) ){
        cat("[", paste(Sys.time()), "] Specimen data acquired from database.\n")
        
      }else{
        
        warning("Missing specimens in the database, using input annotations.")
        
        specimenData <- dplyr::select(
            intsiteSpecimenMetadata, 
            Trial, Patient_ID, Abv_Cell_Type, 
            Sorting_Parameters, Timepoint, GTSP) %>%
          dplyr::rename(
            "trial" = Trial, "patient" = Patient_ID, 
            "timepoint" = Timepoint, "specimenaccnum" = GTSP) %>%       
          mutate(
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
          arrange(specimenaccnum)
        
      }
      
    }else{
      
      stop("Specimen data not acquired from Specimen Database.")
      
    }
    
    # Save specimen data.
    saveRDS(
      specimenData,
      file = file.path(
        outputDir,
        paste0(analysisDate, "_specimen_data.rds")
      )
    )
    
  }
  
  
  ## Acquire all integration sites from patients ----
  if( any(grepl("unprocessed_intsites.rds", list.files(outputDir))) ){
    
    cat("[", paste(Sys.time()), "] Unprocessed integration site data found, loading...\n")
    dbIntsitesData <- readRDS(file.path(
      outputDir, 
      grep("unprocessed_intsites.rds", list.files(outputDir), value = TRUE)
    ))
    
  }else{
    
    cat("[", paste(Sys.time()), "] Connecting to Integration Site Database...\n")
    dbConn <- dbConnect(MySQL(), group = connections$intsites_db_group)
    stopifnot(dbGetQuery(dbConn, "SELECT 1") == 1)
    
    queryUniq <- sprintf(
      "SELECT * FROM %1$s.samples
       JOIN %1$s.sites ON %1$s.samples.sampleID = %1$s.sites.sampleID
       JOIN %1$s.pcrbreakpoints ON %1$s.pcrbreakpoints.siteID = %1$s.sites.siteID",
       connections$intsites_db_name, 
       trial)
    
    uniqSitesData <- dbGetQuery(dbConn, queryUniq)
    uniqSitesData <- uniqSitesData[, !duplicated(colnames(uniqSitesData))]
    
    cat(
      "[", paste(Sys.time()),
      "] Acquired unique integration site data from", 
      format(sum(uniqSitesData$count), big.mark = ","),
      "sequence reads.\n"
    )
    
    queryMulti <- sprintf(
      "SELECT * FROM %1$s.samples 
       JOIN %1$s.multihitpositions 
       ON %1$s.samples.sampleID = %1$s.multihitpositions.sampleID 
       JOIN %1$s.multihitlengths 
       ON %1$s.multihitpositions.multihitID = %1$s.multihitlengths.multihitID",
      connections$intsites_db_name, 
      trial)
    
    multihitData <- dbGetQuery(dbConn, queryMulti)
    multihitData <- multihitData[, !duplicated(colnames(multihitData))]
    multihitSummary <- distinct(multihitData, multihitID, length, count)
    
    cat(
      "[", paste(Sys.time()), "] Acquired multihit integration site data from",
      format(sum(multihitSummary$count), big.mark = ","),
      "sequence reads.\n"
    )
    
    rm(multihitSummary)
    
    dbIntsitesData <- list(
      "uniqSites" = uniqSitesData, 
      "multihits" = multihitData)
    
    null <- dbDisconnect(dbConn)
    rm(uniqSitesData, multihitData, queryUniq, queryMulti, dbConn)
    
    if( exists("dbIntsitesData") ){
      if( nrow(dbIntsitesData$uniqSites) > 0 & 
          nrow(dbIntsitesData$multihits) > 0 
        ){
        
        cat("[", paste(Sys.time()), "] Integration site data acquired from database.\n")
        
        acq_specimens <- unique(
          stringr::str_extract(dbIntsitesData$uniqSites$sampleName, "[\\w]+")
        )
        
        cat(
          "[", paste(Sys.time()), "] Integrations sites acquired for ", 
          sum(acq_specimens %in% specimenData$specimenaccnum), 
          " out of ", nrow(specimenData), "."
        )
        
      }else{
        
        stop("Partial or no acquisition of intSite data from database.")
        
      }
      
    }else{
      
      stop("IntSite data not acquired from IntSite Database.")
      
    }
  
    ## Save unprocessed data
    saveRDS(
      dbIntsitesData, 
      file = file.path(
        outputDir, 
        paste0(analysisDate, "_unprocessed_intsites.rds")
      )
    )
    
  }
    
  
  ## Join GTSP, Patient, CellType, and Timepoint info with sites ----
  if( any(grepl("[0-9]_unique_intsites", list.files(outputDir), perl = TRUE)) &
      any(grepl("multihits.rds", list.files(outputDir))) ){
    
    cat("[", paste(Sys.time()), "] Processed integration site data found, loading...\n")
    
    intsitesData <- list(
      "uniqSites" = readRDS(file.path(
        outputDir, 
        grep(
          "[0-9]_unique_intsites", list.files(outputDir), 
          perl = TRUE, value = TRUE
        )
      )),
      "multihits" = readRDS(file.path(
        outputDir, 
        grep(
          "multihits.rds", list.files(outputDir), 
          value = TRUE
        )
      ))
    )
    
    stopifnot(all(sapply(intsitesData, function(x) class(x) == "GRanges")))
    
    save.image(file.path(outputDir, "temp.image.RData"))
    
  }else{
    
    cat("[", paste(Sys.time()), "] Processing integration site data...\n")
    intsitesData <- lapply(
      seq_along(dbIntsitesData), 
      function(i){
        
        sites <- dbIntsitesData[[i]]
        sites <- sites[, !duplicated(colnames(sites))]
        sites$gtsp <- stringr::str_extract(as.character(sites$sampleName), "[\\w]+")
        
        sites <- left_join(
          sites,
          specimenData[,c("patient", "celltype", "timepoint", "specimenaccnum")],
          by = c("gtsp" = "specimenaccnum")
        )
        
        sites <- db_to_granges(sites)
        sites <- split(sites, sites$patient)
        
        cluster <- makeCluster(numCores)
        clusterExport(
          cl = cluster, 
          varlist = c(),
          envir = environment())
        
        cat(
          "[", paste(Sys.time()), "] Refining integration template breakpoints (", 
          i, " of ", length(dbIntsitesData), ").\n")
        
        sites <- unlist(GenomicRanges::GRangesList(
            parLapply(
              cluster,
              sites,
              function(x){
                library(gintools)
                message(unique(x$patient))
                grl <- split(x, x$samplename)
                gr <- unlist(GenomicRanges::GRangesList(lapply(grl, function(x){
                  message(unique(x$samplename))
                  refine_breakpoints(x, counts = "count")})))
                gr$siteid <- NULL
                unique_granges(gr, sum.cols = "count")
              })))
        
        sites <- split(sites, sites$patient)
        
        cat(
          "[", paste(Sys.time()), "] Standardizing integration template positions (", 
          i, " of ", length(dbIntsitesData), ").\n")
        
        sites <- unlist(GenomicRanges::GRangesList(
            parLapply(
              cluster, 
              sites,
              function(x){
                library(gintools)
                gr <- standardize_sites(x)
                unique_granges(gr, sum.cols = "count")
              })))
        
        names(sites) <- NULL
        stopCluster(cl = cluster)
        sites
        
      }
    )
  
    cat("[", paste(Sys.time()), "] Breakpoint refinement and positional standardizing completed.\n")
    
    names(intsitesData) <- c("uniqSites", "multihits")
    
    stopifnot(all(sapply(intsitesData, function(x) class(x) == "GRanges")))
    
    cat("[", paste(Sys.time()), "] Writing unique integration site data to output directory.\n")
    
    saveRDS(
      intsitesData$uniqSites, 
      file = file.path(
        outputDir, 
        paste0(analysisDate, "_unique_intsites.rds")
      )
    )
  
    save.image(file.path(outputDir, "temp.image.RData"))
    
    cat("[", paste(Sys.time()), "] Normalizing multihit clusters...\n")
    
    intsitesData$multihits$pos.clus <- as.integer(factor(
      generate_posid(intsitesData$multihits)
    ))
    
    intsitesData$multihits <- gintools:::normalize_multihit_clusters(
      intsitesData$multihits,
      grouping = "patient", 
      cores = numCores
    )
    
    cat("[", paste(Sys.time()), "] Multihit cluster IDs normalized.\n")
    
    intsites_data <- within(intsitesData, {
      multihits$called.pos <- multihits$called.bp <- multihits$norm.pos <- NULL
    })
    
    saveRDS(
      intsitesData$multihits,
      file = file.path(
        outputDir, 
        paste0(analysisDate, "_multihits.rds")
      )
    )
  
  }
  
  ## Load and format VCN data ----
  if( any(grepl("vcn_data.rds", list.files(outputDir))) ){
    
    cat("[", paste(Sys.time()), "] VCN data found, loading...\n")
    
    vcnData <- readRDS(file.path(
      outputDir, 
      grep("vcn_data.rds", list.files(outputDir), value = TRUE)
    ))
    
  }else{
    
    cat("[", paste(Sys.time()), "] Loading and consolidating VCN data...\n")
    vcnSets <- grep("cur.csv", list.files(vcnDir), value = TRUE)
    patients <- sapply(strsplit(vcnSets, "-vcn"), "[[", 1)
    
    vcnData <- bind_rows(lapply(
      patients, 
      function(x){
        
        data <- read.csv(
          file.path(vcnDir, paste0(x, "-vcn_cur.csv")),
          header = TRUE,
          stringsAsFactors = FALSE
        )
        colnames(data) <- tolower(colnames(data))
        data
        
      }
    ))
    
    rm(patients, vcnSets)
    
    cat("[", paste(Sys.time()), "] Writing VCN data to output directory.\n")
    
    saveRDS(
      vcnData,
      file = file.path(
        outputDir,
        paste0(analysisDate, "_vcn_data.rds")
      )
    )
    
  }
  
  # Contamination filtering of unique integration sites ----
  if( any(grepl("filtered_unique_intsites.rds", list.files(outputDir))) ){
    
    cat("[", paste(Sys.time()), "] Crossover filtered data found, loading...\n")
    
    filteredUniqSites <- readRDS(file.path(
      outputDir, 
      grep("filtered_unique_intsites.rds", list.files(outputDir), value = TRUE)
    ))
    
  }else{
  
    cat("[", paste(Sys.time()), "] Conducting contamination filtering of unique integration sites...\n")
    
    uniqSites <- intsitesData$uniqSites
    uniqSites$usid <- seq_along(uniqSites)
    
    ## Filter probably contaminating sites out from samples
    testSites <- split(uniqSites, uniqSites$patient)
    crossoverSites <- track_clones(testSites, gap = 0L)
    
    cat("[", paste(Sys.time()), "] Number of crossover sites identified: ", length(crossoverSites), "\n")
    
    ## Identify contaminates by shared break points
    bpsDuplicated <- sapply(
      crossoverSites, function(siteGrp) any(duplicated(siteGrp))
    )
    
    contamSites <- crossoverSites[bpsDuplicated]
    
    cat(
      "[", paste(Sys.time()), "] Number of crossover sites sharing breakpoints: ", 
      length(contamSites), "\n"
    )
    
    ## Correct the duplicated sites to the correct sample based on read counts
    # Only correct if duplicated sites have the same breakpoint and are in 
    # different patients. Requiring correction to occur within the same 
    # sequencing run. Requiring read counts to dictate original sample, if no 
    # difference then reads remain untouched. Correcting duplicated sites within
    # the same specimen will skew sonic abundance in expanded cases. Correcting 
    # duplicated sites across specimens of the same patient will also skew 
    # abundance.
    
    contamSites <- unlist(GRangesList(contamSites))
    
    contamSites$bpid <- generate_posid(
      seqnames = seqnames(contamSites),
      strand = strand(contamSites),
      start = end(contamSites),
      end = start(contamSites))
    
    contamSites <- split(contamSites, contamSites$bpid)
    
    contamSites <- contamSites[
      sapply(contamSites, function(x) any(duplicated(x)))
    ]
    
    contamSites <- contamSites[
      sapply(contamSites, function(x) length(unique(x$patient)) > 1)
    ]
    
    contamSites <- contamSites[
      sapply(contamSites, function(x) length(unique(x$miseqid)) == 1)
    ]
    
    contamSites <- contamSites[
      sapply(contamSites, function(x) length(unique(x$count)) > 1)
    ]
    
    cat("[", paste(Sys.time()), "] Final count of sites needing correction: ", length(contamSites), "\n")
    cat("[", paste(Sys.time()), "] Reassigning reads to correct patient...\n")
    
    modifiedSites <- lapply(
      contamSites, 
      function(siteGrp){
        
        modifiedUsids <- siteGrp$usid
        
        df <- data.frame(
            'patient' = siteGrp$patient,
            'readCounts' = siteGrp$count,
            stringsAsFactors = FALSE
          ) %>%
          dplyr::group_by(patient) %>%
          dplyr::summarise(readCounts = sum(readCounts))
        
        oriPat <- df[grep(max(df$readCounts), df$readCounts),]$patient
        oriPatReads <- siteGrp[siteGrp$patient == oriPat]
        topOriRead <- oriPatReads[oriPatReads$count == max(oriPatReads$count)]
        
        otherOriReads <- oriPatReads[
          oriPatReads$count != max(oriPatReads$count)
        ]
        
        false_reads <- siteGrp[siteGrp$patient != oriPat]
        topOriRead$count <- topOriRead$count + sum(false_reads$count)
        modifiedReads <- c(topOriRead, otherOriReads)
        
        return(
          list("modifiedReads" = modifiedReads, "modifiedUsids" = modifiedUsids)
        )
        
    })
    
    modifiedReads <- unlist(GRangesList(
      lapply(modifiedSites, "[[", 1)), use.names = FALSE)
    
    names(modifiedReads) <- NULL
    modifiedUsids <- suppressWarnings(
      do.call(c, lapply(modifiedSites, "[[", 2))
    )
    
    ## Removed modified reads from uniqSites and add back the modified reads
    modifiedReads$origin <- modifiedReads$posid <- modifiedReads$bpid <- NULL
    filteredUniqSites <- uniqSites[!uniqSites$usid %in% modifiedUsids]
    filteredUniqSites <- suppressWarnings(c(filteredUniqSites, modifiedReads))
    filteredUniqSites$usid <- NULL
    
    cat("[", paste(Sys.time()), "] Writing filtered unique integration sites to output directory.\n")
    
    saveRDS(
      filteredUniqSites, 
      file = file.path(
        outputDir, 
        paste0(analysisDate, "_filtered_unique_intsites.rds")
      )
    )
  
  }
  
  # Completed ----
  cat(
    "[ ", paste(Sys.time()), " ] Initial processing has competed.\n", 
    "[ ", paste(Sys.time()), " ] Output files can be found at: ", outputDir, "\n",
    sep = ""
  )
  
  if(all(c(
    any(grepl("specimen_data.rds", list.files(outputDir))),
    any(grepl("unprocessed_intsites.rds", list.files(outputDir))),
    file.exists(file.path(
      outputDir, paste0(analysisDate, "_unique_intsites.rds"))),
    file.exists(file.path(
      outputDir, paste0(analysisDate, "_multihits.rds"))),
    file.exists(file.path(
      outputDir, paste0(analysisDate, "_vcn_data.rds"))),
    file.exists(file.path(
      outputDir, paste0(analysisDate, "_filtered_unique_intsites.rds")))
    ))){
    
    system(paste0("rm ", file.path(outputDir, "temp.image.RData"))) 
    
  }

}