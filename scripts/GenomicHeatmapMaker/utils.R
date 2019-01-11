#    This source code file is a component of the larger INSPIIRED genomic analysis software package.
#    Copyright (C) 2016 Frederic Bushman
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

library(hotROCs)

getRefSeq_genes <- function(reference_genome) {
  # Identify table and track names
  bsession <- makeUCSCsession(reference_genome)
  trk_names <- names(trackNames(ucscTableQuery(bsession)))
  tbl_names <- tableNames(ucscTableQuery(bsession))

  trk <- grep("RefSeq", trk_names, value = TRUE)
  trk <- trk[trk %in% c("NCBI RefSeq", "RefSeq Genes")]
  stopifnot(length(trk) == 1)

  tbl <- grep("refGene", tbl_names, value = TRUE)
  stopifnot(length(tbl) == 1)

  refSeq <- makeGRanges(
    getUCSCtable(tbl, trk, freeze=reference_genome),
    freeze=reference_genome
  )
}

getCpG_islands <- function(reference_genome) {
  cpg <- getUCSCtable("cpgIslandExt", "CpG Islands", freeze=reference_genome)
  cpg$strand <- "*" # either strand
  makeGRanges(cpg, freeze=reference_genome, chromCol='chrom')
}

getDNaseI <- function(reference_genome) {
  DNaseI <- getUCSCtable("wgEncodeRegDnaseClustered", 
                         "DNase Clusters", freeze=reference_genome)
  DNaseI$strand <- "*" # either strand
  makeGRanges(DNaseI, freeze=reference_genome, chromCol='chrom')
}

#' for given samples pull sites from database and construct MRCs
get_sites_controls_from_db <- function(sampleName_GTSP, referenceGenome, connection) {
    if ( ! "label" %in% colnames(sampleName_GTSP)) {
        sampleName_GTSP$label <- sampleName_GTSP$GTSP
    }
    sampleName_GTSP <- select(sampleName_GTSP, sampleName, GTSP, label)

    # should have at least two samples
    stopifnot(length(unique(sampleName_GTSP$GTSP)) != 1)

    sampleName_GTSP$refGenome <- rep(referenceGenome, nrow(sampleName_GTSP))

    # samples should have sites
    stopifnot(nrow(getUniqueSiteCounts(sampleName_GTSP, connection)) > 1)
    # also we need at least several sites per sample/replicate
    stopifnot(is_enough_sites(sampleName_GTSP, connection))

    # check that all samples processed with the same reference genome
    is_in_db <- setNameExists(sampleName_GTSP, connection)
    if ( ! all(is_in_db)) {
        print("The following samples are NOT in the database")
        print(sampleName_GTSP[ ! is_in_db, ])
        stop()
    }
    #stopifnot(all(setNameExists(sampleName_GTSP, connection)))

    reference_genome_sequence <- get_reference_genome(referenceGenome)
    get_integration_sites_with_mrcs(sampleName_GTSP, reference_genome_sequence, connection)
}

get_sites_controls_from_file <- function(sampleName_GTSP, referenceGenome, sites){
    if ( ! "label" %in% colnames(sampleName_GTSP)) {
        sampleName_GTSP$label <- sampleName_GTSP$GTSP
    }
    sampleName_GTSP <- select(sampleName_GTSP, sampleName, GTSP, label)

    # should have at least two samples
    stopifnot(length(unique(sampleName_GTSP$GTSP)) != 1)

    sampleName_GTSP$refGenome <- rep(referenceGenome, nrow(sampleName_GTSP))

    # samples should have sites
    stopifnot(all(sapply(split(sites, sites$sampleName), nrow) > 0))
    # also we need at least several sites per sample/replicate

    sample_table <- table(as.character(sites$sampleName))
    message("\nSites observed per sample:")
    print(sample_table)
   
    if(any(sample_table < 3) | any(!sampleName_GTSP$label %in% names(sample_table))){
      print("Not enough sites (minimum 3) where found for each sample.")
      stop()
    }
    reference_genome_sequence <- get_reference_genome(referenceGenome)
    get_integration_sites_with_mrcs(sampleName_GTSP, reference_genome_sequence, sites)    
}

add_label <- function(sites, sampleName_GTSP) {
    sites_GTSP <- merge(sites, sampleName_GTSP)
    sites_GTSP$sampleName <- sites_GTSP$label
    sites_GTSP$refGenome <- NULL # not needed downstream
    sites_GTSP$GTSP <- NULL # not needed downstream
    sites_GTSP$label <- NULL
    sites_GTSP
}

get_integration_sites_with_mrcs <- function(
    sampleName_GTSP, refGenomeSeq, connection
) {
    if(class(connection) == "character"){
      sites <- getUniqueSites(sampleName_GTSP, connection)
    }else{
      sites <- connection
    }
    sites$type <- "insertion"
    
    if(class(connection) == "character"){
      sites <- add_label(sites, sampleName_GTSP)
      mrcs <- getMRCs(sampleName_GTSP, connection)
    }else{
      sites.metadata <- select(sites, siteID, gender, sampleName, refGenome)
      mrcs <- get_N_MRCs(sites.metadata[,c("siteID", "gender")], refGenomeSeq)
      mrcs <- merge(mrcs, sites.metadata[,c("siteID", "sampleName", "refGenome")])
    }

    mrcs$type <- "match"
    mrcs <- add_label(mrcs, sampleName_GTSP)
    sites <- select(sites, sampleName, siteID, chr, strand, position, type)
    sites_mrcs <- rbind(sites, mrcs)

    sites_mrcs <- makeGRanges(sites_mrcs, soloStart=TRUE,
        chromCol='chr', strandCol='strand', startCol='position')

    #seqinfo needs to be exact here or trimming will be wrong
    newSeqInfo <- seqinfo(refGenomeSeq)
    seqInfo.new2old <- match(seqnames(newSeqInfo),
        seqnames(seqinfo(sites_mrcs)))
    seqinfo(sites_mrcs, new2old=seqInfo.new2old) <- newSeqInfo

    sites_mrcs
}

get_annotation_columns <- function(sites) {
  granges_column_names <- c("seqnames", "start", "end", "width", "strand")
  int_site_column_names <- c("siteID", "sampleName", "chr", "strand", "position")
  required_columns <- unique(c(
    granges_column_names, int_site_column_names, "type"))
  stopifnot(all(required_columns %in% names(sites)))
  setdiff(names(sites), required_columns)
}

from_counts_to_density <- function(sites, column_prefix, window_size) {
  metadata <- mcols(sites)
  sapply(seq(window_size), function(i) {
    val <- window_size[i]
    name <- names(window_size)[i]
    column_name <- paste0(column_prefix, ".", name)
    metadata[[column_name]] <<- metadata[[column_name]]/val
  })
  mcols(sites) <- metadata
  sites
}

getPositionalValuesOfFeature <- function(sites, genomicData) {
  #### Boundary Distances #### Nirav Malani code TODO: refactor into several functions
  ## (refSeq boundary.dist), Start (refSeq start.dist), non-width (), General (general.width)
  ## when inGene is FALSE then set following: ref.left.pos, ref.right.pos, ref.left.strand, ref.right.strand
  ## when inGene is TRUE then set following: ref.start.pos, ref.end.pos, ref.gene.strand
  
  ## prepare the new columns ##
  colnam <- paste("ref", c("left.pos", "right.pos", "left.strand", "right.strand", 
                           "start.pos", "end.pos", "gene.strand"), sep=".") 
  mcols(sites)[colnam] <- NA
  
  ## add the respective columns as needed ##
  ## beware: precede returns range which is following the query and
  ## follow returns the range which is preceding the query!
  ## so do a switcheroo in terms of extracting the start & stop ##
  left <- follow(sites, genomicData, ignore.strand=TRUE)
  left[is.na(left) | sites$within_refSeq_gene] <- NA
  rows <- na.omit(left)
  sites$ref.left.pos[!is.na(left)] <- end(genomicData[rows])
  sites$ref.left.strand[!is.na(left)] <- as.character(strand(genomicData[rows]))
  
  right <- precede(sites, genomicData, ignore.strand=TRUE)
  right[is.na(right) | sites$within_refSeq_gene] <- NA
  rows <- na.omit(right)
  sites$ref.right.pos[!is.na(right)] <- start(genomicData[rows])
  sites$ref.right.strand[!is.na(right)] <- as.character(strand(genomicData[rows]))
  
  inIt <- findOverlaps(sites, genomicData, ignore.strand=TRUE, select="arbitrary")
  inIt[is.na(inIt) | !sites$within_refSeq_gene] <- NA
  rows <- na.omit(inIt)
  sites$ref.start.pos[!is.na(inIt)] <- start(genomicData[rows])
  sites$ref.end.pos[!is.na(inIt)] <- end(genomicData[rows])
  sites$ref.gene.strand[!is.na(inIt)] <- as.character(strand(genomicData[rows]))
  
  sites$boundary.dist <-
    eval(expression(pmin((ref.end.pos-position)/(ref.end.pos-ref.start.pos),
                         (position-ref.start.pos)/(ref.end.pos-ref.start.pos),
                         (ref.right.pos-position)/(ref.right.pos-ref.left.pos),
                         (position-ref.left.pos)/(ref.right.pos-ref.left.pos),
                         na.rm=T)), mcols(sites))
  
  sites$start.dist <-
    eval(expression(pmin(ifelse(ref.gene.strand=="-",
                                (ref.end.pos-position)/(ref.end.pos-ref.start.pos),
                                (position-ref.start.pos)/(ref.end.pos-ref.start.pos)),
                         ifelse(ref.right.strand=="-",
                                (ref.right.pos-position)/(ref.right.pos-ref.left.pos),
                                NA),
                         ifelse(ref.left.strand=="+",
                                (position-ref.left.pos)/(ref.right.pos-ref.left.pos),
                                NA),na.rm=T)), mcols(sites))
  
  sites$general.width <- eval(expression(pmin(ref.end.pos-ref.start.pos, 
                                              ref.right.pos-ref.left.pos,na.rm=T)),
                              mcols(sites))
  sites$gene.width <- eval(expression(ref.end.pos-ref.start.pos ), mcols(sites))
  
  meta <- mcols(sites)
  meta <- meta[ , ! (names(meta) %in% colnam)]
  mcols(sites) <- meta
  
  sites 
}

#' ROC.stata does not work with too few sites
is_enough_sites <- function(sampleName_GTSP, connection) {
     MIN_NUMBER_OF_SITES <- 3
     n_sites <- getUniqueSiteCounts(sampleName_GTSP, connection) 
     n_sites$enough_sites <- n_sites$uniqueSites >= MIN_NUMBER_OF_SITES
     if (all(n_sites$enough_sites)) {
        return(TRUE) 
     }
     message("****************************************")
     message("The following samples have too few sites to generate heatmap:")
     print(filter(n_sites, enough_sites == FALSE)) 
     message("****************************************")
     FALSE 
}

#' create a folder containing the SVG outputs and p-value calculations for the (clasically matched) random controls used in ROC calculations
#'
#' @param sites_mrcs Granges with sites, controls and features
sites_to_ROC_old <- function(sites_mrcs, output_dir) {
    sites_mrcs <- as.data.frame(sites_mrcs)
    annotation_columns <- get_annotation_columns(sites_mrcs)
    rset <- with(sites_mrcs, ROC.setup(
      rep(TRUE, nrow(sites_mrcs)), type, siteID, sampleName))
    roc.res <- ROC.strata(annotation_columns, rset, add.var=TRUE, sites_mrcs)
    ROCSVG(roc.res, output_dir)
}

##' Substitute Median for \code{NA}
##'
##' When there is little mising data, a rough-and-ready fill-in method
##' may be preferred to computationally intensive method of handling
##' missingness.  For non-parametric ROC curves based on ranks of the
##' data, using the median (of the non-missing data) as the fill-in is
##' a fairly innocuous choice.  If there is much missing data, this
##' method is not advised as it tends to bias the ROC curve area
##' towards 0.5.
##' @title na.median
##' @param x \code{matrix} or \code{data.frame} possibly containing
##'     \code{NA} values.
##' @return An object like \code{x}, but with the medians of the
##'     columns used in place of the \code{NA} values in the
##'     corresponding columns.
##' @author Charles Berry
na.median <-
  function(x)
  {
    if (!is.matrix(x)) x <- as.matrix(x)
    na.count <- colSums(is.na(x))
    if (any(na.count != 0 ))
    {
      for (i in 1:ncol(x))
        if (na.count[i]>0){
          med <- median(x[,i],na.rm=TRUE)
          x[is.na(x[,i]),i] <- med
        }
    }
    x
  }

# ROC calculation for Matched Random Controls 
# This is similar to the way we did this classically for restriction sites.
sites_to_ROC_matched <- function(sites_mrcs, sampleName_GTSP, output_dir) {
  sites_mrcs <- as.data.frame(sites_mrcs)
  annotation_columns <- get_annotation_columns(sites_mrcs)
    roc.res <- ROC.MRC(
    sites_mrcs[,"type"],
    sites_mrcs[,"siteID"],
    na.median(sites_mrcs[,annotation_columns]),
    sites_mrcs[,"sampleName"],
    origin.levels = unique(as.character(sampleName_GTSP$GTSP)))
  ROCSVG(roc.res, output_dir)
  saveRDS(roc.res, file = file.path(output_dir, "roc.res.rds"))
}


# ROC calculation for Ordinary (unmatched) Random Controls
# Provides the corrected version for unmatched controls.
# This should roughly match. p-values will differ but should not be
# qualitatively different.

sites_to_ROC_ordinary <- function(sites_mrcs, sampleName_GTSP, output_dir) {
  sites_mrcs <- as.data.frame(sites_mrcs)

  write.table(sites_mrcs, file='sites_mrcs.gen')

  annotation_columns <- get_annotation_columns(sites_mrcs)

  #i <- sites_mrcs[c("type", annotation_columns, "sampleName")] 
  #message('sites_mrcs data frame:')
  #write.table(i)

  roc.res <- ROC.ORC(
    sites_mrcs[,"type"],
    na.median(sites_mrcs[,annotation_columns]),
    sites_mrcs[,"sampleName"],
    origin.levels = unique(as.character(sampleName_GTSP$GTSP)))
  ROCSVG(roc.res, output_dir)
  saveRDS(roc.res, file = file.path(output_dir, "roc.res.rds"))
}
