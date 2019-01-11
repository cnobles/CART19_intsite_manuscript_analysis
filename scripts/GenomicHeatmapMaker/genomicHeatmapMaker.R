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

packs <- c("colorspace", "hiAnnotator", "intSiteRetriever", "GCcontent", "BSgenome")
null <- suppressMessages(sapply(packs, require, character.only = TRUE)) 

make_heatmap <- function(sampleName_GTSP, referenceGenome, output_dir, connection) {
    if(class(connection) == "character"){
      sites_mrcs <- get_sites_controls_from_db(
          sampleName_GTSP, referenceGenome, connection)
    }else{
      sites_mrcs <- get_sites_controls_from_file(
          sampleName_GTSP, referenceGenome, connection)
    }
    sites_to_heatmap(sites_mrcs, referenceGenome, sampleName_GTSP, output_dir)
}

#' generate heatmap
#'
#' @param sites_mrcs real sites and mrcs GRanges object with metacolumns: 
#'      sitID, type, label
#' @param referenceGenome reference genome name, e.g. "hg18"
#' @param sampleName_GTSP input sample table, preserves order.
#' @param output_dir name of the directory
sites_to_heatmap <- function(sites_mrcs, referenceGenome, 
                             sampleName_GTSP, output_dir) {
    reference_genome_sequence <- get_reference_genome(referenceGenome)
    # TODO: populate from local database, at present pulled from UCSC web-site
    refSeq_genes <- getRefSeq_genes(referenceGenome)
    CpG_islands <- getCpG_islands(referenceGenome)

    ### oncogene_file <- "allonco_no_pipes.csv"
    oncogene_file <- file.path(codeDir, "allonco_no_pipes.csv")

    if (grepl("^mm", referenceGenome)) {
        ### oncogene_file <- "allonco_no_pipes.mm.csv"
        oncogene_file <- file.path(codeDir, "allonco_no_pipes.mm.csv")
    }

    # @return vector of gene symbols
    get_oncogene_from_file <- function(filename) {
        onco <- read.csv(filename, header=FALSE, stringsAsFactors=FALSE)
        as.character(onco$V1)
    }
    oncogenes <- get_oncogene_from_file(oncogene_file)
    
    if(referenceGenome == "hg38"){
      stability_regions <- readRDS(
        file.path(codeDir, "fragileRegions.hg38.rds"))
      fragile_regions <- stability_regions[
        stability_regions$regionType == "aCFR"]
      non_fragile_regions <- stability_regions[
        stability_regions$regionType == "NFR"]
      
      window_size_fragReg <- c("1M"=1e6)
      
      sites_mrcs <- getFeatureCounts(
        sites_mrcs, non_fragile_regions, "NFR", width = window_size_fragReg)
      sites_mrcs <- getFeatureCounts(
        sites_mrcs, fragile_regions, "aCFR", width = window_size_fragReg)
    }
    
    # is there oncogene closer than 50k
    refSeq_gene_symbols <- refSeq_genes$name2
    #' check if gene is onco gene list(curated by Bushman's lab)
    #' @return TRUE if onco-gene FALSE if not
    is_onco_gene <- function(gene_symbol_sites, oncogenes) {
        toupper(gene_symbol_sites) %in% toupper(oncogenes)
    }
    is_refSeq_oncogene <- is_onco_gene(refSeq_gene_symbols, oncogenes)
    refSeq_oncogene <- refSeq_genes[is_refSeq_oncogene]
    sites_mrcs <- getNearestFeature(
      sites_mrcs, refSeq_oncogene, dists.only=TRUE, colnam="onco")
    #sites_mrcs, refSeq_oncogene, dists.only=TRUE, colnam="onco.100k")
    sites_mrcs$onco.100k <- abs(sites_mrcs$oncoDist) <= 50000
    sites_mrcs$oncoDist <- NULL
    # end oncogene

    # need within_refSeq_gene for getPositionalValuesOfFeature
    sites_mrcs <- getSitesInFeature(
      sites_mrcs, refSeq_genes, "within_refSeq_gene", asBool=TRUE)

    sites_mrcs <- getPositionalValuesOfFeature(sites_mrcs, refSeq_genes)
   
    # redo the work to move it after positional values 
    mcols(sites_mrcs)$within_refSeq_gene <- NULL
    sites_mrcs <- getSitesInFeature(
      sites_mrcs, refSeq_genes, "within_refSeq_gene", asBool=TRUE)

    window_size_refSeq <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
    sites_mrcs <- getFeatureCounts(sites_mrcs, refSeq_genes, "refSeq_counts", 
                                   width=window_size_refSeq)

    window_size_GC <- c("100"=100, 
        "1k"=1000, "10k"=1e4, "100k"=1e5, "1M"=1e6)
    sites_mrcs <- getGCpercentage(
      sites_mrcs, "GC", window_size_GC, reference_genome_sequence)

    window_size_CpG_counts <- c("1k"=1e3, "10k"=1e4)
    sites_mrcs <- getFeatureCounts(sites_mrcs, CpG_islands, "CpG_counts", 
                                   width=window_size_CpG_counts)

    window_size_CpG_density <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
    sites_mrcs <- getFeatureCounts(sites_mrcs, CpG_islands, "CpG_density", 
                                   width=window_size_CpG_density)
    sites_mrcs <- from_counts_to_density(sites_mrcs, 
                                         "CpG_density", window_size_CpG_density)

    if( ! grepl("^mm", referenceGenome)) { # mouse does not have DNaseI track
        DNaseI <- getDNaseI(referenceGenome)
        window_size_DNaseI <- c("1k"=1e3, "10k"=1e4, "100k"=1e5, "1M"=1e6)
        sites_mrcs <- getFeatureCounts(sites_mrcs, DNaseI, "DNaseI_count", 
                                       width=window_size_DNaseI)
    }
    
    if (config$rocControls == 'unmatched')
    {
       sites_to_ROC_ordinary(sites_mrcs, sampleName_GTSP, output_dir)
    } else if (config$rocControls == 'matched') {
       sites_to_ROC_matched(sites_mrcs, sampleName_GTSP, output_dir)
    } else {
       stop('Error, rocControls is not properly defined in the configuration file.') 
    }
}
