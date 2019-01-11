#### supporting functions for the cart19_intsite_analysis.Rmd scripts
# Functions
scan_sites_count <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, sites))
    length(unique(hits))
  })
}

scan_sites_abundance <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, sites))
    sum(sites[hits]$estAbund)
  })
}

scan_patients <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, sites))
    length(unique(sites[hits]$patient))
  })
}

scan_genes <- function(clus_gr, refGenes){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    hits <- subjectHits(findOverlaps(clus, refGenes))
    paste(unique(refGenes[hits]$name2), collapse = ", ")
  })
}

scan_orientation <- function(clus_gr, sites){
  sapply(1:length(clus_gr), function(i){
    clus <- clus_gr[i]
    site_hits <- findOverlaps(clus, sites)
    site_ort <- strand(sites[subjectHits(site_hits)])
    gene_ort <- sites[subjectHits(site_hits)]$in_geneOrt
    if(length(sites[subjectHits(site_hits)]) > 0){
      df <- data.frame(
        "posid" = generate_posid(sites[subjectHits(site_hits)]),
        "site_ort" = as.character(site_ort),
        "gene_ort" = as.character(gene_ort),
        "patient" = sites[subjectHits(site_hits)]$patient,
        stringsAsFactors = FALSE
      )
      df <- distinct(df)
      df$same_orientation <- df$site_ort == df$gene_ort
      score <- paste0(
        "T", length(grep("TRUE", df$same_orientation)), ":",
        "F", length(grep("FALSE", df$same_orientation)), ":",
        "N", length(grep("TRUE", is.na(df$same_orientation)))
      )
    }else{
      score <- "T0:F0:N0"
    }
    score
  })
}

scan_fisher_test_ort <- function(score1, score2){
  sapply(1:length(score1), function(i){
    grp1 <- unlist(strsplit(score1[i], ":"))
    grp2 <- unlist(strsplit(score2[i], ":"))
    x <- matrix(c(
      as.integer(substr(grp1[1], 2, 3)),
      as.integer(substr(grp1[2], 2, 3)),
      as.integer(substr(grp2[1], 2, 3)),
      as.integer(substr(grp2[2], 2, 3))),
      ncol = 2
    )
    fisher.test(x)$p.value
  })
}

annotate_scan_clusters <- function(scanned_ranges, tdn_sites, 
                                   timepoint_sites, refGenes){
  scanned_ranges$n_sites_tdn <- scan_sites_count(
    scanned_ranges, tdn_sites)
  scanned_ranges$n_sites_tp <- scan_sites_count(
    scanned_ranges, timepoint_sites)
  scanned_ranges$sum_abund_tdn <- scan_sites_abundance(
    scanned_ranges, tdn_sites)
  scanned_ranges$sum_abund_tp <- scan_sites_abundance(
    scanned_ranges, timepoint_sites)
  scanned_ranges$n_patients_tdn <- scan_patients(
    scanned_ranges, tdn_sites)
  scanned_ranges$n_patients_tp <- scan_patients(
    scanned_ranges, timepoint_sites)
  scanned_ranges$genes_in_cluster <- scan_genes(
    scanned_ranges, refGenes)
  scanned_ranges$in_gene_ort_tdn <- scan_orientation(
    scanned_ranges, tdn_sites)
  scanned_ranges$in_gene_ort_tp <- scan_orientation(
    scanned_ranges, timepoint_sites)
  scanned_ranges$ort_fisher_test <- scan_fisher_test_ort(
    scanned_ranges$in_gene_ort_tdn, scanned_ranges$in_gene_ort_tp)
  scanned_ranges
}

get_top_sites <- function(sites, percent, rank_by){
  sites <- sites[order(mcols(sites)[,rank_by], decreasing = TRUE)]
  first_sites <- sites[!duplicated(sites$posid)]
  ranks <- rank(-(mcols(first_sites)[,rank_by]), ties.method = "min")
  cutoff <- round(length(first_sites)*percent/100)
  if(cutoff == 0){ cutoff <- 1 }
  first_sites[ranks <= cutoff]$posid
}

format_sites <- function(sites){
  df <- GenomicRanges::as.data.frame(sites, row.names = NULL) %>%
    mutate(geneid = gene_id_wo_annot) %>%
    dplyr::select(
      seqnames, start, strand, patient, timepoint, 
      celltype, posid, estAbund, relAbund) %>%
    as.data.frame()
  ranges <- IRanges(start = df$start, width = rep(1, nrow(df)))
  gr <- GRanges(
    seqnames = df$seqnames,
    ranges = ranges,
    strand = df$strand,
    seqlengths = seqlengths(sites),
    seqinfo = seqinfo(sites)
  )
  mcols(gr) <- dplyr::select(df, -seqnames, -start, -strand)
  unique_granges(gr)
}

get_top_genes <- function(sites, genes, percent, rank_by){
  sites <- sites[order(mcols(sites)[,rank_by], decreasing = TRUE)]
  first_sites <- sites[!duplicated(sites$posid)]
  ranks <- rank(-(mcols(first_sites)[,rank_by]), ties.method = "min")
  cutoff <- round(length(first_sites)*percent/100)
  if(cutoff == 0){ cutoff <- 1 }
  posids <- first_sites[ranks <= cutoff]$posid
  top_sites <- sites[sites$posid %in% posids]
  unique(genes[subjectHits(findOverlaps(
    top_sites, genes, ignore.strand = TRUE))]$name2)
}

compare_format <- function(sites, clusters){
  hits <- findOverlaps(clusters, sites, ignore.strand = TRUE)
  sites$in_cluster <- seq_len(length(sites)) %in% subjectHits(hits)
  df <- GenomicRanges::as.data.frame(sites, row.names = NULL) %>%
    mutate(geneid = gene_id_wo_annot) %>%
    dplyr::select(patient, posid, strand, in_gene, in_geneOrt, nearest_gene, nearest_geneOrt, geneid, in_cluster) %>%
    group_by(patient) %>%
    distinct(posid, strand, in_gene, in_geneOrt, nearest_gene, nearest_geneOrt, geneid, in_cluster) %>% 
    ungroup() %>% as.data.frame()
  df
}

same_ort <- function(intSiteStrand, in_geneOrt){
  mapply(
    function(strand, in_gene_ort){
      orts <- unlist(strsplit(in_gene_ort, ","))
      strand %in% orts
    },
    as.character(intSiteStrand),
    as.character(in_geneOrt)
  )}

get_loci_id <- function(gr){
  paste0(seqnames(gr), strand(gr), start(gr), ":", end(gr))
}

# Generate random sites across the reference genome.
selectRandomSites <- function(num, refGenome, drop_extra_seqs = TRUE, 
                              seqNames = NULL, setSeed = NULL){
  require(GenomicRanges)
  if(!is.null(setSeed)) set.seed(setSeed)
  if(is.null(seqinfo(refGenome))) stop("Ref genome does not have seqinfo.")
  
  if(!is.null(seqNames)){
    seqLengths <- refGenome@seqinfo@seqlengths[
      match(seqNames, seqnames(refGenome))]
    names(seqLengths) <- seqNames
  }else{
    if(drop_extra_seqs){
      seqNames <- grep("chr[0-9XY]+$", seqnames(refGenome), value = TRUE)
      seqLengths <- refGenome@seqinfo@seqlengths[
        match(seqNames, seqnames(refGenome))]
      names(seqLengths) <- seqNames
    }else{
      seqNames <- seqnames(refGenome)
      seqLengths <- refGenome@seqinfo@seqlengths
      names(seqLengths) <- seqNames
    }
  }
  
  chrs <- sort(factor(
    sample(seqNames, num, replace = TRUE, prob = seqLengths),
    levels = seqNames))
  strands <- sample(c("+", "-"), num, replace = TRUE)
  spltChrs <- split(chrs, chrs)
  spltChrs <- spltChrs[sapply(spltChrs, length) > 0]
  positions <- unname(unlist(sapply(
    spltChrs, function(seq, seqLen){
      sample.int(
        n = seqLen[unique(as.character(seq))], size = length(seq), 
        replace = FALSE, useHash = FALSE)
    }, seqLen = seqLengths)))
  
  GRanges(
    seqnames = chrs,
    ranges = IRanges(start = positions, width = 1L),
    strand = strands,
    seqinfo = seqinfo(refGenome))
}

#' Identify intersection between multiple vectors
#' 
#' @usage serial_intersect(..., limit = NULL)
#' 
#' @param ... a series of vectors to compare against eachother.
#' @param limit integer The number of minimum observations of a single value 
#' before it is included in the output. Defaults to the number of input vectors,
#' but can be reduce to increase the number of returned values.
#' 
#' @author Christopher Nobles, Ph.D.
#' 

serial_intersect <- function(..., limit = NULL){
  v <- list(...)
  stopifnot(all(sapply(v, is.vector)))
  if(is.null(limit)) limit <- length(v)
  stopifnot(limit <= length(v))
  cnt <- table(unlist(v))
  cnt <- cnt[cnt >= limit]
  return(names(cnt))
}

#' Convert formats of [dmy]## into numeric days
#' 
#' @usage convert_time(timecode)
#' 
#' @param timecode character vector in the format [dmy]## where ## can be any 
#' number.
#' 
#' @details Output is a numeric vector in the same order and length as the input
#' with values indicating the conversion to days from the input format.
#' 
#' @author Christopher Nobles, Ph.D.
#' 
convert_time <- function(timecode){
  x <- tolower(timecode)
  scale <- ifelse(
    grepl("d", x), 1, 
    ifelse(grepl("m", x), 2, 
           ifelse(grepl("y", x), 3, 0)))
  timepoint <- as.numeric(stringr::str_extract(x, "[0-9.]+"))
  if(any(scale == 0)){
    stop("Not all timecodes in acceptable format.")
  }
  time_convert <- c("d" = 1, "m" = 30, "y" = 365)
  unname(timepoint * time_convert[scale])
}
