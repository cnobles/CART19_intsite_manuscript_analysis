# supporting functions for gene ontology analysis setup

#' Cluster by similarity between gene_lists
#' @param keys character vector Keys corresponding to the names of the `keyList`
#' object. Keys will be the values returned clustered.
#' @param keyList list List of character vectors dictating the contents of each
#' key value. Names need to correspond to `keys` object.
#' @param limitValues character vector Values to consider when comparing 
#' keyLists for similarity. 
#' @param cores number of cores to use during distance matrix calculations. Uses
#' `parallel` R-package.
#' @param returnGraph returned object is changed from a numeric vector of 
#' groupings to a listed object with grouping vector and igraph object.
#' @param type character. Supported types of list similarity are currently 
#' Jaccard Index ("jaccard") and maximal similarity ("maxSim"). Where Jaccard 
#' Index is calculated by the intersection / union of two lists, the maximum
#' similarity is calculated by the intersection / minimum length of the two 
#' lists.
#' @param ... additional argument to pass to 
#' 
cluster_by_list_similarity <- function(keys, keyList, 
                                       limitValues = NULL, cores = 1, 
                                       returnGraph = FALSE, type = "jaccard",
                                       ...){
  packs <- c("igraph")
  packsLoaded <- suppressMessages(sapply(packs, require, character.only = TRUE))
  stopifnot(all(packsLoaded))
  if(cores > 1){
    stopifnot(require("parallel"))
    buster <- makeCluster(min(cores, detectCores()))
  }
  
  # Isolate working list
  keyList <- keyList[keys]
  
  if(!is.null(limitValues)){
    keyList <- lapply(keyList, function(x) x[x %in% limitValues])
  }
  
  # Determine similarity index of each key from each other
  if(type == "jaccard"){
    if(cores > 1){
      dist_mat <- do.call(cbind, parLapply(buster, keyList, function(x, kL){
        sapply(kL, function(y){
          length(which(x %in% y)) / length(unique(c(x, y)))})
      }, kL = keyList))
      stopCluster(buster)
    }else{
      dist_mat <- do.call(cbind, lapply(keyList, function(x, kL){
        sapply(kL, function(y){
          length(which(x %in% y)) / length(unique(c(x, y)))})
      }, kL = keyList))
    }
  }else if(type == "maxSim"){
    if(cores > 1){
      dist_mat <- do.call(cbind, parLapply(buster, keyList, function(x, kL){
        sapply(kL, function(y){
          length(which(x %in% y)) / min(length(unique(x)), length(unique(y)))})
      }, kL = keyList))
      stopCluster(buster)
    }else{
      dist_mat <- do.call(cbind, lapply(keyList, function(x, kL){
        sapply(kL, function(y){
          length(which(x %in% y)) / min(length(unique(x)), length(unique(y)))})
      }, kL = keyList))
    }
  }else{
    stop("Please set 'type' to either 'jaccard' or 'maxSim'.")
  }

  # Generate weighted (distance based) graph and cluster
  g <- graph_from_adjacency_matrix(
    adjmatrix = dist_mat, mode = "lower", weighted = TRUE, diag = FALSE)
  clus <- cluster_louvain(g)
  
  if(returnGraph){
    return(
      list("membership" = membership(clus), "graph" = g, "community" = clus))
  }else{
    return(membership(clus))
  }
}

# Load Gene Ontology GAF v2.1 files downloaded from UniProt
load_go_gaf <- function(file, ..., data.table = FALSE){
  stopifnot(require(data.table))
  # GAF column names
  gaf_cols <- c(
    "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", 
    "DB_Reference", "Evidence_Code", "With_From", "Aspect", "DB_Object_Name", 
    "DB_Object_Synonym", "DB_Object_Type", "Taxon_and_Interacting_taxon", 
    "Date", "Assigned_by", "Annotation_Extension", "Gene_Product_Form_ID")
  # GAF column classes
  gaf_classes <- rep("character", length(gaf_cols))
  # Read table
  df <- fread(file, ..., data.table = data.table, colClasses = gaf_classes)
  names(df) <- gaf_cols
  return(df)
}

#' Fisher non-central hypergeometric test for gene lists against specified Gene
#' Ontology terms.
#' @param gene_list character vector of genes to test for GO enrichment.
#' @param GO_list list of genes mapping to GO terms. Names of the list should be
#' the GO IDs while the contents of each object in the list should matche the 
#' input gene list.
#' @param N numeric, the size of all possible genes for the dataset.
#' @param odds logical, if TRUE (default), odds ratio will be calculated from 
#' gene_list size and input size or from pwd. If FALSE, odds ratio is set to 1,
#' effectively making the test a standard Fisher's Test.
#' @param pwd named numeric vector where names match gene_list terms and the 
#' numeric vector is the probability / weight of selecting the gene.
#' @param cutoff numeric, p-value cutoff to consider significant, 
#' default = 0.05.
#' @param p_adjust character, method for adjusting p-values to account for 
#' multiple comparisons. Default: "bonferroni", for others check 
#' p.adjust.methods.
#' @param filter logical, if TRUE (default), comparisons with zero overlapping
#' genes or above the cutoff value are removed from the output. If FALSE, they 
#' are left in the output data.frame.
#' @param overlap_genes logical, if TRUE, a column listing the overlapping genes
#' is returned. Default is FALSE.
#' @param cores integer If greater than 1, the function will require `parallel`,
#' `multidplyr`, and `foreach` packages and use the designated number of cores 
#' (or max cores) to compute the calculations.
#' @param lower_tail logical to be passed onto BiasedUrn::pFNCHypergeo() 
#' lower.tail argument.

fisher_hyper_GO_test <- function(gene_list, GO_list, odds = FALSE, 
                                 pwd = NULL, cutoff = 0.05, 
                                 p_adjust = "bonferroni", filter = FALSE, 
                                 overlap_genes = FALSE, lower_tail = TRUE,
                                 cores = 1){
  packs <- c("BiasedUrn", "dplyr", "magrittr", "GO.db")
  stopifnot(all(sapply(packs, require, character.only = TRUE)))
  if(cores > 1){
    parPacks <- c("parallel")
    stopifnot(all(sapply(parPacks, require, character.only = TRUE)))
    buster <- makeCluster(min(c(cores, detectCores())))
    clusterExport(
      cl = buster, 
      varlist = c("GO_list", "lower_tail"),
      envir = environment())
  }

  # Only assess intersecting genes
  gene_list <- gene_list[gene_list %in% unique(unlist(GO_list))]
  if(!is.null(pwd)){
    pwd <- pwd[names(pwd) %in% unique(unlist(GO_list))]
    stopifnot(all(gene_list %in% names(pwd)))
  }
  
  N <- length(unique(unlist(GO_list)))
  
  #Each test is independent to the GO Term
  go_df <- data.frame(GO_ID = names(GO_list)) %>%
    mutate(
      GO_size = lengths(GO_list),
      gene_list_size = length(gene_list),
      overlap_size = sapply(GO_list, function(x) sum(gene_list %in% x) ),
      overlap_genes = sapply(GO_list, function(x){
        paste(x[which(x %in% gene_list)], collapse = ";")}),
      non_GO_size = N - GO_size)
  
  if(is.null(pwd)){
    go_df <- dplyr::mutate(go_df, odds = GO_size / non_GO_size)
  }else{
    go_df <- dplyr::group_by(go_df, GO_ID) %>%
      dplyr::mutate(
        odds = sapply(GO_ID, function(id){
          sum(pwd[names(pwd) %in% GO_list[[id]]]) / 
            sum(pwd[!names(pwd) %in% GO_list[[id]]]) })) %>%
      dplyr::ungroup()
  }
  
  if(!odds) go_df <- dplyr::mutate(go_df, odds = 1)
  
  # Calculate p-values and adjust for each term
  if(cores > 1){
    go_df <- split(go_df, ceiling(seq_len(nrow(go_df))/(nrow(go_df)/cores)))
    
    go_df <- dplyr::bind_rows(parLapply(buster, go_df, function(df){
      library(magrittr)
      dplyr::group_by(df, GO_ID) %>%
        dplyr::mutate(
          p.value = BiasedUrn::pFNCHypergeo(
            x = overlap_size, m1 = GO_size, m2 = non_GO_size, 
            n = gene_list_size, odds = odds, lower.tail = lower_tail,
            precision = 1E-32)) %>%
        dplyr::ungroup() }))
    
  }else{
    go_df <- group_by(go_df, GO_ID) %>%
      mutate(
        p.value = BiasedUrn::pFNCHypergeo(
          x = overlap_size, m1 = GO_size, m2 = non_GO_size, n = gene_list_size, 
          odds = odds, lower.tail = lower_tail, precision = 1E-32)) %>%
      ungroup()
  }
  
  # Filter output to only significant adj.p.values
  go_df <- mutate(go_df, adj.p.value = p.adjust(p.value, method = p_adjust))
  
  if(filter){
    go_df <- dplyr::filter(go_df, adj.p.value <= cutoff, overlap_size > 0)
  }
  
  # Organize output and return
  if(cores > 1){
    go_df <- split(go_df, ceiling(seq_len(nrow(go_df))/(nrow(go_df)/cores)))
    
    go_df <- bind_rows(parLapply(buster, go_df, function(df){
      library(GO.db)
      dplyr::group_by(df, GO_ID) %>%
        dplyr::mutate(
          GO_Term = ifelse(
            adj.p.value <= cutoff, Term(GOTERM[[GO_ID]]), NA)) %>%
        dplyr::ungroup() }))
    
    stopCluster(buster)
  }else{
    go_df <- dplyr::group_by(go_df, GO_ID) %>%
      dplyr::mutate(
        GO_Term = ifelse(
          adj.p.value <= cutoff, Term(GOTERM[[GO_ID]]), NA)) %>%
      dplyr::ungroup()
  }
  
  go_df <- dplyr::select(
      go_df, GO_ID, GO_Term, GO_size, overlap_size, odds, 
      p.value, adj.p.value, overlap_genes) %>%
    arrange(adj.p.value, desc(overlap_size))
  
  if(!overlap_genes) go_df <- dplyr::select(go_df, -overlap_genes)
  
  go_df
}

fisher_hyper_KEGG_test <- function(gene_list, KEGG_list, odds = FALSE, 
                                   pwd = NULL, cutoff = 0.05, 
                                   p_adjust = "bonferroni", filter = FALSE, 
                                   overlap_genes = FALSE, lower_tail = TRUE,
                                   cores = 1){
  packs <- c("BiasedUrn", "dplyr", "magrittr", "KEGGREST")
  stopifnot(all(sapply(packs, require, character.only = TRUE)))
  kegg_path <- KEGGREST::keggList("pathway", "hsa")
  if(cores > 1){
    buster <- parallel::makeCluster(min(c(cores, parallel::detectCores())))
    parallel::clusterExport(
      cl = buster, 
      varlist = c("KEGG_list", "lower_tail", "kegg_path"),
      envir = environment())
  }
  
  # Only assess intersecting genes
  gene_list <- gene_list[gene_list %in% unique(unlist(KEGG_list))]
  if(!is.null(pwd)){
    pwd <- pwd[names(pwd) %in% unique(unlist(KEGG_list))]
    stopifnot(all(gene_list %in% names(pwd)))
  }
  
  N <- length(unique(unlist(KEGG_list)))
  
  # Each test is independent to the GO Term
  kegg_df <- data.frame(path_id = names(KEGG_list)) %>%
    dplyr::mutate(
      KEGG_size = lengths(KEGG_list),
      gene_list_size = length(gene_list),
      overlap_size = sapply(KEGG_list, function(x) sum(gene_list %in% x) ),
      overlap_genes = sapply(KEGG_list, function(x){
        paste(x[which(x %in% gene_list)], collapse = ";")}),
      non_KEGG_size = N - KEGG_size)
  
  if(is.null(pwd)){
    kegg_df <- dplyr::mutate(kegg_df, odds = KEGG_size / non_KEGG_size)
  }else{
    kegg_df <- dplyr::group_by(kegg_df, path_id) %>%
      mutate(
        odds = sapply(path_id, function(id){
          sum(pwd[names(pwd) %in% KEGG_list[[id]]]) / 
            sum(pwd[!names(pwd) %in% KEGG_list[[id]]]) })) %>%
      ungroup()
  }
  
  if(!odds) kegg_df <- dplyr::mutate(kegg_df, odds = 1)
  
  # Calculate p-values and adjust for each term
  if(cores > 1){
    kegg_df <- split(
      kegg_df, ceiling(seq_len(nrow(kegg_df))/(nrow(kegg_df)/cores)))
    
    kegg_df <- bind_rows(parLapply(buster, kegg_df, function(df){
      library(magrittr)
      dplyr::group_by(df, path_id) %>%
        dplyr::mutate(
        p.value = BiasedUrn::pFNCHypergeo(
          x = overlap_size, m1 = KEGG_size, m2 = non_KEGG_size, 
          n = gene_list_size, odds = odds, lower.tail = lower_tail,
          precision = 1E-32)) %>%
        dplyr::ungroup() }))
  }else{
    kegg_df <- dplyr::group_by(kegg_df, path_id) %>%
      dplyr::mutate(
        p.value = BiasedUrn::pFNCHypergeo(
          x = overlap_size, m1 = KEGG_size, m2 = non_KEGG_size, 
          n = gene_list_size, odds = odds, lower.tail = lower_tail,
          precision = 1E-32)) %>%
      dplyr::ungroup()
  }
  
  # Filter output to only significant adj.p.values
  kegg_df <- dplyr::mutate(
    kegg_df, adj.p.value = p.adjust(p.value, method = p_adjust))
  
  if(filter){
    kegg_df <- dplyr::filter(kegg_df, adj.p.value <= cutoff, overlap_size > 0)
  }
  
  # Organize output and return
  if(cores > 1){
    kegg_df <- split(
      kegg_df, ceiling(seq_len(nrow(kegg_df))/(nrow(kegg_df)/cores)))
    kegg_df <- bind_rows(parLapply(buster, kegg_df, function(df){
      library(magrittr)
      dplyr::group_by(df, path_id) %>%
        dplyr::mutate(KEGG_Term = kegg_path[path_id]) %>%
        dplyr::ungroup() }))
    stopCluster(buster)
  }else{
    kegg_df <- dplyr::group_by(kegg_df, path_id) %>%
      dplyr::mutate(KEGG_Term = kegg_path[path_id]) %>%
      dplyr::ungroup()
  }
      
  kegg_df <- dplyr::select(
    kegg_df, path_id, KEGG_Term, KEGG_size, overlap_size, odds, 
    p.value, adj.p.value, overlap_genes) %>%
    arrange(desc(overlap_size), adj.p.value)
  
  if(!overlap_genes) kegg_df <- dplyr::select(kegg_df, -overlap_genes)
  
  kegg_df
}

fisher_hyper_LIST_test <- function(gene_list, ref_list, odds = FALSE, 
                                 pwd = NULL, cutoff = 0.05, 
                                 p_adjust = "bonferroni", filter = FALSE, 
                                 overlap_genes = FALSE, lower_tail = TRUE,
                                 cores = 1){
  packs <- c("BiasedUrn", "dplyr", "magrittr")
  stopifnot(all(sapply(packs, require, character.only = TRUE)))
  if(cores > 1){
    parPacks <- c("parallel")
    stopifnot(all(sapply(parPacks, require, character.only = TRUE)))
    buster <- makeCluster(min(c(cores, detectCores())))
    clusterExport(
      cl = buster, 
      varlist = c("ref_list", "lower_tail"),
      envir = environment())
  }
  
  # Only assess intersecting genes
  gene_list <- gene_list[gene_list %in% unique(unlist(ref_list))]
  if(!is.null(pwd)){
    pwd <- pwd[names(pwd) %in% unique(unlist(ref_list))]
    stopifnot(all(gene_list %in% names(pwd)))
  }
  
  N <- length(unique(unlist(ref_list)))
  
  #Each test is independent to the Ref. Term
  ref_df <- data.frame(ref_id = names(ref_list)) %>%
    mutate(
      ref_size = lengths(ref_list),
      gene_list_size = length(gene_list),
      overlap_size = sapply(ref_list, function(x) sum(gene_list %in% x) ),
      overlap_genes = sapply(ref_list, function(x){
        paste(x[which(x %in% gene_list)], collapse = ";")}),
      non_ref_size = N - ref_size)
  
  if(is.null(pwd)){
    ref_df <- dplyr::mutate(ref_df, odds = ref_size / non_ref_size)
  }else{
    ref_df <- dplyr::group_by(ref_df, ref_id) %>%
      dplyr::mutate(
        odds = sapply(ref_id, function(id){
          sum(pwd[names(pwd) %in% ref_list[[id]]]) / 
            sum(pwd[!names(pwd) %in% ref_list[[id]]]) })) %>%
      dplyr::ungroup()
  }
  
  if(!odds) ref_df <- dplyr::mutate(ref_df, odds = 1)
  
  # Calculate p-values and adjust for each term
  if(cores > 1){
    ref_df <- split(ref_df, ceiling(seq_len(nrow(ref_df))/(nrow(ref_df)/cores)))
    
    ref_df <- dplyr::bind_rows(parLapply(buster, ref_df, function(df){
      library(magrittr)
      dplyr::group_by(df, ref_id) %>%
        dplyr::mutate(
          p.value = BiasedUrn::pFNCHypergeo(
            x = overlap_size, m1 = ref_size, m2 = non_ref_size, 
            n = gene_list_size, odds = odds, lower.tail = lower_tail,
            precision = 1E-32)) %>%
        dplyr::ungroup() }))
    
  }else{
    ref_df <- group_by(ref_df, ref_id) %>%
      mutate(
        p.value = BiasedUrn::pFNCHypergeo(
          x = overlap_size, m1 = ref_size, m2 = non_ref_size, n = gene_list_size, 
          odds = odds, lower.tail = lower_tail, precision = 1E-32)) %>%
      ungroup()
  }
  
  # Filter output to only significant adj.p.values
  ref_df <- mutate(ref_df, adj.p.value = p.adjust(p.value, method = p_adjust))
  
  if(filter){
    ref_df <- dplyr::filter(ref_df, adj.p.value <= cutoff, overlap_size > 0)
  }
  
  ref_df <- dplyr::select(
      ref_df, ref_id, ref_size, overlap_size, odds, p.value, adj.p.value, 
      overlap_genes) %>%
    arrange(adj.p.value, desc(overlap_size))
  
  if(!overlap_genes) ref_df <- dplyr::select(ref_df, -overlap_genes)
  
  ref_df
}


gene_rep_binomial_test <- function(nc, N, k, K){
  Pc <- nc / N
  rep <- ifelse(k >= ( Pc * K ), "+", "-")
  if(rep == "+"){
    p <- sum(exp(dbinom(k:K, K, Pc, log = TRUE)))
  }else{
    p <- sum(exp(dbinom(0:k, K, Pc, log = TRUE)))  
  }
  return(list(p.value = p, over_under = rep))
}

# Binomial test to estimate the probability of observing two lists.
binomial_compare_test <- function(gene_list, GO_list, N, cutoff = 0.05, 
                                  filter = TRUE, p_adjust = "fdr", 
                                  overlapping_genes = FALSE, ...){
  packs <- c("dplyr", "magrittr", "stats")
  stopifnot(sapply(packs, require, character.only = TRUE))
  #Each test is independent to the GO Term
  go_df <- data.frame(GO_ID = names(GO_list)) %>%
    dplyr::mutate(
      nc = sapply(GO_list, length),
      Pc = nc / N,
      exp_k = Pc * length(gene_list),
      k = sapply(GO_list, function(x){
        length(which(gene_list %in% x))}),
      overlapping_genes = sapply(GO_list, function(x){
        paste(x[which(x %in% gene_list)], collapse = ";")}))
  
  # Calculate p-values and adjust for each term
  sig_data <- lapply(1:nrow(go_df), function(i){
    gene_rep_binomial_test(
      nc = go_df$nc[i], N = N, k = go_df$k[i], K = length(gene_list))})
  go_df$rep <- sapply(sig_data, function(x) x$over_under)
  go_df$p.value <- sapply(sig_data, function(x) x$p.value)
  go_df$adj.p.value <- p.adjust(go_df$p.value, method = p_adjust)
  
  # Filter output to only significant adj.p.values
  if(filter){go_df <- dplyr::filter(go_df, adj.p.value <= cutoff)}
  
  # Organize output and return
  go_df <- dplyr::select(
      go_df, GO_ID, nc, exp_k, k, rep, 
      p.value, adj.p.value, overlapping_genes) %>%
    rename(
      "GO_size" = nc, "Exp_overlap" = exp_k, "Overlap_size" = k, 
      "Over_Under_Rep" = rep, "Overlap_genes" = overlapping_genes)
  
  if(!overlapping_genes) go_df <- dplyr::select(go_df, -Overlap_genes)
  return(go_df)
}

# Create a wordcloud from GO terms given a set of GO_IDs
# @param GO_IDs character vector of GO identifiers.
# @param GO_terms list of GO terms, as generated by as.list(GO.db::GOTERM). If
# left `NULL`, function will generate on the fly, but this will significantly 
# increase the time the function takes to execute.
# @param ... options to be passed into `wordcloud()`.
plot_GO_cloud <- function(GO_IDs, GO_terms = NULL, ...){
  packs <- c("wordcloud", "GO.db", "tm")
  stopifnot(sapply(packs, require, character.only = TRUE))
  ignored_words <- c("process", "function", "component")
  
  if(is.null(GO_terms)){
    GO_terms <- as.list(GO.db::GOTERM)
  }
  
  # Collect GO_terms
  terms <- sapply(GO_IDs, function(id) GO_terms[[id]]@Term)
  
  # Filter specific words
  terms <- tolower(terms)
  terms <- gsub("t cell", "t-cell", terms)
  terms <- gsub("b cell", "b-cell", terms)
  term_freq <- table(unlist(strsplit(terms, " ")))
  term_freq <- term_freq[!names(term_freq) %in% stopwords()]
  term_freq <- term_freq[!names(term_freq) %in% ignored_words]
  
  # Construct wordcloud
  wordcloud(words = names(term_freq), freq = term_freq, ...)
}