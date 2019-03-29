# CART19 Integration Site Manuscript Analysis

[![DOI](https://zenodo.org/badge/164499573.svg)](https://zenodo.org/badge/latestdoi/164499573)

Analysis provided in the "Linking outcome of CD19-directed CAR T cell therapy with genome modification by vector integration" manuscript.

This repository is meant to provide the source code for the analysis provided in the above manuscript by Christopher L. Nobles, Scott Sherrill-Mix, John Everett, Shantan Reddy, Joseph A. Fraietta, David Porter, Noelle Frey, Saar Gill, Steven Grupp, Shannon Maude, Donald Siegel, Bruce Levine, Carl H. June, Simon F. Lacey, J. Joseph Melenhorst, and Frederic D. Bushman.

Abstract:
    Chimeric antigen receptor-engineered T-cells targeting CD19 (CART19) provide a highly effective treatment for pediatric acute lymphoblastic leukemia (ALL), but are less effective for chronic lymphocytic leukemia (CLL) and other indications, focusing attention on approaches to improve efficacy.  CART19 cells are prepared by integration of the engineered receptor gene into the host T cell chromosome using a lentiviral vector.  Vector integration thus marks T cell lineages uniquely, and modifies the cellular genome by insertional mutagenesis.  Previously we reported that vector integration in the host gene TET2 resulted in gene inactivation and consequently therapeutic expansion of a cell clone associated with remission of CLL.  Here we investigate clonal population structure and therapeutic outcome for an additional 39 subjects by high throughput sequencing of sites of vector integration.  Expansion of cell clones with integration sites in specific genes in responders suggests possible insertional mutagenesis promoting therapeutic proliferation. XXX These data thus establish that a signal is present in integration site data associated with clinical response, and provide tools to help optimize T cell engineering.

To run the analysis, run `Rscript initiate_analysis.R`, which will subsequenty process the input data, populate the directory with informational files, and generate a report with the manuscript figures. Additionally, the user can specify a number of cores to use for parallelization with a `-c` option. i.e. `Rscript initiate_analysis.R -c 20`

This repository also contains a final INSPIIRED input dataset from the manuscript. Anyone can use this data rather than reprocessing with the INSPIIRED software from the SRA submission (BioProject: PRJNA510570). These data were processed with the [intSiteCaller v1.2.0 software](https://github.com/BushmanLab/intSiteCaller), within the [INSPIIRED suite](https://github.com/BushmanLab/INSPIIRED).

This source code was executed under the following session information:

```
R version 3.4.0 (2017-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.5 LTS

Matrix products: default
BLAS: /home/opt/R-3.4.0/lib/R/lib/libRblas.so
LAPACK: /home/opt/R-3.4.0/lib/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
 [1] splines   stats4    grid      parallel  methods   stats     graphics
 [8] grDevices utils     datasets  base

other attached packages:
 [1] caret_6.0-80                glmnet_2.0-13
 [3] doMC_1.3.5                  iterators_1.0.9
 [5] kableExtra_0.8.0            dnet_1.1.2
 [7] supraHex_1.14.0             hexbin_1.27.2
 [9] sna_2.4                     statnet.common_4.0.0
[11] network_1.13.0              intergraph_2.0-2
[13] BSgenome_1.44.2             rtracklayer_1.36.6
[15] GCcontent_0.0.0.9000        intSiteRetriever_0.0.0.9000
[17] hotROCs_0.99.5              bindrcpp_0.2.2
[19] spraphal_0.1.0              ggnet_0.1.0
[21] gintools_0.1.2              GO.db_3.4.1
[23] AnnotationDbi_1.38.2        Biobase_2.36.2
[25] KEGGREST_1.16.1             sonicLength_1.4.4
[27] hiAnnotator_1.11.1          geneRxCluster_1.12.0
[29] GenomicRanges_1.28.6        GenomeInfoDb_1.12.3
[31] Biostrings_2.44.2           XVector_0.16.0
[33] IRanges_2.10.5              S4Vectors_0.14.7
[35] BiocGenerics_0.22.1         forcats_0.3.0
[37] stringr_1.3.1               dplyr_0.7.6
[39] purrr_0.2.5                 readr_1.1.1
[41] tidyr_0.8.1                 tibble_1.4.2
[43] tidyverse_1.2.1             UpSetR_1.3.3
[45] foreach_1.4.4               BiasedUrn_1.07
[47] pander_0.6.1                magrittr_1.5
[49] knitr_1.19                  RColorBrewer_1.1-2
[51] gridExtra_2.3               scales_0.5.0
[53] ggrepel_0.7.0               ggplot2_3.0.0
[55] lubridate_1.7.4             vegan_2.4-6
[57] lattice_0.20-35             permute_0.9-4
[59] data.table_1.10.4-3         RMySQL_0.10.14
[61] DBI_0.8                     reshape2_1.4.3
[63] reldist_1.6-6               Matrix_1.2-9
[65] igraph_1.2.2

loaded via a namespace (and not attached):
  [1] tidyselect_0.2.3                  RSQLite_2.0
  [3] htmlwidgets_1.0                   BiocParallel_1.10.1
  [5] munsell_0.4.3                     codetools_0.2-15
  [7] withr_2.1.1                       colorspace_1.3-2
  [9] rstudioapi_0.7                    geometry_0.3-6
 [11] robustbase_0.93-3                 dimRed_0.1.0
 [13] labeling_0.3                      GenomeInfoDbData_0.99.0
 [15] bit64_0.9-7                       rprojroot_1.3-2
 [17] ipred_0.9-7                       R6_2.2.2
 [19] DRR_0.0.3                         bitops_1.0-6
 [21] DelayedArray_0.2.7                assertthat_0.2.0
 [23] nnet_7.3-12                       gtable_0.2.0
 [25] ddalpha_1.3.4                     timeDate_3043.102
 [27] rlang_0.2.2                       CVST_0.2-2
 [29] RcppRoll_0.3.0                    lazyeval_0.2.1
 [31] ModelMetrics_1.2.0                acepack_1.4.1
 [33] broom_0.5.0                       checkmate_1.8.5
 [35] yaml_2.1.16                       abind_1.4-5
 [37] modelr_0.1.2                      backports_1.1.2
 [39] Hmisc_4.1-1                       tools_3.4.0
 [41] lava_1.6.3                        Rcpp_0.12.17
 [43] plyr_1.8.4                        base64enc_0.1-3
 [45] zlibbioc_1.22.0                   RCurl_1.95-4.10
 [47] rpart_4.1-11                      sfsmisc_1.1-2
 [49] SummarizedExperiment_1.6.5        haven_1.1.2
 [51] cluster_2.0.6                     RSVGTipsDevice_1.0-7
 [53] matrixStats_0.53.0                hms_0.4.2
 [55] evaluate_0.10.1                   XML_3.98-1.9
 [57] readxl_1.1.0                      compiler_3.4.0
 [59] crayon_1.3.4                      htmltools_0.3.6
 [61] mgcv_1.8-17                       Formula_1.2-2
 [63] magic_1.5-9                       MASS_7.3-47
 [65] cli_1.0.1                         bindr_0.1.1
 [67] gower_0.1.2                       pkgconfig_2.0.1
 [69] GenomicAlignments_1.12.2          foreign_0.8-67
 [71] recipes_0.1.3                     xml2_1.2.0
 [73] prodlim_2018.04.18                rvest_0.3.2
 [75] BSgenome.Hsapiens.UCSC.hg38_1.4.1 digest_0.6.15
 [77] pls_2.7-0                         graph_1.54.0
 [79] rmarkdown_1.8                     cellranger_1.1.0
 [81] htmlTable_1.11.2                  curl_3.1
 [83] kernlab_0.9-25                    Rsamtools_1.28.0
 [85] nlme_3.1-131                      jsonlite_1.5
 [87] viridisLite_0.2.0                 pillar_1.1.0
 [89] httr_1.3.1                        DEoptimR_1.0-8
 [91] survival_2.41-3                   glue_1.2.0
 [93] png_0.1-7                         bit_1.1-12
 [95] Rgraphviz_2.20.0                  class_7.3-14
 [97] stringi_1.2.4                     blob_1.1.0
 [99] latticeExtra_0.6-28               memoise_1.1.0
[101] ape_5.0
```
