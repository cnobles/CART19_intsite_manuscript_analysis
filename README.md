# CART19 Integration Site Manuscript Analysis

Analysis provided in the "Linking outcome of CD19-directed CAR T cell therapy with genome modification by vector integration" manuscript.

This repository is meant to provide the source code for the analysis provided in the above manuscript by Christopher L. Nobles, Scott Sherrill-Mix, John Everett, Shantan Reddy, Joseph A. Fraietta, David Porter, Noelle Frey, Saar Gill, Steven Grupp, Shannon Maude, Donald Siegel, Bruce Levine, Carl H. June, Simon F. Lacey, J. Joseph Melenhorst, and Frederic D. Bushman.

Abstract:
    Chimeric antigen receptor-engineered T-cells targeting CD19 (CART19) provide a highly effective treatment for pediatric acute lymphoblastic leukemia (ALL), but are less effective for chronic lymphocytic leukemia (CLL) and other indications, focusing attention on approaches to improve efficacy.  CART19 cells are prepared by integration of the engineered receptor gene into the host T cell chromosome using a lentiviral vector.  Vector integration thus marks T cell lineages uniquely, and modifies the cellular genome by insertional mutagenesis.  Previously we reported that vector integration in the host gene TET2 resulted in gene inactivation and consequently therapeutic expansion of a cell clone associated with remission of CLL.  Here we investigate clonal population structure and therapeutic outcome for an additional 39 subjects by high throughput sequencing of sites of vector integration.  Expansion of cell clones with integration sites in specific genes in responders suggests possible insertional mutagenesis promoting therapeutic proliferation. XXX These data thus establish that a signal is present in integration site data associated with clinical response, and provide tools to help optimize T cell engineering.

To run the analysis, run the R-script 'XXX.R', which will subsequenty process the input data, populate the directory with informational files, and generate a report with the manuscript figures. Additionally, the user can run the 'generate_goi_report.R' script after the intial script to generate the supplemental report included in the manuscript.

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
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] splines   stats4    grid      parallel  stats     graphics  grDevices utils     datasets  methods  
[11] base     

other attached packages:
 [1] spraphal_0.1.0       ggnet_0.1.0          gintools_0.1.1       GO.db_3.4.1         
 [5] AnnotationDbi_1.38.2 Biobase_2.36.2       KEGGREST_1.16.1      sonicLength_1.4.4   
 [9] hiAnnotator_1.11.1   geneRxCluster_1.12.0 GenomicRanges_1.28.6 GenomeInfoDb_1.12.3 
[13] Biostrings_2.44.2    XVector_0.16.0       IRanges_2.10.5       S4Vectors_0.14.7    
[17] BiocGenerics_0.22.1  forcats_0.3.0        stringr_1.3.1        dplyr_0.7.6         
[21] purrr_0.2.5          readr_1.1.1          tidyr_0.8.1          tibble_1.4.2        
[25] tidyverse_1.2.1      UpSetR_1.3.3         foreach_1.4.4        BiasedUrn_1.07      
[29] pander_0.6.1         magrittr_1.5         knitr_1.19           RColorBrewer_1.1-2  
[33] gridExtra_2.3        scales_0.5.0         ggrepel_0.7.0        ggplot2_3.0.0       
[37] lubridate_1.7.4      vegan_2.4-6          lattice_0.20-35      permute_0.9-4       
[41] data.table_1.10.4-3  RMySQL_0.10.14       DBI_0.8              reshape2_1.4.3      
[45] reldist_1.6-6        Matrix_1.2-9         igraph_1.2.2        

loaded via a namespace (and not attached):
 [1] colorspace_1.3-2           htmlTable_1.11.2           base64enc_0.1-3           
 [4] rstudioapi_0.7             bit64_0.9-7                xml2_1.2.0                
 [7] codetools_0.2-15           Formula_1.2-2              jsonlite_1.5              
[10] Rsamtools_1.28.0           broom_0.5.0                cluster_2.0.6             
[13] png_0.1-7                  compiler_3.4.0             httr_1.3.1                
[16] backports_1.1.2            assertthat_0.2.0           lazyeval_0.2.1            
[19] cli_1.0.1                  acepack_1.4.1              htmltools_0.3.6           
[22] tools_3.4.0                bindrcpp_0.2.2             gtable_0.2.0              
[25] glue_1.2.0                 GenomeInfoDbData_0.99.0    Rcpp_0.12.17              
[28] cellranger_1.1.0           nlme_3.1-131               rtracklayer_1.36.6        
[31] iterators_1.0.9            rvest_0.3.2                devtools_1.13.4           
[34] XML_3.98-1.9               zlibbioc_1.22.0            MASS_7.3-47               
[37] BSgenome_1.44.2            hms_0.4.2                  SummarizedExperiment_1.6.5
[40] yaml_2.1.16                memoise_1.1.0              rpart_4.1-11              
[43] RSQLite_2.0                latticeExtra_0.6-28        stringi_1.2.4             
[46] checkmate_1.8.5            BiocParallel_1.10.1        rlang_0.2.2               
[49] pkgconfig_2.0.1            bitops_1.0-6               matrixStats_0.53.0        
[52] bindr_0.1.1                GenomicAlignments_1.12.2   htmlwidgets_1.0           
[55] bit_1.1-12                 tidyselect_0.2.3           plyr_1.8.4                
[58] R6_2.2.2                   Hmisc_4.1-1                DelayedArray_0.2.7        
[61] pillar_1.1.0               haven_1.1.2                foreign_0.8-67            
[64] withr_2.1.1                mgcv_1.8-17                survival_2.41-3           
[67] RCurl_1.95-4.10            nnet_7.3-12                modelr_0.1.2              
[70] crayon_1.3.4               readxl_1.1.0               blob_1.1.0                
[73] digest_0.6.15              munsell_0.4.3          
```
