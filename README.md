# CART19 Integration Site Manuscript Analysis

Analysis provided in the "Linking outcome of CD19-directed CAR T cell therapy with genome modification by vector integration" manuscript.

This repository is meant to provide the source code for the analysis provided in the above manuscript by Christopher L. Nobles, Scott Sherrill-Mix, John Everett, Shantan Reddy, Joseph A. Fraietta, David Porter, Noelle Frey, Saar Gill, Steven Grupp, Shannon Maude, Donald Siegel, Bruce Levine, Carl H. June, Simon F. Lacey, J. Joseph Melenhorst, and Frederic D. Bushman.

Abstract:
    Chimeric antigen receptor-engineered T-cells targeting CD19 (CART19) provide a highly effective treatment for pediatric acute lymphoblastic leukemia (ALL), but are less effective for chronic lymphocytic leukemia (CLL) and other indications, focusing attention on approaches to improve efficacy.  CART19 cells are prepared by integration of the engineered receptor gene into the host T cell chromosome using a lentiviral vector.  Vector integration thus marks T cell lineages uniquely, and modifies the cellular genome by insertional mutagenesis.  Previously we reported that vector integration in the host gene TET2 resulted in gene inactivation and consequently therapeutic expansion of a cell clone associated with remission of CLL.  Here we investigate clonal population structure and therapeutic outcome for an additional 39 subjects by high throughput sequencing of sites of vector integration.  Expansion of cell clones with integration sites in specific genes in responders suggests possible insertional mutagenesis promoting therapeutic proliferation. XXX These data thus establish that a signal is present in integration site data associated with clinical response, and provide tools to help optimize T cell engineering.

To run the analysis, run the R-script 'XXX.R', which will subsequenty process the input data, populate the directory with informational files, and generate a report with the manuscript figures. Additionally, the user can run the 'generate_goi_report.R' script after the intial script to generate the supplemental report included in the manuscript.

This repository also contains a final INSPIIRED input dataset from the manuscript. Anyone can use this data rather than reprocessing with the INSPIIRED software from the SRA submission (BioProject: PRJNA510570). These data were processed with the [intSiteCaller v1.2.0 software](https://github.com/BushmanLab/intSiteCaller), within the [INSPIIRED suite](https://github.com/BushmanLab/INSPIIRED).
