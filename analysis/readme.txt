A number of analyses were performed on the count matrices derived by the read alignment and gene expression quantification steps. This folder includes the R script used to perform each analysis ("analysis.Rnw"), the input data sets ("data/") and the expected results ("output/"). Where the "analysis.Rnw" script takes as input the count matrices and gene annotations stored in the "data/" directory, and store the results of each analysis in the "output/" folder.

The "analysis_launch.Rnw" can be run using R (v.3.6.1) software (available at "https://cran.r-project.org/") in order to reproduce every analysis highlighted in the manuscript. The above R script downloads the count matrices from the GEO repository GSE151009 (available at "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151009"), and stores the manuscript's figures and additional R objects in the "output/" directory. 

The following R libraries are required (in brackets the library version used for data analyses):

- biomaRt (v.2.40.5)
- data.table (v.1.12.8)
- doParallel (v.1.0.15)
- dplyr (v.1.0)
- edgeR (v.3.26.8)
- GEOquery (v.2.52.0)
- ggalluvial (v.0.11.3)
- ggplot2 (v.3.3.2)
- ggrepel (v.0.8.2)
- gridExtra (v.2.3)
- lsr (v.0.5)
- MAST (v.1.10)
- matrixStats (v.0.56.0)
- monocle (v.2.12.0)
- multidplyr (v.0.0.0.9000)
- pcaMethods (v.1.76.0)
- pheatmap (v.1.0.12)
- plyr (v.1.8.6)
- reshape2 (v.1.4.4)
- scran (v.1.12.1)
- tidyr (v.1.1.0)
- umap (v.0.2.6.0)
- velocyto.R (v.0.6)
