A number of analyses were performed on the count matrices derived by the read alignment and gene expression quantification steps. This folder includes the R script used to perform each analysis ("analysis_launch.Rnw"), the input data ("data/") and the expected results ("output/"). 

The "analysis_launch.Rnw" can be run using R (v.3.6.1) software (available at "https://cran.r-project.org/") to reproduce every analysis highlighted in the manuscript. The above R script downloads the count matrices from the GEO repository (available at "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151009") and stores them in the "GSE151009/" directory, launches the scripts stored in the "scripts/" folder, and stores the manuscript's figures and R objects in the "output/" directory. 

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
- readxl (v.1.3.1)
- reshape2 (v.1.4.4)
- scran (v.1.12.1)
- tidyr (v.1.1.0)
- umap (v.0.2.6.0)
- UpSetR (v.1.4.0)
- velocyto.R (v.0.6)
