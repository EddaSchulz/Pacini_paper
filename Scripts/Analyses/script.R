### R script

### 1) Download data from GEO

# 1.1) define output folder
# path <- getwd()
path <- "C:/Users/pacini/Desktop/final_XXnew/NatureCommunications"
datapath <- paste0(path, "/GEOdata/"); dir.create(path = datapath, showWarnings = FALSE, recursive = TRUE)
outpath <- paste0(path, "/Output/"); dir.create(path = outpath, showWarnings = FALSE, recursive = TRUE)
scriptpath <- paste0(path, "/Scripts/"); dir.create(path = outpath, showWarnings = FALSE, recursive = TRUE)
dgepath <- paste0(outpath, "/ProcessedData/"); dir.create(path = dgepath, showWarnings = FALSE, recursive = TRUE)


# 1.2) Install and Load R packages
library(GEOquery)

# 1.3) Download data from GEO
geoID <- "GSE151009"
start <- Sys.time()
filePaths <- getGEOSuppFiles(GEO = geoID, baseDir = datapath)
end <- Sys.time()
time <- end - start; time


### 2) Data loading and pre-processing

start <- Sys.time()
source(paste0(scriptpath, "R/1_Data_PreProcessing.R"))
end <- Sys.time()
time <- end - start; time


### 3) Figure 1: 

## 2.1) Load gene and sample information
library(data.table); library(plyr); library(dplyr); library("readxl")
cell <- read_excel(path = paste0(datapath, "ToAdd/S1_Cell&Gene_Filtering.xlsx"), sheet = "Cell Filtering", skip = 16) %>% as.data.frame()
cell_features <- c("ID", "Day", "Empty", "Multiple", "Dead Stain [intensity]", "Total Reads")
gene <- read_excel(path = paste0(datapath, "ToAdd/S1_Cell&Gene_Filtering.xlsx"), sheet = "Gene Filtering", skip = 15) %>% as.data.frame()
gene_features <- c("Chromosome", "Ensembl ID", "Gene Symbol", "Strand")

# 2.2) Define DGE list object for notAS/AS, exonic/spliced/unspliced gene expression quantification
library(edgeR)
for(allele in c("", "_B6", "_Cast")){
  for(w in c("", "_Spliced", "_Unspliced")){
    x <- read.table(file = paste0(datapath, geoID, "/", geoID, allele, w, "_UMICountMatrix.txt.gz"), row.names = 1)
    genes <- gene[match(rownames(x), gene$`Gene Symbol`), gene_features]; colnames(genes) <- c("chromosome", "ensembl", "symbol", "strand")
    samples <- cell[match(colnames(x), cell$ID), cell_features]; colnames(samples) <- c("id", "day", "empty", "multiple", "dead_stain", "seqdepth")
    dge <- DGEList(counts = x, samples = samples, genes = genes)
    save(dge, file = paste0(outpath, "UF_DGE", allele, w, ".RData"))
  }
}

# 2.3) Data filtering
MADoutlier <- function(value, nmads = 3, type, log = FALSE){
  if(log == TRUE){
    value <- log10(value)
  }
  cur.med <- median(value, na.rm = TRUE)
  cur.mad <- mad(value, center = cur.med, na.rm = TRUE)
  upper.limit <- cur.med + nmads * cur.mad
  lower.limit <- cur.med - nmads * cur.mad
  if(type == "lower"){
    outliers <- value < lower.limit
  }
  if(type == "upper"){
    outliers <- value > upper.limit
  }
  if(type == "both"){
    outliers <- value < lower.limit | value > upper.limit
  }
  return(outliers)
}

# 2.3.1) Cell filtering

## 2.3.1.1) Image and cell feature filtering
load(paste0(outpath, "UF_DGE.RData"))
dge$genes$isERCC <- grepl(dge$genes$symbol, pattern = "^ERCC")
dge$genes$isMT <- grepl(dge$genes$symbol, pattern = "^mt")
dge$samples$perc_exprgenes <- colMeans(dge$counts>0)*100
dge$samples$perc_ERCC <- colSums(dge$counts[dge$genes$isERCC,])/colSums(dge$counts)*100
dge$samples$perc_mtDNA <- colSums(dge$counts[dge$genes$isMT,])/colSums(dge$counts)*100
dge$samples$Red <- as.numeric(as.character(dge$samples$dead_stain))

filter_matrix <- plyr::ddply(dge$samples, .variables = .(group), transform, 
                             seq_filter = MADoutlier(seqdepth, type = "lower", log = TRUE, nmads = 3),
                             lib_filter = MADoutlier(lib.size, type = "lower", log = TRUE, nmads = 3),
                             expr_filter = MADoutlier(perc_exprgenes, type = "lower", log = TRUE, nmads = 3),
                             ercc_filter = MADoutlier(perc_ERCC, type = "upper", nmads = 3),
                             mtdna_filter = MADoutlier(perc_mtDNA, type = "upper", nmads = 3),
                             red_filter = MADoutlier(Red, type = "upper", nmads = 3))
filter_matrix$image_filter <- rowSums(filter_matrix[, c("empty", "multiple")]==TRUE)>0
filter_matrix$filtered <- rowSums(filter_matrix[, grepl(x = colnames(filter_matrix), pattern = "filter")]) > 0
remove_cells <- filter_matrix$id[filter_matrix$filtered]

## 2.3.1.2) XO cells filtering
load(paste0(outpath, "UF_DGE_B6.RData")); b6 <- dge
load(paste0(outpath, "UF_DGE_Cast.RData")); cast <- dge
r <- data.frame(id = b6$samples$id,
                Xist_AS_UMI = (b6$counts["Xist",] + cast$counts["Xist",]),
                b6Xsum = colSums(b6$counts[b6$genes$chromosome %in% "X",]),
                castXsum = colSums(cast$counts[cast$genes$chromosome %in% "X",]),
                b6Xratio = colSums(b6$counts[b6$genes$chromosome %in% "X",])/(colSums(b6$counts[b6$genes$chromosome %in% "X",]) + colSums(cast$counts[cast$genes$chromosome %in% "X",])))
r[is.na(r$b6Xratio),]
XO_threshold <- 0.8
r$XO_filtering <- ifelse(is.na(r$b6Xratio), NA, 
                         (r$Xist_AS_UMI==0)&(abs(r$b6Xratio-0.5)>abs(XO_threshold-0.5)))
remove_XOcells <- as.character(r$id[r$XO_filtering %in% TRUE])

## 2.3.1.3) Remove cells from all DGE objects
for(f in list.files(path = outpath, full.names = TRUE)){
  load(f)
  dge <- dge[, !colnames(dge) %in% c(remove_cells, remove_XOcells)]
  save(dge, file = gsub(x = f, pattern = "UF_DGE", replacement = "CF_DGE"))
}


# 2.3) Gene filtering

## 2.3.1) not-AS dropout rate
load(paste0(outpath, "CF_DGE.RData")); notas <- dge
detrate_notas <- rowMeans(notas$counts>0)

## 2.3.2) AS dropout rate
load(paste0(outpath, "CF_DGE_B6.RData")); b6 <- dge
load(paste0(outpath, "CF_DGE_Cast.RData")); cast <- dge
as <- b6$counts+cast$counts
detrate_as <- rowMeans(as>0)

## 2.3.3) misleading SNP annotation
r <- data.frame(chr = b6$genes$chromosome, 
                symbol = b6$genes$symbol, 
                both_detrate = rowMeans((b6$counts+cast$counts)>0),
                both_detrate_filtering = rowMeans((b6$counts+cast$counts)>0) <= 0.2,
                both_AScount = rowSums(b6$counts+cast$counts), 
                b6_AScount = rowSums(b6$counts),
                B6Total_ASratio = rowSums(b6$counts)/rowSums(b6$counts+cast$counts))
r$misleadingSNP_filtering <- abs(r$B6Total_ASratio - 0.5) > 0.4; table(r$misleadingSNP_filtering)

## 2.3.4) Identify genes to be removed in AS and not-AS DGE lists
GF_notas <- data.frame(Gene = notas$genes$symbol, detrate = detrate_notas, r[match(notas$genes$symbol, r$symbol), -c(1,2)]) 
GF_as <- data.frame(Gene = b6$genes$symbol, detrate = detrate_as, r[match(b6$genes$symbol, r$symbol), -c(1,2)])
remove_notas <- GF_notas$detrate <= 0.2 | (GF_notas$misleadingSNP_filtering & !GF_notas$both_detrate_filtering & !is.na(GF_notas$both_detrate_filtering))
remove_as <- GF_as$detrate <= 0.2 | GF_as$misleadingSNP_filtering

## 2.3.5) Remove genes from AS and not-AS DGE lists
files_DGE <- list.files(path = outpath, full.names = TRUE)[grepl(x = list.files(path = outpath, full.names = TRUE), pattern = "CF_") & !grepl(x = list.files(path = outpath, full.names = TRUE), pattern = "pliced")]
notas_DGE <- files_DGE[!(grepl(x = files_DGE, pattern = "B6") | grepl(x = files_DGE, pattern = "Cast"))]
as_DGE <- files_DGE[grepl(x = files_DGE, pattern = "B6") | grepl(x = files_DGE, pattern = "Cast")]
for(f in files_DGE){
  if(f %in% notas_DGE){
    load(f); dge <- dge[!remove_notas,]; save(dge, file = gsub(x = f, pattern = "CF_", replacement = "GFCF_"))
  }else{
    load(f); dge <- dge[!remove_as,]; save(dge, file = gsub(x = f, pattern = "CF_", replacement = "GFCF_"))
  }
}


# 2.4) Data normalization

## 2.4.1) Compute size factors on not-AS DGE list based on autosomal gene expression
library(scran)
load(paste0(outpath, "GFCF_DGE.RData"))
input <- dge$counts[dge$genes$chromosome %in% c(1:19),]
emp.clusters <- as.numeric(as.factor(strsplit2(colnames(input), split = "\\_")[, 1]))
sizefact <- scran::computeSumFactors(x = input, clusters = emp.clusters)
sf <- data.frame(day = dge$samples$day, ID = dge$samples$id, sf = sizefact)
save(sf, file = paste0(outpath, "PoolClust_sizefactors_autosomal.RData"))

## 2.4.2) Store size factors in DGE lists
files_DGE <- list.files(path = outpath, full.names = TRUE)[grepl(x = list.files(path = outpath, full.names = TRUE), pattern = "GFCF_")]
for(f in files_DGE){
  load(f); dge$samples$sf_notX <- sizefact; dge$samples$eff_libsize_notX <- colSums(dge$counts)*dge$samples$sf_notX; save(dge, file = gsub(x = f, pattern = "GFCF_", replacement = "NGFCF_"))
}


# 2.5) Cell classification

# 2.5.1) Xist notAS expression: Xist UMI
load(paste0(outpath, "NGFCF_DGE.RData")); xist_umi <- dge$counts["Xist",]
df <- data.frame(day = dge$samples$day, ID = dge$samples$id, xist_umi)
df$Xist_class <- ifelse(df$xist_umi == 0, "Undetected", ifelse(df$xist_umi <= 5, "Detected (Xist UMI <= 5)",  "Detected (Xist UMI > 5)"))

# 2.5.2) Xist notAS expression: CPM Xist expression
xist_cpm <- xist_umi/dge$samples$eff_libsize_notX*1e6
k <- 7; kmeans_class <- kmeans(log1p(xist_cpm), centers = k, iter.max = 1000, nstart = 1000)
df <- data.frame(df, xist_cpm, Xist_kmeans_class = kmeans_class$cluster)
subs <- ddply(df, .variables = .(Xist_kmeans_class), summarize, median_Xist = median(xist_cpm))
subs <- subs[order(subs$median_Xist, decreasing = FALSE), ]; subs$ordered <- seq_len(nrow(subs))
newvalue <- subs$ordered; names(newvalue) <- subs$Xist_kmeans_class
df$Xist_kmeans_class <- revalue(as.character(df$Xist_kmeans_class), newvalue); table(day = df$day, Xist_CPM_Kmeans = df$Xist_kmeans_class)

# 2.5.2) Xist AS expression
load(paste0(outpath, "NGFCF_DGE_B6.RData")); b6 <- dge
load(paste0(outpath, "NGFCF_DGE_Cast.RData")); cast <- dge
xist_as <- b6$counts["Xist",] + cast$counts["Xist",]
xist_ratio <- b6$counts["Xist",]/xist_as
xchr_b6 <- colSums(b6$counts[b6$genes$chromosome %in% "X",])
xchr_cast <- colSums(cast$counts[cast$genes$chromosome %in% "X",])
xchr_ratio <- xchr_b6/(xchr_b6 + xchr_cast)
df <- data.frame(df, xist_b6 = b6$counts["Xist",], xist_cast = cast$counts["Xist",], xist_ratio, xchr_b6, xchr_cast, xchr_ratio)
df$Xist_ratio_class <- ifelse(df$xist_ratio == 0, "Cast_MA", ifelse(df$xist_ratio == 1, "BL6_MA", ifelse((df$xist_ratio >= 0.2) & (df$xist_ratio <= 0.8), "Xist_BA", "Middle")))
df$Xist_ratio_class[(df$xist_b6 + df$xist_cast) <= 5] <- "Low-Xist"
df$Xist_ratio_class[df$xist_umi == 0] <- "Undetected"
table(day = df$day, df$Xist_ratio_class)

# 2.5.3) Store cell classifications on DGE lists
files_DGE <- list.files(path = outpath, full.names = TRUE)[grepl(x = list.files(path = outpath, full.names = TRUE), pattern = "NGFCF_")]
for(f in files_DGE){
  load(f); dge$samples <- data.frame(dge$samples, df[, -c(1,2)]); save(dge, file = gsub(x = f, pattern = "NGFCF_", replacement = ""))
}













