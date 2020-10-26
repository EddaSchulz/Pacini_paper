print("1) Download gene annotation")
mart <- useMart("ensembl")
mart <- useDataset("mmusculus_gene_ensembl", mart = mart)
features <- listAttributes(mart)
f <- c('mgi_symbol', 'chromosome_name', 'entrezgene_id', 'strand', 'start_position', 'end_position', 
       'percentage_gene_gc_content', 'ensembl_gene_id', 'gene_biotype', 'go_linkage_type', 'namespace_1003', 'name_1006', 'go_id', 'description')
bm <- getBM(attributes = f, mart = mart)
save(bm, file = paste0(datapath, "biomart.RData"))

print("2) Download count matrices from GEO")
geoID <- "GSE151009"
files <- getGEOSuppFiles(GEO = geoID, baseDir = path, fetch_files = FALSE)
subfiles <- as.list(as.character(files$fname[grepl(x = files$fname, pattern = "CountMatrix")]))
lapply(X = subfiles, FUN = function(x) getGEOSuppFiles(GEO = geoID, baseDir = path, filter_regex = x))