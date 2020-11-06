print("1) Define output path")

outpath <- paste0(path, "output/fig5_XistXchrConcordance/"); dir.create(path = outpath, recursive = T, showWarnings = F)
f3path <- paste0(path, "output/fig3_deXistHighLow/")
f4path <- paste0(path, "output/fig4_deXchrChangeHighLow/")

print("2) A: Concordance between Xist and X-chromosome Change DE and Correlation analyses [d1 & d2]")

print("2.1) Load DE and Correlation analyses results")

results <- c()

print("2.1.1) Xist - DE analysis")

# load and filter data
load(paste0(f3path, "ALLresults.RData"))
fields <- c("day", "chromosome_name", "mgi_symbol", "coef", "fdr")
xist_de <- all_de[(all_de$day %in% c(1,2)), fields]
xist_de$direction <- ifelse(xist_de$coef>0, "Up-regulated", "Down-regulated") 
results <- data.frame(analysis = "DE", test = "Xist", 
                      xist_de[, c("day", "chromosome_name", "mgi_symbol", "direction", "fdr")])

print("2.1.2) Xist - Correlation analysis")

# load and filter data
load(paste0(f3path, "sprmcor.RData"))
colnames(sprmcor) <- c("day", "mgi_symbol", "chromosome_name", "correlation", 
                       "droprate", "droprate_Xist", "n", "pvalue_Xist", "fdr")
fields <- c("day", "chromosome_name", "mgi_symbol", "correlation", "fdr")
xist_cor <- sprmcor[(sprmcor$day %in% c(1,2)), fields]
xist_cor$direction <- ifelse(xist_cor$correlation>0, "Up-regulated", "Down-regulated") 
results <- rbind(results,
                 data.frame(analysis = "Correlation", test = "Xist", 
                            xist_cor[, c("day", "chromosome_name", "mgi_symbol", "direction", "fdr")]))

print("2.1.3) Xchr Change - DE analysis")

# load and filter data
load(paste0(f4path, "XchrChange_HighLow_MAST.RData"))
fields <- c("day", "chromosome_name", "mgi_symbol", "coef", "fdr")
xchr_de <- all_de[(all_de$day %in% c(1,2)), fields]
xchr_de$direction <- ifelse(xchr_de$coef>0, "Up-regulated", "Down-regulated") 
results <- rbind(results,
                 data.frame(analysis = "DE", test = "Xchr",
                            xchr_de[, c("day", "chromosome_name", "mgi_symbol", "direction", "fdr")]))

print("2.1.4) Xchr Change - Correlation analysis")

# load and filter data
load(paste0(f4path, "cor_Xchr.RData"))
colnames(sprmcor) <- c("day", "mgi_symbol", "chromosome_name", "correlation", 
                       "droprate", "n", "ypos", "pvalue", "fdr")
fields <- c("day", "chromosome_name", "mgi_symbol", "correlation", "fdr")
xchr_cor <- sprmcor[(sprmcor$day %in% c(1,2)), fields]
xchr_cor$direction <- ifelse(xchr_cor$correlation>0, "Up-regulated", "Down-regulated") 
results <- rbind(results,
                 data.frame(analysis = "Correlation", test = "Xchr",
                            xchr_cor[, c("day", "chromosome_name", "mgi_symbol", "direction", "fdr")]))


print("2.2) Compute number of significant hits per gene and day: FDR<0.05 & nhits>3")

# number of significant tests per gene
fdr_threshold <- 0.05
results$significant <- results$fdr <= fdr_threshold
sig_gene_day <- results %>% 
  dplyr::group_by(mgi_symbol, day) %>% 
  dplyr::summarise(nhits_pos = sum(significant[direction == "Up-regulated"]),
                   nhits_neg = sum(significant[direction == "Down-regulated"]),
                   nhits = max(c(sum(significant[direction == "Up-regulated"]),
                                 sum(significant[direction == "Down-regulated"]))))
table(sig_gene_day$nhits)
sig_gene <- results %>% dplyr::group_by(mgi_symbol) %>% dplyr::summarise(nhits = sum(significant))
table(sig_gene$nhits)
results$nhits <- sig_gene_day$nhits[match(paste0(results$mgi_symbol, "_", results$day), 
                                          paste0(sig_gene_day$mgi_symbol, "_", sig_gene_day$day))]
results$tnhits <- sig_gene$nhits[match(results$mgi_symbol, sig_gene$mgi_symbol)]

# select genes identified in at least 3 over 8 analyses for each day
n <- 3
tophits <- results[(results$nhits >= n)&(results$significant),] %>%
  arrange(-as.numeric(factor(direction)), -tnhits)
tophits_genes <- unique(as.character(tophits$mgi_symbol)); tophits_genes; length(tophits_genes)
temp <- results[results$mgi_symbol %in% tophits_genes,]
temp$mgi_symbol <- factor(temp$mgi_symbol, levels = tophits_genes)
temp <- temp %>% arrange(mgi_symbol)


print("2.3) Heatmap")

# input matrix
temp <- temp %>% dplyr::group_by(mgi_symbol) %>% dplyr::mutate(n = length(mgi_symbol))
temp <- temp[temp$n == 8,]
temp$lev <- paste0(temp$test, " - ", temp$analysis)
temp$dlev <- paste0("Day ", temp$day, ": ", temp$lev)
col_order <- c("Day 1: Xist - DE", "Day 1: Xist - Correlation", "Day 1: Xchr - DE", "Day 1: Xchr - Correlation",
               "Day 2: Xist - DE", "Day 2: Xist - Correlation", "Day 2: Xchr - DE", "Day 2: Xchr - Correlation")
temp$dlev <- factor(temp$dlev, levels = col_order)
temp <- data.frame(temp) %>% arrange(mgi_symbol, dlev, test)

# Color by direction and -log10(FDR)
temp$lfdr <- -log10(temp$fdr + 1e-6)
temp$col <- ifelse(temp$direction == "Up-regulated", temp$lfdr, -temp$lfdr)
x <- matrix(data = c(temp$col), ncol = 8, byrow = T)
rownames(x) <- unique(temp$mgi_symbol)
colnames(x) <- unique(temp$dlev)
x <- x[, col_order]

# order genes by minimum FDR
y <- apply(x, 1, function(x) max(abs(x))*sign(x[which.max(abs(x))]))
o <- order(y, decreasing = T)
x <- x[o,]

# identify significant ones
thr <- log10(fdr_threshold)
significant_labels <- matrix(ifelse((c(x) <= thr)|(c(x) >= -thr), "*", ""), ncol = ncol(x), byrow = F)

# aesthetics
colorPalette = rev(c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0")); paletteLength = 20;
myColor <- c(colorRampPalette(colorPalette)(paletteLength))
cellheight = 6; fontsize_col = 6; fontsize_row = 6; width = 10; cellwidth = 10; fontsize = 6
height <- cellheight*nrow(x)/30
breaks <- seq(-6, 6)

pheatmap(x,
         cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = TRUE, color = myColor,
         cellwidth = cellwidth, cellheight = cellheight,
         width = width, height = height,
         fontsize_col = fontsize_col, fontsize_row = fontsize_row, fontsize = fontsize,
         gaps_col = 4, gaps_row = 23, border_color = FALSE,
         display_numbers = significant_labels, number_color = "white",
         filename = paste0(outpath, "A_XistXchr_Concordance.pdf"))
