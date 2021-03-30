print("1) Define output path")

outpath <- paste0(path, "output/fig5_XistXchrConcordance/"); dir.create(path = outpath, recursive = T, showWarnings = F)
f4path <- paste0(path, "output/fig4_deXistHighLow/")
suppath <- paste0(path, "output/fig_Supplementary/")
fdr_threshold <- 0.05

print("2) A: Concordance between Xist and X-chromosome Change (dX) DE and Correlation analyses [days 1 and 2]")

print("2.1) Load DE and Correlation analyses results")

results <- c()

print("2.1.1) Xist - DE analysis")

# load data
load(paste0(f4path, "ALLresults.RData"))
fields <- c("day", "chromosome_name", "mgi_symbol", "coef", "fdr")
xist_de <- all_de[(all_de$day %in% c(1,2)), fields]
xist_de$value <- xist_de$coef
xist_de$direction <- ifelse(xist_de$coef>0, "Up-regulated", "Down-regulated") 
results <- data.frame(analysis = "DE", test = "Xist", 
                      xist_de[, c("day", "chromosome_name", "mgi_symbol", "direction", "fdr", "value")])

print("2.1.2) Xist - Correlation analysis")

# load data
load(paste0(f4path, "sprmcor.RData"))
colnames(sprmcor) <- c("day", "mgi_symbol", "chromosome_name", "correlation", 
                       "droprate", "droprate_Xist", "n", "pvalue_Xist", "fdr")
fields <- c("day", "chromosome_name", "mgi_symbol", "correlation", "fdr")
xist_cor <- sprmcor[(sprmcor$day %in% c(1,2)), fields]
xist_cor$value <- xist_cor$correlation
xist_cor$direction <- ifelse(xist_cor$correlation>0, "Up-regulated", "Down-regulated") 
results <- rbind(results,
                 data.frame(analysis = "Correlation", test = "Xist", 
                            xist_cor[, c("day", "chromosome_name", "mgi_symbol", "direction", "fdr", "value")]))

print("2.1.3) Xchr Change - DE analysis")

# load data
load(paste0(suppath, "XchrChange_HighLow_MAST.RData"))
fields <- c("day", "chromosome_name", "mgi_symbol", "coef", "fdr")
xchr_de <- all_de[(all_de$day %in% c(1,2)), fields]
xchr_de$value <- xchr_de$coef
xchr_de$direction <- ifelse(xchr_de$coef>0, "Up-regulated", "Down-regulated") 
results <- rbind(results,
                 data.frame(analysis = "DE", test = "Xchr",
                            xchr_de[, c("day", "chromosome_name", "mgi_symbol", "direction", "fdr", "value")]))

print("2.1.4) Xchr Change - Correlation analysis")

# load and filter data
load(paste0(suppath, "cor_Xchr.RData"))
colnames(sprmcor) <- c("day", "mgi_symbol", "chromosome_name", "correlation", 
                       "droprate", "n", "ypos", "pvalue", "fdr")
fields <- c("day", "chromosome_name", "mgi_symbol", "correlation", "fdr")
xchr_cor <- sprmcor[(sprmcor$day %in% c(1,2)), fields]
xchr_cor$value <- xchr_cor$correlation
xchr_cor$direction <- ifelse(xchr_cor$correlation>0, "Up-regulated", "Down-regulated") 
results <- rbind(results,
                 data.frame(analysis = "Correlation", test = "Xchr",
                            xchr_cor[, c("day", "chromosome_name", "mgi_symbol", "direction", "fdr", "value")]))

print("2.2) Exclude pseudogenes and X-linked genes with significant negative correlation coef. or logFC")

results$gene_biotype <- bm$gene_biotype[match(results$mgi_symbol, bm$mgi_symbol)]
results <- results[!grepl(x = results$gene_biotype, pattern = "pseudogene"),]
results <- results[!(results$chromosome_name %in% "X" & results$value < 0 & results$fdr <= fdr_threshold),]

print("2.3) Compute number of significant hits per gene and day: FDR<=0.05 & nhits>3")

# number of significant tests per gene
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


print("2.4) Heatmap")

# input matrix
temp <- temp %>% dplyr::group_by(mgi_symbol) %>% dplyr::mutate(n = length(mgi_symbol))
temp <- temp[temp$n == 8,]
temp$dlev <- paste0(temp$test, " - ", temp$analysis)
col_order <- c("Xist - DE", "Xist - Correlation", "Xchr - DE", "Xchr - Correlation")
temp$dlev <- factor(temp$dlev, levels = col_order)
temp <- data.frame(temp) %>% arrange(mgi_symbol, dlev, test)

# Color by direction and -log10(FDR)
temp$lfdr <- -log10(temp$fdr + 1e-6)
temp$col <- ifelse(temp$direction == "Up-regulated", temp$lfdr, -temp$lfdr)

# order genes by minimum FDR
x <- temp %>% dplyr::group_by(mgi_symbol) %>% dplyr::summarise(m = max(abs(col))*sign(col[which.max(abs(col))])) %>%
  arrange(-m) %>% as.data.frame()
order_genes <- rev(as.character(x$mgi_symbol))
temp$mgi_symbol <- factor(temp$mgi_symbol, levels = order_genes)

# label significant results
temp$sig_label <- ifelse(temp$fdr <= fdr_threshold, "*", "")

# plot
h <- length(unique(temp$mgi_symbol))*0.225
w <- length(unique(temp$dlev))*0.4
temp$Day <- paste0("Day ", temp$day)

g <- temp %>%
  ggplot() +
  theme_bw() + theme1 + 
  facet_grid(. ~ Day) +
  geom_tile(aes(x = factor(dlev), y = mgi_symbol, fill = col), color = "white") +
  geom_text(aes(x = factor(dlev), y = mgi_symbol, label = sig_label), color = "white", size = geomtext_size) +
  scale_fill_gradient2(low = "#0571b0", high = "#ca0020", mid = "white", 
                       midpoint = 0, limit = c(-6, 6),
                       name= "-log10FDR") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
adjust_size(g = g, panel_width_cm = w, panel_height_cm = h, 
            savefile = paste0(outpath, "A_XistXchr_Concordance.pdf"), 
            height = 10, width = 5)

print("2.5) Heatmap - log2FC: MAST DE analyses")

# subset data
temp <- results[(results$mgi_symbol %in% order_genes)&(results$analysis == "DE"),]
temp$mgi_symbol <- factor(temp$mgi_symbol, levels = order_genes)
temp$dlev <- paste0(temp$test, " - ", temp$analysis)
col_order <- c("Xist - DE", "Xchr - DE")
temp$dlev <- factor(temp$dlev, levels = col_order)
temp$sig_label <- ifelse(temp$fdr <= fdr_threshold, "*", "")

# plot
h <- length(unique(temp$mgi_symbol))*0.225
w <- length(unique(temp$dlev))*0.4
temp$Day <- paste0("Day ", temp$day)

g <- temp %>%
  ggplot() +
  theme_bw() + theme1 + 
  facet_grid(. ~ Day) +
  geom_tile(aes(x = factor(dlev), y = mgi_symbol, fill = value), color = "white") +
  geom_text(aes(x = factor(dlev), y = mgi_symbol, label = sig_label), color = "white", size = geomtext_size) +
  scale_fill_gradient2(low = "#0571b0", high = "#ca0020", mid = "white", 
                       midpoint = 0, limit = c(-3, 3),
                       name= "log2FC") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
adjust_size(g = g, panel_width_cm = w, panel_height_cm = h, 
            savefile = paste0(outpath, "B_XistXchr_DE_log2FC.pdf"), 
            height = 10, width = 5)


print("2.6) Heatmap - Spearman's coefficient: Correlation analyses")

# subset data
temp <- results[(results$mgi_symbol %in% order_genes)&(results$analysis == "Correlation"),]
temp$mgi_symbol <- factor(temp$mgi_symbol, levels = order_genes)
temp$dlev <- paste0(temp$test, " - ", temp$analysis)
col_order <- c("Xist - Correlation", "Xchr - Correlation")
temp$dlev <- factor(temp$dlev, levels = col_order)
temp$sig_label <- ifelse(temp$fdr <= fdr_threshold, "*", "")

# plot
h <- length(unique(temp$mgi_symbol))*0.225
w <- length(unique(temp$dlev))*0.4
temp$Day <- paste0("Day ", temp$day)

g <- temp %>%
  ggplot() +
  theme_bw() + theme1 + 
  facet_grid(. ~ Day) +
  geom_tile(aes(x = factor(dlev), y = mgi_symbol, fill = value), color = "white") +
  geom_text(aes(x = factor(dlev), y = mgi_symbol, label = sig_label), color = "white", size = geomtext_size) +
  scale_fill_gradient2(low = "#0571b0", high = "#ca0020", mid = "white", 
                       midpoint = 0, limit = c(-0.5, 0.5),
                       name= expression(rho)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
adjust_size(g = g, panel_width_cm = w, panel_height_cm = h, 
            savefile = paste0(outpath, "C_XistXchr_Correlation_rho.pdf"), 
            height = 10, width = 5)