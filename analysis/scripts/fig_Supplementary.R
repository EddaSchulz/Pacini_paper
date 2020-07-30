library(UpSetR)

print("1) Define output path")
outpath <- paste0(path, "output/fig_Supplementary/"); dir.create(path = outpath, recursive = T, showWarnings = F)

print("2) Supplementary figure 1: Alignment, Gene counting and Cell filtering")

print("2.1) A: Alignment and gene quantification")

print("2.1.1) Load data")
load(paste0(datapath, "UF_DGE.RData")); notas <- dge
load(paste0(datapath, "UF_DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "UF_DGE_Cast.RData")); cast <- dge

summary_stats <- data.frame(id = notas$samples$id, 
                            day = notas$samples$day,
                            total = notas$samples$seqdepth, 
                            ua = notas$samples$uniquealigned, 
                            exonic = colSums(notas$counts), 
                            b6 = colSums(b6$counts), 
                            cast = colSums(cast$counts))

print("2.1.2) Compute median number and percentage of reads across samples")
sum <- summary_stats %>%
  dplyr::group_by(day) %>% 
  dplyr::summarise(seqdepth = median(total), 
                   ua = median(ua), 
                   exonic = median(exonic),
                   b6 = median(b6), 
                   cast = median(cast))
perc <- data.frame(summary_stats[, c("total", "ua", "exonic", "b6", "cast")]/summary_stats$total*100, day = summary_stats$day)
perc <- perc %>%
  dplyr::group_by(day) %>% 
  dplyr::summarise(seqdepth = median(total), 
                   ua = median(ua), 
                   exonic = median(exonic),
                   b6 = median(b6), 
                   cast = median(cast))

print("2.1.3) Melt variables")
sum_melt <- tidyr::gather(sum, 'variable', 'value', -day)
perc_melt <- tidyr::gather(perc, 'variable', 'value', -day)
perc_melt$value <- paste0(round(perc_melt$value, digits = 1), "%")
sum_melt$perc <- perc_melt$value

newlev <- c("seqdepth" = "Sequencing\nDepth", "ua" = "Uniquely\nAligned", "exonic" = "UMI\nnot-AS", "b6" = "UMI\nB6", "cast" = "UMI\nCast")
sum_melt$variable <- revalue(sum_melt$variable, replace = newlev)
perc_melt$variable <- revalue(perc_melt$variable, replace = newlev)
sum_melt$variable <- factor(sum_melt$variable, levels = newlev)

print("2.1.4) Plot")
g <- sum_melt %>% 
  ggplot() + 
  theme_bw() + theme1 + 
  geom_bar(aes(x = variable, y = value, fill = factor(day)), stat = "identity", position=position_dodge()) + 
  geom_text(aes(x = variable, y = value, group = factor(day), label=perc), hjust = -.1,
            color="black", position = position_dodge(0.9), size = geomtext_size, angle = 90, show.legend = FALSE) +
  scale_fill_grey() +
  scale_y_continuous(label = scientific_label, breaks = seq(0, 6e5, by = 1e5), limits = c(0, 7e5)) + 
  labs(x = "", y = "Number of reads", fill = "Time [days]")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "S1_A_SeqOutput.pdf"))



print("2.2) B: Cell filtering")

print("2.2.1) Load data")
load(paste0(datapath, "UF_DGE.RData"))
dge$genes$isERCC <- grepl(dge$genes$symbol, pattern = "^ERCC")
dge$genes$isMT <- grepl(dge$genes$symbol, pattern = "^mt")
dge$samples$exprgenes <- colSums(dge$counts>0)
dge$samples$perc_ERCC <- colSums(dge$counts[dge$genes$isERCC,])/colSums(dge$counts)*100
dge$samples$perc_mtDNA <- colSums(dge$counts[dge$genes$isMT,])/colSums(dge$counts)*100
dge$samples$Red <- as.numeric(as.character(dge$samples$dead_stain))

print("2.2.2) Define filtering matrix")
filter_matrix <- ddply(dge$samples, .variables = .(group), transform, 
                       seq_filter = MADoutlier(seqdepth, type = "lower", log = TRUE, nmads = 3),
                       lib_filter = MADoutlier(lib.size, type = "lower", log = TRUE, nmads = 3),
                       expr_filter = MADoutlier(exprgenes, type = "lower", log = TRUE, nmads = 3),
                       ercc_filter = MADoutlier(perc_ERCC, type = "upper", nmads = 3),
                       mtdna_filter = MADoutlier(perc_mtDNA, type = "upper", nmads = 3),
                       red_filter = MADoutlier(Red, type = "upper", nmads = 3))
filter_matrix$image_filter <- rowSums(filter_matrix[, c("empty", "multiple")]==1)>0
filter_matrix$filtered <- rowSums(filter_matrix[, grepl(x = colnames(filter_matrix), pattern = "filter")]) > 0
dge$samples$MAD_filtering <- dge$samples$id %in% filter_matrix$id[filter_matrix$filtered == TRUE]
features <- c("id", "day", "seqdepth", "lib.size", "exprgenes", "perc_ERCC", "perc_mtDNA", "Red", 
              "empty", "multiple", 
              "MAD_filtering")
cells <- dge$samples[, features]

print("2.2.3) Melt variables")
features_filter <- c("seqdepth", "lib.size", "exprgenes", "perc_ERCC", "perc_mtDNA", "Red")
temp <- dge$samples[, colnames(dge$samples) %in% c("day", "id", features_filter)]
s_melt <- tidyr::gather(temp, 'variable', 'value', -day, -id)
s_melt <- s_melt[s_melt$variable %in% features_filter,]
s_melt$variable <- revalue(s_melt$variable, c("seqdepth" = "Sequencing\nDepth", 
                                              "lib.size" = "Library\nSize", 
                                              "exprgenes" = "Expressed\nGenes",
                                              "perc_ERCC" = "ERCC (%)", "perc_mtDNA" = "mtDNA (%)", 
                                              "Red" = "Red fluorophore\n(AF594)"))
s_melt$value <- as.numeric(s_melt$value)

print("2.2.4) Compute 3*MAD threshold values")
MAD_threshold <- function(value, nmads = 3, type, log = FALSE){
  if(log == TRUE){
    value <- log10(value)
  }
  cur.med <- median(value, na.rm = TRUE)
  cur.mad <- mad(value, center = cur.med, na.rm = TRUE)
  upper.limit <- cur.med + nmads * cur.mad
  lower.limit <- cur.med - nmads * cur.mad
  return(c(upper.limit, lower.limit))
}
thresholds <- ddply(s_melt, .variables = .(variable), summarize, 
                    threshold_log_low = MAD_threshold(value, log = TRUE, nmads = 3)[2],
                    threshold_log_high = MAD_threshold(value, log = TRUE, nmads = 3)[1],
                    threshold_nolog_low = MAD_threshold(value, log = FALSE, nmads = 3)[2],
                    threshold_nolog_high = MAD_threshold(value, log = FALSE, nmads = 3)[1])

print("2.2.5) Plot - filter low features")
lv <- c("Sequencing\nDepth", 
        "Library\nSize", 
        "Expressed\nGenes")
lv_df <- s_melt[s_melt$variable %in% lv,]
g <- lv_df %>% 
  ggplot() +
  theme_bw() +  theme1 + 
  facet_grid(variable~., scales = "free_y") +
  geom_violin(aes(x = factor(day), y = value), alpha = 1/2, draw_quantiles = c(0.5), size = violin_box_size) + 
  geom_jitter(aes(x = factor(day), y = value), alpha = 1/4, size = outliersize,
              position=position_jitter(width = .1), shape = 21) +
  geom_hline(data = thresholds[thresholds$variable %in% lv,], aes(yintercept = 10^threshold_log_low),
             linetype = "dashed", color = "red", size = linesize) + 
  labs(x="Time [days]", y = "") +
  theme(strip.text.y = element_text(angle = 0))
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, savefile = paste0(outpath, "S1_B_LowValueFilter.pdf"))

print("2.2.5) Plot - filter high features")
hv <- c("ERCC (%)", "mtDNA (%)", "Red fluorophore\n(AF594)")
hv_df <- s_melt[s_melt$variable %in% hv,]
g <- hv_df %>% 
  ggplot() +
  theme_bw() +  theme1 + 
  facet_grid(variable~., scales = "free_y") +
  geom_violin(aes(x = factor(day), y = log10(value + 0.01)), alpha = 1/2, draw_quantiles = c(0.5), size = violin_box_size) + 
  geom_jitter(aes(x = factor(day), y = log10(value + 0.01)), alpha = 1/4, size = outliersize,
              position=position_jitter(width = .1), shape = 21) +
  geom_hline(data = thresholds[thresholds$variable %in% hv,], aes(yintercept = log10(threshold_nolog_high + 0.01)),
             linetype = "dashed", color = "red", size = linesize) + 
  labs(x="Time [days]", y = expression(log[10]*"(value + 0.01)]")) +
  theme(strip.text.y = element_text(angle = 0))
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, savefile = paste0(outpath, "S1_B_HighValueFilter.pdf"))



print("2.3) C: Remove putative XO cells")

print("2.3.1) Load data")
load(paste0(datapath, "UF_DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "UF_DGE_Cast.RData")); cast <- dge
r <- data.frame(id = dge$samples$id,
                day = dge$samples$day,
                Xist_AS_UMI = b6$counts["Xist",] + cast$counts["Xist",],
                libsize_AS = colSums(b6$counts) + colSums(cast$counts),
                b6Xsum = colSums(b6$counts[b6$genes$chromosome %in% "X",]),
                castXsum = colSums(cast$counts[cast$genes$chromosome %in% "X",])
)
r$bothXsum <- r$b6Xsum + r$castXsum
r$b6Xratio <- r$b6Xsum/r$bothXsum

print("2.3.2) Set threshold and highlight putative XO cells")
r[is.na(r$b6Xratio),]
XO_threshold <- 0.8
r$XO_filtering <- ifelse(is.na(r$b6Xratio), NA, 
                         (r$Xist_AS_UMI==0)&(abs(r$b6Xratio-0.5)>abs(XO_threshold-0.5)))
r$Xist_AS_CPM <- (r$Xist_AS_UMI/r$libsize_AS)*1e6
r$Day <- paste0("Day ", r$day)

print("2.3.3) Plot")
g <- r[!is.na(r$b6Xratio),] %>% ggplot() + 
  theme_bw() + theme1 +
  geom_point(aes(x = b6Xratio, y = log10(Xist_AS_CPM + 1)), alpha = 0.5,
             size = small_scattersize, show.legend = FALSE) + 
  geom_vline(xintercept = c(abs(1-XO_threshold), XO_threshold), linetype = "dashed", size = linesize, color = "red", alpha = 1/2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", size = linesize, color = "red", alpha = 1/2) +
  facet_grid(Day ~ .) +  
  scale_x_continuous(breaks = c(0, 0.2, 0.5, 0.8, 1)) +
  labs(x = expression("X Chromosome Ratio [XR = "*B6[chrX]*"/("*B6[chrX]+Cast[chrX]*")"), 
       y = expression(Xist[AS]*" ["*log[10]*"(CPM + 1)]")) +
  theme(strip.text.y = element_text(angle = 0))
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, savefile = paste0(outpath, "S1_C_XOfiltering.pdf"), height = 10)



print("2.4) D: Alignment and gene quantification - spliced/unspliced quantification")

print("2.4.1) Load data")
load(paste0(datapath, "UF_DGE_Spliced.RData")); spliced <- dge
load(paste0(datapath, "UF_DGE_Unspliced.RData")); unspliced <- dge
load(paste0(datapath, "UF_DGE_B6_Spliced.RData")); b6_spliced <- dge
load(paste0(datapath, "UF_DGE_B6_Unspliced.RData")); b6_unspliced <- dge
load(paste0(datapath, "UF_DGE_Cast_Spliced.RData")); cast_spliced <- dge
load(paste0(datapath, "UF_DGE_Cast_Unspliced.RData")); cast_unspliced <- dge

summary_stats <- data.frame(id = b6_spliced$samples$id, 
                            day = b6_spliced$samples$day,
                            total = b6_spliced$samples$seqdepth, 
                            ua = b6_spliced$samples$uniquealigned, 
                            spliced = colSums(spliced$counts),
                            unspliced = colSums(unspliced$counts),
                            spliced_b6 = colSums(b6_spliced$counts), 
                            unspliced_b6 = colSums(b6_unspliced$counts), 
                            spliced_cast = colSums(cast_spliced$counts), 
                            unspliced_cast = colSums(cast_unspliced$counts))

print("2.4.2) Compute median number and percentage of reads across samples")
sum <- summary_stats %>%
  dplyr::group_by(day) %>% 
  dplyr::summarise(seqdepth = median(total), 
                   ua = median(ua), 
                   spliced = median(spliced),
                   spliced_B6 = median(spliced_b6),
                   spliced_Cast = median(spliced_cast),
                   unspliced = median(unspliced),
                   unspliced_B6 = median(unspliced_b6),
                   unspliced_Cast = median(unspliced_cast))
features <- colnames(summary_stats)[!colnames(summary_stats) %in% c("id", "day")]
perc <- data.frame(summary_stats[, features]/summary_stats$total*100, day = summary_stats$day)
perc <- perc %>%
  dplyr::group_by(day) %>% 
  dplyr::summarise(seqdepth = median(total), 
                   ua = median(ua), 
                   spliced = median(spliced),
                   spliced_B6 = median(spliced_b6),
                   spliced_Cast = median(spliced_cast),
                   unspliced = median(unspliced),
                   unspliced_B6 = median(unspliced_b6),
                   unspliced_Cast = median(unspliced_cast))

print("2.4.3) Melt variables")
sum_melt <- tidyr::gather(sum, 'variable', 'value', -day)
perc_melt <- tidyr::gather(perc, 'variable', 'value', -day)
perc_melt$value <- paste0(round(perc_melt$value, digits = 1), "%")
sum_melt$perc <- perc_melt$value

newlev <- c("seqdepth" = "Sequencing\nDepth", 
            "ua" = "Uniquely\nAligned", 
            "spliced" = "UMI\nSpliced", "spliced_B6" = "UMI\nSpliced\nB6", "spliced_Cast" = "UMI\nSpliced\nCast",
            "unspliced" = "UMI\nUnspliced", "unspliced_B6" = "UMI\nUnspliced\nB6", "unspliced_Cast" = "UMI\nUnspliced\nCast")
sum_melt$variable <- revalue(sum_melt$variable, replace = newlev)
perc_melt$variable <- revalue(perc_melt$variable, replace = newlev)
sum_melt$variable <- factor(sum_melt$variable, levels = newlev)

print("2.4.4) Plot")
features <- as.character(unique(sum_melt$variable)[grepl(unique(sum_melt$variable), pattern = "UMI")])
s <- sum_melt[sum_melt$variable %in% features,]
s <- data.frame(s, strsplit2(s$variable, split = "\n")[,-1])
s$X2 <- as.character(s$X2); s$X2[s$X2 == ""] <- "not-AS UMI"
newlev <- c("not-AS UMI" = "not-AS UMI", "B6" = "B6 UMI", "Cast" = "Cast UMI")
s$X2 <- factor(revalue(factor(s$X2), replace = newlev), levels = newlev)

g <- s %>% 
  ggplot() + 
  theme_bw() + theme1 + 
  facet_grid(.~X1) +
  geom_bar(aes(x = X2, y = value, fill = factor(day)), stat = "identity", position=position_dodge()) + 
  geom_text(aes(x = X2, y = value, group = factor(day), label=perc), hjust = -.1,
            color="black", position = position_dodge(0.9), size = geomtext_size, angle = 90, show.legend = FALSE) +
  scale_fill_grey() +
  scale_y_continuous(label = scientific_label, breaks = seq(0, 6e5, by = 5e4), limits = c(0, 1.75e5)) + 
  labs(x = "", y = "Number of UMI counts", fill = "Time [days]")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "S1_D_SeqOutput_SplicedUnspliced.pdf"), width = 10)



print("3) Supplementary figure 3: MAST DE analyses and Correlation results")


print("3.1) UpSetR plots for MAST DE analyses")

print("3.1.1) Load data")
load(paste0(path, "output/fig3_deXistHighLow/ALLresults.RData"))

print("3.1.2) Subset results to significantly DE genes - excluding Xist")
de_threshold <- 0.05
all_de <- all_de[!all_de$mgi_symbol %in% "Xist",]
all_de$isX <- ifelse(all_de$chromosome_name %in% "X", "X-linked", "Autosomal")
all_de$direction <- factor(ifelse(all_de$coef < 0, "Down-regulated", "Up-regulated"), levels = c("Up-regulated", "Down-regulated"))
sig_de <- all_de[(all_de$fdr <= de_threshold),]
de_results <- sig_de[sig_de$test %in% "XistHigh_XistLow",]
de_results$Time <- paste0(de_results$day*24, "h")

print("3.1.3) Define significance of each gene throughout time")
de_upset <- ddply(de_results, .variables = .(mgi_symbol, isX), summarize,
                  de_24h = ifelse(length(fdr[Time == "24h"]) > 0, 
                                  ifelse(fdr[Time == "24h"] < de_threshold, 1, 0), 0),
                  de_48h = ifelse(length(fdr[Time == "48h"]) > 0, 
                                  ifelse(fdr[Time == "48h"] < de_threshold, 1, 0), 0),
                  de_72h = ifelse(length(fdr[Time == "72h"]) > 0, 
                                  ifelse(fdr[Time == "72h"] < de_threshold, 1, 0), 0),
                  de_96h = ifelse(length(fdr[Time == "96h"]) > 0, 
                                  ifelse(fdr[Time == "96h"] < de_threshold, 1, 0), 0),
                  direction = ifelse(all(coef > 0), "Up-regulated",
                                     ifelse(all(coef < 0), "Down-regulated", "Mixed")))

print("3.1.4) Plot - UpSetR for X-linked genes")
maxbar <- 60
x <- de_upset[de_upset$isX == "X-linked", ]
a <- colSums(x[x$direction == "Up-regulated", grepl(x = colnames(x), pattern = "^de_")])
b <- colSums(x[x$direction == "Down-regulated", grepl(x = colnames(x), pattern = "^de_")])
metadata <- data.frame(sets = names(a),
                       Up = as.numeric(a),
                       Down = as.numeric(b))
metadata$sets <- factor(metadata$sets, levels = paste0("de_", seq(24, 96, 24), "h"))
pdf(file = paste0(outpath, "S3_MASTupset_Xlinked.pdf"),
    width = 4, height = 5, onefile = TRUE, useDingbats = FALSE)
upset(x, main.bar.color = "blue",
      sets = paste0("de_", rev(seq(48, 96, 24)), "h"),
      order.by = "degree",
      group.by = "degree",
      keep.order = TRUE,
      decreasing = FALSE,
      mb.ratio = c(0.7, 0.3),
      point.size = scattersize*10, line.size = linesize*2, text.scale = 1.5,
      queries = list(list(query = elements, params = list("direction", c("Up-regulated")),color = "red", active = T)),
      mainbar.y.max = maxbar, mainbar.y.label = "Number of X-linked DE genes\n(Xist: High VS Low cells)\n(Xist excluded)\n"
)
dev.off()

# number of up/down-regulated genes per time point
vars <- c("de_24h", "de_48h", "de_72h", "de_96h")
y <- reshape2::melt(x, id.vars = colnames(x)[!colnames(x) %in% vars], measure.vars = vars)
y <- ddply(y, .variables = .(mgi_symbol), transform, 
           id = paste0(variable[value==1], collapse = "."))
y <- y[y$value ==1,]; ddply(y, .variables = .(id), summarize, 
                            up = length(unique(mgi_symbol[direction == "Up-regulated"])),
                            down = length(unique(mgi_symbol[direction == "Down-regulated"])))

print("3.1.4) Plot - UpSetR for Autosomal genes")
maxbar <- 300
x <- de_upset[de_upset$isX == "Autosomal", ]
a <- colSums(x[x$direction == "Up-regulated", grepl(x = colnames(x), pattern = "^de_")])
b <- colSums(x[x$direction == "Down-regulated", grepl(x = colnames(x), pattern = "^de_")])
metadata <- data.frame(sets = names(a),
                       Up = as.numeric(a),
                       Down = as.numeric(b))
metadata$sets <- factor(metadata$sets, levels = paste0("de_", seq(24, 96, 24), "h"))

pdf(file = paste0(outpath, "S3_MASTupset_Autosomal.pdf"),
    width = 4, height = 5, onefile = TRUE, useDingbats = FALSE)
upset(x, main.bar.color = "blue",
      sets = paste0("de_", rev(seq(24, 96, 24)), "h"),
      order.by = "degree",
      group.by = "degree",
      keep.order = TRUE,
      decreasing = FALSE,
      mb.ratio = c(0.7, 0.3),
      point.size = scattersize*10, line.size = linesize*2, text.scale = 1.5,
      queries = list(
        list(query = elements, params = list("direction", c("Up-regulated")),
             color = "red", active = T)),
      mainbar.y.max = maxbar, mainbar.y.label = "Number of Autosomal DE genes\n(Xist: High VS Low cells)\n"
)
dev.off()

# number of up/down-regulated genes per time point
vars <- c("de_24h", "de_48h", "de_72h", "de_96h")
y <- reshape2::melt(x, id.vars = colnames(x)[!colnames(x) %in% vars], measure.vars = vars)
y <- ddply(y, .variables = .(mgi_symbol), transform, 
           id = paste0(variable[value==1], collapse = "."))
y <- y[y$value ==1,]; ddply(y, .variables = .(id), summarize, 
                            up = length(unique(mgi_symbol[direction == "Up-regulated"])),
                            down = length(unique(mgi_symbol[direction == "Down-regulated"])))


print("3.2) Correlation analysis (Xist CPM) - known Xist regulators")

print("3.2.1) Load data")
load(paste0(path, "output/fig3_deXistHighLow/sprmcor.RData"))

print("3.2.2) Subset to known regulators")
regulators <- c("Nanog", "Klf2", "Klf4", "Prdm14", "Pou5f1", "Sox2", "Myc", "Zfp42", "Ctcf", "Tsix", "Yy1", "Rlim", "Ftx")
temp <- sprmcor[sprmcor$Gene %in% regulators,] %>% as.data.frame()
temp$detrate <- 1-temp$droprate
temp$isX <- ifelse(temp$Chr %in% "X", "X-linked", "Autosomal")
mypalette <- colorRampPalette(c("blue", "white", "red"))(1e3)

print("3.2.3) Gather data frame")
features <- c("day", "Gene", "isX", "detrate", "spr_Xist", "fdr_Xist")
temp_m <- gather(temp[, features], variable, value, 
                 -day, -Gene, -isX, -detrate, -fdr_Xist)
temp_m$group <- revalue(factor(temp_m$variable), replace = c("spr_Xist" = "Xist corr. (all cells)"))
temp_m$group <- factor(temp_m$group, levels = c("Xist corr. (all cells)"))

print("3.2.4) Identify genes with significant correlations")
temp_m$significant <- temp_m$fdr_Xist<0.05
temp_m <- temp_m %>% dplyr::arrange(fdr_Xist)
temp_m$Gene <- factor(temp_m$Gene, levels = rev(regulators))

print("3.2.5) Plot")
# compute plot size by # of genes
ngenes <- length(unique(temp_m$Gene))
ph <- ngenes*0.3
# plot
g <- temp_m[!temp_m$day %in% 0,] %>% 
  ggplot() + 
  theme_bw() + theme1 +
  facet_grid(.~group, scales = "free_y") + 
  geom_point(aes(x = factor(day), y = Gene, fill = value, size = abs(value)), pch=22, stroke = 0) +
  geom_point(data = temp_m[temp_m$significant == TRUE,], aes(x = factor(day), y = Gene), size = small_scattersize, color = "white", pch = 1) +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0, breaks=seq(-0.5, 0.5, .25), limits=c(-0.51,0.51)) + 
  scale_radius(breaks = seq(0, 0.5, by = 0.25), limits = c(0, 0.5), range = c(0.1, 4)) +
  labs(x = "Time [days]", y = "", fill = "Spearman's correlation to Xist CPM", size = "Absolute correlation")
adjust_size(g = g, panel_width_cm = 2.5, panel_height_cm = ph, savefile = paste0(outpath, "S3_SprmCorrToXistCPM_KnownRegulators.pdf"))



print("4) Supplementary figure 4: B6/Cast ratio for Autosomal and X-linked genes over time")

print("4.1) A: B6/Cast ratio for Autosomal (all cells) and X-linked (Xist-Undetected cells) genes over time")

print("4.1.1) Load data")
load(paste0(datapath, "DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "DGE_Cast.RData")); cast <- dge
df_as <- data.frame(day = rep(b6$samples$day, each = nrow(b6)), 
                    Cell = rep(colnames(b6), each = nrow(b6)), 
                    Xist = rep(b6$samples$Xist_ratio_class, each = nrow(b6)),
                    B6_els = rep(b6$samples$eff_libsize_notX, each = nrow(b6)),
                    Cast_els = rep(cast$samples$eff_libsize_notX, each = nrow(cast)),
                    Gene = rep(rownames(b6), times = ncol(b6)),
                    Chromosome = rep(b6$genes$chromosome, times = ncol(b6)),
                    b6 = c(b6$counts), cast = c(cast$counts))
df_as$both <- df_as$b6 + df_as$cast
df_as$isX <- ifelse(df_as$Chromosome %in% "X", "X-linked", "Autosomal")
df_as$Xist_classification <- revalue(factor(df_as$Xist), replace = c("Undetected" = "Xist\nUndetected", "Low-Xist" = "Low",
                                                                     "Middle" = "Skewed", "Xist_BA" = "BA",
                                                                     "BL6_MA" = "Xist-MA\n(Xi=B6)", "Cast_MA" = "Xist-MA\n(Xi=Cast)"))
lev_ref <- c("Xist\nUndetected", "Low", "Skewed", "BA", "Xist-MA\n(Xi=B6)", "Xist-MA\n(Xi=Cast)")
df_as$Xist_classification <- factor(df_as$Xist_classification, levels = lev_ref)

print("4.1.2) Compute gene-wise B6/Cast ratios for each Xist group")
temp <- df_as %>%
  dplyr::group_by(day, Xist_classification, Gene, isX) %>% 
  dplyr::summarise(b6_cast_ratio = (sum(b6) + 0.01)/(sum(cast) + 0.01))
summary(temp$b6_cast_ratio)
temp_median <- temp %>%
  dplyr::group_by(day, Xist_classification, isX) %>% dplyr::summarise(med = median(b6_cast_ratio)) %>%
  dplyr::arrange(Xist_classification)
temp$group <- ifelse(temp$isX == "Autosomal", "Autosomal", paste0("X-linked\n", temp$Xist_classification))

print("4.1.3) Plot")
# setting y-limits to zoom in without changing boxplot -> NA values 
remove_Xistgroups <- c("Low", "Skewed", "BA")
g <- temp[temp$group %in% c("X-linked\nXist\nUndetected", "Autosomal"),] %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = hvalpha, size = linesize) + 
  geom_boxplot(aes(x = group, y = log2(b6_cast_ratio), color = factor(day)), size = violin_box_size,
               outlier.size = outliersize, outlier.alpha = 1/4, outlier.shape = NA,
               position = position_dodge2(preserve = "single")) + 
  scale_color_manual(values=time_colors) +
  labs(y = expression("["*sum[g]*"("*B6[gc]*") + 0.01] / ["*sum[g]*"("*Cast[gc]*") + 0.01] ["*log[2]*"(value)]"), x = "", fill = "",
       title = "Gene-wise log2(total B6/total Cast) expression", color = "Time [days]") +
  scale_y_continuous(breaks = seq(-10, 10, by = 1), limits = c(-3, 3))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "S4_A_B6CastRatio_XistclassDay.pdf"))

print("4.2) B: B6/Cast ratio for X-linked genes (grouping Xist-Undetected, B6-MA or Cast-MA cells) over time")

print("4.2.1) Plot")
# setting y-limits to zoom in without changing boxplot -> NA values 
Xlinked_groups <- c("X-linked\nXist\nUndetected", "X-linked\nXist-MA\n(Xi=B6)", "X-linked\nXist-MA\n(Xi=Cast)")
g <- temp[temp$group %in% Xlinked_groups,] %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = hvalpha, size = linesize) + 
  geom_boxplot(aes(x = group, y = log2(b6_cast_ratio), color = factor(day)), size = violin_box_size,
               outlier.size = outliersize, outlier.alpha = 1/4, outlier.shape = NA,
               position = position_dodge2(preserve = "single")) + 
  scale_color_manual(values=time_colors) +
  labs(y = expression("["*sum[g]*"("*B6[gc]*") + 0.01] / ["*sum[g]*"("*Cast[gc]*") + 0.01] ["*log[2]*"(value)]"), x = "", fill = "",
       title = "Gene-wise log2(total B6/total Cast) expression", color = "Time [days]") +
  scale_y_continuous(breaks = seq(-10, 10, by = 5), limits = c(-7, 12))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "S4_B_B6CastRatio_XistclassDay_Xlinked.pdf"))
