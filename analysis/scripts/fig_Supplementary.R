print("1) Define output path")
outpath <- paste0(path, "output/fig_Supplementary/"); dir.create(path = outpath, recursive = T, showWarnings = F)









print("2) Supplementary figure 1: Alignment, Gene counting and Cell filtering")

print("2.1) A: Alignment and gene quantification - notAS")

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
              position=position_jitter(width = .05), shape = 21) +
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
              position=position_jitter(width = .05), shape = 21) +
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
  scale_y_continuous(label = scientific_label, breaks = seq(0, 6e5, by = 5e4), limits = c(0, 1.9e5)) + 
  labs(x = "", y = "Number of UMI counts", fill = "Time [days]")
adjust_size(g = g, panel_width_cm = 3.5, panel_height_cm = 2.5, savefile = paste0(outpath, "S1_D_SeqOutput_SplicedUnspliced.pdf"), width = 10)


print("2.5) E: Gene quantification - AS")

print("2.5.1) Load AS data")
load(paste0(datapath, "UF_DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "UF_DGE_Cast.RData")); cast <- dge
n_b6 <- colSums(b6$counts>0); n_cast <- colSums(cast$counts>0)
temp <- data.frame(day = b6$samples$day,
                   Cell = colnames(b6),
                   n_b6, n_cast)

print("2.5.2) Melt variables")
temp_melt <- melt(temp, id.vars = c("day", "Cell"), measure.vars = c("n_b6", "n_cast"))
temp_melt$variable <- revalue(factor(temp_melt$variable), replace = c("n_b6" = "B6", 
                                                                      "n_cast" = "Cast"))

print("2.5.3) Plot")

g <- temp_melt %>% 
  ggplot() +
  theme_bw() +  theme1 + 
  facet_grid(.~variable) +
  geom_violin(aes(x = factor(day), y = value), 
              alpha = 0.5, draw_quantiles = c(0.5), size = violin_box_size, show.legend = FALSE) + 
  geom_jitter(aes(x = factor(day), y = value), 
              alpha = 0.5, size = outliersize,
              position=position_jitter(width = .05), shape = 16) +
  labs(x="Time [days]", y = "# Detected genes")
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 2.5, savefile = paste0(outpath, "S1_E_AS_DetectedGenes.pdf"))









print("3) Supplementary figure 2: Marker gene expression, notAS and AS X:A ratios")

print("3.1) A: Compare expression of marker genes over time between Xist+ and Xist- cells")

print("3.1.1) Load data")
load(paste0(datapath, "DGE.RData"))
genes <- c("Nanog", "Esrrb", "Dnmt3a")
dge <- dge[dge$genes$symbol %in% genes,]
df <- data.frame(day = rep(dge$samples$day, each = nrow(dge)), 
                 Cell = rep(colnames(dge), each = nrow(dge)), 
                 Xist = rep(dge$samples$Xist_class, each = nrow(dge)),
                 els = rep(dge$samples$eff_libsize_notX, each = nrow(dge)),
                 Gene = rep(rownames(dge), times = ncol(dge)),
                 Chromosome = rep(dge$genes$chromosome, times = ncol(dge)),
                 umi = c(dge$counts))
df$cpm_value <- (df$umi/df$els)*1e6

print("3.1.2) Test difference between Xist- and Xist+ cells")
test <- df[df$day>0,] %>% dplyr::group_by(day, Gene) %>% 
  dplyr::summarise(w_pvalue = t.test(x = cpm_value[Xist == "Undetected"],
                                          y = cpm_value[Xist == "Detected (Xist UMI > 5)"])$p.value)
test$pvalue <- ifelse(test$w_pvalue >= 0.01, paste0("p = ", round(test$w_pvalue, digits = 2)), 
                    ifelse(test$w_pvalue >= 0.001, paste0("p = ", round(test$w_pvalue, digits = 3)),
                           "p < 0.001"))

print("2.5.3) Plot")
cols <- c("black", "#e7298a")
temp <- df[!df$Xist %in% "Detected (Xist UMI <= 5)",]
temp$Xist <- factor(as.character(temp$Xist), levels = c("Undetected", "Detected (Xist UMI > 5)"))
temp$Xist <- revalue(factor(temp$Xist), replace = c("Undetected" = "Xist- cells",
                                                    "Detected (Xist UMI > 5)" = "Xist+ cells"))

g <- temp %>%  
  ggplot() +
  theme_bw() + theme1 + 
  facet_grid(.~Gene) +
  geom_boxplot(aes(x = factor(day), y = log10(cpm_value + 1), color = Xist), outlier.shape = NA, size = violin_box_size) + 
  geom_jitter(aes(x = factor(day), y = log10(cpm_value + 1), color = Xist), fill = "black",
              alpha = 1/4, position=position_jitterdodge(jitter.width = .25, dodge.width = 0.75), 
              size = outliersize, show.legend = FALSE) +
  scale_color_manual(values = cols) + 
  geom_text(data = test, aes(x = factor(day), y = 3.65, label = pvalue), size = geomtext_size, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 3.75), breaks = seq(0, 3.5, 0.5)) +
  labs(x = "Time [days]", y = expression(log[10]*"(CPM + 1) values"), color = "Xist not-AS classification")
adjust_size(g = g, panel_width_cm = 4.5, panel_height_cm = 3, savefile = paste0(outpath, "S2_A_Xist_Regulators.pdf"), width = 8, height = 2)


print("3.2) B: X:A ratio varying Xist UMI threshold")

print("3.2.1) Load data")
load(paste0(datapath, "DGE.RData"))
nboot <- 1e3
xcounts <- dge$counts[dge$genes$chromosome %in% "X",]
autcounts <- dge$counts[dge$genes$chromosome %in% c(1:19),]
n_xlinked <- nrow(xcounts); n_autlinked <- nrow(autcounts)

print("3.2.2) Define Xist thresholds and compute bootstrapped X:A ratios")
Xist_thr <- c(0, 5, 10, 25); Xist_thr <- Xist_thr[order(Xist_thr, decreasing = F)]
temp_boot <- c()
for (j in seq_len(nboot)) {
  boot <- sample(x = seq_len(n_autlinked), size = n_xlinked, replace = TRUE)
  value <- colSums(xcounts)/colSums(autcounts[boot,])
  temp_boot <- rbind(temp_boot, value)
}
x2a_boot <- apply(temp_boot, 2, median)
x2a <- data.frame(day = dge$samples$day, 
                  id = colnames(dge), 
                  Xist_UMI = dge$counts["Xist",], 
                  x2a = x2a_boot)

print("3.2.3) Store results")
x2a_sub <- data.frame(x2a[x2a$Xist_UMI == 0,], group = "Xist = 0")
for(j in seq_len(length(Xist_thr))){
  thr <- Xist_thr[j]
  x2a_sub <- rbind(x2a_sub,
                   data.frame(x2a[x2a$Xist_UMI > thr,], group = paste0("Xist > ", thr)))
}
summ_cell <- ddply(x2a_sub, .variables = .(day, group), summarize, n = length(day))

print("2.6.4) Plot")
features <- c("id", "day", "group", "x2a")
total <- x2a_sub
total$day <- factor(total$day); summ_cell$day <- factor(summ_cell$day)
cols <- c("black", 
          adjustcolor("#e7298a", alpha.f = 1/5), adjustcolor("#e7298a", alpha.f = 1/3),
          adjustcolor("#e7298a", alpha.f = 1/2), adjustcolor("#e7298a", alpha.f = 1)
          )

g <- total %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = seq(0.5, 2, by = 0.5), linetype = "dashed", alpha = hvalpha, size = linesize) + 
  scale_y_continuous(breaks = seq(0.5, 2, by = 0.5)) + 
  geom_boxplot(aes(x = day, y = x2a, color = group), outlier.shape = NA, size = violin_box_size) + 
  geom_jitter(aes(x = day, y = x2a, color = group), fill = "black", alpha = 1/10,
              position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75), 
              size = outliersize, show.legend = FALSE) +
  geom_text(summ_cell,
            mapping = aes(x = day, y = Inf, label = paste0("n = ", n), fill = group),
            position=position_dodge(.75), angle = 90, alpha = 1, size = geomtext_size, hjust = -0.5) +
  coord_cartesian(clip = "off") +
  labs(x = "Time [days]", y = "X:A ratio", color = "") +
  scale_color_manual(values = cols)
adjust_size(g = g, panel_width_cm = 10, panel_height_cm = 3, 
            savefile = paste0(outpath, "S2_B_x2a_XistnotASthreshold.pdf"), width = 7, height = 3)


print("3.3) C: Allele specific X:A ratio for Xist-MA cells vs Scaled Pseudotime")

print("3.3.1) Load AS data")
load(paste0(datapath, "DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "DGE_Cast.RData")); cast <- dge
xcounts_b6 <- b6[b6$genes$chromosome %in% "X",]$counts
autcounts_b6 <- b6[b6$genes$chromosome %in% c(1:19),]$counts
xcounts_cast <- cast[cast$genes$chromosome %in% "X",]$counts
autcounts_cast <- cast[cast$genes$chromosome %in% c(1:19),]$counts
n_xlinked <- nrow(xcounts_b6); n_autlinked <- nrow(autcounts_b6)

print("3.3.2) Compute AS bootstrapped X:A ratios")
b6_XiAi <- cast_XiAi <- c()
nboot <- 1000
for (j in seq_len(nboot)) {
  boot <- sample(x = seq_len(n_autlinked), size = n_xlinked, replace = TRUE)
  
  # define values
  X_b6 <- colSums(xcounts_b6)
  X_cast <- colSums(xcounts_cast)
  A_b6 <- colSums(autcounts_b6[boot,])
  A_cast <- colSums(autcounts_cast[boot,])
  
  # compute ratios
  b6_XiAi <- rbind(b6_XiAi, X_b6/A_b6)
  cast_XiAi <- rbind(cast_XiAi, X_cast/A_cast)
}
l <- list(b6_XiAi = b6_XiAi, 
          cast_XiAi = cast_XiAi)
x2a_boot <- lapply(l, function(x) apply(x, 2, median))
x2a_boot <- do.call(cbind, x2a_boot)
x2a <- data.frame(day = dge$samples$day,
                  id = colnames(dge),
                  xist_classification = b6$samples$Xist_ratio_class, 
                  x2a_boot)

print("3.3.3) Select only Xist_MA cells and define Xi and Xa accordingly")
temp <- x2a[x2a$xist_classification %in% c("BL6_MA", "Cast_MA"),]
temp$Xi_Ai <- ifelse(temp$xist_classification == "BL6_MA", temp$b6_XiAi, temp$cast_XiAi)
temp$Xa_Aa <- ifelse(temp$xist_classification == "BL6_MA", temp$cast_XiAi, temp$b6_XiAi)
temp$xist_classification <- plyr::revalue(factor(temp$xist_classification, levels = c("BL6_MA", "Cast_MA")), 
                                          replace = c("BL6_MA" = "Xist-MA (Xi=B6)", 
                                                      "Cast_MA" = "Xist-MA (Xi=Cast)"))
temp_melt <- melt(temp, id.vars = c("day", "id", "xist_classification"), measure.vars = c("Xi_Ai", "Xa_Aa"))
temp_melt$variable <- factor(temp_melt$variable, levels = c("Xi_Ai", "Xa_Aa"))
temp_melt$pdt <- pData(XX)$Scaled_PDT[match(temp_melt$id, pData(XX)$id)]

print("3.3.4) Plot")
g <- temp_melt %>%
  ggplot() + 
  theme_bw() + theme1 + 
  facet_grid(variable ~ xist_classification) + 
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 1/2, size = linesize) + 
  geom_point(aes(x = pdt, y = value, color = factor(day)), size = scattersize, alpha = 1/2, shape = 20) + 
  scale_color_manual(values = time_colors) + 
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  guides(color=guide_legend(override.aes = list(size=2, alpha = 1))) + 
  theme(legend.position = "right", 
        strip.text.y = element_text(angle = 0)) +
  labs(x="Scaled Pseudotime", y = "Allele Specific X:A ratio", color = "Time [days]")
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, 
            savefile = paste0(outpath, "S2_C_asX2A_PDT.pdf"))










print("4) Supplementary figure 4: MAST DE analyses and Correlation results")


print("4.1) UpSetR plots for MAST DE analyses")

print("4.1.1) Load data")
load(paste0(path, "output/fig3_deXistHighLow/ALLresults.RData"))

print("4.1.2) Subset results to significantly DE genes - excluding Xist")
de_threshold <- 0.05
all_de <- all_de[!all_de$mgi_symbol %in% "Xist",]
all_de$isX <- ifelse(all_de$chromosome_name %in% "X", "X-linked", "Autosomal")
all_de$direction <- factor(ifelse(all_de$coef < 0, "Down-regulated", "Up-regulated"), levels = c("Up-regulated", "Down-regulated"))
sig_de <- all_de[(all_de$fdr <= de_threshold),]
de_results <- sig_de[sig_de$test %in% "XistHigh_XistLow",]
de_results$Time <- paste0(de_results$day*24, "h")

print("4.1.3) Define significance of each gene throughout time")
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

print("4.1.4) Plot - UpSetR for X-linked genes")
maxbar <- 60
x <- de_upset[de_upset$isX == "X-linked", ]
a <- colSums(x[x$direction == "Up-regulated", grepl(x = colnames(x), pattern = "^de_")])
b <- colSums(x[x$direction == "Down-regulated", grepl(x = colnames(x), pattern = "^de_")])
metadata <- data.frame(sets = names(a),
                       Up = as.numeric(a),
                       Down = as.numeric(b))
metadata$sets <- factor(metadata$sets, levels = paste0("de_", seq(24, 96, 24), "h"))
pdf(file = paste0(outpath, "S4_A_MASTupset_Xlinked.pdf"),
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

print("4.1.4) Plot - UpSetR for Autosomal genes")
maxbar <- 300
x <- de_upset[de_upset$isX == "Autosomal", ]
a <- colSums(x[x$direction == "Up-regulated", grepl(x = colnames(x), pattern = "^de_")])
b <- colSums(x[x$direction == "Down-regulated", grepl(x = colnames(x), pattern = "^de_")])
metadata <- data.frame(sets = names(a),
                       Up = as.numeric(a),
                       Down = as.numeric(b))
metadata$sets <- factor(metadata$sets, levels = paste0("de_", seq(24, 96, 24), "h"))

pdf(file = paste0(outpath, "S4_A_MASTupset_Autosomal.pdf"),
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


print("4.2) Correlation analysis (Xist CPM) - known Xist regulators")

print("4.2.1) Load data")
load(paste0(path, "output/fig3_deXistHighLow/sprmcor.RData"))

print("4.2.2) Subset to known regulators")
regulators <- c("Nanog", "Klf2", "Klf4", "Prdm14", "Pou5f1", "Sox2", "Myc", "Zfp42", "Ctcf", "Tsix", "Yy1", "Rlim", "Ftx")
temp <- sprmcor[sprmcor$Gene %in% regulators,] %>% as.data.frame()
temp$detrate <- 1-temp$droprate
temp$isX <- ifelse(temp$Chr %in% "X", "X-linked", "Autosomal")
mypalette <- colorRampPalette(c("blue", "white", "red"))(1e3)

print("4.2.3) Gather data frame")
features <- c("day", "Gene", "isX", "detrate", "spr_Xist", "fdr_Xist")
temp_m <- gather(temp[, features], variable, value, 
                 -day, -Gene, -isX, -detrate, -fdr_Xist)
temp_m$group <- revalue(factor(temp_m$variable), replace = c("spr_Xist" = "Xist corr. (all cells)"))
temp_m$group <- factor(temp_m$group, levels = c("Xist corr. (all cells)"))

print("4.2.4) Identify genes with significant correlations")
temp_m$significant <- temp_m$fdr_Xist<0.05
temp_m <- temp_m %>% dplyr::arrange(fdr_Xist)
temp_m$Gene <- factor(temp_m$Gene, levels = rev(regulators))

print("4.2.5) Plot")
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
adjust_size(g = g, panel_width_cm = 2.5, panel_height_cm = ph, savefile = paste0(outpath, "S4_B_SprmCorrToXistCPM_KnownRegulators.pdf"))


print("4.3) not-AS RNA-velocity normalized change in chrX expression in Xist+ and Xist- cells")

print("4.3.1) Load data")
f1path <- paste0(path, "output/fig1_NotAS/")
load(paste0(datapath, "NCF_DGE_Spliced.RData")); spliced <- dge
load(paste0(datapath, "DGE.RData")); notas <- dge
load(paste0(f1path, "notAS_vel.RData"))

print("4.3.2) Compute RNA-velocity predicted change in normalized X-chromosome expression")

p <- notAS_vel$projected; c <- notAS_vel$current
colnames(p) <- apply(cbind(as.numeric(gsub(strsplit2(colnames(p), split = "\\_")[,1], pattern = "h", replacement = ""))/24,
                           strsplit2(colnames(p), split = "\\_")[, -1]), 1, function(x) paste0("d", paste0(x, collapse = "_")))
colnames(c) <- apply(cbind(as.numeric(gsub(strsplit2(colnames(c), split = "\\_")[,1], pattern = "h", replacement = ""))/24,
                           strsplit2(colnames(c), split = "\\_")[, -1]), 1, function(x) paste0("d", paste0(x, collapse = "_")))
m <- match(rownames(p), rownames(spliced)); table(is.na(m))
genes <- spliced$genes[m,]; genes_x <- genes[genes$chromosome %in% "X",]
p_xchr <- p[genes$chromosome %in% "X",]; c_xchr <- c[genes$chromosome %in% "X",]; table(genes$chromosome %in% "X")
diff_xchr <- log2(colSums(c_xchr)/colSums(p_xchr))
m <- match(names(diff_xchr), rownames(notas$samples)); table(is.na(m))
samples <- notas$samples[m,]

print("4.3.3) Define data frame")

df <- data.frame(day = samples$day, 
                 id = rownames(samples),
                 Xist_AS = samples$Xist_ratio_class,
                 Xist_notAS = samples$Xist_class,
                 diff = diff_xchr)
lv_notas <- c("Undetected" = "Xist- cells", 
              "Detected (Xist UMI <= 5)" = "Xist+ cells [UMI<=5]",
              "Detected (Xist UMI > 5)" = "Xist+ cells")
lv_as <- c("Undetected" = "Undetected", "Low-Xist" = "Low",
        "Middle" = "Skewed", "Xist_BA" = "BA",
        "BL6_MA" = "Xist-MA (Xi=B6)", "Cast_MA" = "Xist-MA (Xi=Cast)")
df$xist_notas <- factor(revalue(factor(df$Xist_notAS), replace = lv_notas), levels = lv_notas)
df$xist_as <- factor(revalue(factor(df$Xist_AS), replace = lv_as), levels = lv_as)

print("4.3.4) Plot")
g <- df[!df$xist_notas %in% "Xist+ cells [UMI<=5]",] %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = hvalpha, size = linesize) +
  geom_boxplot(aes(x = factor(day), y = diff, color = factor(xist_notas)), outlier.shape = NA, size = violin_box_size) + 
  geom_jitter(aes(x = factor(day), y = diff, color = factor(xist_notas)), fill = "black",
              alpha = 1/4, position=position_jitterdodge(jitter.width = .25, dodge.width = 0.75), 
              size = outliersize, show.legend = FALSE) +
  scale_color_manual(values = rev(comparison_colors)) +
  labs(x = "Time [days]", 
       y = expression(log[2]*"( "*Xchr[Current]*" / "*Xchr[Predicted]*" )"),
       color = "Xist not-AS classification")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "S4_C_XchrChange_XistnotAS.pdf"), 
            width = 5, height = 2)



print("4.4) UMAP embedding colored by Xchr predicted change in Xist- and Xist+ cells")

print("4.4.1) Load data")
load(paste0(f1path, "umap.RData"))
xpred <- colSums(p_xchr); xcurr <- colSums(c_xchr)
m <- match(umap$id, names(xpred)); table(is.na(m))
umap <- data.frame(umap,
                   xpred = xpred[m],
                   xcurr = xcurr[m])

print("4.4.2) Plot")
temp <- umap[!umap$xist_detection %in% "Xist+ cells [UMI<=5]",]
mid <- median(temp$xpred[temp$day == 0]); print(paste0("The median normalized predicted X-chromosome change in d0 cells is ", round(mid, digits = 3)))
g <- temp %>%
  ggplot() + 
  theme_bw() + theme1 +
  facet_grid(.~xist_detection) +
  geom_point(aes(x = umap1, y = umap2, color = xpred), size = outliersize, shape = 20) + 
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "grey", high = "red") +
  labs(x="UMAP 1", y = "UMAP 2", color = "RNA-velocity predicted variation\nin chrX gene expression") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4)) + 
  theme(legend.position = "right",
        strip.text.y = element_text(angle = 0))
adjust_size(g = g, panel_width_cm = 2.5, panel_height_cm = 2.5, savefile = paste0(outpath, "S4_D_UMAP_XchrPredicted.pdf"))


print("4.5) UpSetR plots for MAST DE analyses - Xchr change")

print("4.5.1) Load data")
load(paste0(path, "output/fig4_deXchrChangeHighLow/XchrChange_HighLow_MAST.RData"))

print("4.5.2) Subset results to significantly DE genes")
de_threshold <- 0.05
all_de$isX <- ifelse(all_de$chromosome_name %in% "X", "X-linked", "Autosomal")
all_de$direction <- factor(ifelse(all_de$coef < 0, "Down-regulated", "Up-regulated"), levels = c("Up-regulated", "Down-regulated"))
sig_de <- all_de[(all_de$fdr <= de_threshold),]
de_results <- sig_de[sig_de$test %in% "HighXchrChange_LowXchrChange",]
de_results$Time <- paste0(de_results$day*24, "h")

print("4.5.3) Define significance of each gene throughout time")
de_upset <- ddply(de_results, .variables = .(mgi_symbol, isX), summarize,
                  de_0h = ifelse(length(fdr[Time == "0h"]) > 0, 
                                 ifelse(fdr[Time == "0h"] < de_threshold, 1, 0), 0),
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

print("4.5.4) Plot - UpSetR for X-linked genes")
maxbar <- 20
x <- de_upset[de_upset$isX == "X-linked", ]
a <- colSums(x[x$direction == "Up-regulated", grepl(x = colnames(x), pattern = "^de_")])
b <- colSums(x[x$direction == "Down-regulated", grepl(x = colnames(x), pattern = "^de_")])
metadata <- data.frame(sets = names(a),
                       Up = as.numeric(a),
                       Down = as.numeric(b))
metadata$sets <- factor(metadata$sets, levels = paste0("de_", seq(0, 96, 24), "h"))
pdf(file = paste0(outpath, "S4_E_MASTupsetXchrChange_Xlinked.pdf"),
    width = 4, height = 5, onefile = TRUE, useDingbats = FALSE)
upset(x, main.bar.color = "blue",
      sets = paste0("de_", rev(seq(24, 96, 24)), "h"),
      order.by = "degree",
      group.by = "degree",
      keep.order = TRUE,
      decreasing = FALSE,
      mb.ratio = c(0.7, 0.3),
      point.size = scattersize*10, line.size = linesize*2, text.scale = 1.5,
      queries = list(list(query = elements, params = list("direction", c("Up-regulated")),color = "red", active = T)),
      mainbar.y.max = maxbar, mainbar.y.label = "Number of X-linked DE genes\n(Xchr-Change: High VS Low cells)\n"
)
dev.off()

# number of up/down-regulated genes per time point
vars <- c("de_0h", "de_24h", "de_48h", "de_72h", "de_96h")
y <- reshape2::melt(x, id.vars = colnames(x)[!colnames(x) %in% vars], measure.vars = vars)
y <- ddply(y, .variables = .(mgi_symbol), transform, 
           id = paste0(variable[value==1], collapse = "."))
y <- y[y$value ==1,]; ddply(y, .variables = .(id), summarize, 
                            up = length(unique(mgi_symbol[direction == "Up-regulated"])),
                            down = length(unique(mgi_symbol[direction == "Down-regulated"])))

print("4.1.4) Plot - UpSetR for Autosomal genes")
maxbar <- 300
x <- de_upset[de_upset$isX == "Autosomal", ]
a <- colSums(x[x$direction == "Up-regulated", grepl(x = colnames(x), pattern = "^de_")])
b <- colSums(x[x$direction == "Down-regulated", grepl(x = colnames(x), pattern = "^de_")])
metadata <- data.frame(sets = names(a),
                       Up = as.numeric(a),
                       Down = as.numeric(b))
metadata$sets <- factor(metadata$sets, levels = paste0("de_", seq(0, 96, 24), "h"))

pdf(file = paste0(outpath, "S4_E_MASTupsetXchrChange_Autosomal.pdf"),
    width = 4, height = 5, onefile = TRUE, useDingbats = FALSE)
upset(x, main.bar.color = "blue",
      sets = paste0("de_", rev(seq(48, 96, 24)), "h"),
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


print("4.6) Correlation analysis (Xchr change) - known Xist regulators")

print("4.6.1) Load data")
load(paste0(path, "output/fig4_deXchrChangeHighLow/cor_Xchr.RData"))

print("4.6.2) Subset to known regulators")
regulators <- c("Nanog", "Klf2", "Klf4", "Prdm14", "Pou5f1", "Sox2", "Myc", "Zfp42", "Ctcf", "Tsix", "Yy1", "Rlim", "Ftx")
temp <- sprmcor[sprmcor$Gene %in% regulators,] %>% as.data.frame()
temp$detrate <- 1-temp$droprate
temp$isX <- ifelse(temp$Chr %in% "X", "X-linked", "Autosomal")
mypalette <- colorRampPalette(c("blue", "white", "red"))(1e3)

print("4.6.3) Gather data frame")
features <- c("day", "Gene", "isX", "detrate", "spr", "fdr")
temp_m <- gather(temp[, features], variable, value, 
                 -day, -Gene, -isX, -detrate, -fdr)
temp_m$group <- revalue(factor(temp_m$variable), replace = c("spr" = "Xchr corr."))
temp_m$group <- factor(temp_m$group, levels = c("Xchr corr."))

print("4.6.4) Identify genes with significant correlations")
temp_m$significant <- temp_m$fdr <= de_threshold
temp_m <- temp_m %>% dplyr::arrange(fdr)
temp_m$Gene <- factor(temp_m$Gene, levels = rev(regulators))

print("4.6.5) Plot")
# compute plot size by # of genes
ngenes <- length(unique(temp_m$Gene))
ph <- ngenes*0.3
# plot
g <- temp_m %>% 
  ggplot() + 
  theme_bw() + theme1 +
  facet_grid(.~group, scales = "free_y") + 
  geom_point(aes(x = factor(day), y = Gene, fill = value, size = abs(value)), pch=22, stroke = 0) +
  geom_point(data = temp_m[temp_m$significant == TRUE,], aes(x = factor(day), y = Gene), size = small_scattersize, color = "white", pch = 1) +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0, breaks=seq(-0.5, 0.5, .25), limits=c(-0.51,0.51)) + 
  scale_radius(breaks = seq(0, 0.5, by = 0.25), limits = c(0, 0.5), range = c(0.1, 4)) +
  labs(x = "Time [days]", y = "", fill = "Pearson's gene-wise correlation:\nrho(CPM_gene, X-chromosome change)", size = "Absolute correlation")
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = ph, savefile = paste0(outpath, "S4_F_PearsCorrToXistCPM_KnownRegulators.pdf"))









print("5) Supplementary figure 5: B6/Cast ratio for Autosomal and X-linked genes over time")

print("5.1) A: B6/Cast ratio for Autosomal (all cells) and X-linked (Xist-Undetected cells) genes over time")

print("5.1.1) Load data")
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

print("5.1.2) Compute gene-wise B6/Cast ratios for each Xist group")
temp <- df_as %>%
  dplyr::group_by(day, Xist_classification, Gene, isX) %>% 
  dplyr::summarise(b6_cast_ratio = (sum(b6) + 0.01)/(sum(cast) + 0.01))
summary(temp$b6_cast_ratio)
temp_median <- temp %>%
  dplyr::group_by(day, Xist_classification, isX) %>% dplyr::summarise(med = median(b6_cast_ratio)) %>%
  dplyr::arrange(Xist_classification)
temp$group <- ifelse(temp$isX == "Autosomal", "Autosomal", paste0("X-linked\n", temp$Xist_classification))

print("5.1.3) Plot")
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
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "S5_A_B6CastRatio_XistclassDay.pdf"))

print("5.2) B: B6/Cast ratio for X-linked genes (grouping Xist-Undetected, B6-MA or Cast-MA cells) over time")

print("5.2.1) Plot")
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
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "S5_B_B6CastRatio_XistclassDay_Xlinked.pdf"))


print("5.3) C: Differential silencing analysis - varying XP threshold and number of bins")

print("5.3.1) Load function")
XiXa_metacell_ratio_final <- function(x = df_sk, XistMAgroup = c("BL6_MA", "Cast_MA"), idplot = "B6", useCounts = TRUE,
                                      nbins = 5, mincount = 20, mincells = 5, includeXistUnd_firstbin = TRUE,
                                      offset_XiXa = 0.1, logXiXa = TRUE, #XiXa computing
                                      bsl_count = TRUE, XistUnd_XCRthr_low = 0.4, XistUnd_XCRthr_high = 0.6, XistUndTime = "0hrs", 
                                      outpath = report_path, kseed = 0, kiter = 10^5, minbin = 3,
                                      normalize2baseline = TRUE, mincount_cellfit = 5, mincell_cellfit = 10, mincount_cellfit_bin=5,
                                      percentageXinactive = TRUE, equallysized = TRUE, minsilencing_perc = 10, zeroIntercept = FALSE,
                                      weighted = TRUE){
  
  print(paste0("0) subset to Xist-MA cells only & define inactive counts..."))
  id <- paste0(x$Time, "_", x$Cell)
  temp <- data.frame(id, x)
  data <- temp[temp$Xist %in% XistMAgroup,]
  data <- data[data$both > 0,]
  
  print(paste0("0bis) compute X silencing percentage..."))
  data <- data[!data$Gene %in% c("Xist"),] %>% #don't include Xist in XCR
    dplyr::group_by(id) %>%
    dplyr::mutate(b6_X = sum(b6[Chromosome == "X"]),
                  cast_X = sum(cast[Chromosome == "X"]))
  data$AXCR <- (1-((data$b6_X + offset_XiXa)/(data$cast_X + offset_XiXa)))*100
  data$AXCR[data$Xist == "Cast_MA"] <- (1-((data$cast_X[data$Xist == "Cast_MA"] + offset_XiXa)/(data$b6_X[data$Xist == "Cast_MA"] + offset_XiXa)))*100
  
  print(paste0("0ter) exclude cells with silencing below ", minsilencing_perc, "%..."))
  data <- data[data$AXCR >= minsilencing_perc,]
  
  print(paste0("0quat) breaks for percentage of silencing"))
  if(!equallysized){
    break_XCR <- seq(0, 100, length = nbins+1)
    levels_bin <- levels(cut(seq(0, 100, by = 1), breaks = break_XCR, include.lowest = TRUE))
    data$bin_XCR <- cut(data$AXCR, breaks = break_XCR, include.lowest = TRUE)
  }else{
    data$bin_XCR <- quantileCut(x = data$AXCR, n = nbins)
    levels_bin <- levels(data$bin_XCR)
  }
  
  print(paste0("1) define metacell group based on nbins and AXCR and store # of cells..."))
  data <- ddply(data, .variables = .(bin_XCR), transform, bin_ncells = length(unique(id)))
  d1 <- daply(data, .variables = .(Time, bin_XCR), summarize, n = length(unique(id))); print(d1)
  d2 <- ddply(data, .variables = .(bin_XCR), summarize, n = length(unique(id))); print(d2)
  
  print(paste0("2) Define active and inactive allele..."))
  data$inactive_count <- data$b6*I(data$Xist == "BL6_MA") + data$cast*I(data$Xist == "Cast_MA")
  data$active_count <- data$cast*I(data$Xist == "BL6_MA") + data$b6*I(data$Xist == "Cast_MA")
  
  print(paste0("3) compute bins' ratios and XT values, store # cells, and filter out metacell if too small..."))
  data_cell <- data
  data <- data %>%
    dplyr::group_by(Xist, bin_XCR, Gene) %>%
    dplyr::summarise(chr = unique(Chromosome),
                     ncells = length(unique(id)),
                     sumXi_count = sum(inactive_count),
                     sumXa_count = sum(active_count)) %>%
    dplyr::group_by(Xist, bin_XCR) %>%
    dplyr::mutate(sumXi_X = sum(sumXi_count[chr=="X"]),
                  sumXa_X = sum(sumXa_count[chr=="X"]))
  data$AXCR <- (1-((data$sumXi_X + offset_XiXa)/(data$sumXa_X + offset_XiXa)))*100
  data$XiXa <- (data$sumXi_count + offset_XiXa)/(data$sumXa_count + offset_XiXa)
  
  print(paste0("4) Filter out gene expression defined by few counts, compute bins' Xi/Xa and XCR values..."))
  condition <- ((data$sumXi_count + data$sumXa_count) >= mincount) & (data$ncells >= mincells)
  data_filt <- data[condition,]
  
  print(paste0("5) Store the baseline and compute Xi/Xa gene-skewing..."))
  bsl <- temp[(temp$Xist %in% "Undetected") & 
                (temp$Time %in% XistUndTime) & 
                (temp$XCR >= XistUnd_XCRthr_low) & 
                (temp$XCR <= XistUnd_XCRthr_high),]
  bsl <- bsl[(bsl$b6 + bsl$cast) > 0,]
  
  baseline <- bsl %>% 
    dplyr::group_by(Gene) %>%
    dplyr::summarise(n = length(unique(id)),
                     sum_b6_count = sum(b6), 
                     sum_cast_count = sum(cast))
  baseline$b6cast_ratio <- (baseline$sum_b6_count + offset_XiXa)/(baseline$sum_cast_count + offset_XiXa)
  condition <- ((baseline$sum_b6_count + baseline$sum_cast_count) >= mincount) & (baseline$n >= mincells); table(condition)
  baseline$b6cast_ratio[!condition] <- NA
  
  print(paste0("6) normalize XiXa ratios to baseline..."))
  m <- match(data_filt$Gene, baseline$Gene); table(is.na(m))
  data_filt$baseline_XiXa <- baseline$b6cast_ratio[m]
  data_filt <- data_filt[!is.na(data_filt$baseline_XiXa),]
  data_filt$baseline_XiXa[data_filt$Xist == "Cast_MA"] <- 1/data_filt$baseline_XiXa[data_filt$Xist == "Cast_MA"]
  data_filt$XiXa_ratio_XistUnd <- data_filt$XiXa/data_filt$baseline_XiXa
  
  ### single cell
  print(paste0("6bis) normalize cell-wise XiXa ratios to baseline (only cells with UMI>", mincount_cellfit, "..."))
  data_cell <- data_cell[data_cell$both>mincount_cellfit_bin,]
  data_cell$XiXa <- (data_cell$inactive_count + offset_XiXa)/(data_cell$active_count + offset_XiXa)
  m <- match(data_cell$Gene, baseline$Gene); table(is.na(m))
  data_cell$baseline_XiXa <- baseline$b6cast_ratio[m]
  data_cell <- data_cell[!is.na(data_cell$baseline_XiXa),]
  data_cell$baseline_XiXa[data_cell$Xist == "Cast_MA"] <- 1/data_cell$baseline_XiXa[data_cell$Xist == "Cast_MA"]
  data_cell$XiXa_ratio_XistUnd <- data_cell$XiXa/data_cell$baseline_XiXa
  
  print(paste0("6ter) filter cell-wise values separately"))
  both_genes <- data_cell %>% dplyr::group_by(Gene) %>% dplyr::summarise(n = length(unique(Xist)))
  data_cell_fit <- data_cell[data_cell$Gene %in% as.character(both_genes$Gene[both_genes$n==2]),]
  select_genes <- data_cell_fit %>% 
    dplyr::group_by(Xist, Gene) %>% dplyr::summarise(ncells = sum(both>mincount_cellfit)) %>%
    dplyr::group_by(Gene) %>% dplyr::summarise(n_b6 = ncells[Xist=="BL6_MA"], n_cast = ncells[Xist=="Cast_MA"])
  which_genes <- as.character(select_genes$Gene[(select_genes$n_b6 >= mincell_cellfit)&(select_genes$n_cast >= mincell_cellfit)])
  data_cell_fit <- data_cell_fit[(data_cell_fit$Gene %in% which_genes)&(data_cell_fit$both > mincount_cellfit),]
  ### single cell
  
  print(paste0("7) take log ratios"))
  data_filt$logXiXa <- log2(data_filt$XiXa_ratio_XistUnd)
  data_cell_fit$logXiXa <- log2(data_cell_fit$XiXa_ratio_XistUnd)
  
  print(paste0("7bis) Remove genes that don't have at least ", minbin, " observations per Xist-MA group..."))
  data_filt <- data.frame(data_filt)
  repres <- data_filt %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(n_b6 = length(unique(bin_XCR[Xist == "BL6_MA"])), n_cast = length(unique(bin_XCR[Xist == "Cast_MA"])))
  
  repres$keep_b6 <- repres$n_b6 >= minbin; print(table(minbin_B6 = repres$keep_b6))
  repres$keep_cast <- repres$n_cast >= minbin; print(table(minbin_Cast = repres$keep_cast))
  repres$keep <- repres$keep_b6 & repres$keep_cast; print(table(minbin_Both = repres$keep))
  
  test_b6 <- as.character(repres$Gene[repres$keep_b6]); test_cast <- as.character(repres$Gene[repres$keep_cast])
  test_both <- intersect(test_b6, test_cast)
  res_either <- repres[repres$keep_b6 | repres$keep_cast, ]; res_keep <- repres[repres$keep,]
  data_filt <- data_filt[data_filt$Gene %in% as.character(unique(res_either$Gene)),]
  
  
  print(paste0("8) Fit linear model with free intercept for each gene separating fits for Xist-MA Cast/B6 cells..."))
  allgenes <- unique(c(as.character(data_filt$Gene), as.character(data_cell_fit$Gene)))
  # allgenes <- unique(as.character(unique(res_either$Gene))); common_genes <- unique(as.character(res_keep$Gene))
  
  model_fit <- model_fit_sc <- c()
  if(zeroIntercept){
    
    for(i in seq_len(length(allgenes))){
      gene <- allgenes[i]
      
      temp <- data_filt[data_filt$Gene %in% gene,]
      temp$w <- temp$sumXa_count + temp$sumXi_count
      
      temp_sc <- data_cell_fit[data_cell_fit$Gene %in% gene,]
      temp_sc$w <- temp_sc$inactive_count + temp_sc$active_count
      
      # bin analysis
      if((length(unique(temp$Xist))==2) & (gene %in% test_both)){
        
        if(weighted){
          
          complete <- lm(logXiXa ~ 0 + AXCR + AXCR:Xist, weights = w, data = temp)
          oneslope <- lm(logXiXa ~ 0 + AXCR, weights = w, data = temp)
          
          int_b6 <- 0; slope_b6 <- sum(complete$coefficients[c("AXCR", "AXCR:XistBL6_MA")])
          int_cast <- 0; slope_cast <- complete$coefficients[c("AXCR")]
          
          anv <- anova(oneslope, complete)
          slope_pvalue <- anv$`Pr(>F)`[2]
          
        }else{
          complete <- lm(logXiXa ~ 0 + AXCR + AXCR:Xist, data = temp)
          oneslope <- lm(logXiXa ~ 0 + AXCR, data = temp)
          
          int_b6 <- 0; slope_b6 <- sum(complete$coefficients[c("AXCR", "AXCR:XistBL6_MA")])
          int_cast <- 0; slope_cast <- complete$coefficients[c("AXCR")]
          
          anv <- anova(oneslope, complete)
          slope_pvalue <- anv$`Pr(>F)`[2]
        }
      }else{
        int_cast <- int_b6 <- slope_cast <- slope_b6 <- slope_pvalue <- NA
      }
      
      #single allele analysis
      
      # test b6
      if(gene %in% test_b6){
        if(weighted){
          fit <- lm(logXiXa ~ 0 + AXCR, weights = w, data = temp[temp$Xist == "BL6_MA",])
        }else{
          fit <- lm(logXiXa ~ 0 + AXCR, data = temp[temp$Xist == "BL6_MA",])
        }
        int_b6 <- 0; slope_b6 <- fit$coefficients["AXCR"]; bic_b6 <- BIC(fit)
      }else{
        int_b6 <- slope_b6 <- bic_b6 <- NA
      }
      
      # test cast
      if(gene %in% test_cast){
        if(weighted){
          fit <- lm(logXiXa ~ 0 + AXCR, weights = w, data = temp[temp$Xist == "Cast_MA",])
        }else{
          fit <- lm(logXiXa ~ 0 + AXCR, data = temp[temp$Xist == "Cast_MA",])
        }
        int_cast <- 0; slope_cast <- fit$coefficients["AXCR"]; bic_cast <- BIC(fit)
      }else{
        int_cast <- slope_cast <- bic_cast <- NA
      }
      
      # single cell analysis
      if(length(unique(temp_sc$Xist))==2){
        
        if(weighted){
          
          complete <- lm(logXiXa ~ 0 + AXCR + AXCR:Xist, weights = w, data = temp_sc)
          oneslope <- lm(logXiXa ~ 0 + AXCR, weights = w, data = temp_sc)
          
          int_b6_sc <- 0; slope_b6_sc <- sum(complete$coefficients[c("AXCR", "AXCR:XistBL6_MA")])
          int_cast_sc <- 0; slope_cast_sc <- complete$coefficients[c("AXCR")]
          
          anv <- anova(oneslope, complete)
          slope_pvalue_sc <- anv$`Pr(>F)`[2]
          
        }else{
          complete <- lm(logXiXa ~ 0 + AXCR + AXCR:Xist, data = temp_sc)
          oneslope <- lm(logXiXa ~ 0 + AXCR, data = temp_sc)
          
          int_b6_sc <- 0; slope_b6_sc <- sum(complete$coefficients[c("AXCR", "AXCR:XistBL6_MA")])
          int_cast_sc <- 0; slope_cast_sc <- complete$coefficients[c("AXCR")]
          
          anv <- anova(oneslope, complete)
          slope_pvalue_sc <- anv$`Pr(>F)`[2]
        }
      }else{
        int_cast_sc <- int_b6_sc <- slope_cast_sc <- slope_b6_sc <- slope_pvalue_sc <- NA
      }
      
      # store results
      model_fit <- rbind(model_fit, 
                         data.frame(Gene = gene, 
                                    int_b6 = int_b6, slope_b6 = slope_b6,
                                    int_cast = int_cast, slope_cast = slope_cast,
                                    slope_pvalue = slope_pvalue))
      model_fit_sc <- rbind(model_fit_sc, 
                            data.frame(Gene = gene, 
                                       int_b6 = int_b6_sc, slope_b6 = slope_b6_sc,
                                       int_cast = int_cast_sc, slope_cast = slope_cast_sc,
                                       slope_pvalue = slope_pvalue_sc))
    }
    
  }else{
    
    for(i in seq_len(length(allgenes))){
      gene <- allgenes[i]
      
      temp <- data_filt[data_filt$Gene %in% gene,]
      temp$w <- temp$sumXa_count + temp$sumXi_count
      
      temp_sc <- data_cell_fit[data_cell_fit$Gene %in% gene,]
      temp_sc$w <- temp_sc$inactive_count + temp_sc$active_count
      
      # bin analysis
      if(length(unique(temp$Xist))==2){
        
        if(weighted){
          
          complete <- lm(logXiXa ~ AXCR*Xist, weights = w, data = temp)
          oneslope <- lm(logXiXa ~ AXCR, weights = w, data = temp)
          
          int_b6 <- complete$coefficients["(Intercept)"]
          slope_b6 <- complete$coefficients["AXCR"] 
          int_cast <- sum(complete$coefficients[c("(Intercept)", "XistCast_MA")])
          slope_cast <- sum(complete$coefficients[c("AXCR", "AXCR:XistCast_MA")])
          anv <- anova(oneslope, complete)
          slope_pvalue <- anv$`Pr(>F)`[2]
          
        }else{
          complete <- lm(logXiXa ~ AXCR*Xist, data = temp)
          oneslope <- lm(logXiXa ~ AXCR, data = temp)
          
          int_b6 <- complete$coefficients["(Intercept)"]
          slope_b6 <- complete$coefficients["AXCR"] 
          int_cast <- sum(complete$coefficients[c("(Intercept)", "XistCast_MA")])
          slope_cast <- sum(complete$coefficients[c("AXCR", "AXCR:XistCast_MA")])
          anv <- anova(oneslope, complete)
          slope_pvalue <- anv$`Pr(>F)`[2]
        }
      }else{
        int_cast <- int_b6 <- slope_cast <- slope_b6 <- slope_pvalue <- NA
      }
      
      #single allele analysis
      
      # test b6
      if(gene %in% test_b6){
        if(weighted){
          fit <- lm(logXiXa ~ AXCR, weights = w, data = temp[temp$Xist == "BL6_MA",])
        }else{
          fit <- lm(logXiXa ~ AXCR, data = temp[temp$Xist == "BL6_MA",])
        }
        int_b6 <- fit$coefficient["(Intercept)"]; slope_b6 <- fit$coefficients["AXCR"]; bic_b6 <- BIC(fit)
      }else{
        int_b6 <- slope_b6 <- bic_b6 <- NA
      }
      
      # test cast
      if(gene %in% test_cast){
        if(weighted){
          fit <- lm(logXiXa ~ AXCR, weights = w, data = temp[temp$Xist == "Cast_MA",])
        }else{
          fit <- lm(logXiXa ~ AXCR, data = temp[temp$Xist == "Cast_MA",])
        }
        int_cast <- fit$coefficient["(Intercept)"]; slope_cast <- fit$coefficients["AXCR"]; bic_cast <- BIC(fit)
      }else{
        int_cast <- slope_cast <- bic_cast <- NA
      }
      
      # single cell analysis
      if(length(unique(temp_sc$Xist))==2){
        
        if(weighted){
          
          complete <- lm(logXiXa ~ AXCR*Xist, weights = w, data = temp_sc)
          oneslope <- lm(logXiXa ~ AXCR, weights = w, data = temp_sc)
          
          int_b6_sc <- complete$coefficients["(Intercept)"]
          slope_b6_sc <- complete$coefficients["AXCR"] 
          int_cast_sc <- sum(complete$coefficients[c("(Intercept)", "XistCast_MA")])
          slope_cast_sc <- sum(complete$coefficients[c("AXCR", "AXCR:XistCast_MA")])
          anv <- anova(oneslope, complete)
          slope_pvalue_sc <- anv$`Pr(>F)`[2]
          
        }else{
          complete <- lm(logXiXa ~ AXCR*Xist, data = temp_sc)
          oneslope <- lm(logXiXa ~ AXCR, data = temp_sc)
          
          int_b6_sc <- complete$coefficients["(Intercept)"]
          slope_b6_sc <- complete$coefficients["AXCR"] 
          int_cast_sc <- sum(complete$coefficients[c("(Intercept)", "XistCast_MA")])
          slope_cast_sc <- sum(complete$coefficients[c("AXCR", "AXCR:XistCast_MA")])
          anv <- anova(oneslope, complete)
          slope_pvalue_sc <- anv$`Pr(>F)`[2]
        }
      }else{
        int_cast_sc <- int_b6_sc <- slope_cast_sc <- slope_b6_sc <- slope_pvalue_sc <- NA
      }
      
      # store results
      model_fit <- rbind(model_fit, 
                         data.frame(Gene = gene, 
                                    int_b6 = int_b6, slope_b6 = slope_b6,
                                    int_cast = int_cast, slope_cast = slope_cast,
                                    slope_pvalue = slope_pvalue))
      model_fit_sc <- rbind(model_fit_sc, 
                            data.frame(Gene = gene, 
                                       int_b6 = int_b6_sc, slope_b6 = slope_b6_sc,
                                       int_cast = int_cast_sc, slope_cast = slope_cast_sc,
                                       slope_pvalue = slope_pvalue_sc))
    }
    
  }
  model_fit <- model_fit[order(model_fit$slope_pvalue, decreasing = FALSE),]
  model_fit$FDR <- p.adjust(model_fit$slope_pvalue, method = "BH")
  
  model_fit_sc <- model_fit_sc[order(model_fit_sc$slope_pvalue, decreasing = FALSE),]
  model_fit_sc$FDR <- p.adjust(model_fit_sc$slope_pvalue, method = "BH")
  
  print(paste0("9) Return results..."))
  return <- list(model_fit = model_fit,
                 model_fit_sc = model_fit_sc,
                 baseline = baseline, 
                 data = data_filt,
                 data_cell = data_cell, 
                 data_cell_fit = data_cell_fit)
  return(return)
}


print("5.3.2) Load data")
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
df_as$Time <- paste0(df_as$day*24, "hrs")
df_as$both <- df_as$b6 + df_as$cast
df_sk <- ddply(df_as[df_as$Chromosome == "X",], .variables = .(Time, Cell), transform, 
               Xist = unique(Xist),
               Xist_b6 = b6[Gene == "Xist"],
               Xist_b6_cpm = (b6[Gene == "Xist"]/unique(B6_els))*1e6,
               Xist_cast = cast[Gene == "Xist"],
               Xist_cast_cpm = (cast[Gene == "Xist"]/unique(Cast_els))*1e6,
               b6_X = sum(b6[Gene != "Xist"]),
               cast_X = sum(cast[Gene != "Xist"]))
df_sk$XCR <- df_sk$b6_X/(df_sk$b6_X + df_sk$cast_X)

print("5.3.3) Launch function varying XP threshold")
sigthreshold <- 0.05
mincount <- 25; nbins <- 10
mincells <- 5; minbin <- 5; minsilencing_perc <- 10; mincount_cellfit_bin <- 5
mincount_cellfit <- 5; mincell_cellfit <- 10
XistUndTime <- "0hrs"; XistUnd_XCRthr_low <- 0.4; XistUnd_XCRthr_high <- 0.6
equallysized <- TRUE; zeroIntercept <- TRUE; weighted <- FALSE

results <- c()
for(minsilencing_perc in seq(0, 20, by = 2)){
  binresults <- XiXa_metacell_ratio_final(x = df_sk, mincount = mincount, nbins = nbins, minbin = minbin, mincells = mincells,
                                          minsilencing_perc = minsilencing_perc, mincount_cellfit_bin = mincount_cellfit_bin,
                                          mincount_cellfit = mincount_cellfit, mincell_cellfit = mincell_cellfit,
                                          equallysized = equallysized, zeroIntercept = zeroIntercept, weighted = weighted)
  
  # store results
  x <- binresults$model_fit[!is.na(binresults$model_fit$FDR),]
  variables <- c("Gene", "slope_b6", "slope_cast", "slope_pvalue", "FDR")
  results <- rbind(results,
                   data.frame(x[, variables], silencing_threshold_percentage = minsilencing_perc))
}
results[results$FDR <= sigthreshold,] %>% arrange(silencing_threshold_percentage, FDR)

print("5.3.4) Plot")

# plot genes with 20 smallest FDR values
temp <- results %>% dplyr::group_by(Gene) %>% dplyr::summarise(minFDR = min(FDR[!is.na(FDR)])) %>% arrange(minFDR)
ngenes <- 20
topgenes <- as.character(temp$Gene[seq_len(ngenes)])
temp <- results[results$Gene %in% topgenes,]

# color by -log10(FDR)
temp$lfdr <- -log10(temp$FDR + 1e-6)
temp$significant <- temp$FDR <= sigthreshold

# plot
ph <- length(unique(temp$Gene))*0.3
temp$Gene <- factor(temp$Gene, levels = rev(topgenes))
temp$lev <- paste0("XP >= ", temp$silencing_threshold_percentage)
temp$lev <- factor(temp$lev, levels = unique(temp$lev))

g <- temp %>% 
  ggplot() + 
  theme_bw() + theme1 +
  geom_point(aes(x = lev, y = Gene, color = significant, size = abs(lfdr)), stroke = 0, shape = 20) +
  scale_color_manual(values = rev(c("red", "black"))) +
  scale_radius(breaks = seq(0, 6, by = 1), limits = c(0, 6), range = c(0.5, 6)) +
  labs(x = "", y = "", color = "Significance") +
  theme(axis.text.x = element_text(angle = 90))
adjust_size(g = g, panel_width_cm = 4, panel_height_cm = ph, 
            savefile = paste0(outpath, "S5_C_DSanalysis_varyXPthreshold.pdf"), 
            height = 10, width = 5)

print("5.3.5) Launch function varying number of bins")
sigthreshold <- 0.05
mincount <- 25; nbins <- 10
mincells <- 5; minbin <- 5; minsilencing_perc <- 10; mincount_cellfit_bin <- 5
mincount_cellfit <- 5; mincell_cellfit <- 10
XistUndTime <- "0hrs"; XistUnd_XCRthr_low <- 0.4; XistUnd_XCRthr_high <- 0.6
equallysized <- TRUE; zeroIntercept <- TRUE; weighted <- FALSE

results <- c()
for(nbins in seq(6, 20, by = 2)){
  binresults <- XiXa_metacell_ratio_final(x = df_sk, mincount = mincount, nbins = nbins, minbin = minbin, mincells = mincells,
                                          minsilencing_perc = minsilencing_perc, mincount_cellfit_bin = mincount_cellfit_bin,
                                          mincount_cellfit = mincount_cellfit, mincell_cellfit = mincell_cellfit,
                                          equallysized = equallysized, zeroIntercept = zeroIntercept, weighted = weighted)
  # store results
  x <- binresults$model_fit[!is.na(binresults$model_fit$FDR),]
  variables <- c("Gene", "slope_b6", "slope_cast", "FDR")
  results <- rbind(results,
                   data.frame(x[, variables], nbins = nbins))
}
results[results$FDR <= sigthreshold,] %>% arrange(nbins, FDR)

print("5.3.6) Plot")

# plot genes with 20 smallest FDR values
temp <- results %>% dplyr::group_by(Gene) %>% dplyr::summarise(minFDR = min(FDR[!is.na(FDR)])) %>% arrange(minFDR)
ngenes <- 20
topgenes <- as.character(temp$Gene[seq_len(ngenes)])
temp <- results[results$Gene %in% topgenes,]

# color by -log10(FDR)
temp$lfdr <- -log10(temp$FDR + 1e-6)
temp$significant <- temp$FDR <= sigthreshold

# plot
ph <- length(unique(temp$Gene))*0.3
temp$Gene <- factor(temp$Gene, levels = rev(topgenes))
temp$lev <- paste0("# Bins = ", temp$nbins)
temp$lev <- factor(temp$lev, levels = unique(temp$lev))

g <- temp %>% 
  ggplot() + 
  theme_bw() + theme1 +
  geom_point(aes(x = lev, y = Gene, color = significant, size = abs(lfdr)), stroke = 0, shape = 20) +
  scale_color_manual(values = rev(c("red", "black"))) +
  scale_radius(breaks = seq(0, 6, by = 1), limits = c(0, 6), range = c(0.5, 6)) +
  labs(x = "", y = "", color = "Significance") +
  theme(axis.text.x = element_text(angle = 90))
adjust_size(g = g, panel_width_cm = 4, panel_height_cm = ph, 
            savefile = paste0(outpath, "S5_D_DSanalysis_varyNbins.pdf"), 
            height = 10, width = 5)









print("6) Supplementary figure 7: Experimental validation (Pyrosequencing, qPCR, RNA-FISH) on dXic line")

print("6.1) A: Pyrosequencing - Xist B6-ratio on dXic lines")

print("6.1.1) Load data")
file <- paste0(datapath, "validation_experiments.txt")
data <- data.frame(read.table(file = file, header = T, sep = "\t"))
data$group <- revalue(as.factor(as.character(data$CellLine)), 
                      replace = c("XX" = "WT", 
                                  "TX1072" = "TX1072",
                                  "dXIC_B6" = "Xist-MA (Xi=Cast)",
                                  "dXIC_Cast" = "Xist-MA (Xi=B6)"))
data <- data[!is.na(data$Value),]

pyroseq <- data[(data$Experiment %in% "Pyrosequencing")&(!data$group %in% "WT"),]
hits <- c("Klhl13", "Pir", "Hprt")
pyroseq$is_hit <- ifelse(pyroseq$Gene %in% hits, "DS gene", "Control gene")
pyroseq$B6SNPfreq <- 100-pyroseq$Value

print("6.1.2) Plot")
g <- pyroseq[pyroseq$Gene %in% "Xist",] %>%
  ggplot(aes(x = factor(Day), y = B6SNPfreq, color = group, fill = group)) +
  theme_bw() + theme1 +  
  geom_hline(yintercept = 50, size = linesize, linetype = "dashed", alpha = 1/4) +
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75),
              size = small_scattersize, show.legend = FALSE, alpha = 1/2) +
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = color_alleles) +
  labs(x = "Time [days]", 
       y = "B6 molecules [%]", 
       color = "",
       title = "Xist (Pyrosequencing)") +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100))
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 2, 
            savefile = paste0(outpath, "S7_A_PyrodXic_Xist_B6ratio.pdf"), 
            width = 3, height = 2)


print("6.2) B: qPCR - Xist relative expression")

print("6.2.1) Load data")
qpcr <- data[(data$Experiment %in% "qPCR")&(!data$group %in%  "WT"),]
color_qpcr <- c(rev(color_alleles)[1:2], "WT" = "black")

print("6.2.2) Plot")
g <- qpcr[qpcr$Gene %in% "Xist",] %>%
  ggplot(aes(x = factor(Day), y = log2(Value), color = group, fill = group)) +
  theme_bw() + theme1 +  
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75),
              size = small_scattersize, show.legend = FALSE, alpha = 1/2) +
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = color_qpcr) +
  labs(x = "Time [days]", 
       y = "Relative expression [log2]", 
       color = "",
       title = "Xist (qPCR)") +
  theme(strip.text.y = element_text(angle = 0))
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 2, 
            savefile = paste0(outpath, "S7_B_qPCRdXic_Xist_RelExp.pdf"),
            width = 3, height = 2)



print("6.3) C: FISH - Xist+ cells in dXic lines")

print("6.3.1) Load data")
fish <- data[(data$Experiment %in% "FISH")&(!data$group %in% "TX1072"),]

print("6.3.2) Plot")
g <- fish %>%
  ggplot(aes(x = factor(Day), y = Value, color = group, fill = group)) +
  theme_bw() + theme1 +  
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75),
              size = small_scattersize, show.legend = FALSE, alpha = 1/2) +
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = color_alleles) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Time [days]", 
       y = "Cells [%]",
       color = "",
       title = "Xist (RNA-FISH)")
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 2, 
            savefile = paste0(outpath, "S7_C_FISHdXic.pdf"),
            width = 3, height = 2)


print("6.4) D: Pyrosequencing - Control genes: B6/Total")

print("6.4.1) Plot")
g <- pyroseq[(pyroseq$is_hit != "DS gene")&(!pyroseq$Gene %in% "Xist"),] %>%
  ggplot(aes(x = factor(Day), y = B6SNPfreq, color = group, fill = group)) +
  theme_bw() + theme1 +  
  facet_grid(.~Gene, scales = "free_y") +
  geom_hline(yintercept = 50, size = linesize, linetype = "dashed", alpha = 1/4) +
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75),
              size = small_scattersize, show.legend = FALSE, alpha = 1/2) +
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = color_alleles) +
  labs(x = "Time [days]", 
       y = "B6/Total [%]", 
       color = "") +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100))
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 2, 
            savefile = paste0(outpath, "S7_D_PyrodXic_ControlGenes_B6ratio.pdf"), 
            width = 10, height = 2)


print("6.5) E: Pyrosequencing - Control genes: normalized Xi:Xa values")

print("6.5.1) Load data")
temp <- pyroseq[(pyroseq$is_hit != "DS gene")&(!pyroseq$Gene %in% "Xist"),]
temp$Xiperc <- ifelse(temp$group == "Xist-MA (Xi=B6)", 
                      100-temp$Value, 
                      temp$Value)

print("6.5.2) Compute Xi:Xa ratios, normalize to d0 values, average across replicates")

# compute XiXa ratios
temp$XiXa <- temp$Xiperc/(100-temp$Xiperc)

# normalize to average 0h ratio per cell line
temp <- temp %>%
  dplyr::group_by(Gene, group) %>%
  dplyr::mutate(mean_baseline_0h = mean(XiXa[Day==0]),
                XiXa_scaled = XiXa/mean(XiXa[Day==0]))

# average across replicates
mr <- temp %>%
  dplyr::group_by(Day, group, Gene) %>% 
  dplyr::summarise(meanratio = mean(XiXa_scaled))

# paired test for each day between two cell lines
wt <- ddply(mr[mr$Day != 0,], .variables = .(Day), summarise, 
            wilcox_pvalue = wilcox.test(x = meanratio[group == "Xist-MA (Xi=Cast)"], 
                                        y = meanratio[group == "Xist-MA (Xi=B6)"], 
                                        paired = TRUE)$p.value)
wt$p <- ifelse(wt$wilcox_pvalue >= 0.01, paste0("p = ", round(wt$wilcox_pvalue, digits = 2)), 
               ifelse(wt$wilcox_pvalue >= 0.001, paste0("p = ", round(wt$wilcox_pvalue, digits = 3)),
                      "p < 0.001"))

print("6.5.3) Plot")
g <- mr  %>%
  ggplot(aes(x = factor(Day), y = log2(meanratio), color = group, fill = group)) +
  theme_bw() + theme1 +  
  geom_text(data = data.frame(wt, group = NA), 
            aes(x = factor(Day), y = 1, label = p),
            size = geomtext_size, color = "black") +
  geom_hline(yintercept = 0, size = linesize, linetype = "dashed", alpha = 1/4) +
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75),
              size = small_scattersize, show.legend = FALSE, alpha = 1/2) +
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = color_alleles) +
  labs(x = "Time [days]", 
       y = "norm. Xi/Xa [log2]",
       color = "",
       title = "5 X-linked genes") +
  scale_y_continuous(limits = c(-2.5, 1.5))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 2, 
            savefile = paste0(outpath, "S7_E_PyrodXic_ControlGenes_Avglog2NormRatios.pdf"),
            width = 5, height = 2)


print("6.6) F: Pyrosequencing - DS genes: B6/Total")

print("6.6.1) Plot")
g <- pyroseq[pyroseq$is_hit == "DS gene",] %>%
  ggplot(aes(x = factor(Day), y = B6SNPfreq, color = group, fill = group)) +
  theme_bw() + theme1 +  
  facet_grid(.~Gene, scales = "free_y") +
  geom_hline(yintercept = 50, size = linesize, linetype = "dashed", alpha = 1/4) +
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75),
              size = small_scattersize, show.legend = FALSE, alpha = 1/2) +
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = color_alleles) +
  labs(x = "Time [days]", 
       y = "B6/Total [%]", 
       color = "") +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100))
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 2, 
            savefile = paste0(outpath, "S7_F_PyrodXic_DSgenes_B6ratio.pdf"), 
            width = 10, height = 2)


print("6.7) G: Pyrosequencing - DS genes: normalized Xi:Xa values")

print("6.7.1) Load data")
temp <- pyroseq[(pyroseq$is_hit %in% "DS gene"),]
temp$Xiperc <- ifelse(temp$group == "Xist-MA (Xi=B6)", 
                      100-temp$Value, 
                      temp$Value)

print("6.7.2) Compute Xi:Xa ratios, normalize to d0 values, test for each gene")

# compute XiXa ratios
temp$XiXa <- temp$Xiperc/(100-temp$Xiperc)

# normalize to average 0h ratio per cell line
temp <- temp %>%
  dplyr::group_by(Gene, group) %>%
  dplyr::mutate(mean_baseline_0h = mean(XiXa[Day==0]),
                XiXa_scaled = XiXa/mean(XiXa[Day==0]))

# test
tt <- ddply(temp, .variables = .(Day, Gene), summarise, 
            t_pvalue = t.test(x = XiXa_scaled[group == "Xist-MA (Xi=Cast)"], 
                              y = XiXa_scaled[group == "Xist-MA (Xi=B6)"], 
                              paired = FALSE)$p.value) %>% arrange(Gene)
tt$p <- ifelse(tt$t_pvalue >= 0.01, paste0("p = ", round(tt$t_pvalue, digits = 2)), 
               ifelse(tt$t_pvalue >= 0.001, paste0("p = ", round(tt$t_pvalue, digits = 3)),
                      "p < 0.001"))

print("6.7.3) Plot")

g <- temp %>%
  ggplot(aes(x = factor(Day), y = log2(XiXa_scaled), color = group, fill = group)) +
  theme_bw() + theme1 +  
  facet_grid(.~Gene) +
  geom_text(data = data.frame(tt, group = NA), 
            aes(x = factor(Day), y = 1.75, label = p),
            size = geomtext_size, color = "black") +
  geom_hline(yintercept = 0, size = linesize, linetype = "dashed", alpha = 1/4) +
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75),
              size = small_scattersize, show.legend = FALSE, alpha = 1/2) +
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = color_alleles) +
  labs(x = "Time [days]", 
       y = "norm. Xi/Xa [log2]",
       color = "") +
  scale_y_continuous(limits = c(-2, 2))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 2, 
            savefile = paste0(outpath, "S7_G_PyrodXic_DSgenes_log2NormRatios.pdf"),
            width = 10, height = 2)









print("7) Supplementary figure 8: Bulk RNA-Sequencing of TX1072 and dXic lines")

print("7.1) A: Bulk RNA-Seq: dXic - Percentage of X-linked reads mapping to each allele")

print("7.1.1) Load data")
load(paste0(datapath, "DGE_dXic_B6.RData")); b6 <- dge
load(paste0(datapath, "DGE_dXic_Cast.RData")); cast <- dge
df_dXic <- data.frame(day = rep(b6$samples$day, each = nrow(b6)), 
                      sample = rep(colnames(b6), each = nrow(b6)),
                      cell_line = rep(b6$samples$dXic, each = nrow(b6)),
                      Chr = rep(b6$genes$chromosome, times = ncol(b6)),
                      Gene = rep(b6$genes$symbol, times = ncol(b6)),
                      Ensembl = rep(b6$genes$ensembl, times = ncol(b6)),
                      b6 = c(b6$counts), cast = c(cast$counts))
df_dXic <- df_dXic %>% 
  dplyr::group_by(sample) %>% 
  dplyr::mutate(libsize_b6 = sum(b6[!Gene %in% "Xist_5prime"]), libsize_cast = sum(cast[!Gene %in% "Xist_5prime"]),
                libsize_b6_Xist5prime = sum(b6[!Gene %in% "Xist"]), libsize_cast_Xist5prime = sum(cast[!Gene %in% "Xist"])) %>%
  as.data.frame()


print("7.1.2) Subset to X-linked genes removing the ones on dXic locus")
dXic_locus <- 103182257:103955531
dXic_genes <- bm[(bm$chromosome_name == "X")&(bm$start_position %in% dXic_locus | bm$end_position %in% dXic_locus),]
g <- unique(dXic_genes$mgi_symbol); e <- unique(dXic_genes$ensembl_gene_id)
df_noXic <- df_dXic[!df_dXic$Ensembl %in% e,]

print("7.1.3) Compute percentage of X-linked reads from each allele in two cell lines over time")
temp <- df_noXic %>% 
  dplyr::group_by(day, cell_line, sample) %>% 
  dplyr::summarise(b6_X = sum(b6[(Chr %in% "X")]),
                   b6_Tot = sum(b6),
                   cast_X = sum(cast[(Chr %in% "X")]),
                   cast_Tot = sum(cast)) %>% 
  as.data.frame()
temp$B6 <- (temp$b6_X/(temp$b6_Tot + temp$cast_Tot))*100
temp$Cast <- (temp$cast_X/(temp$b6_Tot + temp$cast_Tot))*100
temp$cl <- revalue(factor(temp$cell_line), replace = c("dB6" = "Xi = Cast [dXic on B6]",
                                                       "dCast" = "Xi = B6 [dXic on Cast]"))
temp_melt <- melt(temp, id.vars = c("day", "cell_line", "cl", "sample"), measure.vars = c("B6", "Cast"))
temp_melt$variable <- factor(temp_melt$variable, levels = c("Cast", "B6"))

print("7.1.4) Plot")
g <- temp_melt %>% 
  ggplot(aes(x = factor(day), y = value, color = variable)) + 
  theme_bw() + theme1 + 
  facet_grid(.~cl) +
  geom_jitter(alpha = 0.5, size = 1,
              position=position_jitterdodge(jitter.width = .1, dodge.width = 0.9), shape = 16) +
  guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.9), 
               size = linesize*2, show.legend = FALSE) +
  labs(x="Time [days]", 
       y = "X-chromosome\nAS expression [%]", 
       color = "Allele") +
  scale_color_manual(values = rev(c("#1b9e77", "#d95f02"))) + 
  scale_y_continuous(breaks = seq(0.5, 2.5, by = 0.5), limits = c(0.5, 2.5))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 2, 
            savefile = paste0(outpath, "S8_A_BulkdXic_XchrASpercs.pdf"), 
            height = 2, width = 7)


print("7.2) B: Bulk RNA-Seq: dXic - Xist gene expression using only 5'-SNP")

print("7.2.1) Identify Xist 5prime expression from allele without deletion")
temp <- df_dXic[df_dXic$Gene %in% "Xist_5prime",]
temp$Xist <- ifelse(temp$cell_line == "dB6", temp$cast, temp$b6)
temp$libsize <- ifelse(temp$cell_line == "dB6", temp$libsize_cast_Xist5prime, temp$libsize_b6_Xist5prime)

print("7.2.2) Test two cell lines for each time point")
temp$Xist_cpm <- (temp$Xist/temp$libsize)*1e6
test <- temp %>% dplyr::group_by(day) %>% 
  dplyr::summarise(w_pvalue = t.test(x = Xist_cpm[cell_line == "dB6"],
                                     y = Xist_cpm[cell_line == "dCast"])$p.value) %>%
  as.data.frame()
test$pvalue <- ifelse(test$w_pvalue >= 0.01, paste0("p = ", round(test$w_pvalue, digits = 2)), 
                      ifelse(test$w_pvalue >= 0.001, paste0("p = ", round(test$w_pvalue, digits = 3)),
                             "p < 0.001"))
temp$cell_line <- revalue(temp$cell_line, replace = c("dB6" = "Xi = Cast [dXic on B6]",
                                                      "dCast" = "Xi = B6 [dXic on Cast]"))
temp$cell_line <- factor(temp$cell_line, levels = c("Xi = Cast [dXic on B6]",
                                                    "Xi = B6 [dXic on Cast]"))

print("7.2.3) Plot")
cols <- c("#1b9e77", "#d95f02")
g <- temp %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_jitter(aes(x = factor(day), y = log10(Xist_cpm + 1), color = cell_line), fill = "black",
              alpha = 1/2, position=position_jitterdodge(jitter.width = .25, dodge.width = 0.75), 
              size = scattersize, show.legend = FALSE) +
  stat_summary(fun=mean, aes(x = factor(day), y = log10(Xist_cpm + 1), ymin=..y.., ymax=..y.., color = cell_line), 
               geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = rev(cols)) + 
  geom_text(data = test, aes(x = factor(day), y = 3.5, label = pvalue), size = geomtext_size, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 3.75), breaks = seq(0, 4, 1)) +
  labs(x = "Time [days]", 
       y = expression("5'-Xist CPM + 1 [ "*log[10]*" ]"), 
       color = "",
       title = "dXic: RNA-Seq - Xist 5' expression")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 2, 
            savefile = paste0(outpath, "S8_B_BulkdXic_Xist5prime.pdf"), 
            height = 2, width = 4)



print("7.3) C: Bulk RNA-Seq: TX1072 - Xist gene expression")

print("7.3.1) Load data")
load(paste0(datapath, "DGE_TX1072_B6.RData")); b6 <- dge
load(paste0(datapath, "DGE_TX1072_Cast.RData")); cast <- dge
df_wt <- data.frame(day = rep(b6$samples$day, each = nrow(b6)), 
                    sample = rep(colnames(b6), each = nrow(b6)),
                    Chr = rep(b6$genes$chromosome, times = ncol(b6)),
                    Gene = rep(b6$genes$symbol, times = ncol(b6)),
                    Ensembl = rep(b6$genes$ensembl, times = ncol(b6)),
                    b6 = c(b6$counts), cast = c(cast$counts))
df_wt <- df_wt %>% 
  dplyr::group_by(sample) %>% 
  dplyr::mutate(libsize_b6 = sum(b6[!Gene %in% "Xist_5prime"]), libsize_cast = sum(cast[!Gene %in% "Xist_5prime"]),
                libsize_b6_Xist5prime = sum(b6[!Gene %in% "Xist"]), libsize_cast_Xist5prime = sum(cast[!Gene %in% "Xist"])) %>%
  as.data.frame()
temp <- df_wt[df_wt$Gene %in% "Xist",]

print("7.3.2) Paired t-test comparing AS Xist expression levels")
temp$b6_cpm <- (temp$b6/temp$libsize_b6)*1e6
temp$cast_cpm <- (temp$cast/temp$libsize_cast)*1e6
test <- temp %>% dplyr::group_by(day) %>% 
  dplyr::summarise(w_pvalue = t.test(x = b6_cpm,
                                     y = cast_cpm, 
                                     paired = T)$p.value) %>%
  as.data.frame()
test$pvalue <- ifelse(test$w_pvalue >= 0.01, paste0("p = ", round(test$w_pvalue, digits = 2)), 
                      ifelse(test$w_pvalue >= 0.001, paste0("p = ", round(test$w_pvalue, digits = 3)),
                             "p < 0.001"))

print("7.3.3) Plot")
cols <- c("#d95f02", "#1b9e77")
temp_melt <- melt(temp, id.vars = c("day", "Gene"), measure.vars = .("b6_cpm", "cast_cpm"))
temp_melt$variable <- revalue(temp_melt$variable, replace = c("b6_cpm" = "B6", "cast_cpm" = "Cast"))
temp_melt$variable <- factor(temp_melt$variable, levels = c("Cast", "B6"))
g <- temp_melt %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_jitter(aes(x = factor(day), y = log10(value + 1), color = variable), fill = "black",
              alpha = 1/2, position=position_jitterdodge(jitter.width = .25, dodge.width = 0.75), 
              size = scattersize, show.legend = FALSE) +
  stat_summary(fun=mean, aes(x = factor(day), y = log10(value + 1), ymin=..y.., ymax=..y.., color = variable), 
               geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = cols) + 
  geom_text(data = test, aes(x = factor(day), y = 3.5, label = pvalue), size = geomtext_size, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 3.75), breaks = seq(0, 4, 1)) +
  labs(x = "Time [days]", 
       y = expression("Xist CPM + 1 [ "*log[10]*" ]"), 
       color = "Allele",
       title = "TX1072: RNA-Seq - Xist expression")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 2, 
            savefile = paste0(outpath, "S8_C_BulkdTX1072_Xist.pdf"), 
            height = 2, width = 4)



print("7.4) D: Bulk RNA-Seq: TX1072 - Xist gene expression using only 5'-SNP")

print("7.4.1) Load data")
temp <- df_wt[df_wt$Gene %in% "Xist_5prime",]
temp$b6_cpm <- (temp$b6/temp$libsize_b6_Xist5prime)*1e6
temp$cast_cpm <- (temp$cast/temp$libsize_cast_Xist5prime)*1e6

print("7.4.2) Paired t-test comparing AS Xist expression levels")
test <- temp %>% dplyr::group_by(day) %>% 
  dplyr::summarise(t_pvalue = t.test(x = b6_cpm,
                                     y = cast_cpm, 
                                     paired = T)$p.value) %>%
  as.data.frame()
test$pvalue <- ifelse(test$t_pvalue >= 0.01, paste0("p = ", round(test$t_pvalue, digits = 2)), 
                      ifelse(test$t_pvalue >= 0.001, paste0("p = ", round(test$t_pvalue, digits = 3)),
                             "p < 0.001"))
print("7.4.3) Plot")
cols <- c("#d95f02", "#1b9e77")
temp_melt <- melt(temp, id.vars = c("day", "Gene"), measure.vars = .("b6_cpm", "cast_cpm"))
temp_melt$variable <- revalue(temp_melt$variable, replace = c("b6_cpm" = "B6", "cast_cpm" = "Cast"))
temp_melt$variable <- factor(temp_melt$variable, levels = c("Cast", "B6"))
g <- temp_melt %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_jitter(aes(x = factor(day), y = log10(value + 1), color = variable), fill = "black",
              alpha = 1/2, position=position_jitterdodge(jitter.width = .25, dodge.width = 0.75), 
              size = scattersize, show.legend = FALSE) +
  stat_summary(fun=mean, aes(x = factor(day), y = log10(value + 1), ymin=..y.., ymax=..y.., color = variable), 
               geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = cols) + 
  geom_text(data = test[!is.na(test$pvalue),], aes(x = factor(day), y = 3.5, label = pvalue), size = geomtext_size, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 3.75), breaks = seq(0, 4, 1)) +
  labs(x = "Time [days]", 
       y = expression("5'-Xist CPM + 1 [ "*log[10]*" ]"), 
       color = "Allele",
       title = "TX1072: RNA-Seq - Xist 5' expression")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 2, 
            savefile = paste0(outpath, "S8_D_BulkdTX1072_Xist5prime.pdf"), 
            height = 2, width = 4)









print("8) Supplementary tables")

print("8.1) ST1 - Cell and Gene filtering")

print("8.2) ST2 - Cell and Gene measures")

print("8.3) ST3 - DE and Correlation analyses")

f3path <- paste0(path, "output/fig3_deXistHighLow/")
f4path <- paste0(path, "output/fig4_deXchrChangeHighLow/")

print("8.3.1) Load data")

# Xist - DE analysis
load(paste0(f3path, "ALLresults.RData"))
fields <- c("day", "chromosome_name", "mgi_symbol", "coef", "fdr")
xist_de <- all_de[, fields]

# Xist - Correlation analysis
load(paste0(f3path, "sprmcor.RData"))
colnames(sprmcor) <- c("day", "mgi_symbol", "chromosome_name", "correlation", 
                       "droprate", "droprate_Xist", "n", "pvalue_Xist", "fdr")
fields <- c("day", "chromosome_name", "mgi_symbol", "correlation", "fdr")
xist_cor <- sprmcor[, fields]

# Xchr Change - DE analysis
load(paste0(f4path, "XchrChange_HighLow_MAST.RData"))
fields <- c("day", "chromosome_name", "mgi_symbol", "coef", "fdr")
xchr_de <- all_de[, fields]

# Xchr Change - Correlation analysis
load(paste0(f4path, "cor_Xchr.RData"))
colnames(sprmcor) <- c("day", "mgi_symbol", "chromosome_name", "correlation", 
                       "droprate", "n", "ypos", "pvalue", "fdr")
fields <- c("day", "chromosome_name", "mgi_symbol", "correlation", "fdr")
xchr_cor <- sprmcor[, fields]

print("8.3.2) Combine results")
names <- c("day", "chromosome", "gene", "value", "fdr") 
colnames(xist_de) <- colnames(xist_cor) <- colnames(xchr_de) <- colnames(xchr_cor) <- names
res <- rbind(data.frame(test = "Xist-DE", xist_de),
             data.frame(test = "Xist-COR", xist_cor),
             data.frame(test = "Xchr-DE", xchr_de),
             data.frame(test = "Xchr-COR", xchr_cor))
genes <- unique(as.character(res$gene))
genes <- genes[!grepl(genes, pattern = "^ERCC")]
genes <- genes[order(genes)]
ensembl <- gene_features$Ensembl.ID[match(genes, gene_features$Gene.Symbol)]
chromosome <- gene_features$Chromosome[match(genes, gene_features$Gene.Symbol)]
load(paste0(datapath, "DGE.RData")); dge$cpm <- t(t(dge$counts)/(dge$samples$eff_libsize_notX))*1e6; cpm <- dge$cpm[match(ensembl, dge$genes$ensembl),]

for(d in 1:4){
  temp <- res[res$day %in% d,]
  x1 <- temp[temp$test %in% "Xist-DE",]
  x2 <- temp[temp$test %in% "Xist-COR",]
  x3 <- temp[temp$test %in% "Xchr-DE",]
  x4 <- temp[temp$test %in% "Xchr-COR",]
  avgcpm <- rowMeans(cpm[, dge$samples$day %in% d])
  x <- data.frame(d, 
                  chromosome, 
                  genes,
                  avgcpm,
                  x1$value[match(genes, x1$gene)],
                  x1$fdr[match(genes, x1$gene)],
                  x2$value[match(genes, x2$gene)],
                  x2$fdr[match(genes, x2$gene)],
                  x3$value[match(genes, x3$gene)],
                  x3$fdr[match(genes, x3$gene)],
                  x4$value[match(genes, x4$gene)],
                  x4$fdr[match(genes, x4$gene)])
  x <- as.matrix(x)
  colnames(x) <- c("Day", "Chromosome", "Gene Symbol", "Average CPM",
                   "log2FC [MAST: Xist High vs Low]", "FDR [MAST: Xist High vs Low]",
                   "Coefficient [Spearman: Gene CPM vs Xist CPM]", "FDR [Spearman: Gene CPM vs Xist CPM]",
                   "log2FC [MAST: Xchr.Change High vs Low]", "FDR [MAST: Xchr.Change High vs Low]",
                   "Coefficient [Pearson: Gene CPM vs Xchr.Change]", "FDR [Pearson: Gene CPM vs Xchr.Change]")
  if(d==1){
    wb <- createWorkbook()
  }
  addWorksheet(wb, sheetName = paste0("Day ", d))
  writeData(wb, sheet = paste0("Day ", d), x = x, rowNames = F, colNames = T, keepNA = T)
}
saveWorkbook(wb, file = paste0(outpath, "ST3_DEA_Correlation.xlsx"), overwrite = T)

print("8.4) ST4 - Experimental Design")

print("8.5) ST5 - Xist not-AS and AS cell classification changing UMI thresholds")

print("8.5.1) Xist not-AS classification changing UMI not-AS threshold")

print("8.5.1.1) Load data")
load(paste0(datapath, "NGFCF_DGE.RData")); notas <- dge
df <- data.frame(day = notas$samples$day, ID = notas$samples$id, xist_umi = notas$counts["Xist",])

print("8.5.1.2) Loop and AS Xist classification")
notas_threshold <- c(0, 5, 10, 25, 50, 100, 200)
notas_classification <- c()
lv <- c("Xist-", "Low", "Xist+")
for(thr in notas_threshold){
  temp <- df
  temp$Xist_notAS <- ifelse(df$xist_umi == 0, "Xist-", 
                            ifelse(df$xist_umi <= thr, "Low",  "Xist+"))
  temp$Xist_notAS <- factor(temp$Xist_notAS, levels = lv)
  temp <- temp %>% dplyr::group_by(day, Xist_notAS) %>% dplyr::summarise(n = length(ID))
  notas_classification <- rbind(notas_classification,
                                data.frame(temp, notAS_threshold = thr))
}

print("8.5.1.3) Write table")
x <- notas_classification %>% arrange(day, Xist_notAS, notAS_threshold)
colnames(x) <- c("Time [Days]", "Xist notAS class", "# Cells", "Xist notAS UMI threshold")
wb <- createWorkbook()
sheetname <- "Xist notAS classification"
addWorksheet(wb, sheetName = sheetname)
writeData(wb, sheet = sheetname, x = x, rowNames = F, colNames = T, keepNA = T)

print("8.5.2) Xist AS classification changing UMI AS threshold")

print("8.5.2.1) Load data")
load(paste0(datapath, "NGFCF_DGE.RData")); notas <- dge
load(paste0(datapath, "NGFCF_DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "NGFCF_DGE_Cast.RData")); cast <- dge
df <- data.frame(day = b6$samples$day, ID = b6$samples$id, 
                 xist_umi = notas$counts["Xist",],
                 xist_b6 = b6$counts["Xist",], 
                 xist_cast = cast$counts["Xist",], 
                 xist_ratio = b6$counts["Xist",]/(b6$counts["Xist",] + cast$counts["Xist",]))

print("8.5.2.2) Loop and AS Xist classification")
as_threshold <- 5:20
as_classification <- c()
lv <- c("Undetected", "Low", "Xist-MA (Xi=Cast)", "Xist-MA (Xi=B6)", "Skewed", "BA")
for(thr in as_threshold){
  temp <- df
  temp$Xist_AS <- ifelse(temp$xist_ratio == 0, "Xist-MA (Xi=Cast)", 
                         ifelse(temp$xist_ratio == 1, "Xist-MA (Xi=B6)", 
                                ifelse((temp$xist_ratio >= 0.2) & (temp$xist_ratio <= 0.8), 
                                       "BA", "Skewed")))
  temp$Xist_AS[(temp$xist_b6 + temp$xist_cast) <= thr] <- "Low"
  temp$Xist_AS[temp$xist_umi == 0] <- "Undetected"
  temp$Xist_AS <- factor(temp$Xist_AS, levels = lv)
  temp <- temp %>% dplyr::group_by(day, Xist_AS) %>% dplyr::summarise(n = length(ID))
  as_classification <- rbind(as_classification,
                             data.frame(temp, AS_threshold = thr))
}

print("8.5.2.3) Write table")
x <- as_classification %>% arrange(day, Xist_AS, AS_threshold)
colnames(x) <- c("Time [Days]", "Xist AS class", "# Cells", "Xist AS UMI threshold")
sheetname <- "Xist AS classification"
addWorksheet(wb, sheetName = sheetname)
writeData(wb, sheet = sheetname, x = x, rowNames = F, colNames = T, keepNA = T)

print("8.5.3) Store supplementary table")
saveWorkbook(wb, file = paste0(outpath, "ST5_XistClassification_ChangingThresholds.xlsx"), overwrite = T)