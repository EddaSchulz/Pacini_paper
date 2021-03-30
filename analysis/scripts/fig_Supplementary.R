print("1) Define output path")
outpath <- paste0(path, "output/fig_Supplementary/"); dir.create(path = outpath, recursive = T, showWarnings = F)









print("2) Supplementary figure 1: Mapping statistics and data pre-processing")

print("2.1) A: Alignment and gene quantification - notAS")

print("2.1.1) Load data")
load(paste0(datapath, "UF_DGE.RData")); notas <- dge
load(paste0(datapath, "UF_DGE_Spliced.RData")); spliced <- dge
load(paste0(datapath, "UF_DGE_Unspliced.RData")); unspliced <- dge

summary_stats <- data.frame(id = notas$samples$id, 
                            day = notas$samples$day,
                            total = notas$samples$seqdepth, 
                            ua = notas$samples$uniquealigned, 
                            exonic = colSums(notas$counts), 
                            spliced = colSums(spliced$counts),
                            unspliced = colSums(unspliced$counts))

print("2.1.2) Compute median number and percentage of reads across samples")
sum <- summary_stats %>%
  dplyr::group_by(day) %>% 
  dplyr::summarise(seqdepth = median(total), 
                   ua = median(ua), 
                   exonic = median(exonic),
                   spliced = median(spliced), 
                   unspliced = median(unspliced))
perc <- data.frame(summary_stats[, c("total", "ua", "exonic", "spliced", "unspliced")]/summary_stats$total*100, day = summary_stats$day)
perc <- perc %>%
  dplyr::group_by(day) %>% 
  dplyr::summarise(seqdepth = median(total), 
                   ua = median(ua), 
                   exonic = median(exonic),
                   spliced = median(spliced), 
                   unspliced = median(unspliced))

print("2.1.3) Melt variables")
sum_melt <- tidyr::gather(sum, 'variable', 'value', -day)
perc_melt <- tidyr::gather(perc, 'variable', 'value', -day)
perc_melt$value <- paste0(round(perc_melt$value, digits = 1), "%")
sum_melt$perc <- perc_melt$value

newlev <- c("seqdepth" = "Total\nReads", 
            "ua" = "Uniquely\nAligned", 
            "exonic" = "UMI\nnot-AS", 
            "spliced" = "UMI\nSpliced", 
            "unspliced" = "UMI\nUnspliced")
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
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "S1_A_SeqOutput_notAS.pdf"))


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









print("3) Supplementary figure 3: Marker gene expression and X:A expression ratios")

print("3.1) A: X:A ratio varying Xist UMI threshold")

print("3.1.1) Load data")
load(paste0(datapath, "DGE.RData"))
nboot <- 1e3
xcounts <- dge$counts[dge$genes$chromosome %in% "X",]
autcounts <- dge$counts[dge$genes$chromosome %in% c(1:19),]
n_xlinked <- nrow(xcounts); n_autlinked <- nrow(autcounts)

print("3.1.2) Define Xist thresholds and compute bootstrapped X:A ratios")
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

print("3.1.3) Store results")
x2a_sub <- data.frame(x2a[x2a$Xist_UMI == 0,], group = "Xist = 0")
for(j in seq_len(length(Xist_thr))){
  thr <- Xist_thr[j]
  x2a_sub <- rbind(x2a_sub,
                   data.frame(x2a[x2a$Xist_UMI > thr,], group = paste0("Xist > ", thr)))
}
summ_cell <- ddply(x2a_sub, .variables = .(day, group), summarize, n = length(day))

print("3.1.4) Plot")
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
            savefile = paste0(outpath, "S3_A_x2a_XistnotASthreshold.pdf"), width = 7, height = 3)


print("3.2) B: Gene-wise X:A ratios in Xist+ and Xist- cells")

print("3.2.1) Load data")
load(paste0(datapath, "DGE.RData"))
df <- data.frame(day = rep(dge$samples$day, each = nrow(dge)), 
                 Cell = rep(colnames(dge), each = nrow(dge)), 
                 Xist = rep(dge$samples$Xist_class, each = nrow(dge)),
                 els = rep(dge$samples$eff_libsize_notX, each = nrow(dge)),
                 Gene = rep(rownames(dge), times = ncol(dge)),
                 Chromosome = rep(dge$genes$chromosome, times = ncol(dge)),
                 umi = c(dge$counts))
df$cpm_value <- (df$umi/df$els)*1e6
df$isX <- ifelse(df$Chromosome == "X", "X-linked", "Autosomal")
df <- df[!is.na(df$isX),]

# filter out lowly expressed genes at d0
x <- df[df$day == 0,] %>% dplyr::group_by(Chromosome, Gene) %>% dplyr::summarise(avecpm = mean(cpm_value))
table(highexp = x$avecpm > 0, Xlinked = x$Chromosome %in% "X")
keep <- as.character(x$Gene[x$avecpm > 0])
df <- df[df$Gene %in% keep,]

print("3.2.2) Gene-wise X:A ratios Xist+ and Xist- cells")
df <- df %>% dplyr::group_by(day) %>% dplyr::mutate(aut_mean = mean(cpm_value[Chromosome %in% c(1:19)]))
df <- df %>% dplyr::group_by(Gene) %>% dplyr::mutate(d0_ratio = mean(cpm_value[day == 0])/unique(aut_mean[day == 0]))

x2a <- df[(!df$Xist %in% "Detected (Xist UMI <= 5)")&(df$Chromosome %in% "X"), ] %>%
  dplyr::group_by(day, Xist, Chromosome, Gene) %>%
  dplyr::summarise(xgene_meancpm = mean(cpm_value),
                   autosomes_meancpm = unique(aut_mean),
                   ratio = mean(cpm_value)/unique(aut_mean),
                   d0_ratio = unique(d0_ratio)) %>%
  as.data.frame()
x2a$Xist <- revalue(factor(as.character(x2a$Xist)), 
                    replace = c("Undetected" = "Xist- cells",
                                "Detected (Xist UMI > 5)" = "Xist+ cells"))
x2a$Xist <- factor(x2a$Xist, levels = c("Xist- cells", "Xist+ cells"))

print("3.2.3) Normalize ratios to d0 and test if median ratio is equal to 1")
x2a <- x2a %>% dplyr::group_by(Gene) %>% dplyr::mutate(ratio_d0 = (ratio + 0.01)/(unique(d0_ratio) + 0.01))
test <- x2a[x2a$day > 0,] %>% 
  dplyr::group_by(day, Xist) %>% 
  dplyr::summarise(pos = sum(ratio_d0 >= 1),
                   n = length(ratio_d0),
                   bt_pvalue = binom.test(x = sum(ratio_d0 >= 1), n = length(ratio_d0), p = 0.5)$p.value) %>%
  as.data.frame()
test$pvalue <- ifelse(test$bt_pvalue >= 0.01, paste0("p = ", round(test$bt_pvalue, digits = 2)), 
                      ifelse(test$bt_pvalue >= 0.001, paste0("p = ", round(test$bt_pvalue, digits = 3)),
                             "p < 0.001"))
test$y <- ifelse(test$Xist == "Xist+ cells", 2.05, 1.95)

print("3.2.4) Plot")
cols <- c("black", "#e7298a")
x <- x2a %>% dplyr::group_by(day, Xist) %>% dplyr::summarize(m = median(ratio_d0),
                                                             q1 = quantile(ratio_d0, probs = 1/4),
                                                             q3 = quantile(ratio_d0, probs = 3/4))

g <- x[x$day>0,] %>%  
  ggplot(aes(x = factor(day), y = m, color = Xist)) +
  theme_bw() + theme1 + 
  geom_hline(yintercept = seq(0.5, 1.75, 0.25), linetype = "dashed", size = linesize, alpha = 1/2) +
  geom_point(size = 1, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin=q1, ymax=q3), width=.2,
                position=position_dodge(0.9), size = linesize) +
  geom_text(data = test,
            aes(x = factor(day), y = y, label = pvalue, color = Xist), 
            size = geomtext_size, show.legend = FALSE) +
  scale_color_manual(values = cols) + 
  scale_y_continuous(breaks = seq(0.5, 1.75, 0.25), limits = c(0.5, 2.1)) +
  labs(x = "Time [days]", 
       y = "(ratio+0.01)/(d0_ratio + 0.01) [q1, median, q3]",
       color = "")
adjust_size(g = g, panel_width_cm = 4.5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S3_B_GeneWise_x2a_d0_mediantest.pdf"), width = 8, height = 2)


print("3.3) C: Xist regulators over time")

print("3.3.1) Load data")
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

print("3.3.2) Test difference between Xist- and Xist+ cells")
test <- df[df$day>0,] %>% dplyr::group_by(day, Gene) %>% 
  dplyr::summarise(w_pvalue = t.test(x = cpm_value[Xist == "Undetected"],
                                     y = cpm_value[Xist == "Detected (Xist UMI > 5)"])$p.value)
test$pvalue <- ifelse(test$w_pvalue >= 0.01, paste0("p = ", round(test$w_pvalue, digits = 2)), 
                      ifelse(test$w_pvalue >= 0.001, paste0("p = ", round(test$w_pvalue, digits = 3)),
                             "p < 0.001"))

print("3.3.3) Plot")
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
adjust_size(g = g, panel_width_cm = 4.5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S3_C_Xist_Regulators.pdf"), width = 8, height = 2)









print("4) Supplementary figure 4: Allele-specific mapping statistics")

print("4.1) A: Alignment and gene quantification - spliced/unspliced quantification")

print("4.1.1) Load data")
load(paste0(datapath, "UF_DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "UF_DGE_Cast.RData")); cast <- dge
load(paste0(datapath, "UF_DGE_B6_Spliced.RData")); b6_spliced <- dge
load(paste0(datapath, "UF_DGE_B6_Unspliced.RData")); b6_unspliced <- dge
load(paste0(datapath, "UF_DGE_Cast_Spliced.RData")); cast_spliced <- dge
load(paste0(datapath, "UF_DGE_Cast_Unspliced.RData")); cast_unspliced <- dge

summary_stats <- data.frame(id = b6_spliced$samples$id, 
                            day = b6_spliced$samples$day,
                            total = b6_spliced$samples$seqdepth, 
                            ua = b6_spliced$samples$uniquealigned, 
                            b6 = colSums(b6$counts), 
                            cast = colSums(cast$counts),
                            spliced_b6 = colSums(b6_spliced$counts), 
                            unspliced_b6 = colSums(b6_unspliced$counts), 
                            spliced_cast = colSums(cast_spliced$counts), 
                            unspliced_cast = colSums(cast_unspliced$counts))

print("4.1.2) Compute median number and percentage of reads across samples")
sum <- summary_stats %>%
  dplyr::group_by(day) %>% 
  dplyr::summarise(seqdepth = median(total), 
                   ua = median(ua), 
                   B6 = median(b6),
                   Cast = median(cast),
                   spliced_B6 = median(spliced_b6),
                   spliced_Cast = median(spliced_cast),
                   unspliced_B6 = median(unspliced_b6),
                   unspliced_Cast = median(unspliced_cast))
features <- colnames(summary_stats)[!colnames(summary_stats) %in% c("id", "day")]
perc <- data.frame(summary_stats[, features]/summary_stats$total*100, day = summary_stats$day)
perc <- perc %>%
  dplyr::group_by(day) %>% 
  dplyr::summarise(seqdepth = median(total), 
                   ua = median(ua), 
                   B6 = median(b6),
                   Cast = median(cast),
                   spliced_B6 = median(spliced_b6),
                   spliced_Cast = median(spliced_cast),
                   unspliced_B6 = median(unspliced_b6),
                   unspliced_Cast = median(unspliced_cast))

print("4.1.3) Melt variables")
sum_melt <- tidyr::gather(sum, 'variable', 'value', -day)
perc_melt <- tidyr::gather(perc, 'variable', 'value', -day)
perc_melt$value <- paste0(round(perc_melt$value, digits = 1), "%")
sum_melt$perc <- perc_melt$value

newlev <- c("seqdepth" = "Total\nReads", 
            "ua" = "Uniquely\nAligned", 
            "B6" = "UMI\nB6", 
            "spliced_B6" = "UMI\nSpliced\nB6", 
            "unspliced_B6" = "UMI\nUnspliced\nB6",
            "Cast" = "UMI\nCast",
            "spliced_Cast" = "UMI\nSpliced\nCast",
            "unspliced_Cast" = "UMI\nUnspliced\nCast")
sum_melt$variable <- revalue(sum_melt$variable, replace = newlev)
perc_melt$variable <- revalue(perc_melt$variable, replace = newlev)
sum_melt$variable <- factor(sum_melt$variable, levels = newlev)

print("4.1.4) Plot")
features <- as.character(unique(sum_melt$variable)[grepl(unique(sum_melt$variable), pattern = "UMI")])
s <- sum_melt[sum_melt$variable %in% features,]
g <- s %>% 
  ggplot() + 
  theme_bw() + theme1 + 
  geom_bar(aes(x = variable, y = value, fill = factor(day)), stat = "identity", position=position_dodge()) + 
  geom_text(aes(x = variable, y = value, group = factor(day), label=perc), hjust = -.1,
            color="black", position = position_dodge(0.9), size = geomtext_size, angle = 90, show.legend = FALSE) +
  scale_fill_grey() +
  scale_y_continuous(label = scientific_label, breaks = seq(0, 1.25e4, by = 2e3), limits = c(0, 1.35e4)) + 
  labs(x = "", y = "Number of UMI counts", fill = "Time [days]")
adjust_size(g = g, panel_width_cm = 7, panel_height_cm = 3, 
            savefile = paste0(outpath, "S4_A_SeqOutput_AS.pdf"), width = 10)



print("4.2) B: Gene quantification - AS")

print("4.2.1) Load AS data")
load(paste0(datapath, "UF_DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "UF_DGE_Cast.RData")); cast <- dge
n_b6 <- colSums(b6$counts>0); n_cast <- colSums(cast$counts>0)
temp <- data.frame(day = b6$samples$day,
                   Cell = colnames(b6),
                   n_b6, n_cast)

print("4.2.2) Melt variables")
temp_melt <- melt(temp, id.vars = c("day", "Cell"), measure.vars = c("n_b6", "n_cast"))
temp_melt$variable <- revalue(factor(temp_melt$variable), replace = c("n_b6" = "B6", 
                                                                      "n_cast" = "Cast"))

print("4.2.3) Plot")
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
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, 
            savefile = paste0(outpath, "S4_B_AS_DetectedGenes.pdf"))









print("5) Supplementary figure 5: Xist classification and allelic expression analysis")


print("5.1) A: Xist AS classification changing AS UMI threshold")

print("5.1.1) Load data")
load(paste0(datapath, "NGFCF_DGE.RData")); notas <- dge
load(paste0(datapath, "NGFCF_DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "NGFCF_DGE_Cast.RData")); cast <- dge
df <- data.frame(day = b6$samples$day, ID = b6$samples$id, 
                 xist_umi = notas$counts["Xist",],
                 xist_b6 = b6$counts["Xist",], 
                 xist_cast = cast$counts["Xist",], 
                 xist_ratio = b6$counts["Xist",]/(b6$counts["Xist",] + cast$counts["Xist",]))

print("5.1.2) Loop and AS Xist classification")
as_threshold <- seq(10, 25, by = 5)
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

print("5.1.3) Compute percentage of cells assigned to each Xist AS class for each day and threshold")
# compute percentages
as_classification <- as_classification %>% 
  dplyr::group_by(day, AS_threshold) %>% 
  dplyr::mutate(tot = sum(n)) %>% as.data.frame()
perc <- as_classification %>%
  dplyr::group_by(AS_threshold, day, Xist_AS) %>%
  dplyr::summarise(p = (n/tot)*100) %>% as.data.frame()
perc$Xist_AS <- factor(perc$Xist_AS, 
                       levels = c("Undetected", "Low", "Xist-MA (Xi=Cast)", "Xist-MA (Xi=B6)", "Skewed", "BA"))

# include missing levels
missing <- perc %>% 
  dplyr::group_by(AS_threshold, day) %>%
  dplyr::summarise(Xist_AS = c("Undetected", "Low", "Xist-MA (Xi=Cast)", "Xist-MA (Xi=B6)", "Skewed", "BA")[!c("Undetected", "Low", "Xist-MA (Xi=Cast)", "Xist-MA (Xi=B6)", "Skewed", "BA") %in% unique(as.character(Xist_AS))],
                   p = 0) %>%
  as.data.frame()
perc <- rbind(perc, missing)

print("5.1.4) Plot")
perc$Xist_AS_threshold <- paste0("Xist AS threshold = ", perc$AS_threshold)
g <- ggplot(perc, aes(x = day, y = p, fill = Xist_AS)) + 
  theme_bw() + theme1 +
  facet_grid(. ~ Xist_AS_threshold) +
  geom_area(size=0, colour="white") +
  scale_fill_manual(values=color_alleles) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("Time [days]") + ylab("Cells [%]") +
  guides(fill = guide_legend(title = "", override.aes = list(alpha = 1)))
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, 
            savefile = paste0(outpath, "S5_A_XistAS_VaryUMIthreshold.pdf"),
            width = 10, heigh = 3)


print("5.2) A: Bulk RNA-Seq: TX1072 - Xist overall and 5' gene expression")

print("5.2.1) Load data")
load(paste0(datapath, "NGF_DGE_TX1072_B6.RData")); b6 <- dge
load(paste0(datapath, "NGF_DGE_TX1072_Cast.RData")); cast <- dge
df_wt <- data.frame(day = rep(b6$samples$day, each = nrow(b6)), 
                    sample = rep(colnames(b6), each = nrow(b6)),
                    sf = rep(b6$samples$sf_notX, each = nrow(b6)),
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

print("5.2.2) Paired t-test comparing AS Xist expression levels")
temp$b6_cpm <- (temp$b6/(temp$libsize_b6*temp$sf))*1e6
temp$cast_cpm <- (temp$cast/(temp$libsize_cast*temp$sf))*1e6
test <- temp %>% dplyr::group_by(day) %>% 
  dplyr::summarise(t_pvalue = t.test(x = b6_cpm,
                                     y = cast_cpm, 
                                     paired = T)$p.value) %>%
  as.data.frame()
test$pvalue <- ifelse(test$t_pvalue >= 0.01, paste0("p = ", round(test$t_pvalue, digits = 2)), 
                      ifelse(test$t_pvalue >= 0.001, paste0("p = ", round(test$t_pvalue, digits = 3)),
                             "p < 0.001"))

print("5.2.3) Plot")
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
  geom_text(data = test, aes(x = factor(day), y = 2.85, label = pvalue), size = geomtext_size, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 4, 1)) +
  labs(x = "Time [days]", 
       y = expression("Xist CPM + 1 [ "*log[10]*" ]"), 
       color = "Allele",
       title = "TX1072: RNA-Seq - Xist expression")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S5_B_BulkTX1072_Xist.pdf"), 
            height = 3, width = 4)


print("5.2.4) Xist 5' expression")
temp <- df_wt[df_wt$Gene %in% "Xist_5prime",]
temp$b6_cpm <- (temp$b6/(temp$libsize_b6_Xist5prime*temp$sf))*1e6
temp$cast_cpm <- (temp$cast/(temp$libsize_cast_Xist5prime*temp$sf))*1e6
test <- temp %>% dplyr::group_by(day) %>% 
  dplyr::summarise(t_pvalue = t.test(x = b6_cpm,
                                     y = cast_cpm, 
                                     paired = T)$p.value) %>%
  as.data.frame()
test$pvalue <- ifelse(test$t_pvalue >= 0.01, paste0("p = ", round(test$t_pvalue, digits = 2)), 
                      ifelse(test$t_pvalue >= 0.001, paste0("p = ", round(test$t_pvalue, digits = 3)),
                             "p < 0.001"))

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
  geom_text(data = test, aes(x = factor(day), y = 2.85, label = pvalue), size = geomtext_size, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 4, 1)) +
  labs(x = "Time [days]", 
       y = expression("Xist 5' CPM + 1 [ "*log[10]*" ]"), 
       color = "Allele",
       title = "TX1072: RNA-Seq - Xist expression")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S5_B_BulkTX1072_Xist5prime.pdf"), 
            height = 3, width = 4)









print("6) Supplementary figure 6: Scaled pseudotime, AS X:A ratio and dX")

print("6.1) A: Allele specific X:A ratio for Xist-MA cells vs Scaled Pseudotime")

print("6.1.1) Load AS data")
load(paste0(datapath, "DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "DGE_Cast.RData")); cast <- dge
xcounts_b6 <- b6[(b6$genes$chromosome %in% "X") & (!b6$genes$symbol %in% "Xist"),]$counts
autcounts_b6 <- b6[b6$genes$chromosome %in% c(1:19),]$counts
xcounts_cast <- cast[(cast$genes$chromosome %in% "X") & (!cast$genes$symbol %in% "Xist"),]$counts
autcounts_cast <- cast[cast$genes$chromosome %in% c(1:19),]$counts
n_xlinked <- nrow(xcounts_b6); n_autlinked <- nrow(autcounts_b6)

print("6.1.2) Compute AS bootstrapped X:A ratios")
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

print("6.1.3) Select only Xist_MA cells and define Xi and Xa accordingly")
temp <- x2a
temp_melt <- melt(temp, id.vars = c("day", "id", "xist_classification"), 
                  measure.vars = c("b6_XiAi", "cast_XiAi"))
temp_melt$allele <- revalue(factor(temp_melt$variable),
                            replace = c("b6_XiAi" = "B6",
                                        "cast_XiAi" = "Cast"))
temp_melt$xist_classification <- factor(revalue(temp_melt$xist_classification,
                                                replace = c("Undetected" = "Undetected",
                                                            "Low-Xist" = "Low",
                                                            "Middle" = "Skewed",
                                                            "Xist_BA" = "BA",
                                                            "BL6_MA" = "Xist-MA (Xi=B6)",
                                                            "Cast_MA" = "Xist-MA (Xi=Cast)")),
                                        levels = names(color_alleles))
temp_melt$pdt <- pData(XX)$Scaled_PDT[match(temp_melt$id, pData(XX)$id)]

print("6.1.4) Plot")
g <- temp_melt[!temp_melt$xist_classification %in% c("BL6_MA", "Cast_MA"),] %>%
  ggplot() + 
  theme_bw() + theme1 + 
  facet_grid(xist_classification ~ allele) + 
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 1/2, size = linesize) + 
  geom_point(aes(x = pdt, y = value, color = factor(day)), size = scattersize, alpha = 1/2, shape = 20) + 
  scale_color_manual(values = time_colors) + 
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  guides(color=guide_legend(override.aes = list(size=2, alpha = 1))) + 
  theme(legend.position = "right", 
        strip.text.y = element_text(angle = 0)) +
  labs(x="Scaled Pseudotime", y = "Allele Specific X:A ratio", color = "Time [days]")
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, 
            height = 7,
            savefile = paste0(outpath, "S6_A_asX2A_PDT.pdf"))



print("6.2) B: Scaled PDT vs Xist expression - coloring by time point")

g <- pData(XX) %>% 
  ggplot() + 
  theme_bw() + theme1 +
  geom_point(aes(x = Scaled_PDT, y = log10(xist_cpm + 1), color = factor(day)),
             size = small_scattersize) + 
  scale_color_manual(values = time_colors) +
  labs(x="Scaled Pseudotime", 
       y = expression("Xist ["*log[10]*"(CPM + 1)]"), 
       color = "Time [days]") +
  guides(color = guide_legend(override.aes = list(size=2)))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S6_B_XistCPM_ScaledPDT.pdf"))


print("6.3) C: dX in Xist+ and Xist- cells")

print("6.3.1) Load data")
f1path <- paste0(path, "output/fig1_NotAS/")
load(paste0(datapath, "NCF_DGE_Spliced.RData")); spliced <- dge
load(paste0(datapath, "DGE.RData")); notas <- dge
load(paste0(f1path, "notAS_vel.RData"))

print("6.3.2) Compute RNA-velocity predicted change in normalized X-chromosome expression")

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

print("6.3.3) Define data frame")

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

print("6.3.4) Plot")
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
       y = expression(log[2]*"( "*Delta*"X )"),
       color = "Xist not-AS classification")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S6_C_XchrChange_XistnotAS.pdf"), 
            width = 5, height = 2)










print("7) Supplementary figure 7: Identification of putative Xist regulators through velocity-based dX analysis")


f1path <- paste0(path, "output/fig1_NotAS/")

print("7.1) A: X-chromosome change K-means classification")

print("7.1.1) Load data")

load(paste0(datapath, "NCF_DGE_Spliced.RData")); spliced <- dge
load(paste0(datapath, "DGE.RData")); notas <- dge
load(paste0(f1path, "notAS_vel.RData"))

print("7.1.2) Compute RNA-velocity predicted change in normalized X-chromosome expression")

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

print("7.1.3) Define data frame")

df <- data.frame(day = samples$day, 
                 id = rownames(samples),
                 Xist_AS = samples$Xist_ratio_class,
                 Xist_notAS = samples$Xist_class,
                 diff = diff_xchr)
lv <- c("Undetected" = "Undetected", "Low-Xist" = "Low",
        "Middle" = "Skewed", "Xist_BA" = "BA",
        "BL6_MA" = "Xist-MA (Xi=B6)", "Cast_MA" = "Xist-MA (Xi=Cast)")
df$xist_as <- factor(revalue(factor(df$Xist_AS), replace = lv), levels = lv)

print("7.1.4) Compute and order K-means clusters, separately for each time point")

# k-means classification
k <- 3
df_de <- df %>%
  dplyr::group_by(day) %>%
  dplyr::mutate(km = kmeans(diff, centers = k, iter.max = 1000, nstart = 1000)$cluster)

# order k-means clusters
df_de <- df_de %>%
  dplyr::group_by(day, km) %>% dplyr::mutate(kmax = max(diff)) %>%
  dplyr::group_by(day) %>% dplyr::mutate(km = ifelse(kmax == max(kmax), "high",
                                                     ifelse(kmax == min(kmax), "low", "medium"))) %>% as.data.frame()
df_de$km <- factor(df_de$km, levels = c("low", "medium", "high"))
table(K = df_de$km, day = df_de$day)

# include variable in DGE list
dge_ext <- notas
m <- match(colnames(dge_ext), df_de$id)
dge_ext$samples$km_XchrChange <- df_de$km[m]

print("7.1.5) Plot")

summ_cell <- ddply(df_de, .variables = .(day, km), summarize, n = length(day))
g <- df_de %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 1/3, size = linesize) +
  geom_jitter(aes(x = day, y = diff, group = km, color = km), 
              size = outliersize, show.legend = FALSE,
              position=position_jitterdodge(jitter.width = .5, dodge.width = 0.75)) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_text(summ_cell, mapping = aes(x = day, y = Inf, label = paste0("n = ", n), fill = km),
            position=position_dodge(0.75), angle = 90, alpha = 1, size = geomtext_size, hjust = -0.5) + 
  coord_cartesian(clip = "off") +
  labs(x="Time [days]", y = expression(log[2]*"( "*Delta*"X )"))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 2, 
            savefile = paste0(outpath, "S7_A_XchrChange_Kmeans.pdf"))


print("7.2) B: X-chromosome change High vs Xist Low MAST DE analysis per time point")

print("7.2.1) Define contrast and launch DE analyses")

# define contrast and launch MAST
contrasts <- list(HighXchrChange_LowXchrChange = c("km_XchrChange", list("high", "low")))
min_sample_test <- 10; de_threshold <- 5e-2; days <- 1:4
comparison <- "HighXchrChange_LowXchrChange"
variable <- "km_XchrChange"
variable_levels <- list("high", "low")

if(!dir.exists(paste0(outpath, comparison))){
  dir.create(path = paste0(outpath, comparison), showWarnings = FALSE, recursive = TRUE)
  
  for (i in 1:length(days)) {
    time <- days[i]; cat('Processing day', time, '..')
    dge <- dge_ext[, dge_ext$samples$day == time]
    dge$samples$sample <- rownames(dge$samples)
    
    # subset to grouping conditions only
    temp <- dge[!(grepl(dge$genes$symbol, pattern = "^ERCC") | duplicated(dge$genes$ensembl)), dge$samples[,variable] %in% unlist(variable_levels)]
    grp <- temp$samples[,variable]
    grp <- ifelse(grp %in% variable_levels[[1]], "case", "control"); table(grp)
    names(grp) <- temp$samples$sample
    grp <- factor(grp, levels = c("control", "case"))
    
    # store the number of cases and controls per time point
    if(!file.exists(paste0(outpath, comparison, "/Nsamples_perTime.txt"))){
      m <- matrix(c(time, sum(grp == "case"), sum(grp == "control")), nrow = 1)
      colnames(m) <- c("Day", c(strsplit2(comparison, split = "\\_")))
      write.table(m, quote =FALSE, row.names = FALSE, col.names = TRUE, 
                  file = paste0(outpath, comparison, "/Nsamples_perTime.txt"))
    }else{
      m <- matrix(c(time, sum(grp == "case"), sum(grp == "control")), nrow = 1)
      write.table(m, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE,
                  file = paste0(outpath, comparison, "/Nsamples_perTime.txt"))
    }
    
    if(min(table(grp)) >= min_sample_test){
      # compute detection rate and cpm matrix
      cdr <- scale(colMeans(temp$counts > 0))
      els <- colSums(temp$counts)*(temp$samples$sf_notX/mean(temp$samples$sf_notX))
      cpm <- t(t(temp$counts)/els)*1e6; rownames(cpm) <- temp$genes$ensembl
      sca <- FromMatrix(exprsArray = log2(cpm + 1), 
                        cData = data.frame(wellKey = names(grp), grp = grp, cdr = cdr),
                        fData = data.frame(chromosome = temp$genes$chromosome, symbol = temp$genes$symbol, ensembl = temp$genes$ensembl, primerid = temp$genes$ensembl))
      zlmdata <- zlm(~ cdr + grp, sca)
      summaryCond <- summary(zlmdata, doLRT='grpcase') ; summaryDt <- summaryCond$datatable
      mast <- merge(summaryDt[contrast=='grpcase' & component=='H',.(primerid, `Pr(>Chisq)`)],
                    summaryDt[contrast=='grpcase' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
      mast$AveCPM <- rowMeans(cpm)
      mast$fdr <- p.adjust(mast$`Pr(>Chisq)`, method = "BH")
      mast <- data.frame(gene_features[match(mast$primerid, gene_features$Ensembl.ID), c("Chromosome", "Gene.Symbol", "Ensembl.ID")], mast)
      colnames(mast) <- c("chromosome_name", "mgi_symbol", "ensembl_gene_id", "primerid", "Pr..Chisq.", "coef", "ci.hi", "ci.lo", "AveCPM", "fdr")
      mast <- mast[order(mast$fdr, decreasing = FALSE),]
      
      # store results
      save(mast, file = paste0(outpath, comparison, "/", time, ".RData"))
      save(zlmdata, file = paste0(outpath, comparison, "/zlm_", time, ".RData"))
      save(sca, file = paste0(outpath, comparison, "/sca_", time, ".RData"))
    }
  }
}

# store results
all_de <- c()
for (i in 1:length(days)) {
  time <- days[i]; cat('Processing day', time, '..')
  
  for(cont in 1:length(contrasts)){
    comparison <- names(contrasts)[cont]
    
    # de results
    file <- paste0(outpath, comparison, "/", time, ".RData")
    if(file.exists(file)){
      load(file); all <- mast; all_de <- rbind(all_de, data.frame(test = comparison, day = time, all))
    }
  }
}
all_de <- all_de[(!is.na(all_de$mgi_symbol))&(!is.na(all_de$coef)),]
save(all_de, file = paste0(outpath, "XchrChange_HighLow_MAST.RData"))

print("7.2.2) Load MAST results excluding Xist from DE gene list")
all_de <- all_de[!all_de$mgi_symbol %in% "Xist",]
all_de$isX <- ifelse(all_de$chromosome_name %in% "X", "X-linked", "Autosomal")
all_de$direction <- ifelse(all_de$coef < 0, "Down-regulated", "Up-regulated")
all_de$direction <- factor(all_de$direction, levels = c("Up-regulated", "Down-regulated"))
sig_de <- all_de[(all_de$fdr <= de_threshold),]
de_results <- sig_de[sig_de$test %in% "HighXchrChange_LowXchrChange",]
de_barplot <- ddply(sig_de[sig_de$test %in% "HighXchrChange_LowXchrChange",], 
                    .variables = .(day, isX, direction), summarize, 
                    n = length(day))

all <- expand.grid(unique(all_de$day), unique(de_barplot$isX), unique(de_barplot$direction)) 
all$id <- apply(all, 1, function(x) paste(x, collapse = "_"))
present <- apply(de_barplot[,1:3], 1, function(x) paste(x, collapse = "_"))
missing <- cbind(all[!all$id %in% present, !colnames(all) %in% "id"], 0)
colnames(missing) <- colnames(de_barplot)
de_barplot_ext <- rbind(de_barplot, missing)
de_barplot_ext$direction <- factor(de_barplot_ext$direction, levels = c("Up-regulated", "Down-regulated"))
de_barplot_ext$higherexp <- revalue(de_barplot_ext$direction, replace = c("Up-regulated" = "Xchr. down-regulated cells",
                                                                          "Down-regulated" = "Xchr. up-regulated cells"))
de_barplot_ext$alpha_numb <- ifelse(de_barplot_ext$n>0, "yes", "no")

print("7.2.3) Plot - Summarize MAST DE results with barplot")
g <- de_barplot_ext %>% 
  ggplot() + 
  theme_bw() + theme1 + 
  facet_grid(.~isX) + 
  geom_bar(aes(x = factor(day), y = n, fill = higherexp), stat = "identity", position=position_dodge()) + 
  scale_fill_manual(values = c("red", "blue")) +
  geom_text(aes(x = factor(day), y = n+1, group = higherexp, label = n), 
            hjust = -.1,
            color="black", position = position_dodge(0.9), size = geomtext_size, angle = 90, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 200, by = 50)) + 
  ylab("# Differentially expressed genes") + xlab("Time [days]") +
  guides(fill = guide_legend(title="Higher expression in:", nrow=2, byrow=TRUE), alpha = 'none')
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 3, 
            savefile = paste0(outpath, "S7_B_DE_barplot.pdf"))


print("7.3) C: Represent DE genes through heatmap")

print("7.3.1) Define plotting function")
plot_heatmap <- function(all_de, 
                         minFDR = 0.05, 
                         minAbsFC = 1.5,
                         outpath){
  
  # plot settings
  cellheight <- 6; fontsize_col <- 6; fontsize_row <- 6; width <- 10; cellwidth  <- .5; fontsize <- 6
  paletteLength <- 500; colorPalette <- c("black", "gold", "red"); abslog2FC_break <- 7
  
  # define heatmap coloring
  myColor <- c(colorRampPalette(colorPalette)(paletteLength))
  
  # select levels of two groups being compared
  variable <- unlist(contrasts[[comparison]][1]); variable_levels <- contrasts[[comparison]][2:3]
  
  # subset results
  all_de$log2FC <- all_de$coef
  sig_all <- all_de[(all_de$test == comparison)&
                      (abs(all_de$log2FC) >= log2(minAbsFC))&
                      (all_de$fdr <= minFDR),]
  
  for(addlabel in c("all")){
    if(addlabel == "all"){sig <- sig_all[!(sig_all$chromosome_name %in% "X" & sig_all$coef<0),]}
    heat_path <- paste0(outpath, comparison, "/Heatmap/PerCell/")
    if(nrow(sig)>1){
      for(d in days){
        dir.create(heat_path, showWarnings = FALSE, recursive = TRUE)
        s <- sig[sig$day == d,]; s$log2FC <- s$coef
        
        # order by FC
        s <- s[order(s$log2FC, decreasing = TRUE),]
        orderlabel <- "_byFC"
        
        if(nrow(s)>1){
          # store DE genes
          degenes_s <- s$ensembl_gene_id
          
          # load normalized matrices
          load(paste0(datapath, "DGE.RData"))
          dge <- dge[, dge$samples$day == d]
          
          if(comparison == "XistHigh_XistLow"){
            levels <- 1:7
            dge$samples[[variable]] <- as.numeric(dge$samples[[variable]])
          }
          if(comparison == "HighXchrChange_LowXchrChange"){
            dge$samples[[variable]] <- as.numeric(df_de$km[match(dge$samples$id, df_de$id)])
            levels <- seq_len(length(unique(dge$samples[[variable]])))
          }
          dge <- dge[, dge$samples[[variable]] %in% levels]
          counts <- log10(t(t(dge$counts)/(dge$samples$sf_notX*colSums(dge$counts)))*1e6 + 1)
          
          # subset matrix to cells in contrast and DE genes only --> compute average gene level per variable group
          avgexp <- counts[match(degenes_s, dge$genes$ensembl), ]
          
          # define column and row annotation
          annot_col <- data.frame(Levels = dge$samples[[variable]]); rownames(annot_col) <- colnames(avgexp)
          if(comparison == "XistHigh_XistLow"){
            colors_col <- c(colorRampPalette(c("lightblue", "blue"))(3), "grey",
                            colorRampPalette(c("orange", "red"))(3))
            names(colors_col) <- 1:7
          }
          if(comparison == "HighXchrChange_LowXchrChange"){
            colors_col <- c("blue", "grey", "red")
            names(colors_col) <- levels
          }
          
          fc_breaks <- seq(-abslog2FC_break, abslog2FC_break, by = 0.5)
          ct <- cut(s$log2FC, breaks = fc_breaks)
          levels <- levels(ct)
          levels[levels %in% c("(-1,0]", "(0,1]")] <- "|log2FC|<1"
          log2fc <- cut(s$log2FC, breaks = fc_breaks, labels = levels)
          isX <- ifelse(s$chromosome_name %in% "X", "X-linked", "Autosomal")
          annot_row <- data.frame(log2FC = log2fc, isX = isX); rownames(annot_row) <- s$mgi_symbol
          colors_row <- colorRampPalette(c("blue", "white", "red"))(length(unique(levels))); names(colors_row) <- unique(levels)
          colors_row <- colors_row[names(colors_row) %in% unique(ct)]
          colors_row_isX <- c("black", "white"); names(colors_row_isX) <- c("X-linked", "Autosomal")
          color_annot <- list(Levels = colors_col, log2FC = colors_row, isX = colors_row_isX)
          
          # produce heatmap
          title_label <- ifelse(is_log10CPMo1, "Avelog10CPM", "AveUMI")
          height <- cellheight*nrow(s)/30
          avg_aut <- avgexp[!log2fc %in% "|log2FC|<1", ]
          
          # order by Xist expression
          annot_col$Xist <- counts["Xist", match(rownames(annot_col), colnames(counts))]
          o <- order(annot_col$Levels, annot_col$Xist)
          avg_aut <- avg_aut[, o]
          m <- match(colnames(avg_aut), rownames(dge$samples))
          annot_col <- data.frame(Levels = dge$samples[[variable]][m]); rownames(annot_col) <- colnames(avg_aut)
          
          # gap between up and down regulated
          gprw <- which(!duplicated(s$log2FC>0))[2] - 1
          if(is.na(gprw)) gprw <- NULL
          
          pheatmap(avg_aut, 
                   gaps_row = gprw,
                   annotation_col = annot_col,
                   # annotation_row = annot_row,
                   annotation_colors = color_annot,
                   cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, color = myColor, 
                   cellwidth = cellwidth, cellheight = cellheight, 
                   fontsize_col = fontsize_col, fontsize_row = fontsize_row, fontsize = fontsize,
                   width = width, height = height,
                   gaps_col = which(!duplicated(annot_col$Levels))[-1] - 1,
                   filename = paste0(heat_path, "day", d, "_", title_label, "_", minFDR, "_", addlabel, orderlabel, ".pdf"))
        }
      }
    }
  }
}

print("7.3.2) Produce heatmaps")
minFDR <- 0.01
minAbsFC <- 1.5
load(paste0(outpath, "XchrChange_HighLow_MAST.RData"))

plot_heatmap(all_de = all_de, minFDR = minFDR, minAbsFC = minAbsFC, outpath = outpath)

print("7.3.3) Move and rename plots in figures folder")
f <- paste0(outpath, names(contrasts), "/Heatmap/PerCell/", 
            c("day1_Avelog10CPM_0.01_all_byFC.pdf", "day2_Avelog10CPM_0.01_all_byFC.pdf"))
file.copy(from = f, to = outpath)
f1 <- paste0(outpath, c("day1_Avelog10CPM_0.01_all_byFC.pdf", "day2_Avelog10CPM_0.01_all_byFC.pdf"))
f2 <- paste0(outpath, c("S7_C_Heatmap_24h_byFC.pdf", "S7_C_Heatmap_48h_byFC.pdf"))
file.rename(from = f1, to = f2)



print("7.4) D: Correlation analysis")

print("7.4.1) Load data and launch gene-wise correlation analysis to X-chromosome change")

### compute correlation between difference in Xchr relative expression and not-AS gene expression (excluding d0 cells)
load(paste0(datapath, "DGE.RData"))
dge$cpm <- t(t(dge$counts)/(colSums(dge$counts)*dge$samples$sf_notX))*1e6
m <- match(df$id, dge$samples$id); table(is.na(m))
dge <- dge[,m]

# compute and test correlations
x <- data.frame(id = rep(dge$samples$day, each = nrow(dge)),
                day = rep(dge$samples$day, each = nrow(dge)), 
                Xchr_notASvelo_diff = rep(df$diff, each = nrow(dge)),
                Gene = rep(dge$genes$symbol, times = ncol(dge)),
                Chr = rep(dge$genes$chromosome, times = ncol(dge)),
                cpm = c(dge$cpm)
)

ncores <- min(c(detectCores(), 10))
cluster <- new_cluster(ncores)
sprmcor <- x[x$day != 0,] %>% 
  dplyr::group_by(day, Gene) %>% 
  partition(cluster) %>%
  dplyr::summarise(Chr = unique(Chr),
                   spr = cor(cpm, Xchr_notASvelo_diff, method = "spearman"),
                   droprate = mean(cpm == 0),
                   n = length(cpm),
                   ypos = quantile(cpm, probs = 0.9),
                   pvalue = cor.test(cpm, Xchr_notASvelo_diff, method = "spearman", use = "complete.obs")$p.value) %>% 
  dplyr::arrange(Gene, day) %>%
  collect()
sprmcor$fdr <- p.adjust(sprmcor$pvalue, method = "BH")
sprmcor$day <- as.numeric(as.character(sprmcor$day))
sprmcor <- sprmcor[!is.na(sprmcor$fdr),] %>% arrange(fdr)
save(sprmcor, file = paste0(outpath, "cor_Xchr.RData"))

print("7.4.2) Select genes significantly correlated to X-chromosome change for each time point")
temp <- sprmcor %>% as.data.frame()
temp$isX <- ifelse(temp$Chr %in% "X", "X-linked", "Autosomal")
temp$direction <- ifelse(temp$spr < 0, "Negative", "Positive")
temp$direction <- factor(temp$direction, levels = c("Positive", "Negative"))
sig_de <- temp[(temp$fdr <= de_threshold),]

print("7.4.3) Plot - Summarize correlation analysis results with barplot")
de_barplot <- ddply(sig_de, .variables = .(day, isX, direction), summarize, n = length(day))
all <- expand.grid(unique(sprmcor$day), unique(de_barplot$isX), unique(de_barplot$direction)) 
all$id <- apply(all, 1, function(x) paste(x, collapse = "_"))
present <- apply(de_barplot[,1:3], 1, function(x) paste(x, collapse = "_"))
missing <- cbind(all[!all$id %in% present, !colnames(all) %in% "id"], 0)
colnames(missing) <- colnames(de_barplot)
de_barplot_ext <- rbind(de_barplot, missing)
de_barplot_ext$alpha_numb <- ifelse(de_barplot_ext$n>0, "yes", "no")

g <- de_barplot_ext %>% 
  ggplot() + 
  theme_bw() + theme1 + 
  facet_grid(.~isX) + 
  geom_bar(aes(x = factor(day), y = n, fill = direction), stat = "identity", position=position_dodge()) + 
  scale_fill_manual(values = c("red", "blue")) +
  geom_text(aes(x = factor(day), y = n+1, group = direction, label = n), 
            hjust = -.1,
            color="black", position = position_dodge(0.9), size = geomtext_size, angle = 90, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 295), breaks = seq(0, 250, by = 50)) +
  labs(x = "Time [days]", y = "# Significantly correlated genes", 
       fill = "Correlation with\nX-chromosome change")
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 3, 
            savefile = paste0(outpath, "S7_D_DE_barplot_corrXchr.pdf"))


print("7.5) E/F: Heatmap - positive/negative regulators over time at day 1/2")

fdr_threshold <- 0.05
maxn <- 20

print(paste0("7.5.1) Subset to genes with FDR<=", fdr_threshold, " and |rho|>=", cor_thr, " at d1 or d2"))

sub <- sprmcor[(sprmcor$day %in% c(1, 2)),]
sub_sig <- sub[sub$fdr <= fdr_threshold,]
table(sub_sig$day, sub_sig$spr>0)
summary(abs(sub_sig$spr[sub_sig$day == 1]))
summary(abs(sub_sig$spr[sub_sig$day == 2]))

print("7.5.2) Remove pseudogenes and X-linked genes with negative correlation")
sub_sig$gene_biotype <- bm$gene_biotype[match(sub_sig$Gene, bm$mgi_symbol)]
sub_sig <- sub_sig[!grepl(x = sub_sig$gene_biotype, pattern = "pseudogene"),]
sub_sig <- sub_sig[!(sub_sig$Chr %in% "X" & sub_sig$spr<0),]

print("7.5.2) Identify positive and negative regulators")

positive <- sub_sig[sub_sig$spr > 0,] %>% dplyr::arrange(-abs(spr))
pos_genes <- unique(as.character(positive$Gene))
negative <- sub_sig[sub_sig$spr<0,] %>% dplyr::arrange(-abs(spr))
neg_genes <- unique(as.character(negative$Gene))

if(length(pos_genes)>maxn) pos_genes <- pos_genes[seq_len(maxn)]
if(length(neg_genes)>maxn) neg_genes <- neg_genes[seq_len(maxn)]


print("7.5.3) Plot results through heatmap")

# positive regulators
ph <- length(pos_genes)*0.3
temp <- sprmcor[sprmcor$Gene %in% pos_genes,]; temp$Gene <- factor(temp$Gene, levels = rev(pos_genes))
temp$significant <- temp$fdr <= fdr_threshold

g <- temp %>% 
  ggplot() + 
  theme_bw() + theme1 +
  geom_point(aes(x = factor(day), y = Gene, fill = spr, size = abs(spr)), pch=22, stroke = 0) +
  geom_point(data = temp[temp$significant == TRUE,], aes(x = factor(day), y = Gene), size = small_scattersize, color = "white", pch = 1) +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0, breaks=seq(-0.5, 0.5, .25), limits=c(-0.51,0.51)) + 
  scale_radius(breaks = seq(0, 0.5, by = 0.25), limits = c(0, 0.5), range = c(0.1, 4)) +
  labs(x = "Time [days]", y = "", 
       fill = "Spearman's gene-wise correlation:\nrho(CPM_gene, X-chromosome change)", 
       size = "Absolute correlation",
       title = "Positively correlated genes")
adjust_size(g = g, panel_width_cm = 2.5, panel_height_cm = ph, 
            savefile = paste0(outpath, "S7_E_corrXchrChange_PositiveRegulators.pdf"), 
            height = 3, width = 5)

# negative regulators
ph <- length(neg_genes)*0.3
temp <- sprmcor[sprmcor$Gene %in% neg_genes,]; temp$Gene <- factor(temp$Gene, levels = rev(neg_genes))
temp$significant <- temp$fdr <= fdr_threshold

g <- temp %>% 
  ggplot() + 
  theme_bw() + theme1 +
  geom_point(aes(x = factor(day), y = Gene, fill = spr, size = abs(spr)), pch=22, stroke = 0) +
  geom_point(data = temp[temp$significant == TRUE,], aes(x = factor(day), y = Gene), size = small_scattersize, color = "white", pch = 1) +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0, breaks=seq(-0.5, 0.5, .25), limits=c(-0.51,0.51)) + 
  scale_radius(breaks = seq(0, 0.5, by = 0.25), limits = c(0, 0.5), range = c(0.1, 4)) +
  labs(x = "Time [days]", y = "", 
       fill = "Spearman's gene-wise correlation:\nrho(CPM_gene, X-chromosome change)", 
       size = "Absolute correlation",
       title = "Negatively correlated genes")
adjust_size(g = g, panel_width_cm = 2.5, panel_height_cm = ph, 
            savefile = paste0(outpath, "S7_F_corrXchrChange_NegativeRegulators.pdf"), 
            height = 3, width = 5)










print("8) Supplementary figure 8: Identification of putative XCI Xist regulators")


print("8.1) A/B: UpSetR plots - concordance between Xist and dX analyses [d1 & d2]")

f4path <- paste0(path, "output/fig4_deXistHighLow/")
suppath <- paste0(path, "output/fig_Supplementary/")
fdr_threshold <- 0.05

print("8.1.1) Load DE and Correlation analyses results")

results <- c()

print("8.1.1.1) Xist - DE analysis")

# load data
load(paste0(f4path, "ALLresults.RData"))
fields <- c("day", "chromosome_name", "mgi_symbol", "coef", "fdr")
xist_de <- all_de[(all_de$day %in% c(1,2)), fields]
xist_de$value <- xist_de$coef
xist_de$direction <- ifelse(xist_de$coef>0, "Up-regulated", "Down-regulated") 
results <- data.frame(analysis = "DE", test = "Xist", 
                      xist_de[, c("day", "chromosome_name", "mgi_symbol", "direction", "fdr", "value")])

print("8.1.1.2) Xist - Correlation analysis")

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

print("8.1.1.3) Xchr Change - DE analysis")

# load data
load(paste0(suppath, "XchrChange_HighLow_MAST.RData"))
fields <- c("day", "chromosome_name", "mgi_symbol", "coef", "fdr")
xchr_de <- all_de[(all_de$day %in% c(1,2)), fields]
xchr_de$value <- xchr_de$coef
xchr_de$direction <- ifelse(xchr_de$coef>0, "Up-regulated", "Down-regulated") 
results <- rbind(results,
                 data.frame(analysis = "DE", test = "Xchr",
                            xchr_de[, c("day", "chromosome_name", "mgi_symbol", "direction", "fdr", "value")]))

print("8.1.1.4) Xchr Change - Correlation analysis")

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

print("8.1.1.5) Exclude pseudogenes and X-linked genes with significant negative correlation coef. or logFC")

results$gene_biotype <- bm$gene_biotype[match(results$mgi_symbol, bm$mgi_symbol)]
results <- results[!grepl(x = results$gene_biotype, pattern = "pseudogene"),]
results <- results[!(results$chromosome_name %in% "X" & results$value < 0 & results$fdr <= fdr_threshold),]

print("8.1.2) Define UpSetR data frame")

upset <- results[!grepl(results$mgi_symbol, pattern = "^ERCC"),] %>%
  dplyr::group_by(chromosome_name, mgi_symbol, day) %>%
  dplyr::summarise(de_xist = as.numeric(fdr[analysis == "DE" & test == "Xist"] <= fdr_threshold),
                   de_Xchr = as.numeric(fdr[analysis == "DE" & test == "Xchr"] <= fdr_threshold),
                   cor_xist = as.numeric(fdr[analysis == "Correlation" & test == "Xist"] <= fdr_threshold),
                   cor_Xchr = as.numeric(fdr[analysis == "Correlation" & test == "Xchr"] <= fdr_threshold)
  )

print("8.1.3) Plot")

sets <- colnames(upset)[-seq_len(3)]

# day 1
x <- data.frame(upset[upset$day == 1,])
pdf(file = paste0(outpath, "S8_A_UpSetR_d1.pdf"),
    width = 4, height = 5, onefile = TRUE, useDingbats = FALSE)
upset(x, 
      main.bar.color = "black",
      sets = sets,
      order.by = "degree",
      group.by = "degree",
      keep.order = TRUE,
      decreasing = FALSE,
      mb.ratio = c(0.7, 0.3),
      point.size = scattersize*10, line.size = linesize*2, text.scale = 1.5,
      # mainbar.y.max = maxbar, 
      mainbar.y.label = "Number of significant genes at d1\nacross the four analyses"
)
dev.off()

# day 2
x <- data.frame(upset[upset$day == 2,])
pdf(file = paste0(outpath, "S8_B_UpSetR_d2.pdf"),
    width = 4, height = 5, onefile = TRUE, useDingbats = FALSE)
upset(x, 
      main.bar.color = "black",
      sets = sets,
      order.by = "degree",
      group.by = "degree",
      keep.order = TRUE,
      decreasing = FALSE,
      mb.ratio = c(0.7, 0.3),
      point.size = scattersize*10, line.size = linesize*2, text.scale = 1.5,
      # mainbar.y.max = maxbar, 
      mainbar.y.label = "Number of significant genes at d1\nacross the four analyses"
)
dev.off()


print("8.2) C: DE analyses - known Xist regulators")

regulators <- c("Nanog", "Klf2", "Klf4", "Prdm14", "Pou5f1", "Sox2", "Myc", "Zfp42", "Ctcf", "Tsix", "Yy1", "Rlim", "Ftx")

print("8.2.1) Load Xist DE results")

load(paste0(path, "output/fig4_deXistHighLow/ALLresults.RData"))
xist <- all_de[(all_de$mgi_symbol %in% regulators)&(all_de$day %in% 1:2),] %>% as.data.frame()
features <- c("day", "mgi_symbol", "coef", "fdr")
xist <- xist[, features]
colnames(xist) <- c("day", "Gene", "value", "fdr")

print("8.2.2) Load dX DE results")

load(paste0(path, "output/fig_Supplementary/XchrChange_HighLow_MAST.RData"))
dx <- all_de[(all_de$mgi_symbol %in% regulators)&(all_de$day %in% 1:2),] %>% as.data.frame()
features <- c("day", "mgi_symbol", "coef", "fdr")
dx <- dx[, features]
colnames(dx) <- c("day", "Gene", "value", "fdr")

print("8.2.3) Plot")

de <- rbind(data.frame(xist, type = "Xist"),
            data.frame(dx, type = "dX"))
de$Gene <- factor(de$Gene, levels = rev(regulators))
de$sig_label <- ifelse(de$fdr <= 0.05, "*", "")
h <- length(unique(de$Gene))*0.3
w <- length(unique(de$day))*0.5

g <- de %>%
  ggplot() +
  theme_bw() + theme1 + 
  facet_grid(.~type) +
  geom_tile(aes(x = factor(day), y = Gene, fill = value), color = "white") +
  geom_text(aes(x = factor(day), y = Gene, label = sig_label), color = "white", size = geomtext_size) +
  scale_fill_gradient2(low = "#0571b0", high = "#ca0020", mid = "white", 
                       midpoint = 0, limit = c(-3, 3),
                       name= "log2FC") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Time [days]",
       y = "")
adjust_size(g = g, panel_width_cm = w, panel_height_cm = h, 
            savefile = paste0(outpath, "S8_C_DElogFC_KnownRegulators.pdf"), 
            height = 3, width = 5)


print("8.3) D: Correlation analyses - known Xist regulators")

print("8.3.1) Load Xist correlation results")

load(paste0(path, "output/fig4_deXistHighLow/sprmcor.RData"))
xist <- sprmcor[(sprmcor$Gene %in% regulators)&(sprmcor$day %in% 1:2),] %>% as.data.frame()
features <- c("day", "Gene", "spr_Xist", "fdr_Xist")
xist <- xist[, features]
colnames(xist) <- c("day", "Gene", "value", "fdr")


print("8.3.2) Load dX correlation results")

load(paste0(outpath, "cor_Xchr.RData"))
dx <- sprmcor[(sprmcor$Gene %in% regulators)&(sprmcor$day %in% 1:2),] %>% as.data.frame()
features <- c("day", "Gene", "spr", "fdr")
dx <- dx[, features]
colnames(dx) <- c("day", "Gene", "value", "fdr")

print("8.3.3) Plot")

cor <- rbind(data.frame(xist, type = "Xist"),
             data.frame(dx, type = "dX"))
cor$Gene <- factor(cor$Gene, levels = rev(regulators))
cor$sig_label <- ifelse(cor$fdr <= 0.05, "*", "")
h <- length(unique(cor$Gene))*0.3
w <- length(unique(cor$day))*0.5

g <- cor %>%
  ggplot() +
  theme_bw() + theme1 + 
  facet_grid(.~type) +
  geom_tile(aes(x = factor(day), y = Gene, fill = value), color = "white") +
  geom_text(aes(x = factor(day), y = Gene, label = sig_label), color = "white", size = geomtext_size) +
  scale_fill_gradient2(low = "#0571b0", high = "#ca0020", mid = "white", 
                       midpoint = 0, limit = c(-0.5, 0.5),
                       name=expression(rho)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Time [days]",
       y = "")
adjust_size(g = g, panel_width_cm = w, panel_height_cm = h, 
            savefile = paste0(outpath, "S8_D_SpearmanCorrelation_KnownRegulators.pdf"), 
            height = 3, width = 5)










print("9) Supplementary figure 9: Gene Silencing analysis")


print("9.1) A: B6/Cast ratio for Autosomal (all cells) and X-linked (Xist-Undetected cells) genes over time")

print("9.1.1) Load data")
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

print("9.1.2) Compute gene-wise B6/Cast ratios for each Xist group")
temp <- df_as[!df_as$day %in% 0,] %>%
  dplyr::group_by(day, Xist_classification, Gene, isX) %>% 
  dplyr::summarise(b6_cast_ratio = (sum(b6) + 0.01)/(sum(cast) + 0.01))
summary(temp$b6_cast_ratio)
temp_median <- temp %>%
  dplyr::group_by(day, Xist_classification, isX) %>% dplyr::summarise(med = median(b6_cast_ratio)) %>%
  dplyr::arrange(Xist_classification)
temp$group <- ifelse(temp$isX == "Autosomal", "Autosomal", paste0("X-linked\n", temp$Xist_classification))

print("9.1.3) Plot")
remove_Xistgroups <- c("Low", "Skewed", "BA")
g <- temp[temp$group %in% c("X-linked\nXist\nUndetected", "Autosomal"),] %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = hvalpha, size = linesize) + 
  geom_boxplot(aes(x = group, y = log2(b6_cast_ratio), color = factor(day)), size = violin_box_size,
               outlier.size = outliersize, outlier.alpha = 1/4, outlier.shape = NA,
               position = position_dodge2(preserve = "single")) + 
  scale_color_manual(values=time_colors[-1]) +
  labs(y = expression("["*sum[g]*"("*B6[gc]*") + 0.01] / ["*sum[g]*"("*Cast[gc]*") + 0.01] ["*log[2]*"(value)]"), x = "", fill = "",
       title = "Gene-wise log2(total B6/total Cast) expression", color = "Time [days]") +
  scale_y_continuous(breaks = seq(-10, 10, by = 1), limits = c(-3, 3))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S9_A_B6CastRatio_XistclassDay.pdf"))



print("9.2) B: B6/Cast ratio for X-linked genes (grouping Xist-Undetected, B6-MA or Cast-MA cells) over time")

print("9.2.1) Plot")
Xlinked_groups <- c("X-linked\nXist\nUndetected", "X-linked\nXist-MA\n(Xi=B6)", "X-linked\nXist-MA\n(Xi=Cast)")
g <- temp[temp$group %in% Xlinked_groups,] %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = hvalpha, size = linesize) + 
  geom_boxplot(aes(x = group, y = log2(b6_cast_ratio), color = factor(day)), size = violin_box_size,
               outlier.size = outliersize, outlier.alpha = 1/4, outlier.shape = NA,
               position = position_dodge2(preserve = "single")) + 
  scale_color_manual(values=time_colors[-1]) +
  labs(y = expression("["*sum[g]*"("*B6[gc]*") + 0.01] / ["*sum[g]*"("*Cast[gc]*") + 0.01] ["*log[2]*"(value)]"), x = "", fill = "",
       title = "Gene-wise log2(total B6/total Cast) expression", color = "Time [days]") +
  scale_y_continuous(breaks = seq(-10, 10, by = 5), limits = c(-7, 12))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S9_B_B6CastRatio_XistclassDay_Xlinked.pdf"))



print("9.3) C: Differential silencing analysis - varying XP threshold and number of bins")

print("9.3.1) Load function")
diff_silencing <- diff_silencing <- function(x = df_sk, 
                                             minsilencing_perc = 10,
                                             nbins = 10, minbin = 5, mincount = 25, mincells = 5){
  
  print(paste0("0) subset to Xist-MA cells only & define inactive counts..."))
  XistMAgroup <- c("BL6_MA", "Cast_MA")
  id <- paste0(x$Time, "_", x$Cell)
  temp <- data.frame(id, x)
  data <- temp[temp$Xist %in% XistMAgroup,]
  data <- data[data$both > 0,]
  
  print(paste0("0bis) compute X silencing percentage..."))
  offset_XiXa <- 0.1
  data <- data[!data$Gene %in% c("Xist"),] %>% #don't include Xist in XCR
    dplyr::group_by(id) %>%
    dplyr::mutate(b6_X = sum(b6[Chromosome == "X"]),
                  cast_X = sum(cast[Chromosome == "X"]))
  data$AXCR <- (1-((data$b6_X + offset_XiXa)/(data$cast_X + offset_XiXa)))*100
  data$AXCR[data$Xist == "Cast_MA"] <- (1-((data$cast_X[data$Xist == "Cast_MA"] + offset_XiXa)/(data$b6_X[data$Xist == "Cast_MA"] + offset_XiXa)))*100
  
  print(paste0("0ter) exclude cells with silencing below ", minsilencing_perc, "%..."))
  data <- data[data$AXCR >= minsilencing_perc,]
  
  print(paste0("0quat) breaks for percentage of silencing"))
  data$bin_XCR <- quantileCut(x = data$AXCR, n = nbins)
  levels_bin <- levels(data$bin_XCR)
  
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
  XistUnd_XCRthr_low <- 0.4; XistUnd_XCRthr_high <- 0.6; XistUndTime <- "0hrs"
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
  mincount_cellfit_bin <- 5
  data_cell <- data_cell[data_cell$both>mincount_cellfit_bin,]
  data_cell$XiXa <- (data_cell$inactive_count + offset_XiXa)/(data_cell$active_count + offset_XiXa)
  m <- match(data_cell$Gene, baseline$Gene); table(is.na(m))
  data_cell$baseline_XiXa <- baseline$b6cast_ratio[m]
  data_cell <- data_cell[!is.na(data_cell$baseline_XiXa),]
  data_cell$baseline_XiXa[data_cell$Xist == "Cast_MA"] <- 1/data_cell$baseline_XiXa[data_cell$Xist == "Cast_MA"]
  data_cell$XiXa_ratio_XistUnd <- data_cell$XiXa/data_cell$baseline_XiXa
  ### single cell
  
  print(paste0("7) take log ratios"))
  data_filt$logXiXa <- log2(data_filt$XiXa_ratio_XistUnd)
  
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
  allgenes <- unique(as.character(data_filt$Gene))
  
  model_fit <- c()
  
  for(i in seq_len(length(allgenes))){
    gene <- allgenes[i]
    
    temp <- data_filt[data_filt$Gene %in% gene,]
    
    # bin analysis
    if((length(unique(temp$Xist))==2) & (gene %in% test_both)){
      complete <- lm(logXiXa ~ 0 + AXCR + AXCR:Xist, data = temp)
      oneslope <- lm(logXiXa ~ 0 + AXCR, data = temp)
      
      int_b6 <- 0; slope_b6 <- sum(complete$coefficients[c("AXCR", "AXCR:XistBL6_MA")])
      int_cast <- 0; slope_cast <- complete$coefficients[c("AXCR")]
      
      anv <- anova(oneslope, complete)
      slope_pvalue <- anv$`Pr(>F)`[2]
    }else{
      int_cast <- int_b6 <- slope_cast <- slope_b6 <- slope_pvalue <- NA
    }
    
    #single allele analysis
    
    # test b6
    if(gene %in% test_b6){
      fit <- lm(logXiXa ~ 0 + AXCR, data = temp[temp$Xist == "BL6_MA",])
      int_b6 <- 0; slope_b6 <- fit$coefficients["AXCR"]; bic_b6 <- BIC(fit)
    }else{
      int_b6 <- slope_b6 <- bic_b6 <- NA
    }
    
    # test cast
    if(gene %in% test_cast){
      fit <- lm(logXiXa ~ 0 + AXCR, data = temp[temp$Xist == "Cast_MA",])
      int_cast <- 0; slope_cast <- fit$coefficients["AXCR"]; bic_cast <- BIC(fit)
    }else{
      int_cast <- slope_cast <- bic_cast <- NA
    }
    
    # store results
    model_fit <- rbind(model_fit, 
                       data.frame(Gene = gene, 
                                  int_b6 = int_b6, slope_b6 = slope_b6,
                                  int_cast = int_cast, slope_cast = slope_cast,
                                  slope_pvalue = slope_pvalue))
  }
  model_fit <- model_fit[order(model_fit$slope_pvalue, decreasing = FALSE),]
  model_fit$FDR <- p.adjust(model_fit$slope_pvalue, method = "BH")
  
  print(paste0("9) Return results..."))
  return <- list(model_fit = model_fit,
                 baseline = baseline, 
                 data = data_filt,
                 data_cell = data_cell)
  return(return)
}

print("9.3.2) Load data")
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

print("9.3.3) Launch function varying XP threshold")
sigthreshold <- 0.05

results <- c()
for(minsilencing_perc in seq(0, 20, by = 2)){
  # run function
  binresults <- diff_silencing(x = df_sk, minsilencing_perc = minsilencing_perc)
  
  # store results
  x <- binresults$model_fit[!is.na(binresults$model_fit$FDR),]
  variables <- c("Gene", "slope_b6", "slope_cast", "slope_pvalue", "FDR")
  results <- rbind(results,
                   data.frame(x[, variables], silencing_threshold_percentage = minsilencing_perc))
}
results[results$FDR <= sigthreshold,] %>% arrange(silencing_threshold_percentage, FDR)

print("9.3.4) Plot")

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
            savefile = paste0(outpath, "S9_C_DSanalysis_varyXPthreshold.pdf"), 
            height = 10, width = 5)

print("9.3.5) Launch function varying number of bins")

sigthreshold <- 0.05
results <- c()
for(nbins in seq(6, 20, by = 2)){
  # run function
  binresults <- diff_silencing(x = df_sk, nbins = nbins)
  
  # store results
  x <- binresults$model_fit[!is.na(binresults$model_fit$FDR),]
  variables <- c("Gene", "slope_b6", "slope_cast", "FDR")
  results <- rbind(results,
                   data.frame(x[, variables], nbins = nbins))
}
results[results$FDR <= sigthreshold,] %>% arrange(nbins, FDR)

print("9.3.6) Plot")

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
            savefile = paste0(outpath, "S9_D_DSanalysis_varyNbins.pdf"), 
            height = 10, width = 5)










print("10) Supplementary figure 11: Independent validation of differential silencing dynamics")

print("10.1) A: Pyrosequencing - Xist B6-ratio on dXic lines")

print("10.1.1) Load data")
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

print("10.1.2) Plot")
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
            savefile = paste0(outpath, "S11_A_PyrodXic_Xist_B6ratio.pdf"), 
            width = 3, height = 2)


print("10.2) B: qPCR - Xist relative expression")

print("10.2.1) Load data")
qpcr <- data[(data$Experiment %in% "qPCR")&(!data$group %in%  "WT"),]
color_qpcr <- c(rev(color_alleles)[1:2], "WT" = "black")

print("10.2.2) Plot")
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
            savefile = paste0(outpath, "S11_B_qPCRdXic_Xist_RelExp.pdf"),
            width = 3, height = 2)



print("10.3) C: FISH - Xist+ cells in dXic lines")

print("10.3.1) Load data")
fish <- data[(data$Experiment %in% "FISH")&(!data$group %in% "TX1072"),]

print("10.3.2) Plot")
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
            savefile = paste0(outpath, "S11_C_FISHdXic.pdf"),
            width = 3, height = 2)


print("11.4) D: Xi/Xa ratio in two cell lines")

print("11.4.1) Load data")
load(paste0(datapath, "NGF_DGE_dXic_B6.RData")); b6 <- dge[!rownames(dge) %in% "Xist_5prime",]
load(paste0(datapath, "NGF_DGE_dXic_Cast.RData")); cast <- dge[!rownames(dge) %in% "Xist_5prime",]
df <- data.frame(day = rep(b6$samples$day, each = nrow(b6)), 
                 sample = rep(colnames(b6), each = nrow(b6)),
                 cell_line = rep(b6$samples$dXic, each = nrow(b6)),
                 Chr = rep(b6$genes$chromosome, times = ncol(b6)),
                 Gene = rep(b6$genes$symbol, times = ncol(b6)),
                 Ensembl = rep(b6$genes$ensembl, times = ncol(b6)),
                 b6 = c(b6$counts), cast = c(cast$counts))
df <- df %>% 
  dplyr::group_by(sample) %>% 
  dplyr::mutate(libsize_b6 = sum(b6), libsize_cast = sum(cast)) %>%
  as.data.frame()

print("11.4.2) Remove genes on deletion locus")
df$Xi <- ifelse(df$cell_line == "dB6", df$cast, df$b6)
df$Xa <- ifelse(df$cell_line == "dB6", df$b6, df$cast)
dXic_locus <- 103182257:103955531
dXic_genes <- bm[(bm$chromosome_name == "X")&(bm$start_position %in% dXic_locus | bm$end_position %in% dXic_locus),]
g <- unique(dXic_genes$mgi_symbol); e <- unique(dXic_genes$ensembl_gene_id)
df <- df[!df$Ensembl %in% e,]

print("11.4.3) Compute Xi/Xa ratio and test for each time point")
temp <- df %>% 
  dplyr::group_by(day, cell_line, sample) %>% 
  dplyr::summarise(Xi_X = sum(Xi[(Chr %in% "X")]),
                   Xi_Tot = sum(Xi),
                   Xa_X = sum(Xa[(Chr %in% "X")]),
                   Xa_Tot = sum(Xa)) %>% 
  as.data.frame()
temp$XiXa <- temp$Xi_X/temp$Xa_X
temp$cl <- revalue(factor(temp$cell_line), replace = c("dB6" = "Xi = Cast [dXic on B6]",
                                                       "dCast" = "Xi = B6 [dXic on Cast]"))
tt <- ddply(temp, .variables = .(day), summarise, 
            t_pvalue = t.test(x = log2(XiXa[cell_line == "dB6"]), 
                              y = log2(XiXa[cell_line == "dCast"]), 
                              paired = FALSE)$p.value)
tt$pvalue <- ifelse(tt$t_pvalue >= 0.01, paste0("p = ", round(tt$t_pvalue, digits = 2)), 
                    ifelse(tt$t_pvalue >= 0.001, paste0("p = ", round(tt$t_pvalue, digits = 3)),
                           "p < 0.001"))

print("11.4.3) Plot")
g <- temp %>% 
  ggplot() + 
  theme_bw() + theme1 + 
  geom_jitter(aes(x = factor(day), y = log2(XiXa), color = cl),
              alpha = 0.5, size = 1,
              position=position_jitterdodge(jitter.width = .1, dodge.width = 0.9), shape = 16) +
  stat_summary(fun=mean, aes(x = factor(day), y = log2(XiXa), ymin=..y.., ymax=..y.., color = cl), 
               geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2, show.legend = F) +
  guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
  labs(x="Time [days]", y = "Xi/Xa [log2]\nexcluding genes on dXic\n[chrX:103,182,257-103,955,531]", 
       color = "") +
  scale_color_manual(values = c("#d95f02", "#1b9e77")) +
  geom_text(data = tt, 
            aes(x = day+1, y = 0.25, label = pvalue),
            size = geomtext_size, color = "black")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S11_D_XiXaRatio_deltaXic.pdf"), height = 2)


print("11.5) F: Pyrosequencing - Control genes: B6/Total")

print("11.5.1) Plot")
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
            savefile = paste0(outpath, "S11_F_PyrodXic_ControlGenes_B6ratio.pdf"), 
            width = 10, height = 2)


print("11.6) E: Pyrosequencing - Control genes: normalized Xi:Xa values")

print("11.6.1) Load data")
temp <- pyroseq[(pyroseq$is_hit != "DS gene")&(!pyroseq$Gene %in% "Xist"),]
temp$Xiperc <- ifelse(temp$group == "Xist-MA (Xi=B6)", 
                      100-temp$Value, 
                      temp$Value)

print("11.6.2) Compute Xi:Xa ratios, normalize to d0 values, average across replicates")

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

print("11.6.3) Plot")
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
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S11_E_PyrodXic_ControlGenes_Avglog2NormRatios.pdf"),
            width = 5, height = 3)


print("11.7) G: Pyrosequencing - DS genes: B6/Total")

print("11.7.1) Plot")
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
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S11_G_PyrodXic_DSgenes_B6ratio.pdf"), 
            width = 10, height = 3)


print("11.7.2) F: Pyrosequencing - DS genes: normalized Xi:Xa values")

print("11.7.3) Load data")
temp <- pyroseq[(pyroseq$is_hit %in% "DS gene"),]
temp$Xiperc <- ifelse(temp$group == "Xist-MA (Xi=B6)", 
                      100-temp$Value, 
                      temp$Value)

print("11.7.4) Compute Xi:Xa ratios, normalize to d0 values, test for each gene")

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

print("11.7.5) Plot")
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
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "S11_G_PyrodXic_DSgenes_log2NormRatios.pdf"),
            width = 10, height = 3)










print("12) Supplementary tables")

print("12.1) ST1 - Cell and Gene filtering")

print("12.2) ST2 - Cell and Gene measures")

print("12.3) ST3 - DE and Correlation analyses")

f4path <- paste0(path, "output/fig4_deXistHighLow/")
suppath <- paste0(path, "output/fig_Supplementary/")

print("12.3.1) Load data")

# Xist - DE analysis
load(paste0(f4path, "ALLresults.RData"))
fields <- c("day", "chromosome_name", "mgi_symbol", "coef", "Pr..Chisq.", "fdr")
xist_de <- all_de[, fields]

# Xist - Correlation analysis
load(paste0(f4path, "sprmcor.RData"))
colnames(sprmcor) <- c("day", "mgi_symbol", "chromosome_name", "correlation", 
                       "droprate", "droprate_Xist", "n", "pvalue_Xist", "fdr")
fields <- c("day", "chromosome_name", "mgi_symbol", "correlation", "pvalue_Xist", "fdr")
xist_cor <- sprmcor[, fields]

# Xchr Change - DE analysis
load(paste0(suppath, "XchrChange_HighLow_MAST.RData"))
fields <- c("day", "chromosome_name", "mgi_symbol", "coef", "Pr..Chisq.", "fdr")
xchr_de <- all_de[, fields]

# Xchr Change - Correlation analysis
load(paste0(suppath, "cor_Xchr.RData"))
colnames(sprmcor) <- c("day", "mgi_symbol", "chromosome_name", "correlation", 
                       "droprate", "n", "ypos", "pvalue", "fdr")
fields <- c("day", "chromosome_name", "mgi_symbol", "correlation", "pvalue", "fdr")
xchr_cor <- sprmcor[, fields]

print("12.3.2) Combine results")
names <- c("day", "chromosome", "gene", "value", "pvalue", "fdr") 
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
                  x1$pvalue[match(genes, x1$gene)],
                  x1$fdr[match(genes, x1$gene)],
                  x2$value[match(genes, x2$gene)],
                  x2$pvalue[match(genes, x2$gene)],
                  x2$fdr[match(genes, x2$gene)],
                  x3$value[match(genes, x3$gene)],
                  x3$pvalue[match(genes, x3$gene)],
                  x3$fdr[match(genes, x3$gene)],
                  x4$value[match(genes, x4$gene)],
                  x4$pvalue[match(genes, x4$gene)],
                  x4$fdr[match(genes, x4$gene)])
  x <- as.matrix(x)
  colnames(x) <- c("Day", "Chromosome", "Gene Symbol", "Average CPM",
                   "log2FC [MAST: Xist High vs Low]", "P-value [MAST: Xist High vs Low]", "FDR [MAST: Xist High vs Low]",
                   "Coefficient [Spearman: Gene CPM vs Xist CPM]", "P-value [Spearman: Gene CPM vs Xist CPM]", "FDR [Spearman: Gene CPM vs Xist CPM]",
                   "log2FC [MAST: Xchr.Change High vs Low]", "P-value [MAST: Xchr.Change High vs Low]", "FDR [MAST: Xchr.Change High vs Low]",
                   "Coefficient [Pearson: Gene CPM vs Xchr.Change]", "P-value [Pearson: Gene CPM vs Xchr.Change]", "FDR [Pearson: Gene CPM vs Xchr.Change]")
  if(d==1){
    wb <- createWorkbook()
  }
  addWorksheet(wb, sheetName = paste0("Day ", d))
  writeData(wb, sheet = paste0("Day ", d), x = x, rowNames = F, colNames = T, keepNA = T)
}
saveWorkbook(wb, file = paste0(outpath, "ST3_DEA_Correlation.xlsx"), overwrite = T)

print("12.4) ST4 - Experimental Design")