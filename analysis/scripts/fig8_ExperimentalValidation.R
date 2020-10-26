print("1) Define output path")
outpath <- paste0(path, "output/fig8_ExperimentalValidation/"); dir.create(path = outpath, recursive = T, showWarnings = F)



print("2) B: Xist expression in dXic line")

print("2.1) Load data")
load(paste0(datapath, "DGE_dXic.RData")); table(dge$samples$dXic, dge$samples$day)

print("2.2) Compare Xist expression between two cell lines")
temp <- data.frame(id = dge$samples$id,
                   day = dge$samples$day,
                   cell_line = dge$samples$dXic,
                   libsize = colSums(dge$counts),
                   Xist = dge$counts[dge$genes$symbol %in% "Xist",])
temp$Xist_cpm <- (temp$Xist/temp$libsize)*1e6
test <- temp %>% dplyr::group_by(day) %>% 
  dplyr::summarise(t_pvalue = t.test(x = Xist_cpm[cell_line == "dB6"],
                                     y = Xist_cpm[cell_line == "dCast"])$p.value) %>%
  as.data.frame()
test$pvalue <- ifelse(test$t_pvalue >= 0.01, paste0("p = ", round(test$t_pvalue, digits = 2)), 
                      ifelse(test$t_pvalue >= 0.001, paste0("p = ", round(test$t_pvalue, digits = 3)),
                             "p < 0.001"))
temp$cell_line <- revalue(temp$cell_line, replace = c("dB6" = "Xi = Cast [dXic on B6]",
                                                      "dCast" = "Xi = B6 [dXic on Cast]"))

print("2.3) Plot")
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
  geom_text(data = test, aes(x = factor(day), y = 4, label = pvalue), size = geomtext_size, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 4.1), breaks = seq(0, 4, 1)) +
  labs(x = "Time [days]", y = expression("Xist CPM + 1 [ "*log[10]*" ]"), color = "Cell line")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "B_Xist_deltaXic.pdf"), height = 2)



print("3) C: Xi/Xa ratio in two cell lines")

print("3.1) Load data")
load(paste0(datapath, "DGE_dXic_B6.RData")); b6 <- dge[!rownames(dge) %in% "Xist_5prime",]
load(paste0(datapath, "DGE_dXic_Cast.RData")); cast <- dge[!rownames(dge) %in% "Xist_5prime",]
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

print("3.2) Remove genes on deletion locus")
df$Xi <- ifelse(df$cell_line == "dB6", df$cast, df$b6)
df$Xa <- ifelse(df$cell_line == "dB6", df$b6, df$cast)
dXic_locus <- 103182257:103955531
dXic_genes <- bm[(bm$chromosome_name == "X")&(bm$start_position %in% dXic_locus | bm$end_position %in% dXic_locus),]
g <- unique(dXic_genes$mgi_symbol); e <- unique(dXic_genes$ensembl_gene_id)
df <- df[!df$Ensembl %in% e,]

print("3.3) Compute Xi/Xa ratio and test for each time point")
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


print("3.3) Plot")
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
            savefile = paste0(outpath, "C_XiXaRatio_deltaXic.pdf"), height = 2)



print("4) D: Xi/Xa ratio for differentially silenced genes")

print("4.1) Load data")
load(paste0(datapath, "DGE_dXic_B6.RData")); b6 <- dge[!rownames(dge) %in% "Xist_5prime",]
load(paste0(datapath, "DGE_dXic_Cast.RData")); cast <- dge[!rownames(dge) %in% "Xist_5prime",]
df <- data.frame(day = rep(b6$samples$day, each = nrow(b6)), 
                 sample = rep(colnames(b6), each = nrow(b6)),
                 Line = rep(b6$samples$dXic, each = nrow(b6)),
                 Gene = rep(b6$genes$symbol, times = ncol(b6)),
                 b6 = c(b6$counts), cast = c(cast$counts))
genes <- c("Klhl13", "Pir", "Hprt")
temp <- df[df$Gene %in% genes,]

print("4.2) Compute Xi/Xa expression ratio, normalize to d0 and test difference between two cell lines")
temp$Xi <- ifelse(temp$Line == "dB6", temp$cast, temp$b6)
temp$Xa <- ifelse(temp$Line == "dB6", temp$b6, temp$cast)
temp$ratio <- temp$Xi/temp$Xa
temp <- temp %>% 
  dplyr::group_by(Gene, Line) %>% 
  dplyr::mutate(d0_avgratio = mean(ratio[day==0])) %>%
  as.data.frame()
temp$norm_ratio <- temp$ratio/temp$d0_avgratio

tt <- ddply(temp, .variables = .(day, Gene), summarise, 
            t_pvalue = t.test(x = norm_ratio[Line == "dB6"], 
                              y = norm_ratio[Line == "dCast"], 
                              paired = FALSE)$p.value) %>% 
  arrange(Gene)
tt$pvalue <- ifelse(tt$t_pvalue >= 0.01, paste0("p = ", round(tt$t_pvalue, digits = 2)), 
                    ifelse(tt$t_pvalue >= 0.001, paste0("p = ", round(tt$t_pvalue, digits = 3)),
                           "p < 0.001"))
temp$cell_line <- revalue(temp$Line, replace = c("dB6" = "Xi = Cast [dXic on B6]",
                                                 "dCast" = "Xi = B6 [dXic on Cast]"))

print("4.3) Plot")
g <- temp %>% 
  ggplot() + 
  theme_bw() + theme1 + 
  facet_grid(.~Gene) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = hvalpha, size = linesize) +
  geom_jitter(aes(x = factor(day), y = log2(norm_ratio), color = cell_line),
              alpha = 0.5, size = 1,
              position=position_jitterdodge(jitter.width = .1, dodge.width = 0.9), shape = 16) +
  stat_summary(fun=mean, aes(x = factor(day), y = log2(norm_ratio), ymin=..y.., ymax=..y.., color = cell_line), 
               geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2, show.legend = F) +
  guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
  labs(x="Time [days]", y = "Xi/Xa [log2]", color = "Cell line") +
  scale_color_manual(values = c("#d95f02", "#1b9e77")) +
  geom_text(data = tt, 
            aes(x = day+1, y = 1.5, label = pvalue),
            size = geomtext_size, color = "black")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, 
            savefile = paste0(outpath, "D_XiXaRatio_DSgenes_deltaXic.pdf"), width = 10)