print("1) Define output path")
outpath <- paste0(path, "output/fig6_PyroseqValidation/"); dir.create(path = outpath, recursive = T, showWarnings = F)

print("2) Load data")
file <- paste0(datapath, "PyroseqValidation.xlsx")
data <- data.frame(read_excel(path = file)[, 1:4])
df <- separate(data = data, col = "BioGroupID", into = c("experiment", "group", "time"), sep = "_", extra = "merge")
df$CastSNPfreq <- as.numeric(as.character(df$CastSNPfreq))
df$experiment <- revalue(as.factor(df$experiment), replace = c("dXIC" = "Pyroseq", "qPCR" = "qPCR"))
df$group <- revalue(as.factor(df$group), replace = c("XX" = "WT", "gDNA" = "gDNA", "A1" = "dXIC on B6 allele", "A6" = "dXIC on Cast allele"))

print("3) Process data")
remove_genes <- c("Xist")

# store percentages from inactive allele
df_fc <- df[(df$experiment == "Pyroseq")&(!df$group %in% c("gDNA", "WT"))&(!is.na(df$CastSNPfreq))&(!df$Gene %in% remove_genes),]
df_fc$Xist_classification <- revalue(factor(as.character(df_fc$group)), 
                                     replace = c("dXIC on B6 allele" = "Xist-MA (Xi=Cast)",
                                                 "dXIC on Cast allele" = "Xist-MA (Xi=B6)"))
df_fc$Xiperc <- ifelse(df_fc$Xist_classification == "Xist-MA (Xi=B6)", 100-df_fc$CastSNPfreq, df_fc$CastSNPfreq)
df_fc$day <- as.numeric(gsub(df_fc$time, pattern = "h", replacement = ""))/24

# compute XiXa ratios
df_fc$XiXa <- df_fc$Xiperc/(100-df_fc$Xiperc)

# normalize to average 0h ratio per cell line
df_fc_scaled <- df_fc %>%
  dplyr::group_by(Gene, Xist_classification) %>%
  dplyr::mutate(mean_baseline_0h = mean(XiXa[time=="0h"]),
                XiXa_scaled = XiXa/mean(XiXa[time=="0h"]))

# ANOVA test
anova_test <- c()
for(g in 1:length(unique(df_fc_scaled$Gene))){
  gene <- as.character(unique(df_fc_scaled$Gene))[g]
  temp <- df_fc_scaled[df_fc_scaled$Gene == gene,]
  complete <- lm(log2(XiXa_scaled) ~ 0 + day + day:Xist_classification, data = temp)
  oneslope <- lm(log2(XiXa_scaled) ~ 0 + day, data = temp)
  anv <- anova(oneslope, complete)
  pvalue <- anv$`Pr(>F)`[2]
  
  beta2 <- coef(complete)['day:Xist_classificationXist-MA (Xi=Cast)']
  beta2_ci <- confint(complete, level = 0.95)['day:Xist_classificationXist-MA (Xi=Cast)',]
  
  # intercept <- coefficients(complete)["(Intercept)"]
  slope_XiB6 <- coefficients(complete)["day"]
  slope_XiCast <- coefficients(complete)["day:Xist_classificationXist-MA (Xi=Cast)"] + slope_XiB6
  anova_test <- rbind(anova_test, 
                      data.frame(Gene = gene, 
                                 intercept = 0, slope_XiB6 = slope_XiB6, slope_XiCast = slope_XiCast, pvalue = pvalue,
                                 slope_diff = beta2, 
                                 slope_diff_CIlow = beta2_ci[1], slope_diff_CIhigh = beta2_ci[2]))
}
anova_test$fdr <- p.adjust(anova_test$pvalue, method = "BH")
anova_test <- anova_test %>% dplyr::arrange(fdr)
anova_test$Gene <- factor(anova_test$Gene, levels = unique(as.character(anova_test$Gene)))

print("4) B: Plot - control genes")

controlgene <- c('Rnf12','Renbp','Cul4b','Prdx4','Atrx')
mr <- df_fc_scaled[df_fc_scaled$Gene %in% controlgene,] %>%
  dplyr::group_by(day, Xist_classification, Gene) %>% 
  dplyr::summarise(meanratio = mean(XiXa_scaled))

wt <- ddply(mr[mr$day != 0,], .variables = .(day), summarise, 
            wilcox_pvalue = wilcox.test(x = meanratio[Xist_classification == "Xist-MA (Xi=Cast)"], 
                                        y = meanratio[Xist_classification == "Xist-MA (Xi=B6)"], 
                                        paired = TRUE)$p.value)
g <- mr  %>%
  ggplot(aes(x = day, y = log2(meanratio), color = Xist_classification, fill = Xist_classification)) +
  theme_bw() + theme1 +  
  geom_text(data = data.frame(wt, Xist_classification = NA), aes(x = day, y = 1, label = formatC(wilcox_pvalue, digits = 2)),
            size = geomtext_size, color = "black") +
  geom_hline(yintercept = 0, size = linesize, linetype = "dashed", alpha = 1/4) +
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75),
              size = small_scattersize, show.legend = FALSE, alpha = 1/2) +
  stat_summary(fun=median, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = color_alleles) +
  labs(x = "Time [days]", 
       y = expression("Average scaled ratio: "*"("*X[i]*"/"*X[a]*") [log2]"),
       color = "") +
  scale_y_continuous(limits = c(-2.5, 1.5))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "B_Averagelog2Ratios_ControlGenes.pdf"))


print("5) C/D: Plot - differentially silenced genes")

anv <- gather(anova_test[, c("Gene", "intercept", "slope_XiB6", "slope_XiCast", 'pvalue')], 
              'slope', 'value', -Gene, -intercept, -pvalue)
anv$Xist_classification <- revalue(factor(anv$slope), replace = c('slope_XiB6' = 'Xist-MA (Xi=B6)',
                                                                  'slope_XiCast' = 'Xist-MA (Xi=Cast)'))
anv <- anv[order(anv$pvalue, decreasing = FALSE),]
anv$genelab <- factor(paste0(anv$Gene, '\npvalue=', signif(anv$pvalue, digits = 3)))
anv$genelab <- factor(anv$genelab, levels = unique(anv$genelab))
df_fc_scaled$genelab <- anv$genelab[match(df_fc_scaled$Gene, anv$Gene)]
df_fc_scaled$genelab <- factor(df_fc_scaled$genelab, levels = levels(anv$genelab))
hits <- c("Klhl13", "Pir", "Hprt")

tt <- ddply(df_fc_scaled[df_fc_scaled$Gene %in% hits,], .variables = .(day, Gene), summarise, 
            t_pvalue = t.test(x = XiXa_scaled[Xist_classification == "Xist-MA (Xi=Cast)"], 
                       y = XiXa_scaled[Xist_classification == "Xist-MA (Xi=B6)"], 
                       paired = FALSE)$p.value) %>% arrange(Gene)

# Klhl13
g <- df_fc_scaled[df_fc_scaled$Gene %in% "Klhl13",] %>%
  ggplot(aes(x = day, y = log2(XiXa_scaled), color = Xist_classification, fill = Xist_classification)) +
  theme_bw() + theme1 +  
  facet_grid(.~Gene) +
  geom_text(data = data.frame(tt[tt$Gene %in% "Klhl13",], Xist_classification = NA), 
            aes(x = day, y = 1.75, label = formatC(t_pvalue, digits = 1)),
            size = geomtext_size, color = "black") +
  geom_hline(yintercept = 0, size = linesize, linetype = "dashed", alpha = 1/4) +
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75),
              size = small_scattersize, show.legend = FALSE, alpha = 1/2) +
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = color_alleles) +
  labs(x = "Time [days]", 
       y = expression("Scaled ratio: "*"("*X[i]*"/"*X[a]*") [log2]"),
       color = "") +
  scale_y_continuous(limits = c(-2, 2))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "C_XiXa_ratios_0hScaled_ttest_Klhl13.pdf"))

# Pir/Hprt
g <- df_fc_scaled[df_fc_scaled$Gene %in% c("Pir", "Hprt"),] %>%
  ggplot(aes(x = day, y = log2(XiXa_scaled), color = Xist_classification, fill = Xist_classification)) +
  theme_bw() + theme1 +  
  facet_grid(.~Gene) +
  geom_text(data = data.frame(tt[tt$Gene %in% c("Pir", "Hprt"),], Xist_classification = NA), 
            aes(x = day, y = 1.75, label = formatC(t_pvalue, digits = 1)),
            size = geomtext_size, color = "black") +
  geom_hline(yintercept = 0, size = linesize, linetype = "dashed", alpha = 1/4) +
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75),
              size = small_scattersize, show.legend = FALSE, alpha = 1/2) +
  stat_summary(fun=mean, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  scale_color_manual(values = color_alleles) +
  labs(x = "Time [days]", 
       y = expression("Scaled ratio: "*"("*X[i]*"/"*X[a]*") [log2]"),
       color = "") +
  scale_y_continuous(limits = c(-2, 2))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "D_XiXa_ratios_0hScaled_ttest_PirHprt.pdf"), width = 10)
