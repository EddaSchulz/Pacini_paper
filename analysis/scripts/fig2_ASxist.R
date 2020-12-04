print("1) Define output path")
outpath <- paste0(path, "output/fig2_ASXistAnalysis/"); dir.create(path = outpath, recursive = T, showWarnings = F)


print("2) A: Xist allele specific (B6/Cast) UMI expression")

print("2.1) Load data")
load(paste0(datapath, "DGE_B6.RData")); b6 <- dge
load(paste0(datapath, "DGE_Cast.RData")); cast <- dge

res <- data.frame(day = b6$samples$day, sf = b6$samples$sf_notX, els = b6$samples$eff_libsize_notX,
                  Xist_cast = cast$counts["Xist",], Xist_b6 = b6$counts["Xist",], 
                  Xist_ratio = as.numeric(b6$samples$xist_ratio), Xist_classification = b6$samples$Xist_ratio_class,
                  Xchr_ratio = as.numeric(b6$samples$xchr_ratio),
                  Xist_sign = b6$samples$Xist_class,
                  b6_X = colSums(b6$counts[b6$genes$chromosome %in% "X",]),
                  cast_X = colSums(cast$counts[cast$genes$chromosome %in% "X",]),
                  b6_X_noXist = colSums(b6$counts[(b6$genes$chromosome %in% "X")&(!b6$genes$symbol %in% "Xist"),]),
                  cast_X_noXist = colSums(cast$counts[(cast$genes$chromosome %in% "X")&(!cast$genes$symbol %in% "Xist"),]))
res$Xist_cast_cpm <- (res$Xist_cast/res$els)*1e6; res$Xist_cast_norm <- (res$Xist_cast/res$sf)
res$Xist_b6_cpm <- (res$Xist_b6/res$els)*1e6; res$Xist_b6_norm <- (res$Xist_b6/res$sf)
res$Xist <- revalue(res$Xist_sign, replace = c("Detected (Xist UMI > 5)" = "Xist+ cells [UMI>5]",
                                              "Detected (Xist UMI <= 5)" = "Xist+ cells [UMI<=5]",
                                              "Undetected" = "Xist- cells [UMI=0]"))
res$Xist <- factor(res$Xist, levels = c("Xist- cells [UMI=0]", "Xist+ cells [UMI<=5]", "Xist+ cells [UMI>5]"))
res$Xist_classification <- revalue(res$Xist_classification, c("Undetected" = "Undetected",
                                                              "Low-Xist" = "Low",
                                                              "Middle" = "Skewed",
                                                              "Xist_BA" = "BA",
                                                              "BL6_MA" = "Xist-MA (Xi=B6)",
                                                              "Cast_MA" = "Xist-MA (Xi=Cast)"))
res$Xist_classification <- factor(res$Xist_classification, levels = c("Undetected", "Low", "Skewed", "BA", "Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)"))
res$Day <- factor(paste0("Day ", res$day))
res$id <- rownames(res)

print("2.2) Plot")

g <- res %>%
  ggplot() + 
  theme_bw() + theme1 + facet_grid(. ~ Day) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = linesize, alpha = hvalpha) + 
  geom_point(aes(x = log10(Xist_b6 + 1), y = log10(Xist_cast + 1), color = Xist_classification), size = scattersize) + 
  geom_point(data = res[res$Xist_classification == "Xist Undetected",],
             aes(x = log10(Xist_b6 + 1), y = log10(Xist_cast + 1), color = Xist_classification), size = scattersize) +
  scale_color_manual(values = color_alleles) + 
  scale_y_continuous(breaks = seq(0,2,1), limits = c(0, 2.5)) +  scale_x_continuous(breaks = seq(0,2,1), limits = c(0, 2.5)) +
  labs(color = "") +  
  xlab(expression(Xist[B6]*" ["*log[10]*"(UMI + 1)]")) + 
  ylab(expression(Xist[Cast]*" ["*log[10]*"(UMI + 1)]")) + 
  guides(color=guide_legend(nrow=2, byrow=TRUE, override.aes = list(size=2))) + 
  theme(legend.position = "top")
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 2, savefile = paste0(outpath, "A_XistUMI.pdf"))


print("3) B: Xist classification over time - Area plot")

print("3.1) Load data")

res <- ddply(res, .variables = .(day), transform, n = length(day))
perc <- res %>%
  dplyr::group_by(day) %>% dplyr::mutate(n = length(day)) %>%
  dplyr::group_by(day, Xist_classification) %>% dplyr::summarise(p = length(day)/unique(n)*100)
perc$Xist_classification <- factor(perc$Xist_classification, 
                                   levels = c("Undetected", "Low", "Xist-MA (Xi=Cast)", "Xist-MA (Xi=B6)", "Skewed", "BA"))
ddply(perc, .variables = .(day), summarize, t = sum(p))

# include 0 values for missing elements at day 0
missing_d0 <- levels(perc$Xist_classification)[!levels(perc$Xist_classification) %in% perc$Xist_classification[perc$day == 0]]
perc <- rbind(data.frame(perc), data.frame(day = 0, Xist_classification = missing_d0, p = 0))

print("3.2) Plot")

g <- ggplot(perc, aes(x=day, y=p, fill=Xist_classification)) + 
  theme_bw() + theme1 +
  geom_area(size=0, colour="white") +
  scale_fill_manual(values=color_alleles) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("Time [days]") + ylab("Cells [%]") +
  guides(fill = guide_legend(title = "", override.aes = list(alpha = 1)))
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, savefile = paste0(outpath, "B_XistClassif_AreaPlot.pdf"))



print("4) C: FISH - Xist positive and BA cells")

print("4.1) Load data")
file <- paste0(datapath, "validation_experiments.txt")
data <- data.frame(read.table(file = file, header = T, sep = "\t"))
data <- data[(data$Experiment %in% "FISH")&(data$CellLine %in% "TX1072"),]
data$Xist <- factor(revalue(factor(data$Variable), replace = c("XistPos" = "Xist+", "XistBi" = "BA")), 
                    levels = c("Xist+", "BA"))

print("4.2) Define percentage of Xist-MA cells per replicate")
temp <- data %>% 
  dplyr::group_by(Day, Replicate) %>%
  dplyr::summarise(perc = Value[Xist == "Xist+"] - Value[Xist == "BA"]) %>%
  as.data.frame()
data <- rbind(data,
              data.frame(Experiment = "FISH", 
                         CellLine = "TX1072",
                         Day = temp$Day,
                         Gene = NA,
                         Variable = "XistMA",
                         Replicate = temp$Replicate,
                         Value = temp$perc,
                         Value_Description = "Cells [%]",
                         Xist = "Xist-MA"))
  

print("4.3) Bar-Plot")
data <- data[!data$Xist %in% "Xist+",]
data$Xist <- factor(data$Xist, levels = c("Xist-MA", "BA"))
cols <- c("#6F72B5", "#fa9fb5")
data$day <- paste0("Day ", data$Day)
data %>% dplyr::group_by(day, Xist) %>% dplyr::summarise(mean = round(mean(Value), digits = 2)) %>% as.data.frame()

# define one bar per replicate for each time point
data$x <- data$Day + (data$Replicate - 2)*0.25
g <- data %>%
  ggplot(aes(x = x, y = Value, fill = Xist)) +
  theme_bw() + theme1 + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = seq(1,4)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Time [days]", 
       y = "Cells [%]",
       fill = "Image\nClassification")
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, 
            savefile = paste0(outpath, "C_FISH_percentage.pdf"))



print("5) D: Xist AS CPM over time - Xist.MA cells")

print("5.1) Load data")

res_xistplus <- res[res$Xist_classification %in% c("Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)"),]
ncells <- ddply(res_xistplus, .variables = .(day), summarize, 
                b6 = sum(Xist_classification == "Xist-MA (Xi=B6)"),
                cast = sum(Xist_classification == "Xist-MA (Xi=Cast)"))
wmwtest <- ddply(res_xistplus[res_xistplus$day != 0,], .variables = .(day), summarize, 
                 B6_MA = sum(Xist_classification == "Xist-MA (Xi=B6)"),
                 Cast_MA = sum(Xist_classification == "Xist-MA (Xi=Cast)"),
                 pvalue = wilcox.test(x = log10(Xist_b6_cpm[Xist_classification == "Xist-MA (Xi=B6)"]+1), 
                                      y = log10(Xist_cast_cpm[Xist_classification == "Xist-MA (Xi=Cast)"]+1), 
                                      paired = FALSE)$p.value)
wmwtest$p <- ifelse(wmwtest$pvalue >= 0.01, paste0("p = ", round(wmwtest$pvalue, digits = 2)), 
                    ifelse(wmwtest$pvalue >= 0.001, paste0("p = ", round(wmwtest$pvalue, digits = 3)),
                           "p < 0.001"))

temp <- res_xistplus[res_xistplus$day != 0, c("day", "Xist_classification", "Xist_b6_cpm", "Xist_cast_cpm")]
temp$Xist <- NA
temp$Xist[temp$Xist_classification %in% "Xist-MA (Xi=B6)"] <- log10(temp$Xist_b6_cpm[temp$Xist_classification %in% "Xist-MA (Xi=B6)"]+1)
temp$Xist[temp$Xist_classification %in% "Xist-MA (Xi=Cast)"] <- log10(temp$Xist_cast_cpm[temp$Xist_classification %in% "Xist-MA (Xi=Cast)"]+1)

print("5.2) Plot")

g <- temp %>% 
  ggplot() +
  theme_bw() +  theme1 + 
  scale_y_continuous(breaks = seq(0, 5, by = 0.5), limits = c(2.4, 4.8)) +
  geom_violin(aes(x = factor(day), y = Xist, colour = Xist_classification), alpha = 0.5, draw_quantiles = c(0.5), size = violin_box_size, show.legend = FALSE) + 
  geom_jitter(aes(x = factor(day), y = Xist, colour = Xist_classification), alpha = 0.5, size = scattersize,
              position=position_jitterdodge(jitter.width = .2, dodge.width = 0.9), shape = 21) +
  geom_text(data = wmwtest[!wmwtest$day %in% "0",], aes(x = factor(day), y = 4.75, label = p), 
            size = geomtext_size, angle = 0) +
  scale_color_manual(values = color_alleles) +
  labs(x = "Time [days]", y = expression("Xist CPM + 1 ["* log[10]*"]"), color = "") +
  guides(color = guide_legend(override.aes = list(size = 1, shape = 20, alpha = 1)))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "D_XistMA_B6vCast.pdf"))