print("1) Define output path")
outpath <- paste0(path, "output/fig7_GeneAlleleSpecificSilencing/"); dir.create(path = outpath, recursive = T, showWarnings = F)

print("2) A: Xi/Xa ratio for MA-B6 and MA-Cast cells")

print("2.1) Load data and perform test")

# compute Xi/Xa ratios for Xist-MA cells only - excluding Xist reads
ma <- res[res$Xist_classification %in% c("Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)"),]
ma$XiXaratio <- ma$b6_X_noXist/ma$cast_X_noXist
w <- ma$Xist_classification == "Xist-MA (Xi=Cast)"
ma$XiXaratio[w] <- 1/ma$XiXaratio[w]

# ncells
ncells <- ma[!ma$day %in% "0",] %>% 
  dplyr::group_by(day, Xist_classification) %>%
  dplyr::summarise(n = length(day))

# wilcoxon mann whitney test
test <- ma[!ma$day %in% "0",] %>% 
  dplyr::group_by(day) %>%
  dplyr::summarise(n = length(day),
                   ma_b6 = sum(Xist_classification == "Xist-MA (Xi=B6)"),
                   ma_cast = sum(Xist_classification == "Xist-MA (Xi=Cast)"),
                   pvalue_wmw = wilcox.test(x = log2(XiXaratio[Xist_classification == "Xist-MA (Xi=B6)"]), 
                                            y = log2(XiXaratio[Xist_classification == "Xist-MA (Xi=Cast)"]))$p.value)
test$p <- ifelse(test$pvalue_wmw >= 0.01, paste0("p = ", round(test$pvalue_wmw, digits = 2)), 
                 ifelse(test$pvalue_wmw >= 0.001, paste0("p = ", round(test$pvalue_wmw, digits = 3)),
                        "p < 0.001"))

print("2.2) Plot")

g <- ma[!ma$day %in% "0",] %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = hvalpha, size = linesize) + 
  geom_boxplot(aes(x = factor(day), y = log2(XiXaratio), color = Xist_classification), outlier.shape = NA, size = violin_box_size) + 
  geom_jitter(aes(x = factor(day), y = log2(XiXaratio), color = Xist_classification), fill = "black",
              alpha = 1/4, position=position_jitterdodge(jitter.width = .25, dodge.width = 0.75), 
              size = small_scattersize, show.legend = FALSE) +
  scale_color_manual(values = color_alleles) + 
  geom_text(data = test, aes(x = factor(day), y = 2.25, label = p), 
            size = geomtext_size, angle = 0) +
  geom_text(data = ncells,
            mapping = aes(x = factor(day), y = Inf, label = paste0("n = ", n), group = Xist_classification),
            position=position_dodge(.75), angle = 90, alpha = 1, size = geomtext_size, hjust = -0.5) +
  scale_y_continuous(limits = c(-7.5, 3), breaks = seq(-7.5, 2.5, by = 2.5)) +
  coord_cartesian(clip = "off") + 
  labs(x = "Time [days]", y = expression(log[2]*"("*X[i.chrX]*"/"*X[a.chrX]*")"), color = "")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "A_XiXaRatio_XistMA.pdf"))


print("3) C: XCI time vs PDT in MA-B6 and MA-Cast cells, colored by sequencing day")

print("3.1) Load data")
offset_XiXa <- 0.1; minsilencing_perc <- 10
r <- res[res$Xist_classification %in% c("Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)"),]
m <- match(r$id, pData(XX)$id)
r$scaledPDT <- pData(XX)$Scaled_PDT[m]
r$Xi <- ifelse(r$Xist_classification == "Xist Undetected", apply(r[, c("b6_X_noXist", "cast_X_noXist")], 1, min),
               ifelse(r$Xist_classification == "Xist-MA (Xi=B6)", r$b6_X_noXist, r$cast_X_noXist))
r$Xa <- ifelse(r$Xist_classification == "Xist Undetected", apply(r[, c("b6_X_noXist", "cast_X_noXist")], 1, max),
               ifelse(r$Xist_classification == "Xist-MA (Xi=B6)", r$cast_X_noXist, r$b6_X_noXist))
r$XT <- (1-((r$Xi + offset_XiXa)/(r$Xa + offset_XiXa)))*100
r$XT_floored <- r$XT; r$XT_floored[r$XT_floored<0] <- 0
r$AXCR <- abs((r$Xi/(r$Xi+r$Xa))-0.5)

print("3.2) Plot")
temp <- r; temp$XT[temp$XT < 0] <- 0
temp$silencing <- factor(ifelse(temp$XT<minsilencing_perc, "excluded", "included"))

# plot results: XT vs PDT
g <- temp %>%
  ggplot() + 
  theme_bw() + theme1 + 
  facet_grid(.~Xist_classification) + 
  geom_vline(xintercept = minsilencing_perc, linetype = "dashed", size = linesize, alpha = 1/4) +
  geom_point(aes(x = XT, y = scaledPDT, color = factor(day), shape = silencing, alpha = silencing), size = scattersize) + 
  scale_shape_manual(values = c(21, 20)) +
  scale_alpha_manual(values = c(1/4, 1)) +
  scale_color_manual(values = time_colors) + 
  labs(x = "XCI progress [%]",
       y = "Scaled Pseudotime", 
       color = "Time [days]") +
  guides(alpha = FALSE, shape = FALSE) +
  scale_x_continuous(breaks = c(0, 10, 50, 100), limits = c(0, 100))
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, savefile = paste0(outpath, "C_XT_PDT.pdf"))


print("4) D: Allele specific silencing dynamics for an example gene (Eif1ax)")

print("4.1) Load data to perform differential silencing analysis")
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

print("4.2) Define function to perform differential silencing analysis")
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

print("4.3) Launch function")

sigthreshold <- 0.05
mincount <- 25; nbins <- 10
mincells <- 5; minbin <- 5; minsilencing_perc <- 10; mincount_cellfit_bin <- 5
mincount_cellfit <- 5; mincell_cellfit <- 10
XistUndTime <- "0hrs"; XistUnd_XCRthr_low <- 0.4; XistUnd_XCRthr_high <- 0.6
equallysized <- TRUE; zeroIntercept <- TRUE; weighted <- FALSE

binresults <- XiXa_metacell_ratio_final(x = df_sk, mincount = mincount, nbins = nbins, minbin = minbin, mincells = mincells,
                                        minsilencing_perc = minsilencing_perc, mincount_cellfit_bin = mincount_cellfit_bin,
                                        mincount_cellfit = mincount_cellfit, mincell_cellfit = mincell_cellfit,
                                        equallysized = equallysized, zeroIntercept = zeroIntercept, weighted = weighted)
save(binresults, file = paste0(outpath, "diffsilencing_results.RData"))


print("4.4) Load results for Eif1ax gene")

gene <- "Eif1ax"

# cell values
cell_values <- binresults$data_cell[binresults$data_cell$Gene %in% gene,]

# bin values
bin_values <- binresults$data[binresults$data$Gene %in% gene,]

# revalue Xist
revalue_xist <- c("Undetected" = "Undetected", "Low-Xist" = "Low", "Middle" = "Skewed",
                  "Xist_BA" = "BA", "BL6_MA" = "Xist-MA (Xi=B6)", "Cast_MA" = "Xist-MA (Xi=Cast)")
cell_values$XistGroup <- revalue(cell_values$Xist, revalue_xist)
bin_values$XistGroup <- revalue(bin_values$Xist, revalue_xist)

# extract slopes
bin <- binresults$model_fit
signif <- bin[bin$Gene %in% gene,]
slopes <- data.frame(Gene = unique(signif$Gene),
                     slope = c(signif$slope_b6, signif$slope_cast),
                     intercept = c(signif$int_b6, signif$int_cast),
                     XistGroup = rep(c("Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)"), each = nrow(signif)))
slopes$XistGroup <- factor(slopes$XistGroup, levels = c("Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)"))

# compute XT50s
threshold <- 100
bin$T50_b6 <- -(1 + bin$int_b6)/bin$slope_b6; bin$T50_cast <- -(1 + bin$int_cast)/bin$slope_cast
bin$T50_b6_original <- bin$T50_b6; bin$T50_cast_original <- bin$T50_cast; 
bin$T50_b6[(!is.na(bin$T50_b6)) & (bin$T50_b6<0) | bin$T50_b6 > threshold] <- threshold
bin$T50_cast[(!is.na(bin$T50_cast)) & (bin$T50_cast<0) | bin$T50_cast > threshold] <- threshold
s <- bin[bin$Gene %in% gene,]
t50 <- reshape2::melt(s, id.vars = colnames(signif)[!colnames(signif) %in% c("T50_b6", "T50_cast")], 
                      measure.vars = c("T50_b6", "T50_cast"), value.name = "t50")
t50$XistGroup <- revalue(t50$variable, c("T50_b6" = "Xist-MA (Xi=B6)", "T50_cast" = "Xist-MA (Xi=Cast)"))

print("4.5) Plot results for Eif1ax gene")

minvalue_pergene <- cell_values %>% dplyr::group_by(Gene) %>% dplyr::summarise(min = min(log2(XiXa_ratio_XistUnd)))

for(gg in gene){
  
  minvalue <- minvalue_pergene$min[minvalue_pergene$Gene == gg]
  t <- t50[t50$Gene == gg,]
  
  g <- cell_values[cell_values$Gene %in% gg,] %>%
    ggplot() + 
    theme_bw() + theme1 + facet_grid(XistGroup~.) +
    geom_segment(data = t, aes(x = 0, y = log2(1/2), xend = t50, yend = log2(1/2), color = XistGroup), linetype = "dashed",
                 show.legend = FALSE, size = linesize) +
    geom_segment(data = t, aes(x = t50, y = minvalue, xend = t50, yend = log2(1/2), color = XistGroup), linetype = "dashed",
                 show.legend = FALSE, size = linesize) + 
    geom_hline(yintercept = 0, size = linesize, alpha = 1/2, linetype = "dashed") +
    geom_point(aes(x = AXCR, y = log2(XiXa_ratio_XistUnd), color = XistGroup, alpha = XistGroup), 
               size = scattersize, show.legend = FALSE) + 
    geom_point(data = bin_values[bin_values$Gene %in% gg,], 
               aes(x = AXCR, y = log2(XiXa_ratio_XistUnd), color = XistGroup), size = scattersize*binscatter_mult) + 
    geom_abline(data = slopes[slopes$Gene %in% gg,], 
                aes(intercept = intercept, slope = slope, color = XistGroup),
                show.legend = FALSE, size = linesize) + 
    scale_color_manual(values = color_alleles) + 
    scale_alpha_manual(values = c(1/6, 1/4)) +
    labs(x = "XCI progress [%]",
         y = "Xi/Xa [log2]",
         color = "") + 
    scale_x_continuous(breaks = seq(0, 100, 50), limits = c(0, 100))
  adjust_size(g = g, panel_width_cm = 1.5, panel_height_cm = 2, savefile = paste0(outpath, "D_GeneScheme_", gg, ".pdf"), strptext_col = 0)
}


print("5) G: Gene-wise allele specific XP50 values")

print("5.1) Load data")
mod_binned <- binresults$model_fit
mod_binned$SE50_b6 <- (-1 -mod_binned$int_b6)/(mod_binned$slope_b6); mod_binned$SE50_b6[mod_binned$SE50_b6<0] <- 100
mod_binned$SE50_cast <- (-1 -mod_binned$int_cast)/(mod_binned$slope_cast); mod_binned$SE50_cast[mod_binned$SE50_cast<0] <- 100
mod_binned$SE50_b6[mod_binned$SE50_b6>100] <- 100; mod_binned$SE50_cast[mod_binned$SE50_cast>100] <- 100

# cluster SE50 values into 4 clusters
set.seed(0)
kcenters <- 4; kiter <- 1e6
features <- c("Gene", "FDR", "SE50_b6", "SE50_cast")
class_melt_all <- mod_binned[, features] %>%
  tidyr::gather(variable, se50, -Gene, -FDR)
class_melt <- class_melt_all[!is.na(class_melt_all$se50),] %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(k = kmeans(se50, centers = kcenters, iter.max = kiter)$cluster, id = paste0(variable, "_", k))
class_order <- class_melt %>%
  dplyr::group_by(variable, k) %>% dplyr::summarise(m = max(se50), n = length(Gene)) %>%
  dplyr::group_by(variable) %>% dplyr::mutate(r = rank(m), id = paste0(variable, "_", k))
class_melt$k_ordered <- class_order$r[match(class_melt$id, class_order$id)]
class_melt$Xist <- revalue(class_melt$variable, c("SE50_b6" = "BL6_MA", "SE50_cast" = "Cast_MA"))
ddply(class_melt, .variables = .(k_ordered, Xist), summarize, n = length(unique(Gene)), genes = paste(unique(Gene), collapse = ","))
thrs <- class_melt %>% dplyr::group_by(Xist, k_ordered) %>% dplyr::summarise(t = max(se50)); thrs

print("5.2) Plot XP50 values colored by significance")
mod_binned$significant <- mod_binned$FDR<sigthreshold
g <- mod_binned[!is.na(mod_binned$FDR),] %>%
  ggplot() + 
  theme_bw() + theme1 + 
  geom_abline(slope = 1, intercept = 0, size = linesize, alpha = 1/2, linetype = "dashed") + 
  geom_point(aes(x = SE50_b6, y = SE50_cast, color = factor(significant)), size = 0.1, show.legend = FALSE) + 
  geom_text_repel(aes(x = SE50_b6, y = SE50_cast, label = Gene, 
                      color = factor(significant), alpha = factor(significant)),
                  size = geomtext_size,
                  segment.size = linesize, min.segment.length = 0.1, nudge_y = 1, 
                  show.legend = FALSE, force = 10) + 
  scale_alpha_manual(values = c(1/2, 1)) +
  scale_color_manual(values = c("black", "red")) + 
  labs(x = expression(XT[50] * ": "* X[i]*" = B6"),
       y = expression(XT[50] * ": "* X[i]*" = Cast"),
       color = paste0("FDR < ", sigthreshold)) +
  scale_x_continuous(breaks = seq(0, 100, 50), limits = c(0,100)) +
  scale_y_continuous(breaks = seq(0, 100, 50), limits = c(0,100))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 5, savefile = paste0(outpath, "G_SilencingHalfTime.pdf"))




print("6) E: Comparison of estimated XP values to previously determined silencing classes")

print("6.1) Define function to update gene symbols")
ncores <- min(c(detectCores(), 10))
alias2main_par <- function(genesymbol, annot = Mm_syn, cores = ncores){
  
  x <- as.character(genesymbol)
  
  #define synonyms and ensembl 
  symbol <- annot$Symbol
  syn <- strsplit2(annot$Synonyms, split = "\\|")
  ensembl <- strsplit2(strsplit2(annot$dbXrefs, split = "\\|")[,2], split = ":")[,2]
  
  #get results back --> parallel
  # no_cores <- detectCores() - 1  
  no_cores <- cores
  cl <- makeCluster(no_cores)  
  registerDoParallel(cl) 
  
  results <- foreach(i=1:length(x), .combine = rbind) %dopar% {
    
    gene <- x[i]
    gene_lower <- strsplit(gene, split = "")[[1]]
    gene_lower <- paste(c(gene_lower[1], tolower(paste(gene_lower[2:length(gene_lower)], collapse = ""))), collapse = "")
    
    # look in both main and alias symbols
    if(gene %in% symbol | gene_lower %in% symbol){
      w <- unique(c(which(symbol == gene), which(symbol == gene_lower)))
    }else{
      w <- unique(c(which(apply(syn, 1, function(y){any(gene %in% y)})==TRUE),
                    which(apply(syn, 1, function(y){any(gene_lower %in% y)})==TRUE)))
    }
    
    #for some genes a second line is present, but only the first represent the gene itself
    w <- w[1]
    
    # if a gene is not matched print NAs
    if(length(w)==0){
      temp <- c(gene, NA, NA, NA)
    }else{
      temp <- c(gene, symbol[w], ensembl[w], annot$Synonyms[w])
    }
    return(temp)
  }
  results <- data.frame(results)
  colnames(results) <- c("Gene", "UpdatedSymbol", "Ensembl", "Synonyms")
  stopCluster(cl)
  return(results)
}
Mm_syn <- fread(input = paste0(datapath, "Mus_musculus.gene_info"),
                sep = "\t", header = FALSE, skip = 0L)
colnames(Mm_syn) <- as.character(Mm_syn[1,])
Mm_syn <- Mm_syn[-1,]


print("6.2) Load results of previous publications - Sousa et al. 2019")

# Kinetics of Xist-induced gene silencing can be predicted from combinations of epigenetic and genomic features
# Available at: https://www.dropbox.com/sh/6czxpqgx5r1wdhx/AADXp5Bp0FS42H41aWAU7Ni6a?dl=0

file <- paste0(datapath, "Sousa_Supplemental_Table_S4.xlsx")
lisa <- read_excel(path = file, sheet = "PRO-Seq genes", skip = 11)
genes <- alias2main_par(genesymbol = lisa$`gene name`)
genes <- genes[as.character(genes$Gene) != as.character(genes$UpdatedSymbol),]; genes <- genes[!is.na(genes$UpdatedSymbol),]
lisa$gene <- as.character(lisa$`gene name`); m <- match(genes$Gene, lisa$gene); lisa$gene[m] <- as.character(genes$UpdatedSymbol)
lisa$`half-time` <- as.numeric(as.character(lisa$`half-time`))

# gene classification
lisa$t50 <- lisa$`half-time`
lisa$k <- kmeans(lisa$t50, centers = 4, iter.max = 1000, nstart = 1000)$cluster
avet50 <- ddply(lisa, .variables = .(k), summarize, m = mean(t50)); avet50$r <- rank(avet50$m)
lisa$kmeans <- avet50$r[match(lisa$k, avet50$k)]; lisa$kmeans <- factor(lisa$kmeans)
lisa$Classification <- revalue(lisa$kmeans, replace = c("1" = "Early", "2" = "Intermediate", "3" = "Late", "4" = "Escapee"))
lisa$Classification <- factor(lisa$Classification, levels = c("Early", "Intermediate", "Late", "Escapee"))

print("6.2) Load results of previous publications - Marks et al. 2015")

# Dynamics of gene silencing during X inactivation using allele-specific RNA-seq
# Available at: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4546214/

file <- paste0(datapath, "Marks_classif_parsed.xlsx")
marks <- read_excel(path = file, sheet = "Clusters_and_overlap_Lin", skip = 0)
marks$gene <- as.character(alias2main_par(genesymbol = marks$Timepoint)$UpdatedSymbol)
marks <- marks[!is.na(marks$gene),]
marks$Classification <- factor(marks$'Cluster in Marks', levels = c("Early", "Intermediate", "Late", "Escapee"))

print("6.2) Load results of previous publications - Borensztein et al. 2017")

# Xist-dependent imprinted X inactivation and the early developmental consequences of its failure
# Available at: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5337400/

file <- paste0(datapath, "Borensztein_SuppTable1.txt")
bor <- read.table(file = file, header = TRUE)
bor$gene <- as.character(alias2main_par(genesymbol = bor$Genes)$UpdatedSymbol)
bor <- bor[!is.na(bor$gene),]
bor <- gather(bor[,-1], group, class, -gene)
bor <- bor[!bor$class %in% "Bias",]
bor$Classification <- revalue(factor(bor$class), replace = c("Early" = "Early", "Inter" = "Intermediate",
                                                             "Late" = "Late", "Esc" = "Escapee"))

print("6.2) Load results of previous publications - Combine data sets")
features <- c("gene", "Classification")
gc <- rbind(data.frame(source = "Lisa", lisa[, features], compare = "B6"),
            data.frame(source = "Marks", marks[, features], compare = "B6"),
            data.frame(source = "Borensztein", bor[bor$group == "BC", features], compare = "Cast"),
            data.frame(source = "Borensztein", bor[bor$group == "CB", features], compare = "B6"))

print("6.3) Combine with model XP estimates")
class_melt_all$Xist <- revalue(class_melt_all$variable, c("SE50_b6" = "BL6_MA", "SE50_cast" = "Cast_MA"))
cm <- class_melt_all[!is.na(class_melt_all$se50),]
cm$compare <- revalue(factor(cm$Xist), replace = c("BL6_MA" = "B6", "Cast_MA" = "Cast"))
m <- match(paste0(gc$gene, "_", gc$compare), paste0(cm$Gene, "_", cm$compare))
gc$XT50 <- cm$se50[m]
gg <- gc[!is.na(gc$XT50),]
gg$silgroup <- revalue(gg$compare, replace = c("B6" = "Xi = B6", "Cast" = "Xi = Cast"))
gg$id <- paste0(gg$silgroup, "\n", gg$source)
ngenes <- ddply(gg, .variables = .(Classification, id), summarize, n = length(gene))

print("6.4) Plot XP50 values grouped based on previous silencing dynamics classification - Marks & Borensztein")

ngenes_sub <- ngenes[!ngenes$id %in% "Xi = B6\nLisa",]
gg_sub <- gg[!gg$source %in% "Lisa",]
gg_sub$id <- factor(gg_sub$id, levels = c("Xi = B6\nMarks", "Xi = B6\nBorensztein", "Xi = Cast\nBorensztein"))
cols <- c("#fbb4b9", "#f768a1", "#c51b8a", "#7a0177")
g <- gg_sub %>%  
  ggplot(aes(x = id, y = XT50, color = Classification, fill = Classification)) +
  theme_bw() + theme1 +  
  geom_jitter(color = "black",
              alpha = 1/4, position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75), 
              size = small_scattersize, show.legend = FALSE) +
  stat_summary(fun=median, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5,
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  geom_text(data = ngenes_sub, aes(x = id, y = Inf, group = Classification, label = paste0("n=", n)), 
            position = position_dodge(width = 0.75), angle = 90, alpha = 1, size = geomtext_size,
            color = "black", hjust=-0.5) + 
  coord_cartesian(clip = "off") +
  scale_color_manual(values = cols) +
  scale_y_continuous(breaks = seq(0, 100, 50), limits = c(0,100)) +
  labs(x = "", 
       y = expression(XP[50]*"[%]"),
       color = "")
adjust_size(g = g, panel_width_cm = 4, panel_height_cm = 3, savefile = paste0(outpath, "E_Comparison_Marks&Borensztein.pdf"))

print("6.4) Plot XP50 values grouped based on previous silencing dynamics classification - Sousa")
gc_lisa <- data.frame(lisa[, c("gene", "t50")], compare = "B6")
m <- match(paste0(gc_lisa$gene, "_", gc_lisa$compare), paste0(cm$Gene, "_", cm$compare))
gc_lisa$XT50 <- cm$se50[m]
gc_lisa <- gc_lisa[!is.na(gc_lisa$XT50),]

# Classification: according to the XCI-Escape and Silencing-Dynamics models
gc_lisa$XCI_Escape <- NA
gc_lisa$XCI_Escape[gc_lisa$t50<0.9] <- "Silenced"
gc_lisa$XCI_Escape[gc_lisa$t50>1.6] <- "Not Silenced"

gc_lisa$Sil_Dyn <- NA
gc_lisa$Sil_Dyn[gc_lisa$t50<0.5] <- "Early"
gc_lisa$Sil_Dyn[(gc_lisa$t50>0.9)&(gc_lisa$t50<1.3)] <- "Late"

gcm_lisa <- gather(gc_lisa, model, value, -gene, -t50, -compare, -XT50)
gcm_lisa <- gcm_lisa[!is.na(gcm_lisa$value),] # remove un-classified genes

# plot
gcm_lisa$value <- factor(gcm_lisa$value, levels = c("Early", "Late", "Silenced", "Not Silenced"))
gcm_lisa$model <- revalue(factor(gcm_lisa$model), replace = c("XCI_Escape" = "XCI/Escape\nmodel", "Sil_Dyn" = "Silencing Dynamics\nmodel"))
gcm_lisa$model <- factor(gcm_lisa$model, levels = c("Silencing Dynamics\nmodel", "XCI/Escape\nmodel"))
ngenes_lisa <- gcm_lisa %>% dplyr::group_by(model, value) %>% dplyr::summarise(n = length(unique(gene)))
cols <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")

g <- gcm_lisa %>%  
  ggplot(aes(x = model, y = XT50, color = value, fill = value)) +
  theme_bw() + theme1 +  
  geom_jitter(color = "black",
              alpha = 1/4, position=position_jitterdodge(jitter.width = .1, dodge.width = 0.75), 
              size = small_scattersize, show.legend = FALSE) +
  stat_summary(fun=median, aes(ymin=..y.., ymax=..y..), geom='errorbar', width=0.5, color = "black",
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  geom_text(data = ngenes_lisa, aes(x = model, y = Inf, group = value, label = paste0("n=", n)), 
            position = position_dodge(width = 0.75), angle = 90, alpha = 1, size = geomtext_size,
            color = "black", hjust=-0.5) + 
  coord_cartesian(clip = "off") +
  scale_color_manual(values = cols) +
  scale_y_continuous(breaks = seq(0, 100, 50), limits = c(0,100)) +
  labs(x = "", 
       y = expression(XP[50]*"[%]"),
       color = "")
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, savefile = paste0(outpath, "E_Comparison_Sousa.pdf"))

print("7) F: Distance from Xist across K-means silencing classes")

print("7.1) Compute average distance from Xist")
m <- match(class_melt$Gene, bm$mgi_symbol); table(is.na(m))
class_melt <- data.frame(class_melt, bm[m, c("start_position", "end_position")])
class_melt$Xist_start <- bm$start_position[bm$mgi_symbol == "Xist"]
class_melt$Xist_end <- bm$end_position[bm$mgi_symbol == "Xist"]
class_melt$Xist_avgpos <- (class_melt$Xist_start + class_melt$Xist_end)/2
class_melt$Gene_avgpos <- (class_melt$start_position + class_melt$end_position)/2
class_melt$distance_Xist_Mb <- abs(class_melt$Xist_avgpos - class_melt$Gene_avgpos)/1e6
class_melt$Xist <- revalue(class_melt$Xist, replace = c("BL6_MA" = "Xist-MA (Xi=B6)",
                                                        "Cast_MA" = "Xist-MA (Xi=Cast)"))

print("7.2) Test difference between groups for each Xist-MA class")
# B6-ma
temp <- class_melt[class_melt$Xist == "Xist-MA (Xi=B6)",]
aov <- aov(distance_Xist_Mb~factor(k_ordered), data = temp)
summary(aov)
tt_b6ma <- pairwise.t.test(temp$distance_Xist_Mb, factor(temp$k_ordered), p.adjust.method = "none"); tt_b6ma

# Cast-ma
temp <- class_melt[class_melt$Xist == "Xist-MA (Xi=Cast)",]
aov <- aov(distance_Xist_Mb~factor(k_ordered), data = temp)
summary(aov)
tt_castma <- pairwise.t.test(temp$distance_Xist_Mb, factor(temp$k_ordered), p.adjust.method = "none"); tt_castma

print("7.3) Plot")
class_melt$k_ordered_classes <- plyr::revalue(factor(class_melt$k_ordered), replace = c("1" = "Early",
                                                                                        "2" = "Intermediate",
                                                                                        "3" = "Late",
                                                                                        "4" = "Escapee"))
ngenes <- class_melt %>% 
  dplyr::group_by(Xist, k_ordered_classes) %>% 
  dplyr::summarise(n = length(unique(Gene))) %>%
  as.data.frame()
g <- class_melt %>%  
  ggplot(aes(x = k_ordered_classes, y = distance_Xist_Mb, color = k_ordered_classes, fill = k_ordered_classes)) +
  theme_bw() + theme1 + 
  facet_grid(.~Xist) +
  # geom_boxplot(aes(x = k_ordered_classes, y = distance_Xist_Mb), outlier.shape = NA, size = violin_box_size) +
  geom_jitter(aes(x = k_ordered_classes, y = distance_Xist_Mb), fill = "black", color = "black",
              alpha = 1/4, position=position_jitter(width = .1), 
              size = small_scattersize, show.legend = FALSE) +
  stat_summary(fun=median, aes(x = k_ordered_classes, ymin=..y.., ymax=..y..), geom='errorbar', width=0.5, color = "black",
               position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size = linesize*2) +
  geom_text(data = ngenes, aes(x = k_ordered_classes, y = Inf, label = paste0("n = ", n)), 
            angle = 90, alpha = 1, size = geomtext_size, color = "black", hjust=-0.5) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  labs(x = expression(XP[50]*": K-means classification [K=4]"), y = "Average distance from Xist [Mb]") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(color = FALSE)
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 3, savefile = paste0(outpath, "F_XP50_XistDistance.pdf"), 
            height = 3)


print("8) H: Significantly differentially silenced genes")

print("8.1) Load data")
gene <- as.character(binresults$model_fit$Gene[binresults$model_fit$FDR<sigthreshold])
gene <- gene[!is.na(gene)]
cell_values <- binresults$data_cell[binresults$data_cell$Gene %in% gene,]
bin_values <- binresults$data[binresults$data$Gene %in% gene,]
revalue_xist <- c("Undetected" = "Undetected", "Low-Xist" = "Low", "Middle" = "Skewed",
                  "Xist_BA" = "BA", "BL6_MA" = "Xist-MA (Xi=B6)", "Cast_MA" = "Xist-MA (Xi=Cast)")
cell_values$XistGroup <- revalue(cell_values$Xist, revalue_xist)
bin_values$XistGroup <- revalue(bin_values$Xist, revalue_xist)

# extract slopes
bin <- binresults$model_fit
signif <- bin[bin$Gene %in% gene,]
slopes <- data.frame(Gene = unique(signif$Gene),
                     slope = c(signif$slope_b6, signif$slope_cast),
                     intercept = c(signif$int_b6, signif$int_cast),
                     XistGroup = rep(c("Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)"), each = nrow(signif)))
slopes$XistGroup <- factor(slopes$XistGroup, levels = c("Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)"))

# compute XT50s
threshold <- 100
bin$T50_b6 <- -(1 + bin$int_b6)/bin$slope_b6; bin$T50_cast <- -(1 + bin$int_cast)/bin$slope_cast
bin$T50_b6_original <- bin$T50_b6; bin$T50_cast_original <- bin$T50_cast; 
bin$T50_b6[(!is.na(bin$T50_b6)) & (bin$T50_b6<0) | bin$T50_b6 > threshold] <- threshold
bin$T50_cast[(!is.na(bin$T50_cast)) & (bin$T50_cast<0) | bin$T50_cast > threshold] <- threshold
s <- bin[bin$Gene %in% gene,]
t50 <- reshape2::melt(s, id.vars = colnames(signif)[!colnames(signif) %in% c("T50_b6", "T50_cast")], 
                      measure.vars = c("T50_b6", "T50_cast"), value.name = "t50")
t50$XistGroup <- revalue(t50$variable, c("T50_b6" = "Xist-MA (Xi=B6)", "T50_cast" = "Xist-MA (Xi=Cast)"))

# add labels
signif$label <- paste0(signif$Gene, "\nFDR = ", format(signif$FDR, digits = 1))
lev <- signif$label[rank(signif$FDR)]
bin_values$genelab <- factor(signif$label[match(bin_values$Gene, signif$Gene)], levels = lev)
t50$genelab <- factor(signif$label[match(t50$Gene, signif$Gene)], levels = lev)
slopes$genelab <- factor(signif$label[match(slopes$Gene, signif$Gene)], levels = lev)

print("8.2) Plot")
minvalue <- min(log2(bin_values$XiXa_ratio_XistUnd))
g <- bin_values %>%
  ggplot() + 
  theme_bw() + theme1 + facet_grid(genelab~.) +
  geom_hline(yintercept = 0, size = linesize, alpha = 1/2, linetype = "dashed") +
  geom_point(aes(x = AXCR, y = log2(XiXa_ratio_XistUnd), color = XistGroup), size = scattersize*binscatter_mult) + 
  geom_segment(data = t50, aes(x = 0, y = log2(1/2), xend = t50, yend = log2(1/2), color = XistGroup), linetype = "dashed",
               show.legend = FALSE, size = linesize) +
  geom_segment(data = t50, aes(x = t50, y = minvalue, xend = t50, yend = log2(1/2), color = XistGroup), linetype = "dashed",
               show.legend = FALSE, size = linesize) +
  geom_abline(data = slopes, aes(intercept = intercept, slope = slope, color = XistGroup),
              show.legend = FALSE, size = linesize) + 
  scale_color_manual(values = color_alleles) + 
  scale_alpha_manual(values = c(1/6, 1/4)) +
  labs(x = "XCI progress [%]",
       y = "Xi/Xa [log2]",
       color = "") + 
  scale_x_continuous(breaks = seq(0, 100, 50), limits = c(0, 100)) +
  theme(strip.text.y = element_text(angle = 0))
adjust_size(g = g, panel_width_cm = 1.5, panel_height_cm = 2, savefile = paste0(outpath, "H_MostDEgenes.pdf"), height = 10)


print("9) I: Concordance in allelic silencing dynamics gene-wise classification")

print("9.1) Load data - K-means trend")
bin_value <- binresults$data[binresults$data$Gene %in% class_melt$Gene,]

# include kmeans group
bin_value$k <- NA
bin_value$Xist <- revalue(factor(bin_value$Xist), replace = revalue_xist)
m <- match(paste0(bin_value$Xist, "_", bin_value$Gene), paste0(class_melt$Xist, "_", class_melt$Gene))
bin_value$k <- class_melt$k_ordered[m]
bin_value <- bin_value[!is.na(bin_value$k),]

# define coloring group
color_trend <- c("#b2df8a", "black", "#1f78b4")
genes <- c("Klhl13", "Pir", "Hprt")

# define name for bins
bin_value$bin_id <- as.numeric(gsub(x = strsplit2(bin_value$bin_XCR, split = ",")[,1], pattern = "[(]", replacement = ""))

# remove UndTime0h point
bin_value$kgroup <- paste0("Cluster ", bin_value$k)
bin_value$XistGroup <- revalue(bin_value$Xist, c("BL6_MA" = "Xist-MA (Xi=B6)", "Cast_MA" = "Xist-MA (Xi=Cast)"))
ngenes <- ddply(bin_value, .variables = .(XistGroup, kgroup), summarize, n = length(unique(Gene)))
ds_genes <- bin_value[bin_value$Gene %in% genes,]

print("9.2) Plot - K-means trend")
bin_value$XistGroup <- gsub(bin_value$XistGroup, pattern = ' ', replacement = '\n')
ds_genes$XistGroup <- gsub(ds_genes$XistGroup, pattern = ' ', replacement = '\n')
g <- bin_value %>%
  ggplot() + 
  theme_bw() + theme1 +
  facet_grid(kgroup~XistGroup) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 1/4, size = linesize) + 
  geom_line(aes(x = bin_id, y = log2(XiXa_ratio_XistUnd), group = Gene), alpha = 0.2, size = linesize) + 
  scale_x_continuous(breaks = seq(0, 100, 50), limits = c(0,100)) +
  labs(x = "XCI progress [%]",
       y = "Xi/Xa [log2]") +
  theme(legend.position = "right")
adjust_size(g = g, panel_width_cm = 1.5, panel_height_cm = 1.5, savefile = paste0(outpath, "I1_kmeansTrend.pdf"), height = 10)

print("9.3) Load data - Sankey")
silclasses <- c("fast", "medium", "slow", "escape")
nodes <- data.frame(node = 0:7, 
                    name = paste0(rep(c("b6", "cast"), each = 4), "_", rep(1:4, times = 2)),
                    id = paste0(rep(silclasses, times = 2), " ", rep(c("(Xi = B6)", "(Xi = Cast)"), each = 4)))
temp <- class_melt[!is.na(class_melt$FDR),] %>% 
  dplyr::group_by(Gene) %>%
  dplyr::summarise(b6 = paste0("b6_", k_ordered[Xist=="Xist-MA (Xi=B6)"]), 
                   FDR = unique(FDR),
                   cast = paste0("cast_", k_ordered[Xist=="Xist-MA (Xi=Cast)"]),
                   diff = abs(k_ordered[Xist=="Xist-MA (Xi=B6)"] - k_ordered[Xist=="Xist-MA (Xi=Cast)"]))
temp$significant <- factor(ifelse(temp$FDR<sigthreshold, "significant", "not significant"), levels = c("significant", "not significant"))
temp <- temp %>% 
  dplyr::group_by(Gene, b6, cast, significant) %>%
  dplyr::summarise(n = length(b6), 
                   genes = paste0(unique(Gene), collapse = "."),
                   diff = unique(diff))
temp$b6_new <- factor(as.character(nodes$id[match(temp$b6, nodes$name)]), levels = paste0(silclasses, " (Xi = B6)"))
temp$cast_new <- factor(as.character(nodes$id[match(temp$cast, nodes$name)]), levels = paste0(silclasses, " (Xi = Cast)"))

print("9.4) Plot - Sankey")
temp$kmeans_diff <- ifelse(temp$diff == 0, "same cluster\ndiff = 0",
                           ifelse(temp$diff == 1, "neighboring cluster\ndiff = 1", 
                                  "distal cluster\ndiff > 1"))
temp[temp$diff>1,]
g <- temp %>%
  ggplot(aes(axis1 = b6_new, axis2 = cast_new, y = n)) + 
  theme_bw() + theme1 + 
  geom_alluvium(aes(fill = kmeans_diff), size = linesize, colour = "white") +
  scale_fill_manual(values = c("#ca0020", "#f4a582", "grey")) +
  geom_stratum(size = linesize, width = 1/10) +
  geom_text(stat = "stratum", infer.label = TRUE, size = geomtext_size) +
  scale_x_discrete(limits = c("b6_new", "cast_new"), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="", y = "", fill = paste0("FDR<", sigthreshold))
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 5, savefile = paste0(outpath, "I2_SankeyPlot.pdf"), 
            height = 10, width = 10)
