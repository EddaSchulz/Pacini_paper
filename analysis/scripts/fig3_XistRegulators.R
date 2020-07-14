library(MAST); library(pheatmap); library(multidplyr)

print("1) Define output path")
outpath <- paste0(path, "output/fig3_deXistHighLow/"); dir.create(path = outpath, recursive = T, showWarnings = F)


print("2) A: Xist CPM K-means classification")

print("2.1) Define Xist Low/Medium/High cells based on K-means CPM classification")

load(paste0(datapath, "DGE.RData"))
df <- dge$samples
df$xist_cpm <- (dge$counts["Xist",]/(colSums(dge$counts)*dge$samples$sf_notX))*1e6
df$xist_class <- "Low"; df$xist_class[df$Xist_kmeans_class == 4] <- "Medium"; df$xist_class[df$Xist_kmeans_class %in% 5:7] <- "High"
df$xist_class <- factor(df$xist_class, levels = c("Low", "Medium", "High"))
df$Xist_kmeans_class <- factor(df$Xist_kmeans_class)
summ_cell <- ddply(df, .variables = .(day, xist_class), summarize, n = length(day))
table(df$day, df$xist_class)

print("2.2) Plot Xist CPM classification over time")

thresholds <- c(min(df$xist_cpm[df$xist_class %in% "Medium"]), 
                max(df$xist_cpm[df$xist_class %in% "Medium"]))
subthresholds <- ddply(df, .variables = .(Xist_kmeans_class), summarize, m = max(xist_cpm))
cols <- rev(c(rev(colorRampPalette(c("orange", "red"))(3)), 
              "grey", 
              rev(colorRampPalette(c("lightblue", "blue"))(3))))

g <- df %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = log10(subthresholds$m+1), linetype = "dashed", alpha = 1/3, size = linesize, color = "black") + 
  geom_hline(yintercept = log10(thresholds+1), linetype = "dashed", alpha = hvalpha, size = linesize, color = "red") + 
  geom_jitter(aes(x = day, y = log10(xist_cpm+1), group = xist_class, color = Xist_kmeans_class), 
              size = outliersize, show.legend = FALSE,
              position=position_jitterdodge(jitter.width = .75, dodge.width = 0.75)) +
  scale_color_manual(values = cols) +
  scale_y_continuous(limits = c(0, 4), breaks = 0:4) +
  geom_text(summ_cell, mapping = aes(x = day, y = Inf, label = paste0("n = ", n), fill = xist_class),
            position=position_dodge(0.75), angle = 90, alpha = 1, size = geomtext_size, hjust = -0.5) + 
  coord_cartesian(clip = "off") +
  labs(x="Time [days]", y = expression("Xist CPM + 1 ["* log[10]*"]"))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "A_XistCPM_classification.pdf"))


print("3) B: Xist High vs Xist Low MAST DE analysis per time point")

print("3.1) Define contrast and launch DE analyses")

# define contrast(s)
contrasts <- list(XistHigh_XistLow = c("Xist_kmeans_class", list(5:7, 1:3)))
depath <- outpath
min_sample_test <- 10; de_threshold <- 5e-2; days <- 0:4

# launch MAST analysis for each contrast
for(cont in 1:length(contrasts)){
  comparison <- names(contrasts)[cont]
  variable <- unlist(contrasts[[cont]][1]); variable_levels <- contrasts[[cont]][2:3]
  
  if(!dir.exists(paste0(depath, comparison))){
    dir.create(path = paste0(depath, comparison), showWarnings = FALSE, recursive = TRUE)
    
    for (i in 1:length(days)) {
      time <- days[i]; cat('Processing day', time, '..')
      load(paste0(datapath, "DGE.RData"))
      dge <- dge[, dge$samples$day == time]
      dge$samples$sample <- rownames(dge$samples)
      
      # subset to grouping conditions only
      temp <- dge[!(grepl(dge$genes$symbol, pattern = "^ERCC") | duplicated(dge$genes$ensembl)), dge$samples[,variable] %in% unlist(variable_levels)]
      grp <- temp$samples[,variable]
      grp <- ifelse(grp %in% variable_levels[[1]], "case", "control"); table(grp)
      names(grp) <- temp$samples$sample
      grp <- factor(grp, levels = c("control", "case"))
      
      # store the number of cases and controls per time point
      if(!file.exists(paste0(depath, comparison, "/Nsamples_perTime.txt"))){
        m <- matrix(c(time, sum(grp == "case"), sum(grp == "control")), nrow = 1)
        colnames(m) <- c("Day", c(strsplit2(comparison, split = "\\_")))
        write.table(m, quote =FALSE, row.names = FALSE, col.names = TRUE, 
                    file = paste0(depath, comparison, "/Nsamples_perTime.txt"))
      }else{
        m <- matrix(c(time, sum(grp == "case"), sum(grp == "control")), nrow = 1)
        write.table(m, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE,
                    file = paste0(depath, comparison, "/Nsamples_perTime.txt"))
      }
      
      if(min(table(grp)) >= min_sample_test){
        # compute detection rate and cpm matrix
        cdr <- scale(colMeans(temp$counts > 0))
        els <- colSums(temp$counts)*(temp$samples$sf_notX/mean(temp$samples$sf_notX))
        cpm <- t(t(temp$counts)/els)*1e6; rownames(cpm) <- temp$genes$ensembl
        sca <- FromMatrix(exprsArray = log2(cpm + 1), 
                          cData = data.frame(wellKey = names(grp), grp = grp, cdr = cdr),
                          fData = data.frame(chromosome = temp$genes$chromosome, symbol = temp$genes$symbol, ensembl = temp$genes$ensembl, primerid = temp$genes$ensembl))
        zlmdata <- zlm(~cdr + grp, sca)
        summaryCond <- summary(zlmdata, doLRT='grpcase') ; summaryDt <- summaryCond$datatable
        mast <- merge(summaryDt[contrast=='grpcase' & component=='H',.(primerid, `Pr(>Chisq)`)],
                      summaryDt[contrast=='grpcase' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
        mast$AveCPM <- rowMeans(cpm)
        mast$fdr <- p.adjust(mast$`Pr(>Chisq)`, method = "BH")
        mast <- data.frame(bm[match(mast$primerid, bm$ensembl_gene_id), c("chromosome_name", "mgi_symbol", "entrezgene_id", "ensembl_gene_id")], mast)
        mast <- mast[order(mast$fdr, decreasing = FALSE),]
        
        # store results
        save(mast, file = paste0(depath, comparison, "/", time, ".RData"))
        save(zlmdata, file = paste0(depath, comparison, "/zlm_", time, ".RData"))
        save(sca, file = paste0(depath, comparison, "/sca_", time, ".RData"))
      }
    }
  }
}

# store results
all_de <- c()
for (i in 1:length(days)) {
  time <- days[i]; cat('Processing day', time, '..')
  
  for(cont in 1:length(contrasts[-length(contrasts)])){
    comparison <- names(contrasts)[cont]
    
    # de results
    file <- paste0(depath, comparison, "/", time, ".RData")
    if(file.exists(file)){
      load(file); all <- mast; all_de <- rbind(all_de, data.frame(test = comparison, day = time, all))
    }
  }
}
save(all_de, file = paste0(depath, "ALLresults.RData"))

print("3.2) Load MAST results excluding Xist from DE gene list")
all_de <- all_de[!all_de$mgi_symbol %in% "Xist",]
all_de$isX <- ifelse(all_de$chromosome_name %in% "X", "X-linked", "Autosomal")
all_de$direction <- ifelse(all_de$coef < 0, "Down-regulated", "Up-regulated")
all_de$direction <- factor(all_de$direction, levels = c("Up-regulated", "Down-regulated"))
sig_de <- all_de[(all_de$fdr <= de_threshold),]
de_results <- sig_de[sig_de$test %in% "XistHigh_XistLow",]
de_barplot <- ddply(sig_de[sig_de$test %in% "XistHigh_XistLow",], 
                    .variables = .(day, isX, direction), summarize, 
                    n = length(day))

all <- expand.grid(unique(de_barplot$day), unique(de_barplot$isX), unique(de_barplot$direction)) 
all$id <- apply(all, 1, function(x) paste(x, collapse = "_"))
present <- apply(de_barplot[,1:3], 1, function(x) paste(x, collapse = "_"))
missing <- cbind(all[!all$id %in% present, !colnames(all) %in% "id"], 0)
colnames(missing) <- colnames(de_barplot)
de_barplot_ext <- rbind(de_barplot, missing)
de_barplot_ext$direction <- factor(de_barplot_ext$direction, levels = c("Up-regulated", "Down-regulated"))
de_barplot_ext$higherexp <- revalue(de_barplot_ext$direction, replace = c("Up-regulated" = "Xist high cells",
                                                                          "Down-regulated" = "Xist low cells"))
de_barplot_ext$alpha_numb <- ifelse(de_barplot_ext$n>0, "yes", "no")

print("3.3) Plot - Summarize MAST DE results with barplot")
g <- de_barplot_ext %>% 
  ggplot() + 
  theme_bw() + theme1 + 
  facet_grid(.~isX) + 
  geom_bar(aes(x = factor(day), y = n, fill = higherexp), stat = "identity", position=position_dodge()) + 
  scale_fill_manual(values = c("red", "blue")) +
  geom_text(aes(x = factor(day), y = n+1, group = higherexp, label = n), 
            hjust = -.1,
            color="black", position = position_dodge(0.9), size = geomtext_size, angle = 90, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 185), breaks = seq(0, 150, by = 50)) + 
  ylab("# Differentially expressed genes") + xlab("Time [days]") +
  guides(fill = guide_legend(title="Higher expression in:", nrow=2, byrow=TRUE), alpha = 'none')
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 3, savefile = paste0(outpath, "B_DE_barplot.pdf"))


print("4) C: Represent DE genes through heatmap")

print("4.1) Define plotting function")
plot_heatmap <- function(all_de, contrasts, comparison, minFDR = 0.05, minAbsFC = 2,
                         depath, is_log10CPMo1 = TRUE, abslog2FC_break = 7,
                         colorPalette = c("black", "gold", "red"), paletteLength = 500,
                         cellheight = 15, fontsize_col = 15, fontsize_row = 15, width = 15, cellwidth = 50, fontsize = 15,
                         excludeFC = "|log2FC|<1", orderByXist = TRUE){
  
  # define heatmap coloring
  myColor <- c(colorRampPalette(colorPalette)(paletteLength))
  
  # select levels of two groups being compared
  variable <- unlist(contrasts[[comparison]][1]); variable_levels <- contrasts[[comparison]][2:3]
  
  # subset results
  all_de$log2FC <- all_de$coef
  sig_all <- all_de[(all_de$test == comparison)&
                      (abs(all_de$log2FC) >= log2(minAbsFC))&
                      (all_de$fdr <= minFDR),]
  
  # repeat plot for Xlinked only and Autosomal+PosXlinked
  for(addlabel in c("autosomal_posXlinked", "Xlinked")){
    if(addlabel == "autosomal_posXlinked"){sig <- sig_all[!(sig_all$chromosome_name %in% "X" & sig_all$coef<0),]}
    if(addlabel == "Xlinked"){sig <- sig_all[sig_all$chromosome_name %in% "X",]}
    
    # order results by FC
    sig <- sig[order(sig$log2FC, decreasing = TRUE),]
    
    # plot results for each time point --> Per cell heatmap
    heat_path <- paste0(depath, comparison, "/Heatmap/PerCell/")
    if(nrow(sig)>1){
      for(d in days){
        dir.create(heat_path, showWarnings = FALSE, recursive = TRUE)
        s <- sig[sig$day == d,]; s$log2FC <- s$coef; s <- s[order(s$log2FC, decreasing = TRUE),]
        if(nrow(s)>1){
          # store DE genes
          degenes_s <- s$ensembl_gene_id
          
          # load normalized matrices
          load(paste0(datapath, "DGE.RData"))
          dge <- dge[, dge$samples$day == d]
          
          if(comparison == "XistHigh_XistLow"){
            levels <- 1:7
            dge$samples[[variable]] <- as.numeric(dge$samples[[variable]])
          }else{
            levels <- c(variable_levels[[2]], variable_levels[[1]])
          }
          dge <- dge[, dge$samples[[variable]] %in% levels]
          
          if(is_log10CPMo1){
            counts <- log10(t(t(dge$counts)/(dge$samples$sf_notX*colSums(dge$counts)))*1e6 + 1)
          }else{
            counts <- log10(t(t(dge$counts)/dge$samples$sf_notX) + 1)
          }
          
          # subset matrix to cells in contrast and DE genes only --> compute average gene level per variable group
          avgexp <- counts[match(degenes_s, dge$genes$ensembl), ]
          
          # define column and row annotation
          annot_col <- data.frame(Levels = dge$samples[[variable]]); rownames(annot_col) <- colnames(avgexp)
          if(comparison == "XistHigh_XistLow"){
            colors_col <- c(colorRampPalette(c("lightblue", "blue"))(3), "grey",
                            colorRampPalette(c("orange", "red"))(3))
            names(colors_col) <- 1:7
          }else{
            colors_col <- c(colorRampPalette(c("orange", "red"))(length(variable_levels[[1]])),
                            colorRampPalette(c("blue", "lightblue"))(length(variable_levels[[2]])))
            names(colors_col) <- unlist(variable_levels)
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
          avg_aut <- avgexp[!log2fc %in% excludeFC, ]
          
          # order by Xist expression
          if(orderByXist){
            o <- order(counts["Xist",], decreasing = FALSE)
            avg_aut <- avg_aut[, o]
            m <- match(colnames(avg_aut), rownames(dge$samples))
            annot_col <- data.frame(Levels = dge$samples[[variable]][m]); rownames(annot_col) <- colnames(avg_aut)
          }
          
          # gap between up and down regulated
          gprw <- which(!duplicated(s$log2FC>0))[2] - 1
          if(is.na(gprw)) gprw <- NULL
          
          pheatmap(avg_aut, 
                   gaps_row = c(1, gprw),
                   annotation_colors = color_annot,
                   cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, color = myColor, 
                   cellwidth = cellwidth, cellheight = cellheight, 
                   fontsize_col = fontsize_col, fontsize_row = fontsize_row, fontsize = fontsize,
                   width = width, height = height,
                   gaps_col = which(!duplicated(annot_col$Levels))[-1] - 1,
                   filename = paste0(heat_path, "day", d, "_", title_label, "_", minFDR, "_", addlabel, ".pdf"))
        }
      }
    }
  }
}

print("4.2) Produce heatmaps")
comp <- names(contrasts)
minFDR <- de_threshold
minAbsFC <- 1.5
is_log10CPMo1 <- TRUE
abslog2FC_break <- 7
load(paste0(depath, "ALLresults.RData"))

plot_heatmap(all_de = all_de, contrasts = contrasts, comparison = comp, minFDR = minFDR, minAbsFC = minAbsFC, 
             is_log10CPMo1 = is_log10CPMo1, abslog2FC_break = abslog2FC_break, depath = outpath,
             cellheight = 6, fontsize_col = 6, fontsize_row = 6, width = 10, cellwidth = .5, fontsize = 6, orderByXist = TRUE)

print("4.3) Move and rename plots in figures folder")
f <- paste0(outpath, names(contrasts), "/Heatmap/PerCell/", c("day1_Avelog10CPM_0.05_autosomal_posXlinked.pdf", "day2_Avelog10CPM_0.05_autosomal_posXlinked.pdf"))
file.copy(from = f, to = outpath)
f1 <- paste0(outpath, c("day1_Avelog10CPM_0.05_autosomal_posXlinked.pdf", "day2_Avelog10CPM_0.05_autosomal_posXlinked.pdf"))
f2 <- paste0(outpath, c("C_Heatmap_24h_autosomal_posXlinked.pdf", "C_Heatmap_48h_autosomal_posXlinked.pdf"))
file.rename(from = f1, to = f2)


print("5) D: Correlation analysis")

print("5.1) Load data and launch gene-wise correlation analysis to Xist CPM expression")
load(paste0(datapath, "DGE.RData"))
dge$cpm <- t(t(dge$counts)/(colSums(dge$counts)*dge$samples$sf_notX))*1e6
x <- data.frame(day = rep(dge$samples$day, each = nrow(dge)), 
                Xist = rep(dge$cpm["Xist",], each = nrow(dge)),
                Gene = rep(dge$genes$symbol, times = ncol(dge)),
                Chr = rep(dge$genes$chromosome, times = ncol(dge)),
                count = c(dge$counts),
                cpm = c(dge$cpm)
)

ncores <- min(c(detectCores(), 10))
cluster <- new_cluster(ncores)
sprmcor <- x %>% 
  dplyr::group_by(day, Gene) %>% 
  partition(cluster) %>%
  dplyr::summarise(Chr = unique(Chr),
                   spr_Xist = cor(cpm, Xist, method = "spearman"),
                   droprate = mean(cpm == 0),
                   droprate_Xist = mean(Xist == 0),
                   n = length(cpm),
                   pvalue_Xist = cor.test(cpm, Xist, method = "spearman", use = "complete.obs")$p.value) %>% 
  dplyr::arrange(Gene, day) %>%
  collect()
sprmcor$fdr_Xist <- p.adjust(sprmcor$pvalue_Xist, method = "BH")
save(sprmcor, file = paste0(outpath, "sprmcor.RData"))

print("5.2) Select genes significantly correlated to Xist expression for each time point")
temp <- sprmcor[!is.na(sprmcor$fdr_Xist),] %>% as.data.frame()
temp <- temp[!temp$Gene %in% "Xist",] #exclude Xist
temp$isX <- ifelse(temp$Chr %in% "X", "X-linked", "Autosomal")
temp$direction <- ifelse(temp$spr_Xist < 0, "Negative Xist regulator", "Positive Xist regulator")
temp$direction <- factor(temp$direction, levels = c("Positive Xist regulator", "Negative Xist regulator"))
sig_de <- temp[(temp$fdr_Xist <= de_threshold),]

print("5.3) Plot - Summarize correlation analysis results with barplot")
de_barplot <- ddply(sig_de, .variables = .(day, isX, direction), summarize, n = length(day))
all <- expand.grid(unique(de_barplot$day), unique(de_barplot$isX), unique(de_barplot$direction)) 
all$id <- apply(all, 1, function(x) paste(x, collapse = "_"))
present <- apply(de_barplot[,1:3], 1, function(x) paste(x, collapse = "_"))
missing <- cbind(all[!all$id %in% present, !colnames(all) %in% "id"], 0)
colnames(missing) <- colnames(de_barplot)
de_barplot_ext <- rbind(de_barplot, missing)
de_barplot_ext$direction <- factor(de_barplot_ext$direction, levels = c("Positive Xist regulator", "Negative Xist regulator"))
de_barplot_ext$alpha_numb <- ifelse(de_barplot_ext$n>0, "yes", "no")

g <- de_barplot_ext[de_barplot_ext$day != 0,] %>% 
  ggplot() + 
  theme_bw() + theme1 + 
  facet_grid(.~isX) + 
  geom_bar(aes(x = factor(day), y = n, fill = direction), stat = "identity", position=position_dodge()) + 
  scale_fill_manual(values = c("red", "blue")) +
  geom_text(aes(x = factor(day), y = n+1, group = direction, label = n), 
            hjust = -.1,
            color="black", position = position_dodge(0.9), size = geomtext_size, angle = 90, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 280), breaks = seq(0, 250, by = 50)) +
  xlab("Time [days]") +
  ylab(paste0("Significant genes [FDR<=", de_threshold, "]")) + 
  guides(fill = guide_legend(title="", nrow=2, byrow=TRUE), alpha = 'none')
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 3, savefile = paste0(outpath, "D_DE_barplot_corrXist.pdf"))


print("6) E/F: Heatmap - positive/negative regulators over time at day 1/2")

print("6.1) Load data and remove Xist")
sub <- sprmcor[(sprmcor$day %in% c(1, 2)),]
sub <- sub[!sub$Gene %in% "Xist",]

print("6.2) Remove pseudogenes and select significant genes at day 1/2")
# remove pseudogenes
sub$gene_biotype <- bm$gene_biotype[match(sub$Gene, bm$mgi_symbol)]
sub <- sub[!sub$gene_biotype %in% "pseudogene",]

# select genes showing significant correlation to Xist expression at day 1/2
sub_sig <- sub[(sub$fdr_Xist<0.05)&(!is.na(sub$fdr_Xist)),]
summary(abs(sub_sig$spr_Xist[sub_sig$day == 1]))
summary(abs(sub_sig$spr_Xist[sub_sig$day == 2]))
summary(abs(sub_sig$spr_Xist))

print("6.3) Select genes with absolute correlation greater than 0.25")
sub_sig_high <- sub_sig[abs(sub_sig$spr_Xist)>0.25,]
table(sub_sig_high$spr_Xist>0)
positive <- sub_sig_high[sub_sig_high$spr_Xist>0,] %>% dplyr::arrange(-abs(spr_Xist))
pos_genes <- unique(as.character(positive$Gene))
negative <- sub_sig_high[(sub_sig_high$spr_Xist<0)&(!sub_sig_high$Chr %in% "X"),] %>% dplyr::arrange(-abs(spr_Xist)) # exclude X-linked genes for negative Xist regulators
neg_genes <- unique(as.character(negative$Gene))

print("6.4) E: Plot results through heatmap - positive Xist regulators")
ph <- length(pos_genes)*0.3
temp <- sprmcor[sprmcor$Gene %in% pos_genes,]; temp$Gene <- factor(temp$Gene, levels = rev(pos_genes))
temp$significant <- temp$fdr_Xist<0.05

g <- temp[temp$day != 0,] %>% 
  ggplot() + 
  theme_bw() + theme1 +
  geom_point(aes(x = factor(day), y = Gene, fill = spr_Xist, size = abs(spr_Xist)), pch=22, stroke = 0) +
  geom_point(data = temp[temp$significant == TRUE,], aes(x = factor(day), y = Gene), size = small_scattersize, color = "white", pch = 1) +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0, breaks=seq(-0.5, 0.5, .25), limits=c(-0.51,0.51)) + 
  scale_radius(breaks = seq(0, 0.5, by = 0.25), limits = c(0, 0.5), range = c(0.1, 4)) +
  labs(x = "Time [days]", y = "", fill = "Spearman's correlation to Xist CPM", size = "Absolute correlation")
adjust_size(g = g, panel_width_cm = 2.5, panel_height_cm = ph, savefile = paste0(outpath, "E_corrXist_PositiveRegulators.pdf"))

print("6.5) F: Plot results through heatmap - negative autosomal Xist regulators")
ph <- length(neg_genes)*0.3
temp <- sprmcor[sprmcor$Gene %in% neg_genes,]; temp$Gene <- factor(temp$Gene, levels = rev(neg_genes))
temp$significant <- temp$fdr_Xist<0.05

g <- temp[temp$day != 0,] %>% 
  ggplot() + 
  theme_bw() + theme1 +
  geom_point(aes(x = factor(day), y = Gene, fill = spr_Xist, size = abs(spr_Xist)), pch=22, stroke = 0) +
  geom_point(data = temp[temp$significant == TRUE,], aes(x = factor(day), y = Gene), size = small_scattersize, color = "white", pch = 1) +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0, breaks=seq(-0.5, 0.5, .25), limits=c(-0.51,0.51)) + 
  scale_radius(breaks = seq(0, 0.5, by = 0.25), limits = c(0, 0.5), range = c(0.1, 4)) +
  labs(x = "Time [days]", y = "", fill = "Spearman's correlation\nrho(Xist CPM, Gene CPM)", size = "Absolute correlation")
adjust_size(g = g, panel_width_cm = 2.5, panel_height_cm = ph, savefile = paste0(outpath, "F_corrXist_NegativeRegulators.pdf"))
