print("1) Define output path")

outpath <- paste0(path, "output/fig4_deXchrChangeHighLow/"); dir.create(path = outpath, recursive = T, showWarnings = F)
f1path <- paste0(path, "output/fig1_NotAS/")

print("2) A: X-chromosome change K-means classification")

print("2.1) Load data")

load(paste0(datapath, "NCF_DGE_Spliced.RData")); spliced <- dge
load(paste0(datapath, "DGE.RData")); notas <- dge
load(paste0(f1path, "notAS_vel.RData"))

print("2.2) Compute RNA-velocity predicted change in normalized X-chromosome expression")

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

print("2.3) Define data frame")

df <- data.frame(day = samples$day, 
                 id = rownames(samples),
                 Xist_AS = samples$Xist_ratio_class,
                 Xist_notAS = samples$Xist_class,
                 diff = diff_xchr)
lv <- c("Undetected" = "Undetected", "Low-Xist" = "Low",
        "Middle" = "Skewed", "Xist_BA" = "BA",
        "BL6_MA" = "Xist-MA (Xi=B6)", "Cast_MA" = "Xist-MA (Xi=Cast)")
df$xist_as <- factor(revalue(factor(df$Xist_AS), replace = lv), levels = lv)

print("2.4) Compute and order K-means clusters, separately for each time point")

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

print("2.5) Plot")

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
  labs(x="Time [days]", y = expression(log[2]*"( "*Xchr[Current]*" / "*Xchr[Predicted]*" )"))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 2, savefile = paste0(outpath, "A_XchrChange_Kmeans.pdf"))


print("3) B: X-chromosome change High vs Xist Low MAST DE analysis per time point")

print("3.1) Define contrast and launch DE analyses")

# define contrast(s)
contrasts <- list(HighXchrChange_LowXchrChange = c("km_XchrChange", list("high", "low")))
min_sample_test <- 10; de_threshold <- 5e-2; days <- 1:4

# launch MAST analysis for each contrast
for(cont in 1:length(contrasts)){
  comparison <- names(contrasts)[cont]
  variable <- unlist(contrasts[[cont]][1]); variable_levels <- contrasts[[cont]][2:3]
  
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
        zlmdata <- zlm(~cdr + grp, sca)
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

print("3.2) Load MAST results excluding Xist from DE gene list")
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
  scale_y_continuous(limits = c(0, 250), breaks = seq(0, 200, by = 50)) + 
  ylab("# Differentially expressed genes") + xlab("Time [days]") +
  guides(fill = guide_legend(title="Higher expression in:", nrow=2, byrow=TRUE), alpha = 'none')
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 3, savefile = paste0(outpath, "B_DE_barplot.pdf"))


print("4) C: Represent DE genes through heatmap")

print("4.1) Define plotting function")
plot_heatmap <- function(all_de, contrasts, comparison, minFDR = 0.05, minAbsFC = 2,
                         outpath, is_log10CPMo1 = TRUE, abslog2FC_break = 7,
                         colorPalette = c("black", "gold", "red"), paletteLength = 500,
                         cellheight = 15, fontsize_col = 15, fontsize_row = 15, width = 15, cellwidth = 50, fontsize = 15,
                         excludeFC = "|log2FC|<1", orderByXist = TRUE, orderbyFDR = FALSE){
  
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
  # for(addlabel in c("all", "autosomal", "Xlinked")){
  for(addlabel in c("all")){
    if(addlabel == "all"){sig <- sig_all}
    # if(addlabel == "autosomal"){sig <- sig_all[!sig_all$chromosome_name %in% "X",]}
    # if(addlabel == "Xlinked"){sig <- sig_all[sig_all$chromosome_name %in% "X",]}
    
    # plot results for each time point --> Per cell heatmap
    heat_path <- paste0(outpath, comparison, "/Heatmap/PerCell/")
    if(nrow(sig)>1){
      for(d in days){
        dir.create(heat_path, showWarnings = FALSE, recursive = TRUE)
        s <- sig[sig$day == d,]; s$log2FC <- s$coef; 
        
        if(orderbyFDR){
          # order results by FDR-sign
          s$o <- ifelse(s$log2FC>0, s$fdr, -log10(s$fdr))
          s <- s[order(s$o, decreasing = FALSE),]
          orderlabel <- "_byFDR"
        }else{
          # order by FC
          s <- s[order(s$log2FC, decreasing = TRUE),]
          orderlabel <- "_byFC"
        }
        
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
          avg_aut <- avgexp[!log2fc %in% excludeFC, ]
          
          # order by Xist expression
          if(orderByXist){
            if(comparison == "XistHigh_XistLow"){
              o <- order(counts["Xist",], decreasing = FALSE)
              avg_aut <- avg_aut[, o]
              m <- match(colnames(avg_aut), rownames(dge$samples))
              annot_col <- data.frame(Levels = dge$samples[[variable]][m]); rownames(annot_col) <- colnames(avg_aut)
            }
            if(comparison == "HighXchrChange_LowXchrChange"){
              annot_col$Xist <- counts["Xist", match(rownames(annot_col), colnames(counts))]
              o <- order(annot_col$Levels, annot_col$Xist)
              avg_aut <- avg_aut[, o]
              m <- match(colnames(avg_aut), rownames(dge$samples))
              annot_col <- data.frame(Levels = dge$samples[[variable]][m]); rownames(annot_col) <- colnames(avg_aut)
            }
          }
          
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

print("4.2) Produce heatmaps")
comp <- names(contrasts)
minFDR <- 0.01
minAbsFC <- 1.5
is_log10CPMo1 <- TRUE
abslog2FC_break <- 7
load(paste0(outpath, "XchrChange_HighLow_MAST.RData"))

plot_heatmap(all_de = all_de, contrasts = contrasts, comparison = comp, minFDR = minFDR, minAbsFC = minAbsFC, 
             is_log10CPMo1 = is_log10CPMo1, abslog2FC_break = abslog2FC_break, outpath = outpath,
             cellheight = 4.5, fontsize_col = 4.5, fontsize_row = 4.5, width = 10, cellwidth = .5, fontsize = 4.5)

print("4.3) Move and rename plots in figures folder")
f <- paste0(outpath, names(contrasts), "/Heatmap/PerCell/", 
            c("day1_Avelog10CPM_0.01_all_byFC.pdf", "day2_Avelog10CPM_0.01_all_byFC.pdf"))
file.copy(from = f, to = outpath)
f1 <- paste0(outpath, c("day1_Avelog10CPM_0.01_all_byFC.pdf", "day2_Avelog10CPM_0.01_all_byFC.pdf"))
f2 <- paste0(outpath, c("C_Heatmap_24h_byFC.pdf", "C_Heatmap_48h_byFC.pdf"))
file.rename(from = f1, to = f2)



print("5) D: Correlation analysis")

print("5.1) Load data and launch gene-wise correlation analysis to X-chromosome change")

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

print("5.2) Select genes significantly correlated to X-chromosome change for each time point")
temp <- sprmcor %>% as.data.frame()
temp$isX <- ifelse(temp$Chr %in% "X", "X-linked", "Autosomal")
temp$direction <- ifelse(temp$spr < 0, "Negative", "Positive")
temp$direction <- factor(temp$direction, levels = c("Positive", "Negative"))
sig_de <- temp[(temp$fdr <= de_threshold),]

print("5.3) Plot - Summarize correlation analysis results with barplot")
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
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 3, savefile = paste0(outpath, "D_DE_barplot_corrXchr.pdf"))


print("6) E/F: Heatmap - positive/negative regulators over time at day 1/2")

fdr_threshold <- 0.05
maxn <- 20

print(paste0("6.1) Subset to genes with FDR<=", fdr_threshold, " and |rho|>=", cor_thr, " at d1 or d2"))

sub <- sprmcor[(sprmcor$day %in% c(1, 2)),]
sub_sig <- sub[sub$fdr <= fdr_threshold,]
table(sub_sig$day, sub_sig$spr>0)
summary(abs(sub_sig$spr[sub_sig$day == 1]))
summary(abs(sub_sig$spr[sub_sig$day == 2]))

print("6.2) Identify positive and negative regulators")

positive <- sub_sig[sub_sig$spr > 0,] %>% dplyr::arrange(-abs(spr))
pos_genes <- unique(as.character(positive$Gene))
negative <- sub_sig[sub_sig$spr<0,] %>% dplyr::arrange(-abs(spr))
neg_genes <- unique(as.character(negative$Gene))

if(length(pos_genes)>maxn) pos_genes <- pos_genes[seq_len(maxn)]
if(length(neg_genes)>maxn) neg_genes <- neg_genes[seq_len(maxn)]


print("6.3) Plot results through heatmap")

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
            savefile = paste0(outpath, "E_corrXchrChange_PositiveRegulators.pdf"), 
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
            savefile = paste0(outpath, "F_corrXchrChange_NegativeRegulators.pdf"), 
            height = 3, width = 5)
