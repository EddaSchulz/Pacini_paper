print("1) Define output path")
outpath <- paste0(path, "output/fig1_NotAS/"); dir.create(path = outpath, showWarnings = F, recursive = T)


print("2) B: Pseudotime plot")

print("2.1) Compute pseudotimes using Monocle-DDRTree method")

load(paste0(datapath, "DGE.RData"))
dge$genes$gene_short_name <- dge$genes$symbol

# process and identify top 100 MVGs
pd <- new("AnnotatedDataFrame", data = dge$samples)
fd <- new("AnnotatedDataFrame", data = dge$genes)
XX <- newCellDataSet(as(dge$counts, "sparseMatrix"),
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit=1,
                     expressionFamily=negbinomial.size())
XX <- estimateSizeFactors(XX)
XX <- estimateDispersions(XX)

# identify DE genes over time --> select top 500 DE genes
pData(XX)$hour <- pData(XX)$day*24
de <- differentialGeneTest(XX, fullModelFormulaStr="~hour", cores=detectCores())
top_pdt_genes <- 500
de_genes <- de[de$qval < 1e-2,]
de_genes <- de_genes[order(de_genes$qval, decreasing = FALSE),]
order_genes <- rownames(de_genes[order(de_genes$qval, decreasing = FALSE),])[seq_len(top_pdt_genes)]
XX <- setOrderingFilter(XX, order_genes)

# reduce dimensions and compute pseudotime for XX cells
XX <- reduceDimension(XX, method = 'DDRTree')
XX <- orderCells(XX)

# the origin is defined as the state with the maximum number of 0h cells
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$hour)[,"0"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
XX <- orderCells(XX, root_state = GM_state(XX))

# compute scaled PDT
pData(XX)$Scaled_PDT <- pData(XX)$Pseudotime/max(pData(XX)$Pseudotime)*100

print("2.2) Plot embedding and color by time and PDT")
pData(XX)$day <- factor(pData(XX)$day)
g1 <- plot_cell_trajectory(XX, color_by = "day", 
                           cell_size = scattersize, cell_link_size = linesize, 
                           show_branch_points = FALSE) + theme1 +
  scale_color_manual(values = time_colors) + 
  ggtitle(paste0("XX Pseudotime\nDDRTree: ", top_pdt_genes, " MVGs (FDR<0.01) across time")) +
  labs(color = "Time [days]")
adjust_size(g = g1, panel_width_cm = 5, panel_height_cm = 5, savefile = paste0(outpath, "B_pdtXX_time.pdf"))

g2 <- plot_cell_trajectory(XX, color_by = "Scaled_PDT", 
                           cell_size = scattersize, cell_link_size = linesize, 
                           show_branch_points = FALSE) + theme1 +
  scale_color_gradient2(low = "black", mid = "black", high = "gold") + 
  guides(color = guide_colourbar(barwidth = 4, barheight = 0.7)) +
  ggtitle(paste0("XX Pseudotime\nDDRTree: ", top_pdt_genes, " MVGs (FDR<0.01) across time"))
adjust_size(g = g2, panel_width_cm = 5, panel_height_cm = 5, savefile = paste0(outpath, "B_pdtXX_pdt.pdf"))


print("3) C/D: UMAP embedding and RNA-velocity projections")

print("3.1) Compute RNA velocity")

# load data
load(paste0(datapath, "NCF_DGE_Spliced.RData")); spliced <- dge
load(paste0(datapath, "NCF_DGE_Unspliced.RData")); unspliced <- dge

# include marker genes
markergenes <- c("Xist", "Dnmt3a", "Nanog", "Esrrb")
x <- log10(t(t(spliced$counts)/spliced$samples$sf_notX)+1)
cellfeatures <- spliced$samples
cellfeatures <- data.frame(cellfeatures, t(x[markergenes,]))

# add features in cellfeatures
load(paste0(datapath, "DGE.RData"))
m <- match(cellfeatures$id, dge$samples$id)
cellfeatures <- data.frame(cellfeatures,
                           Time = paste0(dge$samples$day[m]*24, "h"),
                           Xchr_ratio = dge$samples$xchr_ratio[m],
                           Xist_ratio_class = dge$samples$Xist_ratio_class[m])

# filter out lowly expressed genes
rm_s <- rowMeans(spliced$counts); rm_u <- rowMeans(unspliced$counts)
summary(rm_s); summary(rm_u); plot(log10(rm_s), log10(rm_u))
keep <- (rm_s >= 1)&(rm_u >= 0.5); table(keep)
emat <- spliced_filt <- spliced$counts[keep,]; nmat <- unspliced_filt <- unspliced$counts[keep,]

# compute velocities
fit.quantile <- 0.025; kcells <- 20; mincor <- 0.05
arrow.scale=5; cell.alpha=0.4; cell.cex=1; fig.height=4; fig.width=4.5;
notAS_vel <- gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells = kcells, 
                                              fit.quantile = fit.quantile, 
                                              min.nmat.emat.slope = mincor,
                                              min.nmat.emat.correlation = mincor, 
                                              n.cores = 10)
save(notAS_vel, file = paste0(outpath, "notAS_vel.RData"))


print("3.2) Compute notAS UMAP embedding")

# settings
pcount <- 1; mvg <- 500; umap_nPcs <- 50; umap_nn <- 20; neighbor_size <- 100; mult <- 1e3
quantile <- 0.025; min.correlation <- min.slope <- 0.05; 
grid_size <- 30; arrow.pca <- 2.5; arrow <- 2.5; arrow.lwd <- 1; grid.mass <- 5
gapT <- 1; k_range <- 30; controlby <- "time"; cor_dist_measure <- "cor"; velocity_scale <- "sqrt"

# load data and define marker genes
load(paste0(datapath, "DGE.RData"))
x0 <- dge$counts
els <- cellfeatures$sf_notX*colSums(x0)
x0.log <- log10(t(t(x0)/els)*1e6+1)
markers <- c("Xist", "Nanog", "Dnmt3a", "Esrrb", "Tsix", "Pou5f1")
expr_markers <- x0.log[match(markers, rownames(x0.log)), ]

# select MVGs and compute UMAP embedding
rv <- rowVars(x0.log)
o <- order(rv, decreasing = TRUE)
feature_set <- head(o, mvg)
x0.log.mvg <- x0.log[feature_set,]
cent <- rowMeans(x0.log.mvg)
epc <- pcaMethods::pca(t(x0.log.mvg - cent), center = F, nPcs = umap_nPcs)
custom.config = umap.defaults; custom.config$n_neighbors <- umap_nn
umap_both <- umap(epc@scores, custom.config)
emb_umap <- umap_both$layout
features <- c("Time", "Xchr_ratio", "Xist_ratio_class")
sinfo <- cellfeatures[match(rownames(emb_umap), cellfeatures$id), features]
sinfo$Xchr_ratio <- as.numeric(sinfo$Xchr_ratio)
umap_ext <- data.frame(umap1 = emb_umap[,1], umap2 = emb_umap[,2], 
                       sinfo, t(expr_markers))

print("3.3) Plot UMAP embedding and RNA-velocity predictions")
emb <- cbind(umap_ext[,1], umap_ext[,2]); rownames(emb) <- rownames(umap_ext)
m <- match(rownames(emb), cellfeatures$id)
cellfeatures$Xchr_ratio <- as.numeric(cellfeatures$Xchr_ratio)
cellfeatures$AXCR <- abs(cellfeatures$Xchr_ratio-0.5)
colors <- data.frame(sinfo, t(expr_markers))
colors$AXCR <- abs(colors$Xchr_ratio - 0.5)

plot_features <- c("Time", "AXCR", markergenes); plot_list <- list()
x <- show.velocity.on.embedding.cor(emb, vel = notAS_vel, n=neighbor_size, scale=velocity_scale, 
                                    cex=1.5, arrow.scale=10, show.grid.flow=TRUE, 
                                    min.grid.cell.mass=5, grid.n=20,
                                    arrow.lwd=arrow.lwd, do.par=F, cell.border.alpha = 0.5,
                                    return.details = TRUE)

# make the arrow width proportional to the velocity extent
arrowthick <- 3e-1; arrowidth <- 5*1e-2; arrowsharp <- 40
temp <- data.frame(emb, colors, Xist_count = spliced$counts["Xist",]); gridvelo <- x$garrows
x$garrows <- data.frame(x$garrows)
x$garrows$dist <- sqrt((x$garrows$x0 - x$garrows$x1)^2 + (x$garrows$y0 - x$garrows$y1)^2)
x$garrows$arrowidth <- arrowidth*(x$garrows$dist/max(x$garrows$dist))

# include day variable
temp$day <- factor(as.numeric(gsub(temp$Time, pattern = "h", replacement = ""))/24)

# color cells based on sequencing time point
g <- temp %>%
  ggplot() + 
  theme_bw() + theme1 +
  geom_point(aes(x = X1, y = X2, color = day), size = outliersize, shape = 20, alpha = 3/4) + 
  xlab("UMAP 1") + ylab("UMAP 2") + labs(color = "Time [days]") +
  geom_segment(data = data.frame(x$garrows), 
               aes(x = x0, y = y0, xend = x1, yend = y1),
               size = arrowthick, alpha = 3/5,
               arrow = arrow(length = unit(x$garrows$arrowidth, "inches"), angle = arrowsharp)) +
  scale_color_manual(values = time_colors) +
  guides(color = guide_legend(override.aes = list(size=2)))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 5, savefile = paste0(outpath, "C_UMAP_time.pdf"))

# color cells based on marker gene expression
markers <- rev(c("Xist", "Dnmt3a", "Esrrb", "Nanog"))
umap_markers <- reshape2::melt(data = temp, id.vars = c("X1", "X2"), 
                               measure.vars = markers)

g <- umap_markers %>%
  ggplot() + 
  facet_grid(.~variable) +
  theme_bw() + theme1 +
  # scale_y_reverse() +
  geom_point(aes(x = X1, y = X2, color = value), size = outliersize, shape = 20) + 
  scale_colour_gradientn(colours = c("black", "gold")) + 
  labs(x="", y = "", color = expression(log[10]*"(CPM + 1) values")) + 
  xlab("") + ylab("") + 
  guides(color = guide_colourbar(barwidth = 4, barheight = 0.7)) + 
  theme(legend.position = "top")
adjust_size(g = g, panel_width_cm = 2.5, panel_height_cm = 2.5, savefile = paste0(outpath, "D_UMAP_markers.pdf"))


print("4) E: Xist UMI/CPM distribution and % of Xist-UMI>5 cells per sequencing time point")

print("4.1) Load and process data")
load(paste0(datapath, "DGE.RData"))
df <- dge$samples[, c("day", "id", "xist_umi", "xist_cpm")]
nXistUnd <- ddply(df, .variables = .(day), summarize, n = mean(xist_umi <= 5))
nXistUnd$label <- paste0(round((1-nXistUnd$n)*100, digits = 0), "%")
df_melt <- reshape2::melt(df, id.vars = c("day", "id"), measure.vars = c("xist_umi", "xist_cpm"))
df_melt$variable <- revalue(df_melt$variable, replace = c("xist_umi" = "Xist UMI", "xist_cpm" = "Xist CPM"))
df_melt$day <- factor(df_melt$day)

print("4.2) Plot")
# barplot
nXistUnd$day <- factor(nXistUnd$day, levels = rev(levels(factor(nXistUnd$day))))
g <- nXistUnd %>%
  ggplot(aes(x = day, y = (1-n)*100)) + 
  theme_bw() + theme1 + 
  ylab(expression('Xist' ^ '+'*' cells [%]')) + xlab("Time [days]") +
  geom_bar(stat = "identity", width = .5) +
  geom_text(aes(label=label),position="stack", hjust = -0.2, color = "black", size = geomtext_size)+
  coord_flip() + 
  scale_y_continuous(limits = c(0, 130), breaks = seq(0,100,50))
adjust_size(g = g, panel_width_cm = 1.5, panel_height_cm = 3, savefile = paste0(outpath, "E_XistPlus_barplot.pdf"))

# xist UMI and CPM - violin plot
df_melt$day <- factor(df_melt$day, levels = rev(levels(factor(df_melt$day))))
g <- df_melt %>% 
  ggplot(aes(x = day, y = log10(value+1))) + 
  theme_bw() + theme1 + 
  geom_violin(adjust = 1, size = violin_box_size, scale = "width") +
  geom_jitter(alpha = 0.2, width = .15, size = outliersize, shape = 21) + 
  facet_grid(.~variable, scales = "free_x") +
  labs(x="Time [days]", y = expression("Xist ["*log[10]*"(value + 1)]")) + 
  coord_flip()
adjust_size(g = g, panel_width_cm = 1.5, panel_height_cm = 3, savefile = paste0(outpath, "E_XistExpression_violin.pdf"))



print("5) F: X/A ratio using bootstrapping autosomal features (n= #x-linked genes)")

print("5.1) Bootstrapped X/A ratio")
load(paste0(datapath, "DGE.RData"))
nboot <- 1e3
xcounts <- dge$counts[dge$genes$chromosome %in% "X",]; autcounts <- dge$counts[dge$genes$chromosome %in% c(1:19),]
n_xlinked <- nrow(xcounts); n_autlinked <- nrow(autcounts)

temp_boot <- c()
for (j in seq_len(nboot)) {
  boot <- sample(x = seq_len(n_autlinked), size = n_xlinked, replace = TRUE)
  value <- colSums(xcounts)/colSums(autcounts[boot,])
  temp_boot <- rbind(temp_boot, value)
}
x2a_boot <- apply(temp_boot, 2, median)
x2a <- data.frame(day = dge$samples$day, id = colnames(dge), 
                  Xist_UMI = dge$counts["Xist",],
                  xist_detection = dge$samples$Xist_class, 
                  x2a = x2a_boot)
x2a$xist_detection <- revalue(factor(x2a$xist_detection), 
                              replace = c("Detected (Xist UMI > 5)" = "Xist+ cells",
                                          "Detected (Xist UMI <= 5)" = "Xist+ cells [UMI<=5]",
                                          "Undetected" = "Xist- cells"))
x2a$xist_detection <- factor(x2a$xist_detection, levels = c("Xist- cells", "Xist+ cells [UMI<=5]", "Xist+ cells"))
table(x2a$xist_detection, x2a$day)
summ_cell <- ddply(x2a, .variables = .(day, xist_detection), summarize, n = length(day))

print("5.2) Plot Xist_UMI=0 and Xist_UMI>5 cells")
features <- c("id", "day", "xist_detection", "x2a")
total <- rbind(data.frame(x2a[, features], Genes = "Per time & Bootstrap"))
colors <- c("black", "#e7298a")
total$Xist <- total$xist_detection
summ_cell$Xist <- summ_cell$xist_detection
total$day <- factor(total$day); summ_cell$day <- factor(summ_cell$day)

g <- total[!total$Xist %in% "Xist+ cells [UMI<=5]",] %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = seq(0.5, 2, by = 0.5), linetype = "dashed", alpha = hvalpha, size = linesize) + 
  scale_y_continuous(breaks = seq(0.5, 2, by = 0.5), limits = c(0.25, 2.2)) + 
  geom_boxplot(aes(x = day, y = x2a, color = Xist), outlier.shape = NA, size = violin_box_size) + 
  geom_jitter(aes(x = day, y = x2a, color = Xist), fill = "black",
              alpha = 1/10, position=position_jitterdodge(jitter.width = .25, dodge.width = 0.75), 
              size = outliersize, show.legend = FALSE) +
  scale_color_manual(values = colors) + 
  geom_text(summ_cell[!summ_cell$Xist %in% "Xist+ cells [UMI<=5]",],
            mapping = aes(x = day, y = Inf, label = paste0("n = ", n), fill = Xist),
            position=position_dodge(.75), angle = 90, alpha = 1, size = geomtext_size, hjust = -0.5) +
  coord_cartesian(clip = "off") +
  labs(x = "Time [days]", y = "X:A ratio", fill = "")
adjust_size(g = g, panel_width_cm = 3, panel_height_cm = 3, savefile = paste0(outpath, "F_X2A.pdf"))


print("5.3) Compare d4 and d0 cells using an unpaired Wilcoxon test")
wilcox.test(x = total$x2a[(total$day == 0)],
            y = total$x2a[(total$day == 4)&(total$xist_detection == "Xist+ cells")])
wilcox.test(x = total$x2a[(total$day == 0)],
            y = total$x2a[(total$day == 4)&(total$xist_detection == "Xist- cells")])


print("6) G: UMAP colored by X/A ratio for Xist-negative and Xist-positive cells")

print("6.1) Load data")
umap_ext$id <- rownames(umap_ext)
m <- match(umap_ext$id, x2a$id)
var <- c("day", "Xist_UMI", "xist_detection", "x2a")
umap <- data.frame(umap_ext, x2a[m, var])
save(umap, file = paste0(outpath, "umap.RData"))

print("6.2) Plot")
temp <- umap[!umap$xist_detection %in% "Xist+ cells [UMI<=5]",]
mid <- median(temp$x2a[temp$day == 0]); print(paste0("The median X:A ratio for d0 cells is ", round(mid, digits = 3)))
g <- temp %>%
  ggplot() + 
  theme_bw() + theme1 +
  facet_grid(.~xist_detection) +
  geom_point(aes(x = umap1, y = umap2, color = x2a), size = outliersize, shape = 20) + 
  scale_color_gradient2(midpoint = mid, low = "blue", mid = "grey", high = "red") +
  labs(x="UMAP 1", y = "UMAP 2", color = "X:A ratio") + 
  guides(color = guide_colourbar(barwidth = 0.7, barheight = 4)) + 
  theme(legend.position = "right",
        strip.text.y = element_text(angle = 0))
adjust_size(g = g, panel_width_cm = 2.5, panel_height_cm = 2.5, savefile = paste0(outpath, "G_UMAP_X2A.pdf"))
