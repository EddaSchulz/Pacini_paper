print("1) Define output path")
outpath <- paste0(path, "output/fig6_GlobalSilencing/"); dir.create(path = outpath, recursive = T, showWarnings = F)

print("2) A: X Chromosome ratio and Xist detection")

print("2.1) Load data and plot XR for cells with Xist_UMI = 0 and Xist_UMI>5")
colors <- rev(comparison_colors)
temp <- res
temp$Xchr_ratio_noXist <- temp$b6_X_noXist/(temp$b6_X_noXist + temp$cast_X_noXist)
temp$Xist <- revalue(factor(temp$Xist), 
                     replace = c("Xist+ cells [UMI>5]" = "Xist+ cells",
                                 "Xist+ cells [UMI<=5]" = "Xist+ cells [UMI<=5]",
                                 "Xist- cells [UMI=0]" = "Xist- cells"))
temp$Xist <- factor(temp$Xist, levels = c("Xist- cells", 
                                          "Xist+ cells [UMI<=5]",
                                          "Xist+ cells"))
g <- temp[!temp$Xist %in% "Xist+ cells [UMI<=5]",] %>% 
  ggplot() +
  theme_bw() +  theme1 + 
  scale_y_continuous(breaks = c(0, 0.2, 0.5, 0.8, 1)) +
  geom_violin(aes(x = factor(day), y = Xchr_ratio_noXist, colour = Xist), alpha = 1/2, draw_quantiles = c(0.5), size = violin_box_size, show.legend = FALSE) + 
  geom_jitter(aes(x = factor(day), y = Xchr_ratio_noXist, colour = Xist), alpha = 1/4, size = small_scattersize,
              position=position_jitterdodge(jitter.width = .5, dodge.width = 0.9), shape = 21) +
  labs(x = "Time [days]", 
       y = "Chr.X [B6/Total]", 
       colour = "") +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(override.aes = list(size = 1, shape = 20, alpha = 1)))
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "A_XCR_XistDetection.pdf"))


print("3) B: Spliced and Unspliced X-reads")

print("3.1) Load allele specific spliced and unspliced data")

df <- res
load(paste0(datapath, "NCF_DGE_B6_Spliced.RData")); spliced_b6 <- dge
load(paste0(datapath, "NCF_DGE_B6_Unspliced.RData")); unspliced_b6 <- dge
load(paste0(datapath, "UF_DGE_B6.RData")); b6 <- dge[, colnames(spliced_b6)]
load(paste0(datapath, "NCF_DGE_Cast_Spliced.RData")); spliced_cast <- dge
load(paste0(datapath, "NCF_DGE_Cast_Unspliced.RData")); unspliced_cast <- dge
load(paste0(datapath, "UF_DGE_Cast.RData")); cast <- dge[, colnames(spliced_cast)]

commongenes <- intersect(spliced_b6$genes$symbol, b6$genes$symbol)
spliced_b6 <- spliced_b6[match(commongenes, spliced_b6$genes$symbol),]
unspliced_b6 <- unspliced_b6[match(commongenes, unspliced_b6$genes$symbol),]
b6 <- b6[match(commongenes, b6$genes$symbol),]
spliced_cast <- spliced_cast[match(commongenes, spliced_cast$genes$symbol),]
unspliced_cast <- unspliced_cast[match(commongenes, unspliced_cast$genes$symbol),]
cast <- cast[match(commongenes, cast$genes$symbol),]

as_df <- data.frame(day = rep(cast$samples$day, each = nrow(cast$counts)),
                    Cell = rep(cast$samples$id, each = nrow(cast$counts)),
                    Gene = rep(cast$genes$symbol, times = ncol(cast$counts)),
                    Chromosome = rep(cast$genes$chromosome, times = ncol(cast$counts)),
                    b6_exon = c(spliced_b6$counts),
                    b6_intron = c(unspliced_b6$counts),
                    b6_AS = c(b6$counts),
                    cast_exon = c(spliced_cast$counts),
                    cast_intron = c(unspliced_cast$counts),
                    cast_AS = c(cast$counts)
)

fields <- c("Xist_ratio", "Xist_classification")
m <- match(as_df$Cell, df$id)
as_df <- data.frame(as_df, df[m, fields])
as_df <- as_df[as_df$Chromosome %in% c(1:19, "X", "Y"),]

print("3.2) Exclude Xist low and Skewed groups and plot total number of spliced and unspliced X-linked counts per cell")

excludegroup <- c("Low", "Skewed")

# spliced counts
sum_spliced <- ddply(as_df[!as_df$Xist_classification %in% excludegroup,], .variables = .(day, Cell),
                     summarize, 
                     Xist_classification = unique(Xist_classification),
                     b6 = sum(b6_exon[Chromosome %in% "X"]),
                     cast = sum(cast_exon[Chromosome %in% "X"]))
sum_spliced$Day <- factor(paste0("Day ", sum_spliced$day))
g <- sum_spliced %>%
  ggplot() + 
  theme_bw() + theme1 + facet_grid(. ~ Day) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = linesize, alpha = hvalpha) + 
  geom_point(aes(x = b6, y = cast, color = Xist_classification), size = outliersize) + 
  scale_color_manual(values = color_alleles) + 
  scale_y_continuous(breaks = seq(0, 800, 200), limits = c(0, 820)) +  
  scale_x_continuous(breaks = seq(0, 800, 200), limits = c(0, 820)) +  
  labs(x = expression(Spliced[chrX] ^ 'B6'),
       y = expression(Spliced[chrX] ^ 'Cast'),
       color = "") + 
  theme(legend.position = "top")
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 2, savefile = paste0(outpath, "B_Xchr_spliced.pdf"))

# unspliced counts
sum_unspliced <- ddply(as_df[!as_df$Xist_classification %in% excludegroup,], .variables = .(day, Cell), summarize, 
                       Xist_classification = unique(Xist_classification),
                       b6 = sum(b6_intron[Chromosome %in% "X"]),
                       cast = sum(cast_intron[Chromosome %in% "X"]))
sum_unspliced$Day <- factor(paste0("Day ", sum_unspliced$day))
g <- sum_unspliced %>%
  ggplot() + 
  theme_bw() + theme1 + facet_grid(. ~ Day) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = linesize, alpha = hvalpha) + 
  geom_point(aes(x = b6, y = cast, color = Xist_classification), size = outliersize) + 
  scale_color_manual(values = color_alleles) + 
  scale_y_continuous(breaks = seq(0, 200, 50), limits = c(0, 230)) +  
  scale_x_continuous(breaks = seq(0, 200, 50), limits = c(0, 230)) +  
  labs(x = expression(Unspliced[chrX] ^ 'B6'),
       y = expression(Unspliced[chrX] ^ 'Cast'),
       color = "") + 
  theme(legend.position = "top")
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 2, savefile = paste0(outpath, "B_Xchr_unspliced.pdf"))


print("4) C: Allelic RNA velocity")

print("4.1) Load data")

# load cell features
load(paste0(datapath, "DGE.RData"))
features_cell <- c("day", "id", "xist_ratio", "Xist_ratio_class", "xchr_ratio", "sf_notX")
cellfeatures <- dge$samples[, features_cell]
cellfeatures$sf <- cellfeatures$sf_notX
cellfeatures$xist_ratio <- as.numeric(cellfeatures$xist_ratio)
cellfeatures$Xchr_ratio <- as.numeric(cellfeatures$xchr_ratio)
cellfeatures$AXCR <- abs(cellfeatures$Xchr_ratio-0.5)

# load spliced and unspliced count matrices
load(paste0(datapath, "NCF_DGE_B6_Spliced.RData")); spliced_b6 <- dge
load(paste0(datapath, "NCF_DGE_B6_Unspliced.RData")); unspliced_b6 <- dge
load(paste0(datapath, "NCF_DGE_Cast_Spliced.RData")); spliced_cast <- dge
load(paste0(datapath, "NCF_DGE_Cast_Unspliced.RData")); unspliced_cast <- dge

spliced <- list(BL6 = spliced_b6$counts, CastEiJ = spliced_cast$counts)
unspliced <- list(BL6 = unspliced_b6$counts, CastEiJ = unspliced_cast$counts)
spliced$Both <- spliced$BL6 + spliced$CastEiJ
unspliced$Both <- unspliced$BL6 + unspliced$CastEiJ
cellfeatures$Xist <- spliced$Both["Xist",]
cellfeatures$Xist_ratio_correct <- cellfeatures$xist_ratio
cellfeatures$Xist_ratio_correct[is.na(cellfeatures$xist_ratio)|cellfeatures$Xist<5] <- 0.5
cellfeatures$Xist_ratio_correct <- abs(cellfeatures$Xist_ratio_correct-0.5)
genefeatures <- spliced_b6$genes

print("4.2) Compute X-linked B6-ratio matrix")

b6ratio <- spliced$BL6/(spliced$BL6 + spliced$CastEiJ)
b6ratio[is.na(b6ratio)] <- 0.5 # 0/0 ratios are set to 0.5
isX <- spliced_b6$genes$chromosome %in% "X"; table(isX)
cell_Xcor <- cor(b6ratio[spliced_b6$genes$chromosome %in% "X",])

print("4.3) Remove genes with low spliced or unspliced gene expression")

rm_s <- rowMeans(spliced$Both); rm_u <- rowMeans(unspliced$Both)
summary(rm_s); summary(rm_u)
keep <- (rm_s >= 0.5)&(rm_u >= 0.25); table(keep)
spliced_filt <- lapply(spliced, function(x) x[keep,]); unspliced_filt <- lapply(unspliced, function(x) x[keep,])
cells <- colnames(spliced_filt$Both)
keep_genes <- intersect(rownames(spliced_filt$Both), genefeatures$symbol)
genefeatures <- genefeatures[keep_genes,]

print("4.4) Define cell-to-cell distance based on B6-ratio X-linked correlation, and restrict kNN choice to cells sequenced at the same time point")

emat <- spliced_filt$Both; nmat <- unspliced_filt$Both
celldist <- cell_Xcor
days <- strsplit2(rownames(celldist), split = "\\_")[,1]; table(days)
for(i in 1:nrow(celldist)){
  celldist[i,][!days %in% days[i]] <- min(celldist[i,])
}
celldist_adj <- as.dist(1-celldist) # adjust correlation to be a distance matrix to feed to velocyto.R::balanced_knn funct.

print("4.5) Compute velocities")

ncores <- min(c(detectCores(), 10))
fit.quantile <- 0.025; kcells <- 20; mincor <- 0.05
arrow.scale=5; cell.alpha=0.4; cell.cex=1; fig.height=4; fig.width=4.5;
AS_vel <- gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells = kcells, cell.dist = celldist_adj,
                                           fit.quantile = fit.quantile, 
                                           min.nmat.emat.slope = mincor,
                                           min.nmat.emat.correlation = mincor, n.cores = ncores)
save(AS_vel, file = paste0(outpath, "AS_vel.RData"))

print("4.6) Compute PCA cell embedding on X-linked features of B6-ratio matrix")

isX <- as.character(genefeatures$chromosome[match(rownames(AS_vel$current), rownames(genefeatures))]) %in% "X"
Xgenes <- rownames(AS_vel$current)[isX]
x0 <- b6ratio[match(Xgenes, rownames(b6ratio)),]
x0.log <- x0
cent <- rowMeans(x0.log)

# PCA embedding
epc <- pcaMethods::pca(t(x0.log - cent), center = F, nPcs = length(cent))
epc@scores <- scale(epc@completeObs, scale = F, center = T) %*% epc@loadings
curr_pos <- epc@scores
neighbor_size <- 100; velocity_scale <- "sqrt"; arrow.scale <- 5; arrow.lwd <- 1
emb <- cbind(curr_pos[,1], curr_pos[,2]); rownames(emb) <- rownames(curr_pos)
m <- match(rownames(emb), cellfeatures$id)
cellfeatures$Xchr_ratio <- as.numeric(cellfeatures$Xchr_ratio)
cellfeatures$AXCR <- abs(cellfeatures$Xchr_ratio-0.5)
features <- c("day", "Xchr_ratio", "Xist_ratio_class")
sinfo <- cellfeatures[match(rownames(emb), cellfeatures$id), features]
pca_ext <- data.frame(pca1 = emb[,1], pca2 = emb[,2], sinfo)

# project velocities on PCA embedding
colors <- data.frame(sinfo)
colors$AXCR <- abs(colors$Xchr_ratio - 0.5)
x <- show.velocity.on.embedding.cor(emb, vel = AS_vel, n=neighbor_size, scale=velocity_scale, 
                                    cex=1.5, arrow.scale=arrow.scale, show.grid.flow=TRUE, 
                                    min.grid.cell.mass=5, grid.n=20,
                                    arrow.lwd=arrow.lwd, do.par=F, cell.border.alpha = 0.5,
                                    return.details = TRUE)

print("4.7) Plot")
arrowthick <- 4e-1; arrowidth <- 1e-1; arrowsharp <- 40
temp <- data.frame(emb, colors); gridvelo <- x$garrows
x$garrows <- data.frame(x$garrows)
x$garrows$dist <- sqrt((x$garrows$x0 - x$garrows$x1)^2 + (x$garrows$y0 - x$garrows$y1)^2)
x$garrows$arrowidth <- arrowidth*(x$garrows$dist/max(x$garrows$dist))
top_x <- max(temp[,1]); top_y <- max(temp[,2])
cols <- c(allele_colors["B6"], "grey", allele_colors["Cast"])

g <- temp %>%
  ggplot() + 
  theme_bw() + theme1 + 
  geom_point(aes(x = X1, y = X2, fill = Xchr_ratio, color = Xchr_ratio), size = scattersize) + 
  geom_segment(data = data.frame(x$garrows), linejoin='mitre',
               aes(x = x0, y = y0, xend = x1, yend = y1), size = arrowthick, alpha = 3/5, color = "black",
               arrow = arrow(length = unit(x$garrows$arrowidth, "inches"), angle = arrowsharp)) + 
  scale_color_gradientn(colours = cols) +
  scale_fill_gradientn(colours = cols, breaks = c(0.2, 0.5, 0.8)) +
  labs(x = paste0("PCA 1 (", round(epc@R2[1]*100, digits = 1), "%)"),
       y = paste0("PCA 2 (", round(epc@R2[2]*100, digits = 1), "%)"),
       fill = expression("X Chromosome Ratio [XCR = "*B6[chrX]*"/("*B6[chrX]+Cast[chrX]*")]")) + 
  guides(color = FALSE) 
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 5, savefile = paste0(outpath, "C_XCI_velocity.pdf"))


print("5) D: Xist and X-chromosome allelic ratio")

exclude_less5XistUMI <- c("Undetected", "Low")
g <- res[!res$Xist_classification %in% exclude_less5XistUMI,] %>%
  ggplot() + 
  theme_bw() + theme1 + 
  facet_grid(.~Day) + 
  geom_point(aes(x = Xist_ratio, y = Xchr_ratio), size = outliersize) +  
  labs(x = "Xist [B6/Total]",
       y = "Chr.X [B6/Total]") + 
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0,1)) + 
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0,1))
adjust_size(g = g, panel_width_cm = 2, panel_height_cm = 2, savefile = paste0(outpath, "D_XchrRatio_XistRatio.pdf"))


print("6) E: X/A ratio dividing cells by Xist AS classification")

print("6.1) Load notAS data")

load(paste0(datapath, "DGE.RData"))
dge$samples$Xist_classification <- revalue(factor(dge$samples$Xist_ratio_class), 
                                           c("Undetected" = "Undetected", "Low-Xist" = "Low",
                                             "Middle" = "Skewed", "Xist_BA" = "BA",
                                             "BL6_MA" = "Xist-MA (Xi=B6)", "Cast_MA" = "Xist-MA (Xi=Cast)"))
lev_ref <- c("Undetected", "Low", "Skewed", "BA", "Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)")
dge$samples$Xist_classification <- factor(dge$samples$Xist_classification, levels = lev_ref)

print("6.2) Compute X/A ratios")

xcounts <- dge$counts[dge$genes$chromosome %in% "X",]; autcounts <- dge$counts[dge$genes$chromosome %in% c(1:19),]
n_xlinked <- nrow(xcounts); n_autlinked <- nrow(autcounts)
temp_boot <- c()
for (j in seq_len(nboot)) {
  boot <- sample(x = seq_len(n_autlinked), size = n_xlinked, replace = TRUE)
  value <- colSums(xcounts)/colSums(autcounts[boot,])
  temp_boot <- rbind(temp_boot, value)
}
x2a_boot <- apply(temp_boot, 2, median)
x2a <- data.frame(day = dge$samples$day,
                  id = colnames(dge),
                  xist_classification = dge$samples$Xist_classification, 
                  x2a = x2a_boot)

print("6.3) Plot")

summ_cell <- ddply(x2a, .variables = .(day, xist_classification), summarize, n = length(day))
features <- c("id", "day", "xist_classification", "x2a")
total <- rbind(data.frame(x2a[, features], Genes = "Per time & Bootstrap"))
total$Xist <- total$xist_classification
summ_cell$Xist <- summ_cell$xist_classification
restr_lev <- c("Undetected", "BA", "Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)")

# remove BA cells at d4 since there is a single cell
total <- total[!(total$day == 4 & total$Xist %in% "BA"),]
summ_cell <- summ_cell[!(summ_cell$day == 4 & summ_cell$Xist %in% "BA"),]

# numbers
g <- total[(total$Genes %in% "Per time & Bootstrap") & (total$Xist %in% restr_lev),] %>%  
  ggplot() +
  theme_bw() + theme1 + 
  geom_hline(yintercept = seq(0.5, 2, by = 0.5), linetype = "dashed", alpha = hvalpha, size = linesize) + 
  scale_y_continuous(breaks = seq(0.5, 2, by = 0.5), limits = c(0.25, 2.1)) + 
  geom_boxplot(aes(x = factor(day), y = x2a, color = Xist), outlier.shape = NA, size = violin_box_size) + 
  geom_jitter(aes(x = factor(day), y = x2a, color = Xist), fill = "black",
              alpha = 1/4, position=position_jitterdodge(jitter.width = .25, dodge.width = 0.75), 
              size = small_scattersize, show.legend = FALSE) +
  scale_color_manual(values = color_alleles) + 
  geom_text(summ_cell[summ_cell$Xist %in% restr_lev,],
            mapping = aes(x = factor(day), y = Inf, label = paste0("n = ", n), group = Xist),
            position=position_dodge(.75), angle = 90, alpha = 1, size = geomtext_size, hjust=-0.5) +
  coord_cartesian(clip = "off") +
  xlab("Time [days]") + ylab("X:A ratio") + labs(fill = "")
adjust_size(g = g, panel_width_cm = 5, panel_height_cm = 3, savefile = paste0(outpath, "E_X2A_XistGroup.pdf"))

print("6.4) Compare X/A ratios between Xist-BA and Xist-MA cells")
g1 <- "BA"; g2 <- c("Xist-MA (Xi=B6)", "Xist-MA (Xi=Cast)")
test <- ddply(total[total$day %in% 1:3,], .variables = .(day), summarize, 
              XistBA = sum(xist_classification %in% g1), XistBA_medX2A = median(x2a[xist_classification %in% g1]),
              XistMA = sum(xist_classification %in% g2), XistMA_medX2A = median(x2a[xist_classification %in% g2]),
              pvalue_wmw = wilcox.test(x = x2a[xist_classification %in% g1], y = x2a[xist_classification %in% g2])$p.value)
test