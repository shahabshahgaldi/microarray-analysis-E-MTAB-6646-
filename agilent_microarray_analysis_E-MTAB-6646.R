# ============================================================
# R script: Agilent one-color microarray analysis (E-MTAB-6646)
# ============================================================

# Workflow:

# 1. Create design matrix and contrasts
# 2. Import raw Agilent .txt files
# 3. Remove control probes
# 4. Background correction and quantile normalization (with optional batch correction)
# 5. Quality control plots
# 6. Filter low-expressed probes (group-aware)
# 7. Collapse duplicate probes
# 8. Fit linear model (limma)
# 9. Extract differential expression results
# 10. Export results to CSV
# 11. Generate volcano plots with top 20 gene labels

# ----------------------------- Install and load dependencies -----------------------------

req.cran <- c("ggplot2", "dplyr")
for (p in req.cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

req.bioc <- c("limma", "ggrepel")
for (p in req.bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)

library(limma)
library(ggplot2)
library(dplyr)
library(ggrepel)

# ----------------------------- Set working directory -----------------------------

setwd("C:/..../E-MTAB-6646")     # Replace it with the location on your pc containg the raw Agilent.txt files

# ----------------------------- Define sample metadata (targets) -----------------------------

targets <- data.frame(
  FileName = c(
    "257480910350_201704071328_S01_GE1_107_Sep09_1_1.txt", "257480910350_201704071328_S01_GE1_107_Sep09_2_3.txt",
    "257480910351_201704071352_S01_GE1_107_Sep09_1_2.txt", "257480910351_201704071352_S01_GE1_107_Sep09_2_4.txt",
    "257480910350_201704071328_S01_GE1_107_Sep09_1_3.txt", "257480910350_201704071328_S01_GE1_107_Sep09_2_1.txt",
    "257480910351_201704071352_S01_GE1_107_Sep09_2_2.txt", "257480910351_201704071352_S01_GE1_107_Sep09_2_3.txt",
    "257480910350_201704071328_S01_GE1_107_Sep09_1_4.txt", "257480910350_201704071328_S01_GE1_107_Sep09_2_2.txt",
    "257480910351_201704071352_S01_GE1_107_Sep09_1_3.txt", "257480910351_201704071352_S01_GE1_107_Sep09_2_1.txt",
    "257480910350_201704071328_S01_GE1_107_Sep09_2_4.txt", "257480910351_201704071352_S01_GE1_107_Sep09_1_1.txt",
    "257480910351_201704071352_S01_GE1_107_Sep09_1_4.txt", "257480910350_201704071328_S01_GE1_107_Sep09_1_2.txt"
  ),
  Sample = c(
    "Sample5", "Sample6", "Sample7", "Sample8",
    "Sample13", "Sample14", "Sample15", "Sample16",
    "Sample1", "Sample2", "Sample3", "Sample4",
    "Sample10", "Sample11", "Sample12", "Sample9"
  ),
  CellType = c(
    "adipocyte", "adipocyte", "adipocyte", "adipocyte",
    "preadipocyte", "preadipocyte", "preadipocyte", "preadipocyte",
    "adipocyte", "adipocyte", "adipocyte", "adipocyte",
    "preadipocyte", "preadipocyte", "preadipocyte", "preadipocyte"
  ),
  Treatment = c(
    "mock", "mock", "mock", "mock",
    "mock", "mock", "mock", "mock",
    "IAV", "IAV", "IAV", "IAV",
    "IAV", "IAV", "IAV", "IAV"
  ),
  stringsAsFactors = FALSE
)

# Attach path to filenames

targets$FileName <- paste0(getwd(), "/", targets$FileName)

# ----------------------------- Create design matrix & contrasts -----------------------------

group <- factor(paste(targets$CellType, targets$Treatment, sep = "."))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

cm <- makeContrasts(
  pread.IAV.vs.mock = preadipocyte.IAV - preadipocyte.mock,
  adip.IAV.vs.mock  = adipocyte.IAV   - adipocyte.mock,
  levels = design
)

# ----------------------------- Import raw Agilent one-color data -----------------------------

RG <- read.maimages(
  files = targets,
  source = "agilent",
  green.only = TRUE,
  columns = list(G = "gMedianSignal", Gb = "gBGMedianSignal"),
  annotation = c("FeatureNum", "ControlType", "ProbeName", "GeneName", "SystematicName", "Description")
)

# ----------------------------- Filter control probes -----------------------------

table(RG$genes$ControlType)
RG <- RG[RG$genes$ControlType == 0, ]
table(RG$genes$ControlType)
nrow(RG)

# ----------------------------- Background correction & normalization -----------------------------

E.bg <- backgroundCorrect(RG, method = "normexp", offset = 50)

# min <= 0, increase the offset
summary(as.vector(E.bg$E))
min(E.bg$E)

# Normalization and log transformation
summary(as.vector(E.bg$E))
E.norm <- normalizeBetweenArrays(E.bg, method = "quantile")
summary(as.vector(E.norm$E))

# Optional batch correction (check QC plots)
batch <- factor(sub("_.+", "", basename(targets$FileName)))
E.norm$E <- removeBatchEffect(E.norm$E, batch = batch, design = design)

exprs <- E.norm$E

# ----------------------------- QC plots -----------------------------

pdf("QC_plots.pdf", width = 10, height = 8)
par(mfrow = c(2,2))

boxplot(E.bg$E, main = "After Background Correction (linear)", col = "lightblue", ylab = "Expression")
boxplot(exprs, las = 2, main = "After Quantile Norm (log2)", col = "lightblue", ylab = "Expression (log2)")
plotDensities(E.bg$E, main = "Density Before Norm", col = "red", legend = FALSE)
plotDensities(exprs, main = "Density After Norm", col = "blue", legend = FALSE)

par(mfrow = c(1,1))
plotMDS(exprs, labels = paste(targets$CellType, targets$Treatment, sep = "_"),
        col = as.numeric(factor(targets$CellType)), main = "MDS Plot")
legend("topright", legend = levels(factor(targets$CellType)),
       col = 1:length(levels(factor(targets$CellType))), pch = 16)

dev.off()

# ----------------------------- Filter low-expressed probes (group-aware) -----------------------------

thresholds <- c(5, 5.5, 6, 6.5, 7)
for(t in thresholds){
  keep <- rowSums(exprs > t, na.rm=TRUE) >= 4
  cat("th=", t, "  kept=", sum(keep), " (", round(sum(keep)/nrow(exprs)*100,2), "% )\n")
}

group.factor <- factor(paste(targets$CellType, targets$Treatment, sep = "."))
group.expr <- split.data.frame(t(exprs), group.factor)
keep.per.group <- sapply(group.expr, function(g) rowSums(t(g) >= 6.5) >= 2)
keep <- rowSums(keep.per.group) >= 1
exprs.filt <- exprs[keep, , drop = FALSE]

probe.names <- RG$genes$ProbeName[keep]
gene.names <- RG$genes$GeneName[keep]
gene.names[is.na(gene.names) | gene.names==""] <- probe.names

# ----------------------------- Collapse duplicate probes -----------------------------

exprs.filt <- aggregate(exprs.filt, by = list(ProbeName = probe.names), FUN = mean)

rownames(exprs.filt) <- exprs.filt$ProbeName
exprs.filt$ProbeName <- NULL

gene.names.collapsed <- tapply(gene.names, probe.names, function(x) unique(x)[1])

probe.annotation <- data.frame(
  ProbeName = names(gene.names.collapsed),
  GeneName  = unname(gene.names.collapsed),
  row.names = names(gene.names.collapsed),
  stringsAsFactors = FALSE
)

colnames(exprs.filt) <- targets$Sample
targets <- targets[match(colnames(exprs.filt), targets$Sample), ]

# ----------------------------- Fit linear model with limma -----------------------------

fit <- lmFit(exprs.filt, design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)

# ----------------------------- Extract DE results -----------------------------

# For preadipocytes (IAV vs mock)
results.preadip <- topTable(fit2, coef="pread.IAV.vs.mock", number=Inf, adjust.method="fdr")
results.preadip <- cbind(ProbeName = rownames(results.preadip),
                         GeneName  = probe.annotation[rownames(results.preadip), "GeneName"],
                         results.preadip)

# For adipocytes (IAV vs mock)
results.adip <- topTable(fit2, coef="adip.IAV.vs.mock", number=Inf, adjust.method="fdr")
results.adip <- cbind(ProbeName = rownames(results.adip),
                      GeneName  = probe.annotation[rownames(results.adip), "GeneName"],
                      results.adip)

# Save CSV
getwd()
write.csv(results.preadip, "DE.results.preadipocytes.IAV.vs.mock.csv", row.names=FALSE)
write.csv(results.adip,   "DE.results.adipocytes.IAV.vs.mock.csv",   row.names=FALSE)

# ----------------------------- Volcano plots with top 20 gene labels -----------------------------

make.volcano.plot <- function(results.df, title,
                              fdr.threshold=0.05, logfc.threshold=1,
                              p.floor=1e-300, topN=10) {
  df <- as.data.frame(results.df)
  df$P.Value.safe <- pmax(df$P.Value, p.floor)
  df$negLog10P <- -log10(df$P.Value.safe)
  
  # Define significance groups
  df$sig <- "NS"
  df$sig[df$adj.P.Val < fdr.threshold & df$logFC >  logfc.threshold] <- "Upregulated"
  df$sig[df$adj.P.Val < fdr.threshold & df$logFC < -logfc.threshold] <- "Downregulated"
  df$sig <- factor(df$sig, levels=c("NS", "Upregulated", "Downregulated"))
  
  # Pick topN genes by P-value, but remove non-informative names (chr*, *_*)
  top.genes <- df %>%
    arrange(P.Value) %>%
    filter(!grepl("^chr|^LOC_", GeneName)) %>%
    head(topN)
  
  # Build volcano plot
  p <- ggplot(df, aes(x=logFC, y=negLog10P, color=sig)) +
    geom_point(alpha=0.6, size=1.6) +
    geom_vline(xintercept=c(-logfc.threshold, logfc.threshold), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr.threshold), linetype="dashed") +
    scale_color_manual(values=c("NS"="grey", "Upregulated"="red", "Downregulated"="blue")) +
    geom_text_repel(data=top.genes, aes(label=GeneName),
                    size=3, color="black", max.overlaps=50, 
                    box.padding=0.4, point.padding=0.3) +
    labs(title=title, x="log2(Fold Change)", y="-log10(P-value)") +
    theme_minimal() +
    theme(legend.title=element_blank())
  
  return(p)
}

# OR

make.volcano.plot <- function(results.df, title,
                              fdr.threshold=0.05, logfc.threshold=1,
                              p.floor=1e-300, topN=10) {
  df <- as.data.frame(results.df)
  df$P.Value.safe <- pmax(df$P.Value, p.floor)
  df$negLog10P <- -log10(df$P.Value.safe)
  
  df$sig <- "NS"
  df$sig[df$adj.P.Val < fdr.threshold & df$logFC >  logfc.threshold] <- "Upregulated"
  df$sig[df$adj.P.Val < fdr.threshold & df$logFC < -logfc.threshold] <- "Downregulated"
  df$sig <- factor(df$sig, levels=c("NS", "Upregulated", "Downregulated"))
  
  top.genes <- df %>%
    arrange(P.Value) %>%
    filter(!grepl("^chr", GeneName, ignore.case=TRUE),
           !grepl("_", GeneName)) %>%
    head(topN)
  
  p <- ggplot(df, aes(x=logFC, y=negLog10P, color=sig)) +
    geom_point(alpha=0.6, size=1.6) +
    geom_vline(xintercept=c(-logfc.threshold, logfc.threshold), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr.threshold), linetype="dashed") +
    scale_color_manual(values=c("NS"="grey", "Upregulated"="red", "Downregulated"="blue")) +
    geom_text_repel(data=top.genes, aes(label=GeneName),
                    size=3, color="black", max.overlaps=50,
                    box.padding=0.4, point.padding=0.3) +
    labs(title=title, x="log2(Fold Change)", y="-log10(P-value)") +
    theme_minimal() +
    theme(legend.title=element_blank())
  
  return(p)
}

# Volcano plot for Preadipocytes (IAV vs mock)
volc.preadip <- make.volcano.plot(results.preadip, "Volcano: Preadipocytes IAV vs Mock")
volc.preadip
ggsave("volcano.preadipocytes.IAV.vs.mock.pdf", volc.preadip, width=8, height=6)
ggsave("volcano.preadipocytes.IAV.vs.mock.png", volc.preadip, width=8, height=6, dpi=600)

# Volcano plot for Adipocytes (IAV vs mock)
volc.adip <- make.volcano.plot(results.adip, "Volcano: Adipocytes IAV vs Mock")
volc.adip
ggsave("volcano.adipocytes.IAV.vs.mock.pdf", volc.adip, width=8, height=6)
ggsave("volcano.adipocytes.IAV.vs.mock.png", volc.adip, width=8, height=6, dpi=600)

# ----------------------------- Session info -----------------------------

sessionInfo()
