# ============================================================================
# Follow Along Session 5: Over-Representation Analysis (ORA) — 30 minutes
# Advanced R-course 2026
#
# Goal: go from DE results to an interpreted enrichment result.
#
# ============================================================================

if (!require("DT")) install.packages("DT")
if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
if (!require("enrichplot")) BiocManager::install("enrichplot")

# Load necessary libraries
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(DT)

# Read the count matrix
mat <- read.table("./GSE234563_Matrix.txt.gz")
colnames(mat) <- tolower(colnames(mat))
mat <- ceiling(mat)
mat <- mat[rowSums(mat) >= 10, ]

# Prepare the coldata
# coldata <- colnames(mat) |>
#     strsplit("_") |>
#     unlist() |>
#     matrix(ncol = 5, byrow = TRUE) |>
#     as.data.frame() |>
#     setNames(c("age", "sex", "genotype", "tissue", "replicate"))
# rownames(coldata) <- colnames(mat)

coldata <- colnames(mat) |> 
  as.data.frame()|> 
  rename("name"=1) |>
  dplyr::mutate(rowname = name)|>
  tidyr::separate_wider_delim(name,
    "_",
    names = c("age", "sex", "genotype", "tissue", "replicate")) |>
  tibble::column_to_rownames(var="rowname")

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = coldata,
    design = ~ age + sex + tissue + genotype
) |>
    DESeq()

# Define Universe
uni_genes <- results(dds) |>
    as.data.frame() |>
    filter(!is.na(padj)) |>
    rownames()
    
# Define DE genes
de_genes <- results(dds) |>
    as.data.frame() |>
    filter(!is.na(padj), padj <= 0.05, abs(log2FoldChange) > 0) |>
    rownames()

# Run ORA
ego_result <- enrichGO(
    gene = de_genes,
    universe = uni_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
)

# View ORA Result
DT::datatable(head(ego_result))


# Visualizations
# x="Count"/"GeneRatio", color="p.adjust"/"qvalue", size="GeneRatio"/"Count"

dotplot(ego_result, showCategory=10) + 
	ggtitle("GO Enrichment Analysis - DE Genes : Dotplot")

barplot(ego_result, showCategory=10) +
	ggtitle("GO Enrichment Analysis - DE Genes : Barplot")

cnetplot(ego_result, categorySize="pvalue", foldChange=NULL) +
	ggtitle("GO Enrichment Analysis - DE Genes : Cnetplot")

ego_result <- enrichplot::pairwise_termsim(ego_result)
emapplot(ego_result, showCategory=10) +
	ggtitle("GO Enrichment Analysis - DE Genes : Emapplot")




result_up <- results(dds) |>
	as.data.frame() |>
	filter(!is.na(padj), padj <= 0.05, abs(log2FoldChange) > 0) |>
	select(log2FoldChange) |>
	tibble::rownames_to_column(var = "gene") |>
	pull(log2FoldChange,gene)

cnetplot(ego_result_up, foldChange = result_up) +
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Cnetplot with log2FC")












# Define DE genes
de_genes_up <- results(dds) |>
    as.data.frame() |>
    filter(!is.na(padj), padj <= 0.05, log2FoldChange > 0) |>
    rownames()
de_genes_dn <- results(dds) |>
    as.data.frame() |>
    filter(!is.na(padj), padj <= 0.05, log2FoldChange < 0) |>
    rownames()

# Run ORA
ego_result_up <- enrichGO(
    gene = de_genes_up,
    universe = uni_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
)

# View ORA Result
DT::datatable(head(ego_result_up))

ego_result_dn <- enrichGO(
    gene = de_genes_dn,
    universe = uni_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
)

DT::datatable(head(ego_result_dn))

# Visualizations
# x="Count"/"GeneRatio", color="p.adjust"/"qvalue", size="GeneRatio"/"Count"

dotplot(ego_result_up, showCategory=10) + 
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Dotplot")
dotplot(ego_result_dn, showCategory=10) + 
	ggtitle("GO Enrichment Analysis - Downregulated Genes : Dotplot")

barplot(ego_result_up, showCategory=10) +
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Barplot")
barplot(ego_result_dn, showCategory=10) +
	ggtitle("GO Enrichment Analysis - Downregulated Genes : Barplot")

cnetplot(ego_result_up, categorySize="pvalue", foldChange=NULL) +
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Cnetplot")
cnetplot(ego_result_dn, categorySize="pvalue", foldChange=NULL) +
	ggtitle("GO Enrichment Analysis - Downregulated Genes : Cnetplot")

ego_result_up <- enrichplot::pairwise_termsim(ego_result_up)
ego_result_dn <- enrichplot::pairwise_termsim(ego_result_dn)
emapplot(ego_result_up, showCategory=10) +
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Emapplot")
emapplot(ego_result_dn, showCategory=10) +
	ggtitle("GO Enrichment Analysis - Downregulated Genes : Emapplot")




result_up <- results(dds) |>
	as.data.frame() |>
	filter(!is.na(padj), padj <= 0.05, log2FoldChange > 0) |>
	select(log2FoldChange) |>
	tibble::rownames_to_column(var = "gene") |>
	pull(log2FoldChange,gene)

cnetplot(ego_result_up, foldChange = result_up) +
	ggtitle("GO Enrichment Analysis - Upregulated Genes : Cnetplot with log2FC")

