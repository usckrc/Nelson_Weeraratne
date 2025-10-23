---
  title: 'hw #2'
author: "Deshan Weeraratne"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  html_document:
  toc: true
toc_depth: 3
toc_float: true
number_sections: true
theme: bootstrap
df_print: paged
code_folding: hide
highlight: pygments
pdf_document:
  toc: true
toc_depth: '3'
word_document:
  toc: true
toc_depth: '3'
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Load libraries

```{r, echo=TRUE, warning=FALSE, error=FALSE, message=FALSE}
library(Seurat)
library(Matrix)
library(here)
library(dplyr)
library(plotly)   # using this for my bonus question
```

# Question 2a: Load dataset

```{r, echo=TRUE, warning=FALSE, error=FALSE, message=FALSE}
# My dataset location:
data_dir <- here("Week 2 Materials", "Homework Dataset", "Subfolder")

# Loading the Cell Ranger outputs (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
kidney_data <- Read10X(data.dir = here("/Users/deshanweeraratne/Library/CloudStorage/OneDrive-UniversityofSouthernCalifornia/Jonathan Nelson's files - Nelson_Weeraratne/Coding/Week 2 - Rmarkdown and the Whole Kidney/Week 2 Materials/Homework Dataset/Subfolder"))

# Creating a Seurat object
kidney <- CreateSeuratObject(counts = kidney_data, project = "Kidney_snRNAseq")
```

# Question 3i: Filtering parameters and reasoning

```{r, echo=TRUE, warning=FALSE, error=FALSE, message=FALSE}

# Calculate mitochondrial percentage
kidney[["percent.mt"]] <- PercentageFeatureSet(kidney, pattern = "^MT-")

# Check how many cells before filtering
before <- ncol(kidney)

# Choose the thresholds:
# - Keep cells with 200–6000 genes
# - Exclude cells with >10% mitochondrial content
kidney <- subset(kidney, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# After filtering
after <- ncol(kidney)
filtered_out <- before - after
cat("Cells filtered out:", filtered_out, "\n")
```

# Question 3b: Clustering

```{r, echo=TRUE, warning=FALSE, error=FALSE, message=FALSE}
kidney <- NormalizeData(kidney)
kidney <- FindVariableFeatures(kidney)
kidney <- ScaleData(kidney)
kidney <- RunPCA(kidney)

# Choose resolution for clustering diversity
resolution_used <- 0.8
kidney <- FindNeighbors(kidney, dims = 1:20)
kidney <- FindClusters(kidney, resolution = resolution_used)
kidney <- RunUMAP(kidney, dims = 1:20)

cat("Resolution used:", resolution_used, "\n")
cat("Number of clusters found:", length(unique(Idents(kidney))), "\n")
```

# Question 3bii & iii: Name clusters and provide evidence

```{r, echo=TRUE, warning=FALSE, error=FALSE, message=FALSE}
markers <- FindAllMarkers(kidney, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_markers <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
print(top_markers)

# Manually name clusters after checking marker genes
cluster_names <- c(
  "0" = "Proximal Tubule",
  "1" = "Podocyte",
  "2" = "Endothelial",
  "3" = "Fibroblast",
  "4" = "Collecting Duct",
  "5" = "Immune"
)
kidney <- RenameIdents(kidney, cluster_names)
```

# Question 3c: Relevel clusters in biological order

```{r, echo=TRUE, warning=FALSE, error=FALSE, message=FALSE}
bio_order <- c("Proximal Tubule", "Podocyte", "Endothelial", "Fibroblast", "Collecting Duct", "Immune")
Idents(kidney) <- factor(Idents(kidney), levels = bio_order)
```

# Question 3d: Identify DEGs for one cell type (using Podocyte)

```{r, echo=TRUE, warning=FALSE, error=FALSE, message=FALSE}
podocyte_markers <- FindMarkers(kidney, ident.1 = "Podocyte", min.pct = 0.25, logfc.threshold = 0.25)
head(podocyte_markers)
```

# Question 3e: Generate a dot plot to visualize marker expression

```{r, echo=TRUE, warning=FALSE, error=FALSE, message=FALSE}
genes_to_plot <- c("NPHS1", "SLC34A1", "AQP1", "UMOD", "AQP2", "PECAM1", "COL1A1", "PTPRC")
genes_to_plot <- intersect(genes_to_plot, rownames(kidney))

# Create a temporary Seurat object with unique cell names
tmp_kidney <- kidney
colnames(tmp_kidney) <- make.unique(colnames(tmp_kidney))

DotPlot(tmp_kidney, features = genes_to_plot) + 
  RotatedAxis() + 
  ggtitle("DotPlot of Key Kidney Cell Markers")
```

# Question 3f: Save object as .rds file

```{r, echo=TRUE, warning=FALSE, error=FALSE, message=FALSE}
saveRDS(kidney, file = here("Week 2 Materials", "Homework Dataset", "outputs", "kidney_filtered.rds"))
```

# Bonus Question: Interactive UMAP

```{r, echo=TRUE, warning=FALSE, error=FALSE, message=FALSE}
umap_plot <- DimPlot(kidney, reduction = "umap", label = TRUE) +
  ggtitle("Interactive UMAP with Cluster Labels")
ggplotly(umap_plot)
```
