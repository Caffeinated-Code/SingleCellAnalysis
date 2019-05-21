---
title: "PBMC_10k_v3_Seurat_R Notebook"
output: 
  html_document:
    keep_md: true
    theme: united
---
 
1) Load UMI count matrix obtained from the CellRanger Pipeline,10x Genomics 

```r
library(dplyr)      #Package for Data manipulation              
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(Seurat)     #Package for Single Cell Data Analysis to find cell clusters          
pbmc.data <- Read10X(data.dir = "/Users/swatz_bachoo/Downloads/Sanofi/Seurat/filtered_feature_bc_matrix") #10X PBMCS from Healthy Donor
```

2) Exploring the UMI count sparse matrix in pbmc.data, where,
    Rows = Genes/Features.
    Columns = Cells.

```r
cat("Number of rows/genes =", nrow(pbmc.data),"\n")
```

```
## Number of rows/genes = 33538
```

```r
cat("\nNumber of columns/cells =", ncol(pbmc.data), "\n")
```

```
## 
## Number of columns/cells = 11769
```

```r
cat ("Sample of pbmc.data\n")
```

```
## Sample of pbmc.data
```

```r
pbmc.data[1:10,30:50] #UMI count of first 10 features from cells 30 to 50
```

```
## 10 x 21 sparse Matrix of class "dgCMatrix"
```

```
##    [[ suppressing 21 column names 'AAACGCTCATAGTCAC', 'AAACGCTCATGGATCT', 'AAACGCTGTCAGACGA' ... ]]
```

```
##                                                      
## MIR1302-2HG . . . . . . . . . . . . . . . . . . . . .
## FAM138A     . . . . . . . . . . . . . . . . . . . . .
## OR4F5       . . . . . . . . . . . . . . . . . . . . .
## AL627309.1  . . . . . . . . . . . . . 1 . . . . . . .
## AL627309.3  . . . . . . . . . . . . . . . . . . . . .
## AL627309.2  . . . . . . . . . . . . . . . . . . . . .
## AL627309.4  . . . . . . . . . . . . . . . . . . . . .
## AL732372.1  . . . . . . . . . . . . . . . . . . . . .
## OR4F29      . . . . . . . . . . . . . . . . . . . . .
## AC114498.1  . . . . . . . . . . . . . . . . . . . . .
```


3) Advantages of using a sparse matrix

```r
#dgCMatrix class is a class of sparse numeric matrices in the compressed, sparse, column-oriented format
dense.size <- object.size(as.matrix(pbmc.data))
sparse.size <- object.size(pbmc.data)
memorySavings <- dense.size/sparse.size
cat("Memory Savings (Fold change) when using a sparse matrix than a dense matrix =",memorySavings, "\n")
```

```
## Memory Savings (Fold change) when using a sparse matrix than a dense matrix = 10.49255
```


4) Create Seurat Object "pbmc" with the raw data.
   Selecting cells which have atleast 200 unique genes and genes which are expressed in atleast 3 cells.

```r
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 200)
pbmc
```

```
## An object of class Seurat 
## 20292 features across 11537 samples within 1 assay 
## Active assay: RNA (20292 features)
```


5) Pre-Processing: Number of total molecules and unique genes are stored in meta.data while creating seurat object. Additional QC metrices can be added as columns here.
nCount_RNA: total number of molecules detected within a cell.
nFeature_RNA : number of unique genes detected in each cell.

```r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 5) #New QC metric column "percent.mt" is added now
```

```
##                  orig.ident nCount_RNA nFeature_RNA percent.mt
## AAACCCAAGCGCCCAT    pbmc10k       2204         1087   2.359347
## AAACCCAAGGTTCCGC    pbmc10k      20090         4200   6.590343
## AAACCCACAGAGTTGG    pbmc10k       5884         1836  10.757988
## AAACCCACAGGTATGG    pbmc10k       5530         2216   7.848101
## AAACCCACATAGTCAC    pbmc10k       5106         1615  10.830396
```
   
   
6) Visualization of QC Metrics:
   (a) FeatureScatter: Used to visualize feature-feature relationships. Total number of molecules detected within a cell correlates strongly with unique genes.

```r
library(ggplot2)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.05) + theme(axis.text = element_text(size = 8))
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.05) + theme(axis.text = element_text(size = 8))
CombinePlots(plots = list(plot1, plot2))
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
  #r: Strength and Directionality of a linear relationship
  #nCount_RNA & percent.mt : Weak Negative Linear Relationship
  #nFeature_RNA & percent.mt : Strong Positive Linear Relationship
```
 
 
   (b) Violin PLots: Black dots represent the values for individual cells. Wider sections represent a higher probability of a given value; the skinnier sections represent a lower probability.

```r
vlnplot1 <- VlnPlot(object=pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.03)
```


7) QC: 
    Unique gene count/cell should neither be too high (cell doublets/multiplets) nor too low (empty droplets/low quality cells) - cutoff set at not over than 2,500 nor less than 200.
    Percentage of reads mapping to mt genome should be less than 5%.

```r
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
require(gridExtra)
```

```
## Loading required package: gridExtra
```

```
## 
## Attaching package: 'gridExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     combine
```

```r
vlnplot2 <- VlnPlot(object=pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
grid.arrange(vlnplot1, vlnplot2, nrow = 2)
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


8) Data Normalization:
    (i) Feature expression/cell divided by total expression (To compare gene expression across cells).
    (ii) Scale factor of 10000 is multiplied to avoid negative values after Ln Transformation. Now, sum of all transcript counts across all genes in each cell will be 10,000.
    (iii) LogNormaliztion used is a natural log transform (Inclusion of wide ranges, account for skewness, more symmetric distribution of data, intutive sense of fold changes).

```r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc[["RNA"]]@data
```

9) Feature Selection: 
   Identifying features with high cell-to-cell variation for biological insights.
   FindVariableFeature(): Identifies features that are outliers on a mean variability plot.
                          (a) Comparison across features.
                          (b) Computes variance of standardized values across all cells in a way similar to a z-transformation based on expected varaince from the expected means of gene expression (Local polynomial regression). 
                          (c) Used directly to rank the features. 

```r
#nFeatures: Number of features to select as top variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) 
```

```
## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
## parametric, : pseudoinverse used at -1.98
```

```
## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
## parametric, : neighborhood radius 0.30103
```

```
## Warning in simpleLoess(y, x, w, span, degree = degree, parametric =
## parametric, : reciprocal condition number 6.0529e-15
```

```r
top15 <- head(VariableFeatures(pbmc), 15)
plot3 <- VariableFeaturePlot(pbmc)
plot4 <- LabelPoints(plot = plot3, points = top15, repel = TRUE, xnudge = 0, ynudge = 0)
plot4
```

```
## Warning: Transformation introduced infinite values in continuous x-axis
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-10-1.png)<!-- -->
   
10) Linear Tansformation/Scaling across cells [(X-Mean)/Variance]: Done so that highly-expressed genes do not dominate in the analysis results.
    Standard Normal Distribution:
    Mean expression of each gene across cells = 0;
    Variance in expression of each gene across cells = 1. 

```r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes) #Scaling all features
```

```
## Centering and scaling data matrix
```

```r
#pbmc <- ScaleData(pbmc) Scaling only the highly variable features
#pbmc[["RNA"]]@scale.data
```
 
11) Removing known unwanted variation: Focussing only on genomic DNA variation

```r
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
```

```
## Regressing out percent.mt
```

```
## Centering and scaling data matrix
```

12) Linear Dimensional Reduction: Principal Component Analysis (PCA).
      (a) Supervised Method for determining relevant sources of heterogeneity.
      (b) Each PC essentially representing a metafeature that combines information across a correlated feature set. 
      (c) The top principal components therefore represent a robust compression of the dataset.

```r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

```
## PC_ 1 
## Positive:  RPS14, RPLP1, RPL13, RPS12, RPS6, RPS18, RPS29, RPS27, LTB, NEAT1 
## 	   CD37, S100A4, S100A6, CD3D, TRBC2, IL32, TRAC, FOS, ARL4C, SELL 
## 	   CALR, IL7R, IRF1, SYNE2, DUSP1, HSP90B1, ATP5IF1, GABPB1-AS1, SLFN5, CD7 
## Negative:  TUBB1, PPBP, GP9, TMEM40, GNG11, CAVIN2, ACRBP, PF4, TSC22D1, PTCRA 
## 	   NRGN, CMTM5, RGS18, MPP1, HIST1H2AC, PRKAR2B, CTTN, LMNA, MTURN, MMD 
## 	   SPARC, PDLIM1, MAP3K7CL, C2orf88, TSPAN33, RAB32, TREML1, ESAM, CLDN5, PGRMC1 
## PC_ 2 
## Positive:  CD3D, IL32, TRAC, LTB, IL7R, TRBC2, SYNE2, SLFN5, NOSIP, ARL4C 
## 	   RORA, TRBC1, LAT, TRABD2A, RPS29, AC058791.1, CD7, ITM2A, AQP3, GABPB1-AS1 
## 	   CCR7, LEF1, PBXIP1, MAL, RPS27, PIM2, CAMK4, ATP6V0E2, INTS6, ODF2L 
## Negative:  FCN1, MNDA, CSTA, TYROBP, S100A8, LYZ, CSF3R, AC020656.1, CLEC7A, CTSS 
## 	   S100A12, PSAP, FGL2, S100A9, CST3, GRN, MS4A6A, AIF1, FPR1, VCAN 
## 	   CYBB, CD14, TGFBI, TYMP, CD300E, CEBPD, RHOB, NCF2, FCER1G, FAM198B 
## PC_ 3 
## Positive:  MS4A1, IGHM, CD79A, IGHD, LINC00926, SNX22, CD22, HLA-DQA1, BANK1, RALGPS2 
## 	   TNFRSF13C, SPIB, FAM129C, TCL1A, CD19, SWAP70, BCL11A, BLK, HVCN1, IGKC 
## 	   CD79B, PLPP5, HLA-DMB, FCRL1, HLA-DRB5, AFF3, HLA-DRB1, HLA-DQB1, PCDH9, HLA-DQA2 
## Negative:  ANXA1, IL32, S100A4, FYB1, TRAC, RORA, LGALS1, AAK1, CD3D, S100A9 
## 	   ID2, GZMM, NEAT1, KLRD1, S100A6, SRGN, CD7, NKG7, SYNE2, PRF1 
## 	   KLRF1, TYROBP, CTSW, ARL4C, FOS, GZMA, FCER1G, CTSD, FGFBP2, S1PR5 
## PC_ 4 
## Positive:  LTB, IL7R, MAL, TRABD2A, CD3D, SELL, CCR7, SERINC5, FYB1, CAMK4 
## 	   TRAC, TSHZ2, LEF1, SLC2A3, AQP3, PBXIP1, FHIT, PASK, RGCC, SH3YL1 
## 	   USP10, FTH1, LINC00513, RGS10, VCAN, LPIN2, BEX3, MNDA, NELL2, ARMH1 
## Negative:  GNLY, KLRF1, GZMB, PRF1, NKG7, MYOM2, KLRD1, S1PR5, CTSW, HOPX 
## 	   GZMA, SPON2, FGFBP2, NMUR1, SH2D1B, TRDC, CLIC3, CCL4, IGFBP7, EOMES 
## 	   CST7, GZMH, CD300A, CX3CR1, TBX21, PTGDR, C1orf21, FCGR3A, NCR1, KIR2DL3 
## PC_ 5 
## Positive:  MYL9, MGAT4B, CCL5, CLU, GAS2L1, ITGB5, TREML1, UBE2E3, TMEM158, ITGA2B 
## 	   PF4, PARVB, CAVIN2, NRGN, MARCH2, PTGS1, SPARC, HIST1H2AC, ACRBP, PPBP 
## 	   GNG11, PADI4, MAP3K7CL, PGRMC1, BEX3, IRX3, PLEK, CMTM5, ZCCHC17, YWHAH 
## Negative:  DNAH12, SLC35D3, FNBP1L, PCP2, MYCT1, AL365361.1, ADAM9, FRMD4B, STON2, PTK2 
## 	   FAM214B, DENND2C, RNF122, FKBP1B, MBTD1, PANX1, LRRC8B, SPSB1, FGF2, KDM7A 
## 	   MID1IP1, AC023157.3, FYB1, MCTP1, SLC2A3, RDX, ING2, TOP1, F2R, AC012615.1
```

```r
#print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) #Top 5 
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") #Visualize top genes associated with reduction components
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
DimPlot(pbmc, reduction = "pca") #Orthogonal Variation
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-13-2.png)<!-- -->

```r
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) 
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-13-3.png)<!-- -->

```r
#Visualization of heterogeneity. 
#Cells and genes are sorted by their principal component scores
#Can see Smaller sub populations with strong separation
#Negatuve and Positive Markers
#DimHeatmap(pbmc, dims = 1:2, cells = 500, balanced = TRUE)
```

13) Determining Dimensionality/Significant PCs: 
    JackStrawPlot: QQ plot to identify significant PCs.If quantiles of PCs came from the same distribution, we should see the points forming a line that’s roughly straight. Distribution of p-values for the genes in the significant PCs lie parallelly close to the uniform distribution (null distribution of p-value feature scores from repeatedly permutated data). 
    ElbowPlot: An ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs

```r
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
```

```
## Warning: Removed 24752 rows containing missing values (geom_point).
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
ElbowPlot(pbmc) 
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-14-2.png)<!-- -->

14) Clustering of Cells:
    K-nearest neighbor (KNN) graph approach to find neighbours, i.e, cells with similar expression patterns.
    Expect Grouping/Clustering of similar cells.

```r
pbmc <- FindNeighbors(pbmc, dims = 1:10, reduction = "pca") #Dimensions of reduction to use as input
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
pbmc <- FindClusters(pbmc, resolution = 2.4) 
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 191
## Number of edges: 5123
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.3438
## Number of communities: 9
## Elapsed time: 0 seconds
```

```r
#Higher resolution, more communities found. Howvever, overlaps also increase.
head(Idents(pbmc), 10) #Idents() used to find cluster IDs of the first five cells
```

```
## AAACCCAAGCGCCCAT AAAGAACTCTCATAGG AAAGGTAGTCTACACA AAATGGATCATTGCGA 
##                1                0                8                0 
## AACCATGGTATGGAAT AACCTTTCACTCCTGT AACGAAACACCCTAAA AAGATAGCAAACGAGC 
##                8                7                3                0 
## AAGTACCCAGGCCCTA AAGTACCGTGCTAGCC 
##                7                0 
## Levels: 0 1 2 3 4 5 6 7 8
```
   
15) Non-Linear Dimensional Reduction (tSNE/UMAP): To place similar cells together in low-dimensional space

```r
library(reticulate)
library(umap)
pbmc <- RunUMAP(object = pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, xnudge = 0, ynudge = 0) #Same as UMAPPlot()
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

16) Finding differentially expressed features/Cluster Biomarkers:
    (a) min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells
    (b) logfc.threshold argument requires a feature to be differentially expressed (on average) by some amount between the two groups

```r
library(magrittr)
#Identifies positive and negative markers of a single cluster (ident.1), compared to all other cells
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
```

```
##               p_val avg_logFC pct.1 pct.2    p_val_adj
## GBP2   1.508703e-12 0.9529751  0.80 0.175 3.061460e-08
## DUSP4  4.505219e-12 0.5975768  0.28 0.000 9.141991e-08
## IL32   3.568986e-11 1.5100984  0.96 0.470 7.242185e-07
## CDC25B 1.392675e-10 0.7157164  0.68 0.133 2.826016e-06
## MIAT   2.534340e-10 0.5900922  0.56 0.072 5.142682e-06
```

```r
#Identifies positive and negative markers of a single cluster compared to cells from other specified clusters
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(6, 8), min.pct = 0.25)
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.41384219541493, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.99612975768415, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 3.17378850285354, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.99612975768415, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.15554767065752, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.15554767065752, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.61892953633207, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.29301136085359, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 3.27495576226752, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.99612975768415, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.29301136085359, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.60967795564303, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 2.0797040859857, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.15554767065752, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.99612975768415, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.99612975768415, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.99612975768415, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.60967795564303, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.41384219541493, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.99612975768415, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.15554767065752, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.29301136085359, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.15554767065752, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.80638715172529, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 3.06122243117072, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.15554767065752, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.38651940331625, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.86433210635525, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.09881232867678, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.29301136085359, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.15554767065752, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 1.94616733600306, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.26531403665036, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT =
## 0.820705740110727, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 2.15554767065752, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.57200773372164, :
## cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 0, ACCTGTCGTTCAAACC
## = 0.693297214317822, : cannot compute exact p-value with ties
```

```
## Warning in wilcox.test.default(x = c(ACATCCCTCATTTGCT = 1.99612975768415, :
## cannot compute exact p-value with ties
```

```r
head(cluster5.markers, n = 5)
```

```
##               p_val  avg_logFC pct.1 pct.2 p_val_adj
## MT-CO3 0.0001328759 -0.3948694 0.952 1.000         1
## TMSB4X 0.0005350231  0.2829971 1.000 1.000         1
## COX8A  0.0009331907  0.7274471 0.762 0.379         1
## SLC9A9 0.0010004823  0.5176830 0.333 0.000         1
## EMP3   0.0012279196  0.5831173 0.905 0.621         1
```

```r
#Identifies positive markers for all clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(2))
```

```
## Calculating cluster 0
```

```
## Calculating cluster 1
```

```
## Calculating cluster 2
```

```
## Calculating cluster 3
```

```
## Calculating cluster 4
```

```
## Calculating cluster 5
```

```
## Calculating cluster 6
```

```
## Calculating cluster 7
```

```
## Calculating cluster 8
```

```r
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
```

```
## # A tibble: 18 x 7
## # Groups:   cluster [9]
##       p_val avg_logFC pct.1 pct.2 p_val_adj cluster gene  
##       <dbl>     <dbl> <dbl> <dbl>     <dbl> <fct>   <chr> 
##  1 1.15e-40     5.17  1     0.012  2.32e-36 0       PF4   
##  2 6.31e-38     5.57  1     0.031  1.28e-33 0       PPBP  
##  3 3.57e-11     1.51  0.96  0.47   7.24e- 7 1       IL32  
##  4 7.41e-10     1.24  0.96  0.476  1.50e- 5 1       TRAC  
##  5 1.04e-17     3.30  0.667 0.06   2.12e-13 2       NKG7  
##  6 6.64e-11     3.33  0.5   0.066  1.35e- 6 2       GNLY  
##  7 1.77e-16     3.86  1     0.389  3.58e-12 3       S100A8
##  8 2.92e-14     3.69  0.958 0.527  5.92e-10 3       S100A9
##  9 8.02e- 5     0.762 0.565 0.226  1.00e+ 0 4       ALKBH7
## 10 5.10e- 3     0.881 0.435 0.208  1.00e+ 0 4       CHD1  
## 11 5.85e-10     1.32  1     0.418  1.19e- 5 5       IL7R  
## 12 1.91e- 5     0.827 0.905 0.494  3.88e- 1 5       TRBC2 
## 13 5.96e-10     0.861 0.412 0.029  1.21e- 5 6       POLR1E
## 14 1.85e- 3     0.829 0.294 0.069  1.00e+ 0 6       NELL2 
## 15 9.49e-35     2.76  1     0.023  1.93e-30 7       MS4A1 
## 16 6.18e-27     3.58  1     0.062  1.25e-22 7       IGHM  
## 17 5.77e- 6     0.822 0.75  0.207  1.17e- 1 8       TRADD 
## 18 5.40e- 5     1.12  1     0.581  1.00e+ 0 8       LTB
```

17) Assigning cell identity to clusters: Based on available literature in this case

```r
new.cluster.ids <- c("Platelets","Memory CD4+ T","NK", "CD14+ Mono", "Naive CD4+ T", "CD8+ T", "DC", "B", "FCGR3A+ Mono")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 1.5, repel = TRUE, xnudge = 0, ynudge = 0) + NoLegend()
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

18) Additional Visualizations for Marker Expression:

```r
#Violin Plots: Expression probability distributions across clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

```r
#Violin Plots of raw counts
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-19-2.png)<!-- -->

```r
#Visualizes feature expression on a tSNE or PCA plot
FeaturePlot(pbmc, features = c("PPBP", "GNLY", "FCN1", "DNAH12"))
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-19-3.png)<!-- -->

```r
#Expression heatmap for given cells and features for top 10 markers
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene,label = TRUE, size = 2.5, angle = 25) + NoLegend() 
```

```
## Warning in DoHeatmap(pbmc, features = top10$gene, label = TRUE, size =
## 2.5, : The following features were omitted as they were not found in the
## scale.data slot for the RNA assay: CUTA, NPM1, HDAC1, TRADD, EEF1B2, RPS13,
## SPOCK2, CHD1, PMPCA, ALKBH7
```

![](Seurat_PBMC_10K_V3_files/figure-html/unnamed-chunk-19-4.png)<!-- -->

