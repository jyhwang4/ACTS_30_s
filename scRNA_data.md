# R Packages
```{r}
library(MuSiC)
library(CARD)
library(Seurat)
library(presto)
library(dbplyr)
library(dplyr)
library(dittoSeq)
library(ConsensusClusterPlus)
library(pheatmap)
library(ComplexHeatmap)
library(ggpubr)
library(NMF)
```

# Data load and QC
```{r}
ACTJY <- readRDS("~/ACTS30_sss/ACTJY.rds")
```


```{r}
aj <- subset(ACTJY, celltype1 %in% c("Endothelial", "C5_Macrophage_CRIP2", "C0_Macrophage_FABP4", "C03_T/NK_ANK3", "C02_T_CD4_TNFRSF25", "C3_Dendritic_CD1C", "C6_Macrophage_PRKN", "C7_Macrophage_Proliferating", "C04_T/NK_FCGR3A", "C1_Monocyte_FCN1", "C01_T_CCR7", "C2_Macrophage_CCL2", "C8_Type1_Dendritic_IDO1", "C00_T_CD8_XCL2", "C4_Myeloid_CCL5", "C9_Macrophage_CCL5", "C08_T/NK_Proliferating", "C1_Fibroblast_CXCR4", "C3_Fibroblast_MGST1", "Epithelial", "Mast", "C05_T_CXCL13", "C06_T_reg", "Plasma", "C2_Fibroblast_MCAM", "C0_Fibroblast_FAP", "C07_T/NK_FCER1G", "C4_Fibroblast_FCER1G", "C5_Fibroblast_CXCL5"))

rm(ACTJY)

aj$celltype1 <- droplevels(aj$celltype1)
 
unique(aj$celltype1)
```
 
```{r}
aj$celltype1 <- NULL
aj$simple <- NULL
aj$simple2 <- NULL
aj$simple3 <- NULL
aj$simple4 <- NULL
aj$RNA_snn_res.0.8 <- NULL
```
 

```{r}
aj <- NormalizeData(aj, normalization.method = "LogNormalize", scale.factor = 10000)
aj <- FindVariableFeatures(aj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(aj)
aj <- ScaleData(aj, features = all.genes)
aj <- RunPCA(aj, features = VariableFeatures(object = aj))
ElbowPlot(aj, reduction = "pca")
```
![image](https://github.com/jyhwang4/ACTS_30_s/assets/59998490/84b0530c-8141-43d8-bf78-3e07612091a6)

```{r}
aj <- FindNeighbors(aj, dims = 1:13)
aj <- FindClusters(aj, resolution = 0.9)
aj <- RunUMAP(aj, reduction = "pca", dims = 1:20)
DimPlot(aj, reduction = "umap", label = T,repel = T, cols = Cols, group.by = "RNA_snn_res.0.9")
```
![image](https://github.com/jyhwang4/ACTS_30_s/assets/59998490/78a12f85-4e20-4cde-a65f-fda8bfdbd04e)

# Top markers
```{r}
markers<- presto::wilcoxauc(aj, 'RNA_snn_res.0.9', assay = 'data')
markers<- top_markers(markers, n = 20, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
markerss <- markers[,-1]
all_markers<- markerss %>%
  #select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
   .[!is.na(.)]
 
markers
```
# Cell annotation
```{r,fig.width=10,fig.height=4}
DotPlot(aj,features = rev(c("CD14", "LYZ","CD3E","CD3D","NKG7","GNLY","MS4A1", "LY9", "COL1A2", "DCN","FLT1", "PECAM1","MZB1", "JCHAIN","KIT","CPA3","EPCAM", "LAG3")),group.by = "RNA_snn_res.0.9") +
     theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()
```
![image](https://github.com/jyhwang4/ACTS_30_s/assets/59998490/3a4d1461-b68a-4e31-93e1-e16b5b65d9bc)
```{r,fig.width=15,fig.height=6}
dittoBarPlot(aj, "predicted.annotation.l2", group.by = "RNA_snn_res.0.9")
```
![image](https://github.com/jyhwang4/ACTS_30_s/assets/59998490/d119d679-d885-4cad-8881-1bbe7795b1b1)

```{r,fig.width=10,fig.height=5.5}
aj <- SetIdent(aj, value=aj$RNA_snn_res.0.9)
levels(aj)
 
new.cluster.ids <- c(
  "T/NK", #0
  "T/NK", #1
  "T/NK", #2
  "T/NK",#3
  "T/NK",#4
  "Myeloid",#5
  "Myeloid",#6
  "Myeloid",#7
  "T/NK",#8
  "Fibroblast",#9
  "Myeloid",#10
  "Myeloid",#11
  "Epithelial",#12
  "Epithelial",#13
  "Epithelial",#14
  "Mast",#15
  "Endothelial",#16
  "Endothelial",#17
  "Fibroblast",#18
  "Myeloid",#19
  "Epithelial",#20
  "T/NK",#21
  "Epithelial",#22
  "T/NK",#23
  "Fibroblast",#24
  "Epithelial",#25
  "Plasma",#26
  "Myeloid",#27
  "Epithelial",#28
  "Epithelial",#29
  "Fibroblast",#30
  "Fibroblast",#31
  "Endothelial",#32
  "Myeloid",#33
  "Fibroblast",#34
  "Plasma",#35
  "Mast",#36
  "Fibroblast",#37
  "Epithelial",#38
  "Epithelial"#39
)
 
 
names(new.cluster.ids) <- levels(aj)
aj <- RenameIdents(aj, new.cluster.ids)
aj$simple1 = aj@active.ident
 
DimPlot(aj,label=T, repel = T, group.by = "simple1")
```

![image](https://github.com/jyhwang4/ACTS_30_s/assets/59998490/24a58c61-9847-4b3d-b80f-a58c0a094358)

```{r,fig.width=10,fig.height=5.5}
aj <- SetIdent(aj, value=aj$RNA_snn_res.0.9)
levels(aj)
 
new.cluster.ids <- c(
  "C01_T_exhausted", #0
  "C02_T/NK_KLRB1", #1
  "C03_T/NK_NLRC3", #2
  "C04_T/NK_CCL18",#3
  "C05_T_Reg",#4
  "C01_Macrophage_C1QB",#5
  "C02_Myeloid_CD86",#6
  "C03_Macrophage_APOC1",#7
  "C06_T/NK_KLRD1",#8
  "C01_Fibroblast_COL1A2",#9
  "C04_Myeloid_CD68",#10
  "C05_Myeloid_CD14",#11
  "C01_Epithelial_WFDC2",#12
  "C02_Epithelial_CXCL17",#13
  "C03_Epithelial_EGFR",#14
  "C01_Mast_CPA3",#15
  "C01_Endothelial_VWF",#16
  "C02_Endothelial_EGFL7",#17
  "C02_Fibroblast_SOD3",#18
  "C06_Myeloid_MSR1",#19
  "C04_Epithelial_CAPS",#20
  "C07_T/NK_Proliferation",#21
  "C05_Epithelial_RSPH1",#22
  "C08_T/NK_HELLS",#23
  "C03_Fibroblast_SPARC",#24
  "C06_Epithelial_EMP2",#25
  "C01_Plasma_IRF8",#26
  "C07_Myeloid_SPI1",#27
  "C07_Epithelial_SFTPB",#28
  "C08_Epithelial_CYP3A5",#29
  "C04_Fibroblast_MYL9",#30
  "C05_Fibroblast_IGFBP4",#31
  "C03_Endothelial_SPARCL1",#32
  "C08_Myeloid_CXCL3",#33
  "C06_Fibroblast_COL3A1",#34
  "C02_Plasma_Immunoglobulin",#35
  "C02_Mast_TPSAB1",#36
  "C07_Fibroblast_PRG4",#37
  "C09_Epithelial_CXCL6",#38
  "C10_Epithelial_CHI3L2"#39
)
 
 
names(new.cluster.ids) <- levels(aj)
aj <- RenameIdents(aj, new.cluster.ids)
aj$simple2 = aj@active.ident
 
DimPlot(aj,label=T, repel = T, group.by = "simple2", cols = color_assignment2)
```
![image](https://github.com/jyhwang4/ACTS_30_s/assets/59998490/b9045012-b17d-460d-b05f-e5a51395b294)


```{r}

color_assignment <- c("Myeloid" = "red",
                      "T/NK" = "blue", 
                      "Epithelial" = "mediumseagreen",
                      "Fibroblast" = "darkorange", "Endothelial" = "orangered",
                      "Plasma" = "purple", "Mast" = "violetred1")

color_assignment2 <- c( "C01_T_exhausted" = "#ADD8E6", "C02_T/NK_KLRB1" = "#87CEEB", "C03_T/NK_NLRC3" = "#00BFFF",
                        "C04_T/NK_CCL18" = "#1E90FF", "C05_T_Reg" = "#4169E1", "C06_T/NK_KLRD1" = "#0000FF", 
                        "C07_T/NK_Proliferation" = "#0000CD", "C08_T/NK_HELLS" = "#00008B",
                        
                        "C01_Macrophage_C1QB" = "#FF9999", "C02_Myeloid_CD86" = "#FF8080", "C03_Macrophage_APOC1" = "#FF6666",
                        "C04_Myeloid_CD68" = "#FF4D4D", "C05_Myeloid_CD14" = "#FF3333", "C06_Myeloid_MSR1" = "#FF1A1A",
                        "C07_Myeloid_SPI1" = "#FF0000", "C08_Myeloid_CXCL3" = "#CC0000",
                        
                        "C01_Epithelial_WFDC2" = "#00FF00", "C02_Epithelial_CXCL17" = "#00E600", "C03_Epithelial_EGFR" = "#00CC00", 
                        "C04_Epithelial_CAPS" = "#00B300", "C05_Epithelial_RSPH1" = "#3CB371", "C06_Epithelial_EMP2" = "#009900",
                        "C07_Epithelial_SFTPB" = "#007F00", "C08_Epithelial_CYP3A5" = "#006600", "C09_Epithelial_CXCL6" = "#013301",
                        "C10_Epithelial_CHI3L2" ="#001000",
                        
                        "C01_Fibroblast_COL1A2" = "#ecfaa7", "C02_Fibroblast_SOD3" = "#e2f584", "C03_Fibroblast_SPARC" = "#dff765",
                        "C04_Fibroblast_MYL9" = "#ddfa4b", "C05_Fibroblast_IGFBP4" = "#d8fa2a", "C06_Fibroblast_COL3A1" = "#d2fa05",
                        "C07_Fibroblast_PRG4" = "#ccfa00",
                        
                        "C01_Endothelial_VWF" = "#FA5D23", "C02_Endothelial_EGFL7" = "orangered", "C03_Endothelial_SPARCL1" = "#B03102",
                        
                        "C01_Plasma_IRF8" = "#BE67F5",  "C02_Plasma_Immunoglobulin" = "#A020F0",
                        
                        "C01_Mast_CPA3" = "violetred1", "C02_Mast_TPSAB1" = "violetred2")
```
