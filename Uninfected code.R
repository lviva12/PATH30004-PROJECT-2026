
# files include counts matrix, list of features and list of barcodes for the human sample

lr_host_adult_uninfected.data = Read10X(data.dir = "files");
lr_host_adult_uninfected = CreateSeuratObject(counts = lr_host_adult_uninfected.data);


# Quality control
lr_host_adult_uninfected[["percent.mt"]] = PercentageFeatureSet(lr_host_adult_uninfected, 
                                                                pattern = "^MT-");

lr_host_adult_uninfected = subset(lr_host_adult_uninfected, subset = percent.mt < 20 & nFeature_RNA <7500 & nFeature_RNA > 200);



# Normalise data
lr_host_adult_uninfected = NormalizeData(lr_host_adult_uninfected);



# Find variables
lr_host_adult_uninfected = FindVariableFeatures(lr_host_adult_uninfected, selection.method = "vst", nfeatures = 3000);



# scale data
all.genes = rownames(lr_host_adult_uninfected);
lr_host_adult_uninfected = ScaleData(lr_host_adult_uninfected, features = all.genes);


# PCA analysis
lr_host_adult_uninfected = RunPCA(lr_host_adult_uninfected, features = VariableFeatures(lr_host_adult_uninfected));


# Find neighbours and clusters
lr_host_adult_uninfected = FindNeighbors(lr_host_adult_uninfected, dims = 1:10);
lr_host_adult_uninfected = FindClusters(lr_host_adult_uninfected, resolution = 0.5);


lr_host_adult_uninfected = RunUMAP(lr_host_adult_uninfected, dims = 1:10);


# Find markers
lr_host_adult_uninfected.markers = FindAllMarkers(lr_host_adult_uninfected, only.pos = TRUE);
lr_host_adult_uninfected.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1);


# make clusters
cluster0_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 0);

cluster1_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 1);

cluster2_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 2);

cluster3_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 3);

cluster4_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 4);

cluster5_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 5);

cluster6_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 6);

cluster7_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 7);

cluster8_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 8);

cluster9_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 9);

cluster10_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 10);

cluster11_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 11);

cluster12_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 12);

cluster13_uninfected.markers = FindMarkers(lr_host_adult_uninfected, ident.1 = 13);



# cells in each cluster

c0_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 0)
c1_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 1)
c2_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 2)
c3_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 3)
c4_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 4)
c5_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 5)
c6_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 6)
c7_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 7)
c8_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 8)
c9_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 9)
c10_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 10)
c11_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 11)
c12_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 12)
c13_uninfected_cells = WhichCells(lr_host_adult_uninfected, idents = 13)


# UMAP visualisation
DimPlot(lr_host_adult_uninfected, reduction = "umap", pt.size = 0.6) # as in Figure 3


