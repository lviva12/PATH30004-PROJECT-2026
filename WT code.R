
# files include counts matrix, list of features and list of barcodes for the human sample

lr_VICO1_adult_72hpi.data = Read10X(data.dir = "files");
lr_VICO1_adult_72hpi = CreateSeuratObject(counts = lr_VICO1_adult_72hpi.data);


####################################################################################
##################################### WT VIRUS #####################################
####################################################################################

# virus files = counts matrix, list of barcodes and list of features for the virus

lr_VICO1_virus_72hpi.data = Read10X(data.dir = "virus files");
lr_VICO1_virus_72hpi = CreateSeuratObject(counts = lr_VICO1_virus_72hpi.data);


lr_VICO1_virus_72hpi@meta.data = subset(lr_VICO1_virus_72hpi@meta.data, nCount_RNA >=10,)


infected_barcodes = rownames(lr_VICO1_virus_72hpi@meta.data)
infected_barcodes = as.list(infected_barcodes)


lr_VICO1_adult_72hpi@meta.data$infection = ifelse(rownames(lr_VICO1_adult_72hpi@meta.data) %in% infected_barcodes, print("Infected"), print("Uninfected"))


# separating cells based on infection severity
lr_VICO1_adult_72hpi_low_meta.data = subset(lr_VICO1_virus_72hpi@meta.data, nCount_RNA >= 10 & nCount_RNA < 100)
lr_VICO1_adult_72hpi_medium_meta.data = subset(lr_VICO1_virus_72hpi@meta.data, nCount_RNA >= 100 & nCount_RNA < 1000)
lr_VICO1_adult_72hpi_high_meta.data = subset(lr_VICO1_virus_72hpi@meta.data, nCount_RNA >= 1000 & nCount_RNA < 10000)
lr_VICO1_adult_72hpi_very_high_meta.data = subset(lr_VICO1_virus_72hpi@meta.data, nCount_RNA >= 10000)

# cells for each infection tier
lr_VICO1_adult_72hpi_low_barcodes = as.list(rownames(lr_VICO1_adult_72hpi_low_meta.data))
lr_VICO1_adult_72hpi_medium_barcodes = as.list(rownames(lr_VICO1_adult_72hpi_medium_meta.data))
lr_VICO1_adult_72hpi_high_barcodes = as.list(rownames(lr_VICO1_adult_72hpi_high_meta.data))
lr_VICO1_adult_72hpi_very_high_barcodes = as.list(rownames(lr_VICO1_adult_72hpi_very_high_meta.data))


lr_VICO1_adult_72hpi@meta.data$Infection_tier <- ifelse (rownames(lr_VICO1_adult_72hpi@meta.data) %in% lr_VICO1_adult_72hpi_low_barcodes, print("Low"), print("Uninfected"))
lr_VICO1_adult_72hpi@meta.data$Infection_tier <- ifelse (rownames(lr_VICO1_adult_72hpi@meta.data) %in% lr_VICO1_adult_72hpi_medium_barcodes & lr_VICO1_adult_72hpi@meta.data$Infection_tier == "Uninfected", print("Medium"), lr_VICO1_adult_72hpi@meta.data$Infection_tier)
lr_VICO1_adult_72hpi@meta.data$Infection_tier <- ifelse (rownames(lr_VICO1_adult_72hpi@meta.data) %in% lr_VICO1_adult_72hpi_high_barcodes & lr_VICO1_adult_72hpi@meta.data$Infection_tier == "Uninfected", print("High"), lr_VICO1_adult_72hpi@meta.data$Infection_tier)
lr_VICO1_adult_72hpi@meta.data$Infection_tier <- ifelse (rownames(lr_VICO1_adult_72hpi@meta.data) %in% lr_VICO1_adult_72hpi_very_high_barcodes & lr_VICO1_adult_72hpi@meta.data$Infection_tier == "Uninfected", print("Very high"), lr_VICO1_adult_72hpi@meta.data$Infection_tier)


print(nrow(lr_VICO1_adult_72hpi@meta.data[lr_VICO1_adult_72hpi@meta.data$Infection_tier == "Low",]))
print(nrow(lr_VICO1_adult_72hpi@meta.data[lr_VICO1_adult_72hpi@meta.data$Infection_tier == "Medium",]))
print(nrow(lr_VICO1_adult_72hpi@meta.data[lr_VICO1_adult_72hpi@meta.data$Infection_tier == "High",]))
print(nrow(lr_VICO1_adult_72hpi@meta.data[lr_VICO1_adult_72hpi@meta.data$Infection_tier == "Very high",]))
print(nrow(lr_VICO1_adult_72hpi@meta.data[lr_VICO1_adult_72hpi@meta.data$Infection_tier == "Uninfected",]))


####################################################################################
###################################### WT CELLS ####################################
####################################################################################


# Quality control
lr_VICO1_adult_72hpi[["percent.mt"]] = PercentageFeatureSet(lr_VICO1_adult_72hpi, 
                                                              pattern = "^MT-");

lr_VICO1_adult_72hpi = subset(lr_VICO1_adult_72hpi, subset = percent.mt < 20 & nFeature_RNA < 8000 & nFeature_RNA > 200);


# Normalise data
lr_VICO1_adult_72hpi = NormalizeData(lr_VICO1_adult_72hpi);


# Find variables
lr_VICO1_adult_72hpi = FindVariableFeatures(lr_VICO1_adult_72hpi, selection.method = "vst", nfeatures = 3000);
top10 = head(VariableFeatures(lr_VICO1_adult_72hpi), 10);


# scale data
all.genes = rownames(lr_VICO1_adult_72hpi);
lr_VICO1_adult_72hpi = ScaleData(lr_VICO1_adult_72hpi, features = all.genes);


# PCA analysis
lr_VICO1_adult_72hpi = RunPCA(lr_VICO1_adult_72hpi, features = VariableFeatures(lr_VICO1_adult_72hpi));


# Find neighbours and clusters
lr_VICO1_adult_72hpi = FindNeighbors(lr_VICO1_adult_72hpi, dims = 1:10);
lr_VICO1_adult_72hpi = FindClusters(lr_VICO1_adult_72hpi, resolution = 0.5);


# Create UMAP
lr_VICO1_adult_72hpi = RunUMAP(lr_VICO1_adult_72hpi, dims = 1:10);


# Find markers
lr_VICO1_adult_72hpi.markers = FindAllMarkers(lr_VICO1_adult_72hpi, only.pos = TRUE);
lr_VICO1_adult_72hpi.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1);


# save cluster markers
c0_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 0);

c1_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 1);

c2_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 2);

c3_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 3);

c4_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 4);

c5_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 5);

c6_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 6);

c7_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 7);

c8_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 8);

c9_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 9);

c10_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 10);

c11_VIC01_72hpi_adult.markers = FindMarkers(lr_VICO1_adult_72hpi, ident.1 = 11);


# cells in each cluster
c0_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 0)
c1_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 1)
c2_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 2)
c3_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 3)
c4_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 4)
c5_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 5)
c6_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 6)
c7_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 7)
c8_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 8)
c9_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 9)
c10_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 10)
c11_VIC01_72hpi_adult_cells = WhichCells(lr_VICO1_adult_72hpi, idents = 11)


# UMAP visualisation
DimPlot(lr_VICO1_adult_72hpi, reduction = "umap", pt.size = 0.6) # as in Figure 3
DimPlot(lr_VICO1_adult_72hpi, reduction = "umap", pt.size = 0.7, group.by = "Infection_tier"); # as in Figure 4

