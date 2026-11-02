
# files include counts matrix, list of features and list of barcodes for the human sample

lr_UK_host_adult_72hpi.data = Read10X(data.dir = "files");
lr_UK_host_adult_72hpi = CreateSeuratObject(counts = lr_UK_host_adult_72hpi.data);


####################################################################################
#################################### ALPHA VIRUS ###################################
####################################################################################

# virus files = counts matrix, list of barcodes and list of features for the virus

lr_alpha_virus_72hpi.data = Read10X(data.dir = "virus files");
lr_alpha_virus_72hpi = CreateSeuratObject(counts = lr_alpha_virus_72hpi.data);


lr_alpha_virus_72hpi@meta.data = subset(lr_alpha_virus_72hpi@meta.data, nCount_RNA >=10,)


infected_barcodes = rownames(lr_alpha_virus_72hpi@meta.data)
infected_barcodes = as.list(infected_barcodes)


lr_UK_host_adult_72hpi@meta.data$infection = ifelse(rownames(lr_UK_host_adult_72hpi@meta.data) %in% infected_barcodes, print("Infected"), print("Uninfected"))


# separating cells based on infection severity
lr_UK_host_adult_72hpi_low_meta.data = subset(lr_alpha_virus_72hpi@meta.data, nCount_RNA >= 10 & nCount_RNA < 100)
lr_UK_host_adult_72hpi_medium_meta.data = subset(lr_alpha_virus_72hpi@meta.data, nCount_RNA >= 100 & nCount_RNA < 1000)
lr_UK_host_adult_72hpi_high_meta.data = subset(lr_alpha_virus_72hpi@meta.data, nCount_RNA >= 1000 & nCount_RNA < 10000)
lr_UK_host_adult_72hpi_very_high_meta.data = subset(lr_alpha_virus_72hpi@meta.data, nCount_RNA >= 10000)

# cells for each infection tier
lr_UK_host_adult_72hpi_low_barcodes = as.list(rownames(lr_UK_host_adult_72hpi_low_meta.data))
lr_UK_host_adult_72hpi_medium_barcodes = as.list(rownames(lr_UK_host_adult_72hpi_medium_meta.data))
lr_UK_host_adult_72hpi_high_barcodes = as.list(rownames(lr_UK_host_adult_72hpi_high_meta.data))
lr_UK_host_adult_72hpi_very_high_barcodes = as.list(rownames(lr_UK_host_adult_72hpi_very_high_meta.data))


lr_UK_host_adult_72hpi@meta.data$Infection_tier <- ifelse (rownames(lr_UK_host_adult_72hpi@meta.data) %in% lr_UK_host_adult_72hpi_low_barcodes, print("Low"), print("Uninfected"))
lr_UK_host_adult_72hpi@meta.data$Infection_tier <- ifelse (rownames(lr_UK_host_adult_72hpi@meta.data) %in% lr_UK_host_adult_72hpi_medium_barcodes & lr_UK_host_adult_72hpi@meta.data$Infection_tier == "Uninfected", print("Medium"), lr_UK_host_adult_72hpi@meta.data$Infection_tier)
lr_UK_host_adult_72hpi@meta.data$Infection_tier <- ifelse (rownames(lr_UK_host_adult_72hpi@meta.data) %in% lr_UK_host_adult_72hpi_high_barcodes & lr_UK_host_adult_72hpi@meta.data$Infection_tier == "Uninfected", print("High"), lr_UK_host_adult_72hpi@meta.data$Infection_tier)
lr_UK_host_adult_72hpi@meta.data$Infection_tier <- ifelse (rownames(lr_UK_host_adult_72hpi@meta.data) %in% lr_UK_host_adult_72hpi_very_high_barcodes & lr_UK_host_adult_72hpi@meta.data$Infection_tier == "Uninfected", print("Very high"), lr_UK_host_adult_72hpi@meta.data$Infection_tier)


print(nrow(lr_UK_host_adult_72hpi@meta.data[lr_UK_host_adult_72hpi@meta.data$Infection_tier == "Low",]))
print(nrow(lr_UK_host_adult_72hpi@meta.data[lr_UK_host_adult_72hpi@meta.data$Infection_tier == "Medium",]))
print(nrow(lr_UK_host_adult_72hpi@meta.data[lr_UK_host_adult_72hpi@meta.data$Infection_tier == "High",]))
print(nrow(lr_UK_host_adult_72hpi@meta.data[lr_UK_host_adult_72hpi@meta.data$Infection_tier == "Very high",]))
print(nrow(lr_UK_host_adult_72hpi@meta.data[lr_UK_host_adult_72hpi@meta.data$Infection_tier == "Uninfected",]))


####################################################################################
#################################### ALPHA CELLS ###################################
####################################################################################


# Quality control
lr_UK_host_adult_72hpi[["percent.mt"]] = PercentageFeatureSet(lr_UK_host_adult_72hpi, 
                                                              pattern = "^MT-");

lr_UK_host_adult_72hpi = subset(lr_UK_host_adult_72hpi, subset = percent.mt < 20 & nFeature_RNA < 9000 & nFeature_RNA > 200);


# Normalise data
lr_UK_host_adult_72hpi = NormalizeData(lr_UK_host_adult_72hpi);


# Find variables
lr_UK_host_adult_72hpi = FindVariableFeatures(lr_UK_host_adult_72hpi, selection.method = "vst", nfeatures = 3000);
top10 = head(VariableFeatures(lr_UK_host_adult_72hpi), 10);


# scale data
all.genes = rownames(lr_UK_host_adult_72hpi);
lr_UK_host_adult_72hpi = ScaleData(lr_UK_host_adult_72hpi, features = all.genes);


# PCA analysis
lr_UK_host_adult_72hpi = RunPCA(lr_UK_host_adult_72hpi, features = VariableFeatures(lr_UK_host_adult_72hpi));


# Find neighbours and clusters
lr_UK_host_adult_72hpi = FindNeighbors(lr_UK_host_adult_72hpi, dims = 1:10);
lr_UK_host_adult_72hpi = FindClusters(lr_UK_host_adult_72hpi, resolution = 0.5);


# Create UMAP
lr_UK_host_adult_72hpi = RunUMAP(lr_UK_host_adult_72hpi, dims = 1:10);


# Find markers
lr_UK_host_adult_72hpi.markers = FindAllMarkers(lr_UK_host_adult_72hpi, only.pos = TRUE);
lr_UK_host_adult_72hpi.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1);


# save cluster markers
c0_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 0);

c1_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 1);

c2_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 2);

c3_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 3);

c4_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 4);

c5_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 5);

c6_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 6);

c7_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 7);

c8_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 8);

c9_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 9);

c10_72hpi_adult.markers = FindMarkers(lr_UK_host_adult_72hpi, ident.1 = 10);


# cell names
c0_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 0)
c1_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 1)
c2_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 2)
c3_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 3)
c4_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 4)
c5_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 5)
c6_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 6)
c7_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 7)
c8_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 8)
c9_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 9)
c10_72hpi_adult_cells = WhichCells(lr_UK_host_adult_72hpi, idents = 10)


# UMAP visualisation
DimPlot(lr_UK_host_adult_72hpi, reduction = "umap", pt.size = 0.6) # as in Figure 3
DimPlot(lr_UK_host_adult_72hpi, reduction = "umap", pt.size = 0.7, group.by = "Infection_tier"); # as in Figure 4

