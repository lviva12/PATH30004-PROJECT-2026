
#files = counts matrix, list of features and list of barcodes

########################################### Alpha-infected ##########################################

# lr_UK_host_adult_72hpi_merge1.data = Read10X(data.dir = "files");
lr_UK_host_adult_72hpi_merge1 = CreateSeuratObject(counts = lr_UK_host_adult_72hpi_merge1.data);

lr_UK_host_adult_72hpi_merge1[["percent.mt"]] = PercentageFeatureSet(lr_UK_host_adult_72hpi_merge1, pattern = c("^MT-"));
lr_UK_host_adult_72hpi_merge1 = subset(lr_UK_host_adult_72hpi_merge1, subset = percent.mt < 20 & nFeature_RNA < 9000 & nFeature_RNA > 200);


# add QC
lr_UK_host_adult_72hpi_merge1[["percent.mt"]] = PercentageFeatureSet(lr_UK_host_adult_72hpi_merge1, 
                                                                    pattern = "^MT-");
# subset data
lr_UK_host_adult_72hpi_merge1 = subset(lr_UK_host_adult_72hpi_merge1, subset = percent.mt < 20 & nFeature_RNA < 9000 & nFeature_RNA > 200);


# Normalise data
lr_UK_host_adult_72hpi_merge1 = NormalizeData(lr_UK_host_adult_72hpi_merge1);


# Find variables
lr_UK_host_adult_72hpi_merge1 = FindVariableFeatures(lr_UK_host_adult_72hpi_merge1, selection.method = "vst", nfeatures = 3000);
top10 = head(VariableFeatures(lr_UK_host_adult_72hpi_merge1), 10);


# scale data
all.genes = rownames(lr_UK_host_adult_72hpi_merge1);
lr_UK_host_adult_72hpi_merge1 = ScaleData(lr_UK_host_adult_72hpi_merge1, features = all.genes);


# PCA analysis
lr_UK_host_adult_72hpi_merge1 = RunPCA(lr_UK_host_adult_72hpi_merge1, features = VariableFeatures(lr_UK_host_adult_72hpi_merge1));


# Find neighbours and clusters
lr_UK_host_adult_72hpi_merge1 = FindNeighbors(lr_UK_host_adult_72hpi_merge1, dims = 1:10);
lr_UK_host_adult_72hpi_merge1 = FindClusters(lr_UK_host_adult_72hpi_merge1, resolution = 0.5);



new_cluster_ids = c("Suprabasal-infected", "Ciliated-infected", "Suprabasal-infected", "Secretory_infected", "Ciliated-infected", "Secretory_infected", "Basal-infected", "Basal-infected", "Cycling Basal-infected", "Ciliated-infected", "Ionocyte-infected")
names(new_cluster_ids) = levels(lr_UK_host_adult_72hpi_merge1)
lr_UK_host_adult_72hpi_merge1 = RenameIdents(lr_UK_host_adult_72hpi_merge1, new_cluster_ids)



basal_infected = WhichCells(lr_UK_host_adult_72hpi_merge1, idents = "Basal-infected")
suprabasal_infected = WhichCells(lr_UK_host_adult_72hpi_merge1, idents = "Suprabasal-infected")
ionocyte_infected = WhichCells(lr_UK_host_adult_72hpi_merge1, idents = "Ionocyte-infected")
ciliated_infected = WhichCells(lr_UK_host_adult_72hpi_merge1, idents = "Ciliated-infected")
cycling_basal_infected = WhichCells(lr_UK_host_adult_72hpi_merge1, idents = "Cycling Basal-infected")
secretory_infected = WhichCells(lr_UK_host_adult_72hpi_merge1, idents = "Secretory_infected")

lr_UK_host_adult_72hpi_merge1@meta.data$Cluster <- ifelse (rownames(lr_UK_host_adult_72hpi_merge1@meta.data) %in% basal_infected, print("Basal-infected"), print("Unknown"))
lr_UK_host_adult_72hpi_merge1@meta.data$Cluster <- ifelse (rownames(lr_UK_host_adult_72hpi_merge1@meta.data) %in% suprabasal_infected & lr_UK_host_adult_72hpi_merge1@meta.data$Cluster == "Unknown", print("Suprabsal-infected"), lr_UK_host_adult_72hpi_merge1@meta.data$Cluster)
lr_UK_host_adult_72hpi_merge1@meta.data$Cluster <- ifelse (rownames(lr_UK_host_adult_72hpi_merge1@meta.data) %in% ciliated_infected & lr_UK_host_adult_72hpi_merge1@meta.data$Cluster == "Unknown", print("Ciliated-infected"), lr_UK_host_adult_72hpi_merge1@meta.data$Cluster)
lr_UK_host_adult_72hpi_merge1@meta.data$Cluster <- ifelse (rownames(lr_UK_host_adult_72hpi_merge1@meta.data) %in% cycling_basal_infected & lr_UK_host_adult_72hpi_merge1@meta.data$Cluster == "Unknown", print("Cycling basal-infected"), lr_UK_host_adult_72hpi_merge1@meta.data$Cluster)
lr_UK_host_adult_72hpi_merge1@meta.data$Cluster <- ifelse (rownames(lr_UK_host_adult_72hpi_merge1@meta.data) %in% ionocyte_infected & lr_UK_host_adult_72hpi_merge1@meta.data$Cluster == "Unknown", print("Ionocyte-infected"), lr_UK_host_adult_72hpi_merge1@meta.data$Cluster)
lr_UK_host_adult_72hpi_merge1@meta.data$Cluster <- ifelse (rownames(lr_UK_host_adult_72hpi_merge1@meta.data) %in% secretory_infected & lr_UK_host_adult_72hpi_merge1@meta.data$Cluster == "Unknown", print("Secretory_infected"), lr_UK_host_adult_72hpi_merge1@meta.data$Cluster)


# create data for each gene in each cell
ciliated_infected.markers = FindMarkers(lr_UK_host_adult_72hpi_merge1, ident.1 = "Ciliated-infected");
head(ciliated_infected.markers,100)
ciliated_infected.markers = subset(ciliated_infected.markers, subset = avg_log2FC >1)

c_data_total = FetchData(object = lr_UK_host_adult_72hpi_merge1, cells = ciliated_infected, vars = rownames(ciliated_infected.markers))
c_data_expression = as.data.frame(sort(colSums(c_data_total), decreasing = TRUE))

ciliated_infected_data_total = c_data_total
ciliated_infected_data_expression = c_data_expression
top_cil_infected.markers = rownames(ciliated_infected_data_expression)[1:200]



########################################### Uninfected ##########################################


lr_host_adult_uninfected_merge.data = Read10X(data.dir = "files");
lr_host_adult_uninfected_merge = CreateSeuratObject(counts = lr_host_adult_uninfected_merge.data);


# add QC
lr_host_adult_uninfected_merge[["percent.mt"]] = PercentageFeatureSet(lr_host_adult_uninfected_merge, 
                                                                      pattern = "^MT-");

# subset data
lr_host_adult_uninfected_merge = subset(lr_host_adult_uninfected_merge, subset = percent.mt < 20 & nFeature_RNA <7500 & nFeature_RNA > 200);


# Normalise data
lr_host_adult_uninfected_merge = NormalizeData(lr_host_adult_uninfected_merge);


# Find variables
lr_host_adult_uninfected_merge = FindVariableFeatures(lr_host_adult_uninfected_merge, selection.method = "vst", nfeatures = 3000);
#top10 = head(VariableFeatures(lr_host_adult_uninfected_merge), 10);


# scale data
all.genes = rownames(lr_host_adult_uninfected_merge);
lr_host_adult_uninfected_merge = ScaleData(lr_host_adult_uninfected_merge, features = all.genes);


# PCA analysis
lr_host_adult_uninfected_merge = RunPCA(lr_host_adult_uninfected_merge, features = VariableFeatures(lr_host_adult_uninfected_merge));


# Find neighbours and clusters
lr_host_adult_uninfected_merge = FindNeighbors(lr_host_adult_uninfected_merge, dims = 1:10);
lr_host_adult_uninfected_merge = FindClusters(lr_host_adult_uninfected_merge, resolution = 0.5);



lr_host_adult_uninfected_merge = RunUMAP(lr_host_adult_uninfected_merge, dims = 1:10);


# label cells
new_cluster_ids = c("Secretory-uninfected", "Secretory-uninfected", "Ciliated-uninfected", "Secretory-uninfected", "Ciliated-uninfected", "Basal-uninfected", "Suprabasal-uninfected", "Cycling Basal-uninfected", "Ciliated-uninfected", "Ionocyte-uninfected", "Ciliated-uninfected", "Deuterosomal-uninfected")
names(new_cluster_ids) = levels(lr_host_adult_uninfected_merge)
lr_host_adult_uninfected_merge = RenameIdents(lr_host_adult_uninfected_merge, new_cluster_ids)


basal_uninfected = WhichCells(lr_host_adult_uninfected_merge, idents = "Basal-uninfected")
suprabasal_uninfected = WhichCells(lr_host_adult_uninfected_merge, idents = "Suprabasal-uninfected")
ionocyte_uninfected = WhichCells(lr_host_adult_uninfected_merge, idents = "Ionocyte-uninfected")
ciliated_uninfected = WhichCells(lr_host_adult_uninfected_merge, idents = "Ciliated-uninfected")
cycling_basal_uninfected = WhichCells(lr_host_adult_uninfected_merge, idents = "Cycling Basal-uninfected")
secretory_uninfected = WhichCells(lr_host_adult_uninfected_merge, idents = "Secretory-uninfected")
deuterosomal_uninfected = WhichCells(lr_host_adult_uninfected_merge, idents = "Deuterosomal-uninfected")


lr_host_adult_uninfected_merge@meta.data$Cluster <- ifelse (rownames(lr_host_adult_uninfected_merge@meta.data) %in% basal_uninfected, print("Basal-uninfected"), print("Unknown"))
lr_host_adult_uninfected_merge@meta.data$Cluster <- ifelse (rownames(lr_host_adult_uninfected_merge@meta.data) %in% suprabasal_uninfected & lr_host_adult_uninfected_merge@meta.data$Cluster == "Unknown", print("Suprabsal-uninfected"), lr_host_adult_uninfected_merge@meta.data$Cluster)
lr_host_adult_uninfected_merge@meta.data$Cluster <- ifelse (rownames(lr_host_adult_uninfected_merge@meta.data) %in% ciliated_uninfected & lr_host_adult_uninfected_merge@meta.data$Cluster == "Unknown", print("Ciliated-uninfected"), lr_host_adult_uninfected_merge@meta.data$Cluster)
lr_host_adult_uninfected_merge@meta.data$Cluster <- ifelse (rownames(lr_host_adult_uninfected_merge@meta.data) %in% deuterosomal_uninfected & lr_host_adult_uninfected_merge@meta.data$Cluster == "Unknown", print("Deuterosomal-uninfected"), lr_host_adult_uninfected_merge@meta.data$Cluster)
lr_host_adult_uninfected_merge@meta.data$Cluster <- ifelse (rownames(lr_host_adult_uninfected_merge@meta.data) %in% cycling_basal_uninfected & lr_host_adult_uninfected_merge@meta.data$Cluster == "Unknown", print("Cycling basal-uninfected"), lr_host_adult_uninfected_merge@meta.data$Cluster)
lr_host_adult_uninfected_merge@meta.data$Cluster <- ifelse (rownames(lr_host_adult_uninfected_merge@meta.data) %in% ionocyte_uninfected & lr_host_adult_uninfected_merge@meta.data$Cluster == "Unknown", print("Ionocyte-uninfected"), lr_host_adult_uninfected_merge@meta.data$Cluster)
lr_host_adult_uninfected_merge@meta.data$Cluster <- ifelse (rownames(lr_host_adult_uninfected_merge@meta.data) %in% secretory_uninfected & lr_host_adult_uninfected_merge@meta.data$Cluster == "Unknown", print("Secretory_uninfected"), lr_host_adult_uninfected_merge@meta.data$Cluster)



# create gene data for each ciliated uninfected cell

c_data_total = FetchData(object = lr_host_adult_uninfected, vars = c("PIFO", "MORN2", "RSPH1"))
c_data_expression = as.data.frame(sort(colSums(c_data_total), decreasing = TRUE))
c_data_frame = as.data.frame(c_data_total);


ciliated_uninfected.markers = FindMarkers(lr_host_adult_uninfected_merge, ident.1 = "Ciliated-uninfected");
head(ciliated_uninfected.markers,100)
ciliated_uninfected.markers = subset(ciliated_uninfected.markers, subset = avg_log2FC >1)

c_data_total = FetchData(object = lr_host_adult_uninfected_merge, cells = ciliated_uninfected, vars = rownames(ciliated_uninfected.markers))
c_data_expression = as.data.frame(sort(colSums(c_data_total), decreasing = TRUE))

ciliated_uninfected_data_total = c_data_total
ciliated_uninfected_data_expression = c_data_expression
top_cil_uninfected.markers = rownames(ciliated_uninfected_data_expression)[1:200]




# merge Seurat objects together for analysis
ciliated_uninfected_object = subset(lr_host_adult_uninfected_merge, subset = Cluster == "Ciliated-uninfected")
ciliated_infected_object = subset(lr_UK_host_adult_72hpi_merge1, subset = Cluster == "Ciliated-infected")

combined_uninfected_infected = merge(ciliated_infected_object, y = ciliated_uninfected_object, add.cell.ids = c("Infected", "Uninfected"), project = "Project")

VlnPlot(combined_uninfected_infected, "SNTN") # as shown in Fig. 5C

###################################################### Statistical analysis #############################################################

genes = unique(rownames(ciliated_uninfected.markers), rownames(ciliated_infected.markers))


cil_uninfected_df = data.frame("Cell_barcode" = "barcode", "Gene_ID" = "gene", "Gene_exp" = 0, "Sample" = "sample");
for (g in genes){
  for (x in ciliated_uninfected.cells){
    g_exp = ifelse(g %in% rownames(ciliated_uninfected.markers), ciliated_uninfected_data_total[x, g], 0);
    cil_uninfected_df = add_row(cil_uninfected_df, "Cell_barcode" = x, "Gene_ID" = g, "Gene_exp" = g_exp, "Sample" = "Uninfected");
  }
}

cil_infected_df = data.frame("Cell_barcode" = "barcode", "Gene_ID" = "gene", "Gene_exp" = 0, "Sample" = "sample");
for (g in genes){
  for (x in ciliated_infected.cells){
    g_exp = ifelse(g %in% rownames(ciliated_infected.markers), ciliated_infected_data_total[x, g], 0);
    cil_infected_df = add_row(cil_infected_df, "Cell_barcode" = x, "Gene_ID" = g, "Gene_exp" = g_exp, "Sample" = "Infected");
  }
}


df1 = c1_sig_df[-1,]
df1$Gene_exp[df1$Gene_ID==genes[1]]
df1$Gene_ID = factor(df1$Gene_ID)

df4 = c4_sig_df[-1,]
df4$Gene_exp[df4$Gene_ID==genes[1]]
df4$Gene_ID = factor(df4$Gene_ID)


# add gene expression of each cell for each gene
exp1 = sapply(genes, function(genes){df1$Gene_exp[df1$Gene_ID==genes]})
exp4 = sapply(genes, function(genes){df4$Gene_exp[df4$Gene_ID==genes]})


# create table with p-value for each gene
table_gene_pval_cil = data.frame("Gene_ID" = "gene", "p_value" = 0)

num = length(genes)
for(x in 1:num){
  out = wilcox.test(exp1[,x], exp4[,x]);
  table_gene_pval_cil = add_row(table_gene_pval_cil, "Gene_ID" = genes[x], "p_value" = out$p.value)
}






