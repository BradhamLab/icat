if (exists("snakemake")) {
  # .libPaths(c("/projectnb/bradham/RPackages", .libPaths()))
  .libPaths(c(snakemake@params[['seurat']],
            .libPaths()))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(rjson)
})


create_seurat <- function(X, obs, label_col) {
  X_data = as.data.frame(t(read.csv(X, header=FALSE)))
  obs_data = read.csv(obs, row.names=1, check.names=FALSE,
                      stringsAsFactors=FALSE)
  # drop columns with NA b/c tibble binding in Seurat breaks
  keep_cols <- which(colSums(is.na(obs_data)) == 0)
  obs_data <- data.frame(obs_data[ , keep_cols],
                         row.names=row.names(obs_data))
  names(obs_data) <- names(keep_cols)
  obs_data[label_col] <- as.factor(obs_data[label_col])
  # force to population
  # obs_data$Population <- as.factor(obs_data$Population)
  colnames(X) <- row.names(obs_data)
  data <- Seurat::CreateSeuratObject(raw.data=X_data, meta.data=obs_data)
  return(data)
}

preprocess_data <- function(data) {
  # data <- Seurat::NormalizeData(data)
  data <- Seurat::ScaleData(data)
  # data <- Seurat::FindVariableGenes(object=data, do.plot=FALSE)
  return(data)
}

rename_cells <- function(seurat_obj, prefix) {
  ids <- sapply(1:ncol(seurat_obj@raw.data), function(i) paste0(prefix, as.character(i)))
  return(Seurat::RenameCells(seurat_obj, new.names=ids))
}

cluster_across_treatments <- function(ctrl, prtb, k, treatment='treatment') {
  # rename cells, add treatment column
  ctrl <- rename_cells(ctrl, 'ctrl')
  ctrl@meta.data$treatment <- rep('ctrl', ncol(ctrl@raw.data))
  prtb <- rename_cells(prtb, 'prtb')
  prtb@meta.data$treatment <- rep('prtb', ncol(prtb@raw.data))
  combined <- Seurat::RunCCA(object=ctrl, object2=prtb)
  combined <- Seurat::RunPCA(combined)
  combined <- Seurat::CalcVarExpRatio(object=combined, reduction.type='pca',
                                      grouping.var=treatment)
  kept <- Seurat::SubsetData(combined, subset.name='var.ratio.pca',
                             accept.low=0.5)
  discarded <- Seurat::SubsetData(combined, subset.name='var.ratio.pca',
                                  accept.high=0.5)
  kept <- Seurat::AlignSubspace(kept, reduction.type='cca',
                                grouping.var=treatment, dims.align=1:20)
  print('Clustering cells...')
  kept <- Seurat::FindClusters(kept, reduction.type='cca.aligned', k.param=k,
                               dims.use=1:20)
  metadata <- combined@meta.data
  metadata$cluster <- NULL
  combined@meta.data[row.names(kept@meta.data), 'cluster'] <- sapply(kept@meta.data$res.0.8,
                                                            function(x) paste0('aligned', as.character(x)))
  combined@meta.data[row.names(discarded@meta.data), 'cluster'] <- rep('unknown', nrow(discarded@meta.data))
  combined <- Seurat::RunUMAP(combined)
  combined@meta.data <- cbind(combined@meta.data, combined@dr$umap@cell.embeddings)
  return(combined@meta.data)
}

main <- function(X, obs, fit_json, out_csv, treatment, control, label_col) {
  data <- create_seurat(X, obs, label_col)
  # separate controls and treated
  control_cells <- row.names(data@meta.data[ , treatment] == control)
  treated_cells <- row.names(data@meta.data[ , treatment] != control)
  ctrl <- Seurat::SubsetData(data, cells=control_cells)
  prtb <- Seurat::SubsetData(data, cells=treated_cells)
  ctrl <- preprocess_data(ctrl)
  prtb <- preprocess_data(prtb)
  k <- fromJSON(file=fit_json)$n_neighbors
  clustered <- cluster_across_treatments(ctrl, prtb, k)
  write.csv(clustered, out_csv)
}

if (exists('snakemake')) {
  main(snakemake@input[['X']],
       snakemake@input[['obs']],
       snakemake@input[['json']],
       snakemake@output[['csv']],
       snakemake@params[['treatment']],
       snakemake@params[['controls']],
       snakemake@params[['label']])
}


