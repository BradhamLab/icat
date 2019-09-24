if (exists("snakemake")) {
  .libPaths(c("/projectnb/bradham/RPackages", .libPaths()))
}
suppressPackageStartupMessages({
  library(Seurat)
  library(rjson)
})


create_seurat <- function(X, obs) {
  X_data = as.data.frame(t(read.csv(X, header=FALSE)))
  obs_data = read.csv(obs, row.names=1, check.names=FALSE,
                      stringsAsFactors=FALSE)
  # hack, fix this
  drop <- c('poor_quality', 'traj', 'cell_line_demuxlet',
            'demuxlet_cls')
  if (sum(drop %in% colnames(obs_data) > 0)) {
    obs_data <- obs_data[, !(colnames(obs_data) %in% drop)]
  }
  # force to population
  # obs_data$Population <- as.factor(obs_data$Population)
  row.names(obs_data) <- colnames(X_data)
  data <- Seurat::CreateSeuratObject(raw.data=X_data, meta.data=obs_data)
  data <- Seurat::NormalizeData(data)
  data <- Seurat::ScaleData(data)
  data <- Seurat::FindVariableGenes(object=data, do.plot=FALSE)
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

main <- function(X_ctrl, obs_ctrl, X_prtb, obs_prtb, fit_json, out_csv,
                 treatment='treatment') {
  ctrl <- create_seurat(X_ctrl, obs_ctrl)
  prtb <- create_seurat(X_prtb, obs_prtb)
  k <- fromJSON(file=fit_json)$n_neighbors
  clustered <- cluster_across_treatments(ctrl, prtb, k)
  write.csv(clustered, out_csv)
}

if (exists('snakemake')) {
  # snakemake likely only being run on scc, add location to user pkgs
  print(sessionInfo())
  print(.libPaths())
  main(snakemake@input[['ctrl_X']],
       snakemake@input[['ctrl_obs']],
       snakemake@input[['prtb_X']],
       snakemake@input[['prtb_obs']],
       snakemake@input[['json']],
       snakemake@output[['csv']],
       snakemake@params[['name']])
}


