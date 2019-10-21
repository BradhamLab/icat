if (exists("snakemake")) {
  .libPaths(c(snakemake@params$seurat],
              .libPaths()))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(rjson)
  library(reticulate)
})

# check conda installation
if (exists("snakemake")) {
  if (~is.null(snakemake@params$python)) {
    reticulate::use_condaenv(snakemake@params$python)
  }
}

create_seurat <- function(X, obs) {
  X_data = as.data.frame(t(read.csv(X, header=FALSE)))
  obs_data = read.csv(obs, row.names=1, check.names=FALSE,
                      stringsAsFactors=FALSE)
  # drop columns with NA b/c tibble binding in Seurat breaks
  keep_cols <- which(colSums(is.na(obs_data)) == 0)
  obs_data <- data.frame(obs_data[ , keep_cols],
                         row.names=row.names(obs_data))
  names(obs_data) <- names(keep_cols)
  # force to population
  # obs_data$Population <- as.factor(obs_data$Population)
  row.names(obs_data) <- colnames(X_data)
  data <- Seurat::CreateSeuratObject(counts=X_data, meta.data=obs_data)
  return(data)
}

split_and_preprocess <- function(data, treatment) {
  by_treatment <- Seurat::SplitObject(data, split.by=treatment)
  by_treatment <- lapply(by_treatment, function(x) {
    x <- Seurat::NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method='vst', nfeatures=2000)
    return(x)
  })
}

integrate_cells <- function(by_treatment, k) {
  anchors <- Seurat::FindIntegrationAnchors(object.list=by_treatment, dims=1:20)
  combined <- Seurat::IntegrateData(anchorset=anchors, dims=1:20)
  Seurat::DefaultAssay(combined) <- "integrated"
  combined <- Seurat::ScaleData(combined, verbose=FALSE)
  combined <- Seurat::RunPCA(combined, npcs=50, verbose=FALSE)
  combined <- Seurat::RunUMAP(combined, reduction='pca', dims=1:50)
  combined <- Seurat::FindNeighbors(combined, reduction='pca', dims=1:50,
                                    k.param=k)
  combined <- Seurat::FindClusters(combined)
  return(combined)
}

if (exists('snakemake')) {
  seurat <- create_seurat(snakemake@input[['X']], snakemake@input[['obs']])
  k <- rjson::fromJSON(file=snakemake@input[['json']])$n_neighbors
  by_treatment <- split_and_preprocess(seurat, snakemake@params[['treatment']])
  integrated <- integrate_cells(by_treatment, k)
  write.table(t(as.matrix(integrated@assays$integrated@data)),
              snakemake@output[['X']], sep=',',
              row.names=FALSE, col.names=FALSE)
  write.csv(integrated@meta.data, snakemake@output[['obs']])
}