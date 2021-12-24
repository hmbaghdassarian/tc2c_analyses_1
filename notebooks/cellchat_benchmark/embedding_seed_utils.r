# Note on reproducibility/seeds in calculating UMAP via netEmbedding. This version (cellchat commit ID b3ccf96664e29e702e0c23d0bb4e4fc566c72034) has some problems with the seed for UMAP. This may be resolved in more current versions, which explicitly allow seeds to be passed, but this may not be the case. The current version of CellChat::runUMAP does not explicitly pass in the random_state or transform_seed parameter to python's umap-learn. Thus these version of the code do do that:


#' Run UMAP2
#' @param data.use input data matrix
#' @param n.neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param n.components The dimension of the space to embed into.
#' @param distance This determines the choice of metric used to measure distance in the input space.
#' @param n.epochs the number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of the input dataset (200 for large datasets, 500 for small).
#' @param learning.rate The initial learning rate for the embedding optimization.
#' @param min.dist This controls how tightly the embedding is allowed compress points together.
#' Larger values ensure embedded points are moreevenly distributed, while smaller values allow the
#' algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#' @param spread he effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.
#' @param set.op.mix.ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets.
#' @param local.connectivity The local connectivity required - i.e. the number of nearest neighbors
#' that should be assumed to be connected at a local level. The higher this value the more connected
#' the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.
#' @param repulsion.strength Weighting applied to negative samples in low dimensional embedding
#' optimization. Values higher than one will result in greater weight being given to negative samples.
#' @param negative.sample.rate The number of negative samples to select per positive sample in the
#' optimization process. Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
#' @param a More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param b More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param seed.use Set a random seed. By default, sets the seed to 42.
#' @param metric.kwds,angular.rp.forest,verbose other parameters used in UMAP
#' @import reticulate
#' @export
#'
runUMAP2 <- function(
  data.use,
  n.neighbors = 30L,
  n.components = 2L,
  distance = "correlation",
  n.epochs = NULL,
  learning.rate = 1.0,
  min.dist = 0.3,
  spread = 1.0,
  set.op.mix.ratio = 1.0,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  a = NULL,
  b = NULL,
  seed.use = 42L,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  verbose = FALSE){
  if (!reticulate::py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn or reticulate::py_install(packages = 'umap-learn')).")
  }
  set.seed(seed.use)
  reticulate::py_set_seed(seed.use)
  umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(n.neighbors),
    n_components = as.integer(n.components),
    metric = distance,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    metric_kwds = metric.kwds,
    angular_rp_forest = angular.rp.forest,
    verbose = verbose, 
    random_state = seed.use,
    transform_seed = seed.use
  )
  Rumap <- umap$fit_transform
  umap_output <- Rumap(t(data.use))

  colnames(umap_output) <- paste0('UMAP', 1:ncol(umap_output))
  rownames(umap_output) <- colnames(data.use)
  return(umap_output)
}

.error_if_no_Seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat installation required for working with Seurat objects")
  }
}

#' Manifold learning of the signaling networks based on their similarity
#'
#' @param object CellChat object
#' @param slot.name the slot name of object that is used to compute centrality measures of signaling networks
#' @param type "functional","structural"
#' @param comparison a numerical vector giving the datasets for comparison. No need to define for a single dataset. Default are all datasets when object is a merged object
#' @param k the number of nearest neighbors in running umap
#' @param pathway.remove a range of the number of patterns
#' @importFrom methods slot
#' @return
#' @export
#'
#' @examples
netEmbedding2 <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, pathway.remove = NULL, k = NULL, 
                        seed.use = 42L) {
  # b3ccf96664e29e702e0c23d0bb4e4fc566c72034 version from 
  if (object@options$mode == "single") {
    comparison <- "single"
    cat("Manifold learning of the signaling networks for a single dataset", '\n')
  } else if (object@options$mode == "merged") {
    if (is.null(comparison)) {
      comparison <- 1:length(unique(object@meta$datasets))
    }
    cat("Manifold learning of the signaling networks for datasets", as.character(comparison), '\n')
  }
  comparison.name <- paste(comparison, collapse = "-")
  Similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
  if (is.null(pathway.remove)) {
    pathway.remove <- rownames(Similarity)[which(colSums(Similarity) == 1)]
  }
  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(rownames(Similarity) %in% pathway.remove)
    Similarity <- Similarity[-pathway.remove.idx, -pathway.remove.idx]
  }
  if (is.null(k)) {
    k <- ceiling(sqrt(dim(Similarity)[1])) + 1
  }
  options(warn = -1)
  #dimension reduction
  Y <- runUMAP2(Similarity, min.dist = 0.3, n.neighbors = k, seed.use = seed.use)
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$dr)) {
    methods::slot(object, slot.name)$similarity[[type]]$dr <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]] <- Y
  return(object)
}