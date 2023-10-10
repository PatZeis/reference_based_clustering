require(RaceID)
require(Rcpp)
require(Matrix)

reference_clustering <- function(seuratobject_ref=NULL,seurat_assay="RNA", raceIDobject_ref=NULL, referencedata=NULL, part=NULL, normref=T,querydata,normquery=T, FSelect=F, features=NULL,featuremethod="RaceID",scalefactor=10000,sample_ref=NULL, clustsize=NULL,minexpr=5,CGenes=NULL,var_feat_len=NULL, RF=T, logisrec=F,rftrees=NULL, PCAprojection=F, nPCs=100, clustering="kmedoids", sampkmedoid = NULL ){
  feature_genes <- NULL
  if ( !is.null(seuratobject_ref)) {
    rawdata <- seuratobject_ref[[seurat_assay]]@counts
    ndata <- t(t(rawdata)/colSums(rawdata))*scalefactor
    part <- as.numeric(seuratobject_ref$seurat_clusters)
    names(part) <- colnames(seuratobject_ref)
    features <- seuratobject_ref[[seurat_assay]]@var.features
    if ( FSelect) {
      if (eval(parse(text = paste(c("length(", paste0(as.character(substitute(seuratobject)), "@assays$", seurat_assay, "@var.features"), ")"), collapse = ""))) == 0 || !is.null(var_feat_len) && length(features) < var_feat_len) { 
        
        if (!inherits(x = seuratobject, "dgMatrix")) {
          seuratobject2 <- seuratobject[[seurat_assay]]@counts
        }
        if (!inherits(x = seuratobject2, what = "dgCMatrix")) {
          seuratobject2 <- Matrix(as.matrix(seuratobject2), sparse = T)
        }
        clip.max <- sqrt(x = ncol(x = seuratobject2))
        hvf.info <- data.frame(mean = rowMeans(x = seuratobject2))
        hvf.info$variance <- SparseRowVar2(mat = seuratobject2, mu = hvf.info$mean, 
                                           display_progress = T)
        hvf.info$variance.expected <- 0
        hvf.info$variance.standardized <- 0
        not.const <- hvf.info$variance > 0
        fit <- loess(formula = log10(x = variance) ~ log10(x = mean), 
                     data = hvf.info[not.const, ], span = 0.3)
        hvf.info$variance.expected[not.const] <- 10^fit$fitted
        hvf.info$variance.standardized <- SparseRowVarStd(mat = seuratobject2, 
                                                          mu = hvf.info$mean, sd = sqrt(hvf.info$variance.expected), 
                                                          vmax = clip.max, display_progress = T)
        colnames(x = hvf.info) <- paste0("vst.", colnames(x = hvf.info))
        hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 
                                     0), ]  ### mean expression >0
        hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, 
                                   decreasing = TRUE), , drop = FALSE]
        features <- head(rownames(hvf.info), var_feat_len)
      }
      feature_genes <- features }
  }
  else if ( !is.null(raceIDobject_ref)) {
    rawdata <- raceIDobject_ref@expdata
    ndata <- raceIDobject_ref@ndata * scalefactor
    part <- raceIDobject_ref@cpart
    if ( length(part) == 0) { part <- raceIDobject_ref@cluster$kpart }
    if ( FSelect) { 
      feature_genes <- raceIDobject_ref@cluster$features }
  }
  else{
    if(is.null(referencedata)) { stop("referencedata has to be set if no Seurat or RaceID object")}
    if( normref) {
      rawdata <- referencedata
      ndata <- t(t(rawdata)/colSums(rawdata))*scalefactor
    }
    else { ndata <- referencedata }
    if(is.null(part)) { stop("as no object is provided e.g. Seurat object please set cluster partion ")}
    else { part <- part }
    if (FSelect) {
      if (!featuremethod %in% c("RaceID","Seurat")) {stop("feature selection method should be either \"RaceID\" or \"Seurat\"")}
      else{
        if(featuremethod=="Seurat") {
          seuratobject2 <- referencedata
          if (!inherits(x = seuratobject2, what = "dgCMatrix")) {
            seuratobject2 <- Matrix(as.matrix(seuratobject2), sparse = T)
          }
          clip.max <- sqrt(x = ncol(x = seuratobject2))
          hvf.info <- data.frame(mean = rowMeans(x = seuratobject2))
          hvf.info$variance <- SparseRowVar2(mat = seuratobject2, mu = hvf.info$mean, 
                                             display_progress = T)
          hvf.info$variance.expected <- 0
          hvf.info$variance.standardized <- 0
          not.const <- hvf.info$variance > 0
          fit <- loess(formula = log10(x = variance) ~ log10(x = mean), 
                       data = hvf.info[not.const, ], span = 0.3)
          hvf.info$variance.expected[not.const] <- 10^fit$fitted
          hvf.info$variance.standardized <- SparseRowVarStd(mat = seuratobject2, 
                                                            mu = hvf.info$mean, sd = sqrt(hvf.info$variance.expected), 
                                                            vmax = clip.max, display_progress = T)
          colnames(x = hvf.info) <- paste0("vst.", colnames(x = hvf.info))
          hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 
                                       0), ]  ### mean expression >0
          hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, 
                                     decreasing = TRUE), , drop = FALSE]
          features <- head(rownames(hvf.info), var_feat_len)
          feature_genes <- features
        }
        if (featuremethod == "RaceID") {
          g <- apply(rawdata, 1, max , na.rm=T) >= minexpr
          if (!is.null(CGenes)) {
            CGenes <- CGenes[CGenes %in% g]
            h <- NULL
            if (length(CGenes) == 1) {
              k <- cor(as.matrix(t(object@ndata[g, ])), 
                       as.matrix(object@ndata[CGenes, ]))
            }
            else {
              k <- cor(as.matrix(t(object@ndata[g, ])), 
                       as.matrix(t(object@ndata[CGenes, ])))
              
            }
            h <- apply(abs(k), 1, max, na.rm = TRUE) < ccor
            h[is.na(h)] <- TRUE
            if (!is.null(h)) 
              genes <- g[h]
          }
          else {
            genes <- g
          }
          bg <- fitbackground(as.matrix(ndata * min(apply(rawdata, 2, sum)))[genes, ] + 0.1)
          feature_genes <- bg$n
        }
      }
    }
    
  }
  cat(paste("reference_build", "\n"))
  if (is.null(feature_genes)) {
    g <- apply(rawdata, 1, max , na.rm=T) >= minexpr
    intersect_rownames <- intersect(rownames(rawdata)[g], rownames(querydata))
    cat(paste("feature selection done", "\n"))
  }
  else { intersect_rownames <- intersect(feature_genes, rownames(querydata))}
  
  if(normquery) {
    if (!inherits(x = querydata, "matrix")) {
      querydata <- as.matrix(querydata)
    }
    normalized_query <- t(t(querydata)/colSums(querydata)) * 10000
  }
  else { normalized_query <- as.matrix(querydata)}
  cat(paste("query build", "\n"))
  if ( !is.null(clustsize)) {uni_part <- as.numeric(names(table(part)[table(part) >= clustsize]))}
  else { uni_part <- sort(unique(part))}
  
  if ( !is.null(sample_ref)) {
    sample <- c()
    for ( i in 1:length(uni_part)){
      sample <- c(sample, sample(names(part[part == uni_part[i]]), round(sample_ref*sum(part == uni_part[[i]])) ))
    }
    cat(paste("sampling done", "\n"))
  }
  else { sample <- names(part)}
  
  
  training <- t(as.matrix(ndata[intersect_rownames,colnames(ndata) %in% sample]))
  xtest <- t(normalized_query[intersect_rownames,])
  if (PCAprojection) {
    if (is.null(nPCs) || nPCs < 10) { stop("# PCs should be set and at least 10")} 
    training <- prcomp(log10(training + 1), center = F, scale. = F, rank. = nPCs, retx = T)
    xtest <- log10(xtest+1) %*% training$rotation
    training <- training$x
    cat(paste("PCA projections done", "\n"))
  }
  pr <- as.factor(part[names(part) %in% sample])
  if (RF ) {
    if (is.null(rftrees)) {stop("set rftrees to finite number")}
    else {
      rf <- randomForest::randomForest(training, pr, xtest, ntree = rftrees, norm.votes = F, importance = TRUE)
    }
    predictions <- rf$test$predicted
    votes <- rf$test$votes
    cat(paste("random forest classification done", "\n"))
  } 
  else if (logisrec) {
    require(nnet)
    if (FSelect==F) {stop("please set FSelect if running classification based on logistic regression")}
    xtraining <- cbind(xtraining, part=part[names(part) %in% sample])
    model <- nnet::multinom(part ~., data = data.frame(xtraining), MaxNWts=24550)
    xtest <- cbind(xtest, part=rep(1, nrow(xtest)))
    colnames(xtest) <- sub("-", ".", colnames(xtest))
    colnames(xtest) <- sub("-", ".", colnames(xtest))
    colnames(xtest)[which(grepl("^\\d", colnames(xtest)))] <- paste("X", colnames(xtest)[which(grepl("^\\d", colnames(xtest)))], sep = "")
    predictions <- model %>% predict(xtest)
    votes <- model %>% predict(xtest, type = "prob")
    cat(paste("logistic regression based classification done", "\n"))
  }
  
  
  ### clustering 
  if (!clustering %in% c("kmedoids","louvain")) {stop("clustering method should be either \"kmedoids\" or \"louvain\"")}
  else { 
    if ( clustering == "kmedoids") {
      distances <- RaceID:::dist.gen(votes, method = "pearson")
      clustexp <- RaceID:::clustfun(distances, clustnr = 30, bootnr = 50, samp = sampkmedoid, sat = T, cln = NULL, rseed=17000,FUNcluster="kmedoids") 
      partition <- clustexp$clb$result$partition
      umap.pars = umap::umap.defaults
      umap.pars$input <- "dist"
      umap <- as.data.frame(umap(distances, config = umap.pars)$layout)
      cat(paste("kmedoids clustering and embedding done", "\n"))
    }
    else if ( clustering == "louvain" ) {
      require(RANN)
      require(igraph)
      knn.info = RANN::nn2(votes, k=10) 
      knn = knn.info$nn.idx
      adj = matrix(0, nrow(knn), nrow(knn))
      rownames(adj) = colnames(adj) = rownames(votes)
      for (i in seq_len(nrow(votes))) {
        adj[i,rownames(votes)[knn[i,]]] = 1
      }
      g = igraph::graph.adjacency(adj, mode='undirected')
      g = simplify(g)
      km = igraph::cluster_louvain(g)
      partition = km$membership
      names(partition) <- km$names
      distances <- RaceID:::dist.gen(votes, method = "pearson")
      umap.pars = umap::umap.defaults
      umap.pars$input <- "dist"
      umap <- as.data.frame(umap(distances, config = umap.pars)$layout) 
      cat(paste("louvain clustering and embedding done", "\n"))
    }
  }
  medoids <- NULL
  compmedoids <- function (distances=NULL, dimred=NULL, part) 
  {
    if(is.null(distances) && is.null(dimred)) { stop("either distances or dimred has to be set")}
    m <- c()
    for (i in sort(unique(part))) {
      f <- names(part)[part == i]
      if (length(f) == 1) {
        m <- append(m, f)
      }
      else if (!is.null(distances)){
        
        y <- apply(distances[f, f], 2, mean)
        m <- append(m, f[which(y == min(y))[1]])
      }
      else{
        g <- apply(as.matrix(dimred[, part == 
                                      i]) - as.vector(pam(t(dimred[, part == 
                                                                     i]), 1)$medoids), 2, sum) == 0
        m <- append(m, names(part)[part == i][g])
      } 
    }
    m
  }
  
  
  medoids <- compmedoids(distances = distances , part = partition) 
  return(list(partition=partition, umap=umap, votes=votes, predictions=predictions, ndata=ndata, normalized_query=normalized_query, clustsize=clustsize, seuratobject_ref=seuratobject_ref, raceIDobject_ref=raceIDobject_ref, clustsize = clustsize, refpart = part, distances=distances, medoids=medoids )) ### object is better run
  
}
