###################################################
### The main function used to identify doublets ###
###################################################
##scrubDoublets(mtx_tranposed, directory=QC_dir, expected_doublet_rate=0.1)
scrubDoublets <- function(exp,
                          directory = NULL,
                          n_neighbors = NULL,
                          doublet_score_threshold = NULL,
                          sim_doublet_ratio = 2.0,
                          expected_doublet_rate = 0.06,
                          stdev_doublet_rate = 0.02,
                          synthetic_doublet_umi_subsampling = 1.0,
                          use_approx_neighbors = TRUE,
                          distance_metric = 'euclidean',
                          min_counts = 2,
                          min_cells = 3,
                          min_gene_variability_pctl = 85,
                          log_transform = FALSE,
                          z_score = TRUE,
                          n_prin_comps = 30,
                          verbose = TRUE){
  
  if (is.null(n_neighbors)) n_neighbors <- round(0.5 * sqrt(ncol(exp)))
  
  if (verbose) message(green("\nPreprocessing..."))
  E_obs <- t(exp) ### E_obs, ncell * ngene
  total_counts_obs <- apply(E_obs, 1, sum)
  
  E_obs_norm <- pipeline_normalize(E_obs, total_counts_obs)
  gene_filter <- pipeline_get_filter(E_obs_norm)
  E_obs <- E_obs[,gene_filter]
  E_obs_norm <- E_obs_norm[,gene_filter]
  
  if (verbose) message(green("Simulating doublets..."))
  simulateDoublets.res <- simulateDoublets(E_obs, total_counts_obs, sim_doublet_ratio, synthetic_doublet_umi_subsampling)
  E_sim <- simulateDoublets.res$E_sim
  total_counts_sim <- simulateDoublets.res$total_counts_sim
  E_obs_norm <- pipeline_normalize(E_obs, total_counts_obs, postnorm_total = 1e6)
  E_sim_norm <- pipeline_normalize(E_sim, total_counts_sim, postnorm_total = 1e6)
  
  if (log_transform) {
    E_obs_norm <- pipeline_log_transform(E_obs_norm)
    E_sim_norm <- pipeline_log_transform(E_sim_norm)
  }
  
  if (z_score) {
    gene_mean <- apply(E_obs_norm, 2, mean)
    gene_std <- apply(E_obs_norm, 2, sd)
    E_obs_norm <- pipeline_zscore(E_obs_norm, gene_mean, gene_std)
    E_sim_norm <- pipeline_zscore(E_sim_norm, gene_mean, gene_std)
  }
  
  pca.res <- pipeline_pca(E_obs_norm, E_sim_norm, n_prin_comps)
  
  if (verbose) message(green("Calculating doublet scores..."))
  doublet_scores <- calculateDoubletScores(pca.res$pca_obs, pca.res$pca_sim, n_neighbors)
  
  if (is.null(doublet_score_threshold)) {
    if (verbose) message(green("Histogram of doublet scores..."))
    predicted_threshold <- histogramDoubletScores(doublet_scores$doublet_scores_obs, doublet_scores$doublet_scores_sim, directory)
    doublet_score_threshold <- predicted_threshold
  }
  
  if (verbose) message(green("Call transcriptomes as doublets..."))
  predicted_doublets <- callDoublets(doublet_scores$doublet_scores_obs, doublet_scores$doublet_scores_sim, expected_doublet_rate, doublet_score_threshold, verbose)
  return(list(scrubDoublets = predicted_doublets, doublet_scores_obs = doublet_scores$doublet_scores_obs, doublet_scores_sim = doublet_scores$doublet_scores_sim))
  
}

##################################################
### manually reset the doublet_score_threshold ###
##################################################

scrubDoublets_resetThreshold <- function(scrubDoublets_res,
                                         doublet_score_threshold = NULL,
                                         verbose = TRUE){
  
  if(is.null(doublet_score_threshold)) message("Please set doublet_score_threshold.")
  else{
    Ld_obs <- scrubDoublets_res$doublet_scores_obs
    Ld_sim <- scrubDoublets_res$doublet_scores_sim
    threshold <- doublet_score_threshold
    
    predicted_doublets <- Ld_obs > threshold
    
    detected_doublet_rate <- sum(Ld_obs > threshold) / length(Ld_obs)
    detectable_doublet_fraction <- sum(Ld_sim>threshold) / length(Ld_sim)
    overall_doublet_rate <- detected_doublet_rate / detectable_doublet_fraction
    
    if (verbose) {
      message(paste0("Set threshold at doublet score = ", round(doublet_score_threshold,2)))
      message(paste0("Detected doublet rate = ", 100*round(detected_doublet_rate,4), "%"))
      message(paste0("Estimated detectable doublet fraction = ", 100*round(detectable_doublet_fraction,4), "%"))
      message("Overall doublet rate:")
      message(paste0("Estimated  = ", 100*round(overall_doublet_rate,4), "%"))
    }
    
    return(list(scrubDoublets = predicted_doublets, doublet_scores_obs = Ld_obs, doublet_scores_sim = Ld_sim))
    
  }
  
}


############################
### pipeline: preprocess ###
############################

pipeline_normalize <- function(E, total_counts, postnorm_total = NULL){
  
  if (is.null(postnorm_total)){
    total_counts_mean <- mean(total_counts)
  } else {
    total_counts_mean <- postnorm_total
  }
  
  ncell <- nrow(E)
  w <- matrix(0, ncell, ncell)
  diag(w) <- total_counts_mean / total_counts
  Enorm <- w %*% E
  
  return(Enorm)
}

pipeline_get_filter <- function(Enorm, min_counts = 3, min_cells = 3, min_gene_variability_pctl = 85){
  
  vscores.res <- get_vscores(Enorm)
  ix2 <- vscores.res$Vscores > 0
  Vscores <- vscores.res$Vscores[ix2]
  gene_ix <- vscores.res$gene_ix[ix2]
  mu_gene <- vscores.res$mu_gene[ix2]
  FF_gene <- vscores.res$FF_gene[ix2]
  min_vscore <- quantile(Vscores, min_gene_variability_pctl/100)
  ix <- (apply(Enorm[,gene_ix]>=min_counts, 2, sum) >= min_cells) & (Vscores >= min_vscore)
  
  return(gene_ix[ix])
}

get_vscores <- function(Enorm, min_mean = 0, nBins = 50, fit_percentile = 0.1, error_wt = 1){
  
  ncell <- nrow(Enorm)
  mu_gene <- apply(Enorm, 2, mean)
  gene_ix <- c(1:ncol(Enorm))[mu_gene > min_mean]
  mu_gene <- mu_gene[gene_ix]
  
  tmp <- Enorm[,gene_ix]
  tmp <- tmp^2
  var_gene <- apply(tmp, 2, mean) - mu_gene^2
  FF_gene <- var_gene / mu_gene
  
  data_x <- log(mu_gene)
  data_y <- log(FF_gene / mu_gene)
  
  tmp <- runningquantile(data_x, data_y, fit_percentile, nBins)
  x <- tmp$xOut[!is.na(tmp$yOut)]
  y <- tmp$yOut[!is.na(tmp$yOut)]
  
  gLog <- function(x0, x1, x2) log(x1 * exp(-x0) + x2)
  tmp <- log(FF_gene[mu_gene>0])
  tmp <- hist(tmp, breaks=seq(min(tmp), max(tmp), l=201), plot=F)
  h <- tmp$counts
  b <- tmp$breaks
  b <- b[-201] + diff(b)/2
  max_ix <- which.max(h)
  c <- max(exp(b[max_ix]), 1)
  errFun <- function(b2) sum(abs(gLog(x, c, b2)-y)^error_wt)
  b0 <- 0.1
  b <- neldermead::fminsearch(errFun, b0)$simplexopt$x[1,]
  a <- c / (1+b) - 1
  
  v_scores <- FF_gene / ((1+a)*(1+b) + b*mu_gene)
  
  return(list(Vscores=v_scores, gene_ix=gene_ix, mu_gene=mu_gene, FF_gene=FF_gene, a=a, b=b))
}

runningquantile <- function(x, y, p, nBins){
  
  ind <- order(x)
  x <- x[ind]
  y <- y[ind]
  
  dx <- (x[length(x)] - x[1]) / nBins
  xOut <- seq(x[1]+dx/2, x[length(x)]-dx/2, len=nBins)
  yOut <- rep(0, length(xOut))
  
  for (i in 1:length(xOut)){
    ind <- (x >= (xOut[i]-dx/2)) & (x < (xOut[i]+dx/2))
    if (sum(ind)>0){
      yOut[i] <- quantile(y[ind], p/100)
    }
    else{
      if (i>1){
        yOut[i] <- yOut[i-1]
      }
      else {
        yOut[i] <- NA
      }
    }
  }
  
  return(list(xOut=xOut, yOut=yOut))
  
}

pipeline_log_transform <- function(E, pseudocount = 1){
  X <- log10(E + pseudocount)
  return(X)
}

pipeline_zscore <- function(E, gene_mean, gene_std){
  
  E <- t(sweep(E, 2, gene_mean))
  
  nrow <- nrow(E)
  w <- matrix(0, nrow, nrow)
  diag(w) <- 1 / gene_std
  X <- w %*% E
  
  return(t(X))
}

pipeline_pca <- function(X_obs, X_sim, n_prin_comps){
  
  pca <- prcomp(X_obs, rank.=n_prin_comps)
  pca_obs <- pca$x
  pca_sim <- scale(X_sim, pca$center, pca$scale) %*% pca$rotation
  
  return(list(pca_obs = pca_obs, pca_sim = pca_sim))
}

#################################################
### Simulate doublets from observed read cout ###
#################################################

simulateDoublets <- function(E_obs,
                             tatal_counts_obs,
                             sim_doublet_ratio = 2.0,
                             synthetic_doublet_umi_subsampling = 1.0){
  
  n_obs <- nrow(E_obs)
  n_sim <- round(n_obs * sim_doublet_ratio)
  
  pair_ix <- matrix(,n_sim,2)
  for(i in 1:n_sim){
    pair_ix[i,] <- sample(1:n_obs,2)
  }
  
  E1 <- E_obs[pair_ix[,1],]
  E2 <- E_obs[pair_ix[,2],]
  tots1 <- tatal_counts_obs[pair_ix[,1]]
  tots2 <- tatal_counts_obs[pair_ix[,2]]
  
  if (synthetic_doublet_umi_subsampling < 1){
    simulateDoublets.tmp <- subsampleCounts(E1 + E2, synthetic_doublet_umi_subsampling, tots1+tots2)
    E_sim <- simulateDoublets.tmp[[1]]
    total_counts_sim <- simulateDoublets.tmp[[2]]
  } else {
    E_sim <- E1 + E2
    total_counts_sim <- tots1 + tots2
  }
  
  return(list(E_sim = E_sim, total_counts_sim = total_counts_sim, pair_ix = pair_ix))
  
}

subsampleCounts <- function(E, rate, original_totals){
  E <- matrix(rbinom(nrow(E) * ncol(E),round(E),rate),nrow(E),ncol(E))
  current_totals <- apply(E, 1, sum)
  unsampled_orig_totals <- original_totals - current_totals
  unsampled_downsamp_totals <- rbinom(length(unsampled_orig_totals), round(unsampled_orig_totals), rate)
  final_downsamp_totals <- current_totals + unsampled_downsamp_totals
  
  return(list(E, final_downsamp_totals))
}

################################
### Calculate doublet scores ###
################################

calculateDoubletScores <- function(pca_obs,
                                   pca_sim,
                                   n_neighbors,
                                   expected_doublet_rate = 0.06,
                                   stdev_doublet_rate = 0.02,
                                   distance_metric = "euclidean"){
  
  n_obs <- nrow(pca_obs)
  n_sim <- nrow(pca_sim)
  manifold <- rbind(pca_obs, pca_sim)
  doub_labels <- c(rep(0, n_obs), rep(1, n_sim))
  
  # Find k_adj nearest neighbors
  k_adj <- round(n_neighbors * (1+n_sim/n_obs))
  #if (distance_metric %in% c("euclidean")) neighbors <- get.knn(manifold, k = k_adj)$nn.index
  if (distance_metric %in% c("euclidean")) neighbors <- RANN::nn2(manifold,k = k_adj)$nn.idx[,-1]
  
  # Calculate doublet score based on ratio of simulated cell neighbors vs. observed cell neighbors
  doub_neigh_mask <- matrix(doub_labels[neighbors] == 1, nrow(neighbors), ncol(neighbors))
  n_sim_neigh <- apply(doub_neigh_mask, 1, sum)
  n_obs_neigh <- k_adj - n_sim_neigh
  
  rho <- expected_doublet_rate
  r <- n_sim / n_obs
  nd <- n_sim_neigh
  ns <- n_obs_neigh
  N <- k_adj
  
  # Bayesian
  q <- (nd+1)/(N+2)
  Ld <- q*rho/r/(1-rho-q*(1-rho-rho/r))
  
  se_q <- sqrt(q*(1-q)/(N+3))
  se_rho <- stdev_doublet_rate
  
  se_Ld <- q*rho/r / (1-rho-q*(1-rho-rho/r))**2 * sqrt((se_q/q*(1-rho))**2 + (se_rho/rho*(1-q))**2)
  
  doublet_scores_obs <- Ld[doub_labels == 0]
  doublet_scores_sim <- Ld[doub_labels == 1]
  doublet_errors_obs <- se_Ld[doub_labels==0]
  doublet_errors_sim <- se_Ld[doub_labels==1]
  
  return(list(doublet_scores_obs = doublet_scores_obs, doublet_scores_sim = doublet_scores_sim, doublet_errors_obs = doublet_errors_obs, doublet_errors_sim = doublet_errors_sim))
  
}

###################################################################
### Plot the histograme for doublet scores and detect threshold ###
###################################################################

histogramDoubletScores <- function(doublet_scores_obs, doublet_scores_sim, directory){
  
  # estimate the threshold based on kmeans cluster
  km <- kmeans(doublet_scores_sim, centers=2)
  clust <- as.factor(km$cluster)
  predicted_threshold <- (max(doublet_scores_sim[clust==1]) + min(doublet_scores_sim[clust==2]))/2
  
  dat_obs <- data.frame(doublet_scores = doublet_scores_obs, clust = rep(1, length(doublet_scores_obs)))
  dat_obs$clust[dat_obs$doublet_scores > predicted_threshold] <- 2
  dat_obs$clust <- factor(dat_obs$clust)
  
  dat_sim <- data.frame(doublet_scores = doublet_scores_sim, clust = rep(1, length(doublet_scores_sim)))
  dat_sim$clust[dat_sim$doublet_scores > predicted_threshold] <- 2
  dat_sim$clust <- factor(dat_sim$clust)
  
  cat(bold(green("\nPlotting histogram of doublet scores \n")))
  p_obs <- ggplot2::ggplot(dat_obs, aes(x = doublet_scores))
  p_obs <- p_obs + geom_histogram(aes(fill = clust), binwidth = 0.02, color = "grey50")
  p_obs <- p_obs + geom_vline(xintercept = predicted_threshold, color = "blue")
  p_obs <- p_obs + labs(x="Doublet scores", y="Counts", title="Observed Cells") + theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
  
  p_sim <- ggplot2::ggplot(dat_sim, aes(x = doublet_scores))
  p_sim <- p_sim + geom_histogram(aes(fill = clust), binwidth = 0.02, color = "grey50")
  p_sim <- p_sim + geom_vline(xintercept = predicted_threshold, color = "blue")
  p_sim <- p_sim + labs(x="Doublet scores", y="Counts", title="Simulated Doublets") + theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
  
  p_obs2 <- ggplot2::ggplot(dat_obs, aes(x = doublet_scores)) + stat_density(geom="line", color="red") + geom_vline(xintercept = predicted_threshold, color = "blue")
  p_obs2 <- p_obs2 + labs(x="Doublet scores", y="Density", title="") + theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
  p_sim2 <- ggplot2::ggplot(dat_sim, aes(x = doublet_scores)) + stat_density(geom="line", color="red") + geom_vline(xintercept = predicted_threshold, color = "blue")
  p_sim2 <- p_sim2 + labs(x="Doublet scores", y="Density", title="") + theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
  #
  pdf(paste0(directory,"/02i_histogram of doublet scores.pdf"),8,8, useDingbats=FALSE)
  gridExtra::grid.arrange(p_obs, p_sim, p_obs2, p_sim2, nrow = 2, ncol = 2)
  dev.off()
  cat(paste0(silver("Plots saved in: ")),bold(silver("02.QC_Plots\\02i_histogram of doublet scores.pdf \n")))
  
  return(predicted_threshold)
  
}

###################################################
### Call transcriptomes as doublets or singlets ###
###################################################

callDoublets <- function(doublet_scores_obs,
                         doublet_scores_sim,
                         expected_doublet_rate,
                         doublet_score_threshold,
                         verbose){
  
  Ld_obs <- doublet_scores_obs
  Ld_sim <- doublet_scores_sim
  threshold <- doublet_score_threshold
  
  predicted_doublets <- Ld_obs > threshold
  
  detected_doublet_rate <- sum(Ld_obs > threshold) / length(Ld_obs)
  detectable_doublet_fraction <- sum(Ld_sim>threshold) / length(Ld_sim)
  overall_doublet_rate <- detected_doublet_rate / detectable_doublet_fraction
  
  if (verbose) {
    message(paste0("Set threshold at doublet score = ", round(doublet_score_threshold,2)))
    message(paste0("Detected doublet rate = ", 100*round(detected_doublet_rate,4), "%"))
    message(paste0("Estimated detectable doublet fraction = ", 100*round(detectable_doublet_fraction,4), "%"))
    message("Overall doublet rate:")
    message(paste0("Expected   = ", 100*round(expected_doublet_rate,4), "%"))
    message(paste0("Estimated  = ", 100*round(overall_doublet_rate,4), "%"))
  }
  
  return(predicted_doublets)
  
}

plotGene <- function(genelist, umi, dir){
  ### Density plot
  suppressWarnings({pdf(file.path(dir, paste0("02d_QC_Hist_Check.pdf")), useDingbats=FALSE)
    cat(bold(green("Plotting QC per gene Histograms \n")))
    for(gene in genelist)
    {
      if(gene %in% row.names(umi@assays$RNA@counts)) {
        expr_gene <- paste0(gene, "_expressed")
        umi@meta.data[, expr_gene] <- ifelse(GetAssayData(object=umi, slot="counts")[gene,]>0, "TRUE", "FALSE")
        plot.title <- "nGene"
        p1 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=get(expr_gene), fill = get(expr_gene))) + ggplot2::geom_density(size=0.5, alpha=0.2) + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
        p1 <- p1 + guides(color=guide_legend("Expressed:"), fill =guide_legend("Expressed:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)
        plot.title <- "nGene zoom"
        p2 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=get(expr_gene), fill = get(expr_gene))) + ggplot2::geom_density(size=0.5, alpha=0.2) + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5)) + ggplot2::xlim(0,2000) + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
        p2 <- p2 + guides(color=guide_legend("Expressed:"), fill =guide_legend("Expressed:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)
        plot.title <- "nUMI"
        p3 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=get(expr_gene), fill = get(expr_gene))) + ggplot2::geom_density(size=0.5, alpha=0.2) + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
        p3 <- p3 + guides(color=guide_legend("Expressed:"), fill =guide_legend("Expressed:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)
        plot.title <- "nUMI zoom"
        p4 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=get(expr_gene), fill = get(expr_gene))) + ggplot2::geom_density(size=0.5, alpha=0.2) + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5)) + ggplot2::xlim(0,5000) + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
        p4 <- p4 + guides(color=guide_legend("Expressed:"), fill =guide_legend("Expressed:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)
        print(patchwork::wrap_plots(p1+p2+p3+p4+plot_layout(guides = 'collect'))
              + plot_annotation(title = gene, theme = theme(plot.title = element_text(hjust = 0.5, face="bold"))))
      }
    }
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("02.QC_Plots\\02d_QC_Hist_Check.pdf \n")))
  ### Scatter Plot
  suppressWarnings({pdf(file.path(dir, paste0("02e_QC_Scatter_Check.pdf")), useDingbats=FALSE)
    cat(bold(green("Plotting QC per gene Scatter plots \n")))
    for(gene in genelist)
    {
      if(gene %in% row.names(umi@assays$RNA@counts)) {
        expr_gene <- paste0(gene, "_expressed")
        umi@meta.data[, expr_gene] <- ifelse(GetAssayData(object=umi, slot="counts")[gene,]>0, "TRUE", "FALSE")
        #genes vs mt%
        vars1 <- table(umi@meta.data[, expr_gene]=="TRUE")["TRUE"][[1]]
        print(FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_mt", pt.size=0.5, group.by=expr_gene) + ggplot2::ggtitle(paste0(as.character(gene),"\n(Expressed in ", vars1, " cells)"))+ theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + guides(colour = guide_legend("Expressed:"))) + theme(plot.title = element_text(hjust = 0.5))
        print(FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_mt", pt.size=0.5, group.by=expr_gene) + ggplot2::xlim(0,2000)+ ggplot2::ggtitle(as.character(gene)) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + guides(colour = guide_legend("Expressed:"))) + theme(plot.title = element_text(hjust = 0.5))
        #print(vars1)
        print(FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_mt", pt.size=0.5, cells=colnames(umi)[umi@meta.data[, expr_gene]=="TRUE"], col="#00BFC4") + ggplot2::labs(colour="Sample ID:") + ggplot2::ylim(0,100) + ggplot2::ggtitle(as.character(gene)) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)) + theme(plot.title = element_text(hjust = 0.5))
        print(FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_mt", pt.size=0.5, cells=colnames(umi)[umi@meta.data[, expr_gene]=="TRUE"], col="#00BFC4") + ggplot2::labs(colour="Sample ID:") + ggplot2::xlim(0,2000) + ggplot2::ylim(0,100) + ggplot2::ggtitle(as.character(gene)) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)) + theme(plot.title = element_text(hjust = 0.5))
        # UMI vs mt%
        print(FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_mt", pt.size=0.5, group.by=expr_gene) + ggplot2::xlim(0,5000) + ggplot2::ggtitle(as.character(gene)) + guides(colour = guide_legend("Expressed:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)) + theme(plot.title = element_text(hjust = 0.5))
        print(FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_mt", pt.size=0.5, cells=colnames(umi)[umi@meta.data[, expr_gene]=="TRUE"], col="#00BFC4") + ggplot2::labs(colour="Sample ID:") + ggplot2::xlim(0,5000) + ggplot2::ylim(0,100) + ggplot2::ggtitle(as.character(gene)) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)) + theme(plot.title = element_text(hjust = 0.5))
        # genes vs CD3 expression
        print(FeatureScatter(umi, feature1="nFeature_RNA", feature2=gene, pt.size=0.5, group.by=expr_gene) + ggplot2::ylab("Counts") + ggplot2::ggtitle(as.character(gene)) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + guides(colour = guide_legend("Expressed:"))) + theme(plot.title = element_text(hjust = 0.5))
        print(FeatureScatter(umi, feature1="nFeature_RNA", feature2=gene, pt.size=0.5, group.by=expr_gene) + ggplot2::ylab("Counts") + ggplot2::xlim(0,2000)+ ggplot2::ggtitle(as.character(gene)) + guides(colour = guide_legend("Expressed:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)) + theme(plot.title = element_text(hjust = 0.5))
      }
    }
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("02.QC_Plots\\02e_QC_Scatter_Check.pdf \n")))
  cat(paste0(cyan("\nNow check the graphs, choose your thresholds and then run")),bold(cyan("FilterPlots \n")))
  
}

algo_plot <- function(dir, data, type, plot_name, res){
  pdf(file.path(paste0(dir, plot_name, "_", type,".pdf")), width=10, height=10, useDingbats=FALSE)
  for (res.i in res){
  print(DimPlot(data, reduction=type, group.by=paste0("RNA_snn_res.",res.i), label=T, pt.size=1, repel=TRUE) +
          ggplot2::xlab(paste0(toupper(type), " 1")) +
          ggplot2::ylab(paste0(toupper(type), " 2")) +
          ggplot2::ggtitle(paste0("Clusters @ res ", res.i)) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }  
  invisible(dev.off())
  
  pdf(file.path(paste0(dir, plot_name, "_", type, "_DimReduction.pdf")), width=14, height=12, useDingbats=FALSE)
  popsicler:::four_plots(data, type)
  invisible(dev.off())
  
  
}

SR_plots <- function(db_name, annot_db, data, directory, cluster_res) { ### aggiungere BBpar per settare diversi core
  sc_anal <- paste0(db_name,".sc.main.labels")
  var.sc <- SingleR(test=data@assays$RNA@data, ref=annot_db, labels=annot_db$label.main, method="single")
  data[[sc_anal]] <- var.sc$labels
  
  pdf(paste0(directory, "/05a_UMAP_", sc_anal, ".pdf"), width=12, height=10, useDingbats=FALSE)
  print(UMAPPlot(data, label=T, pt.size=1, group.by=sc_anal, repel=TRUE) +
          ggplot2::ggtitle(paste0(db_name," single cell annotation")) +
          ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
          ggplot2::guides(col=ggplot2::guide_legend(ncol=1))+ xlab("UMAP 1") + ylab("UMAP 2"))
  invisible(dev.off())
  
  pdf(paste0(directory, "/05a_tSNE_", sc_anal, ".pdf"), width=12, height=10, useDingbats=FALSE)
  print(DimPlot(data, reduction="tsne", group.by=sc_anal, label=T, pt.size=1, repel=T) +
          ggplot2::ggtitle(paste0(db_name," single cell annotation")) +
          ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
          ggplot2::guides(col=ggplot2::guide_legend(ncol=1))+xlab("tSNE 1") + ylab("tSNE 2"))
  invisible(dev.off())  
  
  
  if(is.null(cluster_res)) {
    selected_clusters<-"seurat_clusters"  
    cl_anal <- paste0(db_name,".cl.main.labels")
    var.cl <- SingleR(test=data@assays$RNA@data, ref=annot_db, labels=annot_db$label.main, method="cluster", clusters=data@meta.data[[sel_cluster]])
    data[[cl_anal]] <- paste0(data[[]][[sel_cluster]], ":", var.cl$labels[match(data[[]][[sel_cluster]], rownames(var.cl))])
    data@meta.data[[cl_anal]]<-factor(data@meta.data[[cl_anal]], levels=mixedsort(unique(data@meta.data[[cl_anal]])))

    pdf(paste0(directory, "/05b_UMAP_", paste0(db_name,".cl.main.labels"), ".pdf"), width=12, height=10, useDingbats=FALSE)  

    print(UMAPPlot(data, label=T, pt.size=1, group.by=cl_anal, repel=T) +
              NoLegend() +
              ggplot2::ggtitle(paste0(db_name," cluster annotation")) +
              ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+ xlab("UMAP 1") + ylab("UMAP 2"))

    invisible(dev.off()) 
    
    pdf(paste0(directory, "/05b_tSNE_", paste0(db_name,".cl.main.labels"), ".pdf"), width=12, height=10, useDingbats=FALSE)
    print(DimPlot(data, reduction="tsne", group.by=cl_anal, label=T, pt.size=1, repel=T) +
              NoLegend() +
              ggplot2::ggtitle(paste0(db_name," cluster annotation")) +
              ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+xlab("tSNE 1") + ylab("tSNE 2"))
    invisible(dev.off()) 
    
  } else {
  cl_analyses<-c()
  for (res.i in cluster_res) {
    sel_cluster<-paste0("RNA_snn_res.",res.i)
    var.cl <- SingleR(test=data@assays$RNA@data, ref=annot_db, labels=annot_db$label.main, method="cluster", clusters=data@meta.data[[sel_cluster]])
    cl_anal <- paste0(db_name,".cl.main.labels.res",res.i)
    data[[cl_anal]] <- paste0(data[[]][[sel_cluster]], ":", var.cl$labels[match(data[[]][[sel_cluster]], rownames(var.cl))])
    data@meta.data[[cl_anal]]<-factor(data@meta.data[[cl_anal]], levels=mixedsort(unique(data@meta.data[[cl_anal]])))
    cl_analyses<-c(cl_analyses,cl_anal)
  } #end for
 
  pdf(paste0(directory, "/05b_UMAP_", paste0(db_name,".cl.main.labels"), ".pdf"), width=12, height=10, useDingbats=FALSE)  
  for (cl_anal in cl_analyses) {
  print(UMAPPlot(data, label=T, pt.size=1, group.by=cl_anal, repel=T) +
          NoLegend() +
          ggplot2::ggtitle(paste0(db_name," cluster ",gsub(paste0(db_name,".cl.main.labels."),"",cl_anal)," annotation")) +
          ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+ xlab("UMAP 1") + ylab("UMAP 2"))
  }
  invisible(dev.off()) 
  
  pdf(paste0(directory, "/05b_tSNE_", paste0(db_name,".cl.main.labels"), ".pdf"), width=12, height=10, useDingbats=FALSE)
  for (cl_anal in cl_analyses) {  
  print(DimPlot(data, reduction="tsne", group.by=cl_anal, label=T, pt.size=1, repel=T) +
          NoLegend() +
          ggplot2::ggtitle(paste0(db_name," cluster ",gsub(paste0(db_name,".cl.main.labels."),"",cl_anal)," annotation")) +
          ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+xlab("tSNE 1") + ylab("tSNE 2"))
  }
  invisible(dev.off())
  } #end else  
  
  return(data)
}

firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

FTP <- function(data, directory, graph_value, dimensional_redux, H, to_be_plotted, filename){
  pdf(paste0(directory, "/", graph_value, dimensional_redux, filename, ".pdf"), width=18, height=H, useDingbats=FALSE)
  print(FeaturePlot(data, features = to_be_plotted, reduction=dimensional_redux) +
          ggplot2::xlab(paste0(toupper(dimensional_redux), " 1")) +
          ggplot2::ylab(paste0(toupper(dimensional_redux), " 2")))
  invisible(dev.off())
}

VLN <- function(data, H, feats, colours=NULL, point=FALSE){
  VlnPlot(data, features=feats, cols=colours, pt.size=point) + ggplot2::geom_boxplot(width=0.1, outlier.shape=NA)
}

DTP <- function(data, markers, annotation){
  intercepts<-(cumsum(sapply(markers,length))+0.5)[-length(markers)]
  mrk_x <- c(1,(cumsum(sapply(markers,length))+1)[-length(markers)])-0.3
  mrk_y <- length(table(data@meta.data[,annotation]))+1
  markers.to.plot<-as.character(unlist(markers))
  #head.angle <- 0
  if(length(markers.to.plot)>10)
  {
    head.angle <- 45
  }else
  {
    head.angle <- 0
  }
  if (compareVersion(as.character(packageVersion("Seurat")), '3.2.0')==-1)
  {
    markers.to.plot <- rev(x=markers.to.plot)
  }
  return(DotPlot(object=data, features=markers.to.plot, cols=c("white", "red"),
                dot.scale=8, assay="RNA", group.by=annotation) + RotatedAxis() +
          ggplot2::geom_vline(xintercept=intercepts, linetype="dashed", color="darkgrey") +
          ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.5,3.5))) + # gglot2 3.3
          ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.5,2))) + # gglot2 3.3
          RotatedAxis() + ggplot2::xlab("marker") + ggplot2::ylab(annotation) +
          ggplot2::annotate("text", x = mrk_x, y = mrk_y, label=names(markers), angle=head.angle, hjust = 0))
}

four_plots <- function(data, redux) {
  PreA <- FeaturePlot(data, pt.size=1, features="nCount_RNA", reduction=redux) + ggplot2::xlab(paste0(toupper(redux), " 1")) + ggplot2::ylab(paste0(toupper(redux), " 2"))
  PreB <- FeaturePlot(data, pt.size=1, features="nFeature_RNA", reduction=redux) + ggplot2::xlab(paste0(toupper(redux), " 1")) + ggplot2::ylab(paste0(toupper(redux), ' 2'))
  PreC <- FeaturePlot(data, pt.size=1, features="percent_mt", reduction=redux) + ggplot2::xlab(paste0(toupper(redux), " 1")) + ggplot2::ylab(paste0(toupper(redux), " 2"))
  PreD <- DimPlot(data, pt.size=1, group.by="Phase", reduction=redux) + ggplot2::xlab(paste0(toupper(redux), " 1")) + ggplot2::ylab(paste0(toupper(redux), " 2")) +
    ggplot2::ggtitle("Phase") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(patchwork::wrap_plots(PreA, PreB, PreC, PreD, ncol=2))
}

annotation_plot <- function(directory, graph_value, data, redux, annotation, filename){
  pdf(paste0(directory, "/", graph_value, redux, "_", filename, "_", annotation,".pdf"), useDingbats=FALSE)
  cell_type <- names(table(data[[annotation]]))
  for (l in 1:length(cell_type)){
    print(DimPlot(data, reduction=tolower(redux), label=F, pt.size=0.5, group.by=annotation,
                  cells.highlight=row.names(data[[annotation]])[data[[annotation]][,annotation] %in% cell_type[l]],
                  cols.highlight="red", sizes.highlight=0.7, repel=T) + ggplot2::ggtitle(cell_type[l]) +
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) + NoLegend()) +
      ggplot2::xlab(paste0(redux, " 1")) + ggplot2::ylab(paste0(redux, " 2"))
  }
  invisible(dev.off())
}

PrePlots <- function(sample, input_data, genelist, percentage=0.1, gene_filter=200, cellranger=TRUE, input_matrix=NULL, organism="human", out_folder=getwd()){
  ### Root directory
  root <- out_folder
  if (!file.exists(root)){dir.create(root, recursive=T)}
  ### QC Plots directory
  QC_dir <- paste0(root,"/02.QC_Plots/")
  if (!file.exists(QC_dir)){dir.create(QC_dir, recursive=T)}
  ### PreProcessing directory
  ###########################
  ###   dataset loading   ###
  ###########################
  if (cellranger == TRUE){
    umi_data <- Read10X(data.dir=input_data) #check if 10X data or not
  } else {
    umi_data <- read.table(file=input_matrix, header=T, sep="\t", row.names=1) #check if 10X data or not
  }
  ### creation of a Seurat object
  # Keep all genes expressed in >= ~0.1% of cells
  # Keep all cells with at least 200 detected genes
  n_cells <- round((ncol(umi_data)*percentage)/100)
  n_genes <- gene_filter
  # min.cells: Include features detected in at least this many cells
  # min.features: Include cells where at least this many features are detected
  umi <- CreateSeuratObject(counts=umi_data, min.cells=n_cells, min.features=n_genes, project=sample)
  Starting_Cells <- ncol(umi)
  #####################################
  ###	  calculation of QC metrics   ###
  #####################################
  ### calculate the percentage of mitochondria, ribosomal and dissociation genes
  if(organism == 'human'){
    dissociation_genes <- c("ACTG1","ANKRD1","ARID5A","ATF3","ATF4","BAG3","BHLHE40","BRD2","BTG1","BTG2","CCNL1","CCRN4L","CEBPB","CEBPD","CEBPG","CSRNP1","CXCL1","CYR61","DCN","DDX3X","DDX5","DES","DNAJA1","DNAJB1","DNAJB4","DUSP1","DUSP8","EGR1","EGR2","EIF1","EIF5","ERF","ERRFI1","FAM132B","FOS","FOSB","FOSL2","GADD45A","GADD45G","GCC1","GEM","H3F3B","HIPK3","HSP90AA1","HSP90AB1","HSPA1A","HSPA1B","HSPA5","HSPA8","HSPB1","HSPE1","HSPH1","ID3","IDI1","IER2","IER3","IER5","IFRD1","IL6","IRF1","IRF8","ITPKC","JUN","JUNB","JUND","KCNE4","KLF2","KLF4","KLF6","KLF9","LITAF","LMNA","MAFF","MAFK","MCL1","MIDN","MIR22HG","MT1","MT2","MYADM","MYC","MYD88","NCKAP5L","NCOA7","NFKBIA","NFKBIZ","NOP58","NPPC","NR4A1","ODC1","OSGIN1","OXNAD1","PCF11","PDE4B","PER1","PHLDA1","PNP","PNRC1","PPP1CC","PPP1R15A","PXDC1","RAP1B","RASSF1","RHOB","RHOH","RIPK1","SAT1","SBNO2","SDC4","SERPINE1","SKIL","SLC10A6","SLC38A2","SLC41A1","SOCS3","SQSTM1","SRF","SRSF5","SRSF7","STAT3","TAGLN2","TIPARP","TNFAIP3","TNFAIP6","TPM3","TPPP3","TRA2A","TRA2B","TRIB1","TUBB4B","TUBB6","UBC","USP2","WAC","ZC3H12A","ZFAND5","ZFP36","ZFP36L1","ZFP36L2","ZYX")
    dissociation_genes <- paste(dissociation_genes, collapse='|')
  } else if(organism == 'mouse') {
    dissociation_genes <- c("ACTG1","ANKRD1","ARID5A","ATF3","ATF4","BAG3","BHLHE40","BRD2","BTG1","BTG2","CCNL1","CCRN4L","CEBPB","CEBPD","CEBPG","CSRNP1","CXCL1","CYR61","DCN","DDX3X","DDX5","DES","DNAJA1","DNAJB1","DNAJB4","DUSP1","DUSP8","EGR1","EGR2","EIF1","EIF5","ERF","ERRFI1","FAM132B","FOS","FOSB","FOSL2","GADD45A","GADD45G","GCC1","GEM","H3F3B","HIPK3","HSP90AA1","HSP90AB1","HSPA1A","HSPA1B","HSPA5","HSPA8","HSPB1","HSPE1","HSPH1","ID3","IDI1","IER2","IER3","IER5","IFRD1","IL6","IRF1","IRF8","ITPKC","JUN","JUNB","JUND","KCNE4","KLF2","KLF4","KLF6","KLF9","LITAF","LMNA","MAFF","MAFK","MCL1","MIDN","MIR22HG","MT1","MT2","MYADM","MYC","MYD88","NCKAP5L","NCOA7","NFKBIA","NFKBIZ","NOP58","NPPC","NR4A1","ODC1","OSGIN1","OXNAD1","PCF11","PDE4B","PER1","PHLDA1","PNP","PNRC1","PPP1CC","PPP1R15A","PXDC1","RAP1B","RASSF1","RHOB","RHOH","RIPK1","SAT1","SBNO2","SDC4","SERPINE1","SKIL","SLC10A6","SLC38A2","SLC41A1","SOCS3","SQSTM1","SRF","SRSF5","SRSF7","STAT3","TAGLN2","TIPARP","TNFAIP3","TNFAIP6","TPM3","TPPP3","TRA2A","TRA2B","TRIB1","TUBB4B","TUBB6","UBC","USP2","WAC","ZC3H12A","ZFAND5","ZFP36","ZFP36L1","ZFP36L2","ZYX")
    dissociation_genes <- tools::toTitleCase(tolower(dissociation_genes))
    dissociation_genes <- paste(dissociation_genes, collapse='|')
  }
  #
  if(organism == 'human'){
    umi[["percent_mt"]] <- PercentageFeatureSet(umi, pattern = "^MT-")
    umi[["percent_ribo"]] <- PercentageFeatureSet(umi, pattern = '^RPL|^RPS|^MRPL|^MRPS')
    umi[["percent_disso"]] <- PercentageFeatureSet(umi, pattern = dissociation_genes)
  } else if(organism == 'mouse') {
    umi[["percent_mt"]] <- PercentageFeatureSet(umi, pattern = "^mt-")
    umi[["percent_ribo"]] <- PercentageFeatureSet(umi, pattern ='^Rpl|^Rps|^Mrpl|^Mrps')
    umi[["percent_disso"]] <- PercentageFeatureSet(umi, pattern = dissociation_genes)
  }
  ### Violin Plot on number of genes, number of UMI and fraction of mitochondrial genes
  cat(bold(green("\nPlotting QC Violin plots \n")))
  suppressWarnings({pdf(paste0(QC_dir, "02a_QC_violin_plots.pdf"), width=24, useDingbats=FALSE)
    vln1 <- popsicler:::VLN(umi, 10, feats="nFeature_RNA", colours= "tomato", 0.01)+ NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggplot2::xlab("")
    vln2 <- popsicler:::VLN(umi, 10, feats="nCount_RNA", colours= "tomato", 0.01)+ NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggplot2::xlab("")
    vln3 <- popsicler:::VLN(umi, 10, feats="percent_mt", colours= "dodgerblue", 0.01)+ NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggplot2::xlab("")
    vln4 <- popsicler:::VLN(umi, 10, feats="percent_ribo", colours= "yellow", 0.01)+ NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggplot2::xlab("")
    vln5 <- popsicler:::VLN(umi, 10, feats="percent_disso", colours= "firebrick", 0.01)+ NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggplot2::xlab("")
    print(patchwork::wrap_plots(vln1 | vln2 | vln3 | vln4 | vln5 + plot_layout(guides = 'collect') + NoLegend()) +
            plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))))
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("02.QC_Plots\\02a_violin_plots.pdf \n")))
  ### Density plot
  cat(bold(green("Plotting QC Density plots \n")))
  suppressWarnings({pdf(file.path(QC_dir, "02b_QC_Hist_nGene_nUMI_MTf_Ribo.pdf"), useDingbats=FALSE)
    plot.title <- "nGene"
    nGene <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    plot.title <- "nGene zoom"
    nGene_zoom <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title) + ggplot2::xlim(0,2000) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    plot.title <- "nUMI"
    nUMI <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    plot.title <- "nUMI zoom"
    nUMI_zoom <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title) + ggplot2::xlim(0,5000) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    print(patchwork::wrap_plots(nGene + nGene_zoom + nUMI + nUMI_zoom + plot_layout(guides = 'collect')+ NoLegend()) +
            plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))) + NoLegend())
    plot.title <- "Mitochondrial Fraction"
    mt_fraction <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_mt, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="blue", fill="dodgerblue") + ggplot2::ggtitle(plot.title) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    plot.title <- "Ribosomal Fraction"
    ribosomal_fraction <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_ribo, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="darkgoldenrod3", fill="yellow") + ggplot2::ggtitle(plot.title) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    plot.title <- "Dissociation Fraction"
    dissociation_fraction <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_disso, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="firebrick") + ggplot2::ggtitle(plot.title) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    print(patchwork::wrap_plots((mt_fraction + ribosomal_fraction) / dissociation_fraction + plot_layout(guides = 'collect') + NoLegend()) +
            plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))) + NoLegend())
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("02.QC_Plots\\02b_QC_Hist_nGene_nUMI_MTf_Ribo.pdf \n")))
  ### Scatter Plot
  cat(bold(green("Plotting QC Scatter plots \n")))
  suppressWarnings({pdf(file.path(QC_dir, "02c_QC_Scatter_nGene_nUMI_MTf.pdf"), width=18, height=12, useDingbats=FALSE)
    plotA <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="nCount_RNA", pt.size=0.5, cols="tomato") +  theme(plot.title=element_blank()) + NoLegend()
    plotB <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_mt", pt.size=0.5, cols="tomato") +  theme(plot.title=element_blank()) + NoLegend()
    plotC <- FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_mt", pt.size=0.5, cols="dodgerblue") +  theme(plot.title=element_blank()) + NoLegend()
    plotD <- FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_ribo", pt.size=0.5, cols="yellow") +  theme(plot.title=element_blank())+ NoLegend()
    plotE <- FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_disso", pt.size=0.5, cols="firebrick") +  theme(plot.title=element_blank())+ NoLegend()
    print(patchwork::wrap_plots(plotA + plotB + ggExtra::ggMarginal(plotC, type="density", color="blue", fill="dodgerblue") + ggExtra::ggMarginal(plotD, type="density", color="darkgoldenrod3", fill="yellow") + ggExtra::ggMarginal(plotE, type="density", color="firebrick", fill="red") + plot_layout(guides = 'collect')+ NoLegend()) +
            plot_annotation(title = unique(as.character(umi@meta.data$orig.ident)), theme=theme(plot.title = element_text(hjust = 0.5, face="bold")), tag_levels = 'A') + NoLegend())
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("02.QC_Plots\\02c_QC_Scatter_nGene_nUMI_MTf.pdf \n")))
  #
  ###########################################################
  ###   check for presence and number of selected genes   ###
  ###########################################################
  #
  #umi@meta.data$[, paste0(opt$gene, '_expressed')]
  if(length(genelist)!=0){
    #for (gene in genelist)
    #{
    popsicler:::plotGene(genelist, umi, QC_dir)
    # }
  }
  return(umi)
}

FilterPlots <- function(UMI, G_RNA_low = 0, G_RNA_hi = Inf, U_RNA_low = 0, U_RNA_hi = Inf, percent_mt_hi = 100, percent_ribo_hi = 100, percent_disso_hi = 100, out_folder=getwd()){
  #######################################
  ###   define and apply thresholds   ###
  #######################################
  QC_dir <- paste0(out_folder,"/02.QC_Plots/")
  suppressWarnings(if (!file.exists(QC_dir)){dir.create(QC_dir, recursive=T)})
  # nFeature: the number of genes detected in each cell
  # nCount_RNA: the total number of molecules detected within a cell
  # percent_mt: percentage of mitochondrial genes (previously calculated)
  # percent_ribo: percentage of ribosomal genes (previously calculated)
  # percent_disso: percentage of dissociation genes (previously calculated)
  umi <- UMI
  umi$filtered <- umi$nFeature_RNA <= G_RNA_low | umi$nFeature_RNA >= G_RNA_hi | umi$percent_mt >= percent_mt_hi | umi$nCount_RNA >=  U_RNA_hi | umi$nCount_RNA <=  U_RNA_low | umi$percent_ribo >= percent_ribo_hi| umi$percent_disso >= percent_disso_hi
  cat(bold(green("The selected thresholds will filter",sum(umi$filtered) ,"cells\n")))

    ### set thresholds vlines and hlines for ggplot
  genes.dn.lim <- ggplot2::geom_vline(xintercept=G_RNA_low, linetype="dashed", color="darkgrey")
  genes.up.lim <- ggplot2::geom_vline(xintercept=G_RNA_hi, linetype="dashed", color="darkgrey")
  #nUMI -> nCount_RNA
  umi.dn.lim <- ggplot2::geom_vline(xintercept=U_RNA_low, linetype="dashed", color="darkgrey")
  umi.up.lim <- ggplot2::geom_vline(xintercept=U_RNA_hi, linetype="dashed", color="darkgrey")
  #umi lines on y axis
  umi.dn.lim.y <- ggplot2::geom_hline(yintercept=U_RNA_low, linetype="dashed", color="darkgrey")
  umi.up.lim.y <- ggplot2::geom_hline(yintercept=U_RNA_hi, linetype="dashed", color="darkgrey")
  #Mt
  mt.dn.lim <- ggplot2::geom_vline(xintercept=percent_mt_hi, linetype="dashed", color="darkgrey")
  mt.dn.lim.y <- ggplot2::geom_hline(yintercept=percent_mt_hi, linetype="dashed", color="darkgrey")
  #Ribo
  ribo.dn.lim <- ggplot2::geom_vline(xintercept=percent_ribo_hi, linetype="dashed", color="darkgrey")
  ribo.dn.lim.y <- ggplot2::geom_hline(yintercept=percent_ribo_hi, linetype="dashed", color="darkgrey")
  #Disso
  disso.dn.lim <- ggplot2::geom_vline(xintercept=percent_disso_hi, linetype="dashed", color="darkgrey")
  disso.dn.lim.y <- ggplot2::geom_hline(yintercept=percent_disso_hi, linetype="dashed", color="darkgrey")
  ### distribution of total number of gene detected per cell
  cat(bold(green("\nPlotting QC final plots")))
  suppressWarnings({pdf(paste0(QC_dir, "/02f_final_Hist_plots.pdf"),width=20, height=10, useDingbats=FALSE)
    plot.title <- "density total genes"
    final <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) +
      genes.dn.lim + genes.up.lim + NoLegend()
    plot.title <- "density total genes zoom"
    final1 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title)  + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + ggplot2::xlim(0,2000) + theme(axis.title.y = element_blank()) +
      genes.dn.lim + genes.up.lim + NoLegend()
    plot.title <- "density total UMI"
    final2 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title)  + theme(plot.title = element_text(hjust = 0.5, face ="bold")) +
      umi.dn.lim + umi.up.lim + NoLegend()
    plot.title <- "density total UMI zoom"
    final3 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title)  + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + ggplot2::xlim(0,5000) + theme(axis.title.y = element_blank()) +
      umi.dn.lim + umi.up.lim + NoLegend()
    print(patchwork::wrap_plots(final + final1 + final2 + final3 + plot_layout(guides = 'collect')))
    plot.title <- "density MT fraction"
    final4 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_mt, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="blue", fill="dodgerblue") + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + NoLegend() + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + mt.dn.lim
    plot.title <- "density Ribosomal fraction"
    final5 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_ribo, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="darkgoldenrod3", fill="yellow") + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + NoLegend() + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + ribo.dn.lim
    plot.title <- "density Dissociation fraction"
    final6 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_disso, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="firebrick") + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + NoLegend() + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + disso.dn.lim
    print(patchwork::wrap_plots(final4 + final5 + final6 + plot_layout(guides = 'collect')))})
  invisible(dev.off())
  cat(paste0(silver("\nPlots saved in: ")),bold(silver("02.QC_Plots\\02f_final_Hist_plots.pdf \n")))

  cat(bold(green("Plotting QC final scatter plots \n")))
  suppressWarnings({pdf(paste0(QC_dir, "/02g_final_Scatter_plots.pdf"),width=20, height=10, useDingbats=FALSE)
    plot1 <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="nCount_RNA", pt.size = 0.3, group.by = "filtered")+ scale_color_manual(values = c("tomato", "#666666")) +  theme(plot.title=element_blank()) +
      genes.dn.lim + genes.up.lim + umi.dn.lim.y + umi.up.lim.y + NoLegend()
    plot2 <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_mt", pt.size = 0.3, group.by = "filtered")+ scale_color_manual(values = c("dodgerblue", "#666666")) + theme(plot.title=element_blank()) +
      genes.dn.lim + genes.up.lim + mt.dn.lim.y + NoLegend()
    plot3 <- FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_mt", pt.size = 0.3, group.by = "filtered")+ scale_color_manual(values = c("dodgerblue", "#666666")) + theme(plot.title=element_blank()) +
      umi.dn.lim + umi.up.lim + mt.dn.lim.y + NoLegend()
    plot4 <- FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_ribo", pt.size = 0.3, group.by = "filtered")+ scale_color_manual(values = c("yellow", "#666666")) + theme(plot.title=element_blank()) +
      umi.dn.lim + umi.up.lim + ribo.dn.lim.y + ggplot2::xlim(0,5000) + NoLegend()
    plot5 <- FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_disso", pt.size = 0.3, group.by = "filtered")+ scale_color_manual(values = c("firebrick", "#666666")) + theme(plot.title=element_blank()) +
      umi.dn.lim + umi.up.lim + disso.dn.lim.y + ggplot2::xlim(0,5000) + NoLegend()
    print(patchwork::wrap_plots(plot1 | plot2 + plot3 + plot4 + plot5+ plot_layout(guides = 'collect') + NoLegend()))
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("02.QC_Plots\\02g_final_Scatter_plots \n")))
  umi <- subset(umi, subset = nFeature_RNA > G_RNA_low &
                  nFeature_RNA < G_RNA_hi &
                  nCount_RNA   > U_RNA_low &
                  nCount_RNA   < U_RNA_hi &
                  percent_mt   < percent_mt_hi &
                  percent_ribo < percent_ribo_hi &
                  percent_disso < percent_disso_hi)
  
  cat(paste0(cyan("\nNext suggested step is Doublets Calculation, run")),bold(cyan("Calculate Doublets \n")))
  cat(paste0(bold(cyan("\nWARNING: \n")),cyan("It is strictly recommended to previously run that step setting "), bold(cyan("dbs_thr ='none' ")),cyan("and "), bold(cyan("dbs_remove= FALSE. \n")), cyan("Once checked the graphs it is possible to re-run this step specifying a custom threshold \nthrough the 'dbs_thr' parameter and removing all the cells identified as doublets setting 'dbs_remove' parameter as TRUE. \n")))
  return(umi)
}

CalculateDoublets <- function(UMI, dbs_thr='none', dbs_remove=TRUE, out_folder=getwd()){
  QC_dir <- paste0(out_folder,"/02.QC_Plots/")
  suppressWarnings(if (!file.exists(QC_dir)){dir.create(QC_dir, recursive=T)})
  if(dbs_thr == 'none' | !"doublets_score" %in% colnames(UMI@meta.data)){
    ### calculate doublets
    doublets <- scrubDoublets(as.matrix(UMI@assays$RNA@counts), directory=QC_dir, expected_doublet_rate=0.1)
    #return(list(scrubDoublets = predicted_doublets, doublet_scores_obs = doublet_scores$doublet_scores_obs, doublet_scores_sim = doublet_scores$doublet_scores_sim))
    names(doublets) <- c("predicted", "score_predicted", "score_simulated")
    ### if doublets exist, calculate, set to 0 otherwise
    ifelse(TRUE %in% doublets$predicted, dbs_found <- table(doublets$predicted)[[2]], dbs_found <- 0)
    thr <- round(mean(c(min(doublets$score_predicted[doublets$predicted]),max(doublets$score_predicted[!doublets$predicted]))), 3)
    ### Threshold definition
    doublets$predicted <- c(doublets$score_predicted > thr)
    ### add features to umi object
    UMI$doublets_score <- doublets$score_predicted
    UMI$doublets <- doublets$predicted
    
    UMI2 <- NormalizeData(UMI, verbose=FALSE)
    UMI2 <- FindVariableFeatures(UMI2, selection.method="vst", nfeatures=2000, verbose=FALSE)
    UMI2 <- ScaleData(UMI2, vars.to.regress=NULL, verbose=FALSE)
    UMI2 <- RunPCA(UMI2, n_pcs=30, verbose=FALSE)
    UMI2 <- suppressWarnings(RunUMAP(UMI2, dims=1:10, verbose=FALSE))
    Idents(UMI2) <- "doublets"
    highlight_labels <- list("doublet"= WhichCells(UMI2, idents = TRUE), "singlet"= WhichCells(UMI2, idents = FALSE))
    cat(bold(green("Plotting doublets UMAP \n")))
    pdf(paste0(QC_dir,"/02h_doublets_umap.pdf"),15,8, useDingbats=FALSE)
    p1 <- DimPlot(UMI2, reduction="umap", group.by = "doublets", pt.size=0.5, cols=c("lightgrey"), cells.highlight = highlight_labels, cols.highlight = "black")
    p2 <- FeaturePlot(UMI2, reduction="umap", features="doublets_score", pt.size=0.5) +scale_colour_gradientn(colours=c("lightgrey", "red", "darkred", "black"))
    print(patchwork::wrap_plots(p1 | p2 + plot_layout(guides = 'collect')))
    dev.off()
    cat(paste0(silver("Plots saved in: ")),bold(silver("02.QC_Plots\\02h_doublets_umap.pdf \n")))
    if(dbs_remove == TRUE){
      cat(paste0(green("Removing", sum(UMI$doublets) ,"Doublets \n")))
      UMI <- UMI[,!UMI$doublets]
    }
    cat(paste0(cyan("Once checked the graphs it is possible to re-run this step specifying a custom threshold through the 'dbs_thr' parameter and removing al the cells identified as doublets setting 'dbs_remove' parameter as TRUE. \n")))
  }else{
    UMI$doublets <- ifelse(UMI$doublets_score > dbs_thr, TRUE, FALSE)
    if(dbs_remove == TRUE){
      cat(paste0(green("Removing", sum(UMI$doublets) ,"Doublets \n")))
      UMI <- UMI[,!UMI$doublets]
      cat(paste0(cyan("Next suggested step is data normalization, run")),bold(cyan("Normalize \n")))
    }
  }
  return(UMI)
}

Normalize <- function(UMI, variable_genes=2000, out_folder=getwd()){
  PP_dir <- paste0(out_folder,"/03.PreProcessing/")
  if (!file.exists(PP_dir)){dir.create(PP_dir, recursive=T)}
  suppressWarnings({umi <- NormalizeData(object=UMI, normalization.method="LogNormalize", scale.factor=1e4)
  cat(bold(green("\nPlotting Normalization graphs \n")))
  pdf(paste0(PP_dir, "/03a_total_expression_after_before_norm.pdf"), useDingbats=FALSE)
  par(mfrow = c(2,1))
  hist(colSums(as.matrix(umi@assays$RNA@counts)), breaks=100, main="Total expression before normalization", xlab="Sum of expression")
  hist(colSums(as.matrix(umi@assays$RNA@data)), breaks=100, main="Total expression after normalization", xlab="Sum of expression")
  invisible(dev.off())
  cat(paste0(silver("Plots saved in: ")),bold(silver("03.PreProcessing\\03a_total_expression_after_before_norm.pdf \n")))
  umi <- FindVariableFeatures(umi, selection.method = "vst", nfeatures = variable_genes)
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(umi), 10)
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(umi)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # When using repel, set xnudge and ynudge to 0 for optimal results
  cat(bold(green("Plotting High Variables Genes \n")))
  pdf(paste0(PP_dir, "/03b_plot_FindVariableGenes.pdf"), useDingbats=FALSE)
  print(patchwork::wrap_plots(plot1 / plot2 + plot_layout(guides = 'collect')), ncol=1, nrow=2)
  cat(paste0(silver("Plots saved in: ")),bold(silver("03.PreProcessing\\03b_plot_FindVariableGenes.pdf \n")))
  invisible(dev.off())})
  cat(paste0(cyan("Next suggested step is calculate Cell Cycle Score, run ")),bold(cyan("CCScore \n")))
  return(umi)
}

CCScore <- function(UMI, organism='human'){
  if(organism == 'human') {
    ### Assign Cell-Cycle Scores
    cc.genes$s.genes <- intersect(cc.genes$s.genes, row.names(UMI@assays$RNA@counts))
    cc.genes$g2m.genes <- intersect(cc.genes$g2m.genes, row.names(UMI@assays$RNA@counts))
    umi <- CellCycleScoring(UMI, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=T)
  } else if(organism == 'mouse'){
    m.s.genes <- c("Mcm4", "Exo1", "Slbp", "Gmnn", "Cdc45", "Msh2", "Mcm6", "Rrm2", "Pold3", "Blm", "Ubr7", "Mcm5", "Clspn", "Hells", "Nasp", "Rpa2", "Rad51ap1", "Tyms", "Rrm1", "Rfc2", "Prim1", "Brip1", "Usp1", "Ung", "Pola1", "Mcm2", "Fen1", "Tipin", "Pcna", "Cdca7", "Uhrf1", "Casp8ap2", "Cdc6", "Dscc1", "Wdr76", "E2f8", "Dtl", "Ccne2", "Atad2", "Gins2", "Chaf1b", "Pcna-ps2")
    m.g2m.genes <- c("Nuf2", "Psrc1", "Ncapd2", "Ccnb2", "Smc4", "Lbr", "Tacc3", "Cenpa", "Kif23", "Cdca2", "Anp32e", "G2e3", "Cdca3", "Anln", "Cenpe", "Gas2l3", "Tubb4b", "Cenpf", "Dlgap5", "Hjurp", "Cks1brt", "Gtse1", "Bub1", "Birc5", "Ube2c", "Rangap1", "Hmmr", "Ect2", "Tpx2", "Ckap5", "Cbx5", "Nek2", "Ttk", "Cdca8", "Nusap1", "Ctcf", "Cdc20", "Cks2", "Mki67", "Tmpo", "Ckap2l", "Aurkb", "Kif2c", "Cdk1", "Kif20b", "Top2a", "Aurka", "Ckap2", "Hmgb2", "Cdc25c", "Ndc80", "Kif11")
    m.s.genes <- intersect(m.s.genes, row.names(UMI@assays$RNA@counts))
    m.g2m.genes <- intersect(m.g2m.genes, row.names(UMI@assays$RNA@counts))
    cat(bold(green("Calculating Cell Cycle Score \n")))
    umi <- CellCycleScoring(UMI, s.features=m.s.genes, g2m.features=m.g2m.genes, set.ident=T)
  }
  cat(paste0(cyan("Next suggested step is regression, run ")),bold(cyan("ApplyRegression")),cyan("\nIt is suggested to apply regression without providing any variables to regress and without exploring PCs (explore_PC=FALSE) to save time and computational effort.\nRegression variables can be chosen through visual inspection of this function outputs \nOnce selected, set explore_PC as TRUE and explore graphs to identify the PCs number to use in your analysis \n"))
  return(umi)
}


ApplyRegression <- function(UMI, variables='none', explore_PC=FALSE,  out_folder=getwd()){
  PP_dir <- file.path(out_folder,"03.PreProcessing")
  suppressWarnings(if (!file.exists(PP_dir)){dir.create(PP_dir, recursive=T)})
  cycle.dir <- file.path(PP_dir,paste0("No_Regression"))
  suppressWarnings(if (!file.exists(cycle.dir)){dir.create(cycle.dir)})
  all.genes <- rownames(UMI)
  UMI <- ScaleData(object=UMI, features=all.genes)
  ### perform PCA on the scaled data.
  UMI <- suppressMessages(RunPCA(UMI, features=UMI[["RNA"]]@var.features, npcs=30, do.print=F, verbose=FALSE))
  UMI2 <- suppressMessages(RunTSNE(UMI, dims = 1:20))
  UMI2 <- suppressWarnings(RunUMAP(UMI2, dims = 1:20, verbose=FALSE))
  ### Single pdf with 3 pages, 4 plots per page
  cat(bold(green("Plotting dimensional reduction graphs with no regression \n")))
  pdf(file.path(cycle.dir, "/03c_DimReduction_NoRegression.pdf"), width=14, height=12, useDingbats=FALSE)
  popsicler:::four_plots(UMI2, "pca")
  popsicler:::four_plots(UMI2, "tsne")
  popsicler:::four_plots(UMI2, "umap")
  invisible(dev.off())
  cat(paste0(silver("Plots saved in: ")),bold(silver("03.PreProcessing dedicated subfolder \n")))
  ### inizio non fare se NULL ###
  if (variables != 'none'){
    cycle.dir <- file.path(PP_dir,paste0("Regression_on_",paste(unlist(variables), collapse='_')))
    suppressWarnings(if (!file.exists(cycle.dir)){dir.create(cycle.dir)})
    UMI <- ScaleData(object=UMI, vars.to.regress=variables, features=all.genes)
    ### perform PCA on the scaled data.
    UMI <- RunPCA(UMI, features=UMI[["RNA"]]@var.features, npcs=30, do.print=F, verbose=FALSE)
    UMI2 <- suppressMessages(RunTSNE(UMI, dims = 1:20))
    UMI2 <- suppressWarnings(RunUMAP(UMI2, dims = 1:20, verbose=FALSE))
    ### PCA plots
    cat(bold(green("Plotting dimensional reduction graphs after regression \n")))
    pdf(file.path(cycle.dir, "/03c_DimReduction_PostRegression.pdf"), width=14, height=12, useDingbats=FALSE)
    popsicler:::four_plots(UMI2, "pca")
    popsicler:::four_plots(UMI2, "tsne")
    popsicler:::four_plots(UMI2, "umap")
    invisible(dev.off())
    cat(paste0(silver("Plots saved in: ")),bold(silver("03.PreProcessing dedicated subfolder \n")))
  }
  ### fine non fare se null ####
  ### NEW (old VizPCA)
  ### se explore_PC = true ###
  
  if (explore_PC == TRUE){
    cat(bold(green("Plotting graphs to explore PCs \n")))
    pdf(paste0(cycle.dir, "/03d_VizPCA_HVG.pdf"), useDingbats=FALSE)
    print(VizDimLoadings(UMI, dims = 1:2, reduction = "pca"))
    invisible(dev.off())
    ### Determine statistically significant principal components ###
    # 1. exploring PCs to determine relevant sources of heterogeneity (could be used in conjunction with GSEA for example)
    ### NEW function (OLD PCHeatmap)
    pdf(paste0(cycle.dir, "/03e_PCHeatmap.pdf"), width=9, height=12, useDingbats=FALSE)
    DimHeatmap(UMI, dims=1:9, cells = 500, balanced = TRUE)
    DimHeatmap(UMI, dims=10:18, cells = 500, balanced = TRUE)
    DimHeatmap(UMI, dims=19:27, cells = 500, balanced = TRUE)
    invisible(dev.off())
    # 2. Significant PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line).
    UMI2 <- suppressWarnings(JackStraw(object=UMI, dims=30, num.replicate=100))
    UMI2 <- suppressWarnings(ScoreJackStraw(UMI2, dims = 1:30))
    pdf(paste0(cycle.dir, "/03f_JackStrawPlot.pdf"), width=15, height=15, useDingbats=FALSE)
    print(JackStrawPlot(object=UMI2, dims=1:30))
    invisible(dev.off())
    # 3. plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph
    ### NEW function (OLD PCElbowPlot)
    pdf(paste0(cycle.dir, "/03g_PCElbowPlot.pdf"), useDingbats=FALSE)
    print(ElbowPlot(object=UMI, ndims=30))
    invisible(dev.off())
    cat(paste0(silver("Plots saved in: ")),bold(silver("03.PreProcessing dedicated subfolder \n")))
    cat(paste0(cyan("\nOnce identified the PCs number to use in your analysis, perform clustering running "),bold(cyan("Calculate Cluster \n"))))
    }
  ### fine explore_pc ###
  return(UMI)
}

CalculateCluster <- function(UMI, dim_pca, organism="human", marker.list='none', PCA=TRUE, cluster_res=0.8, out_folder=getwd()){
  #dimpca --> numero di principal components
  #cluster_res --> cluster resolution (default 0.8)
  

  Cluster_dir <- paste0(out_folder,"/04.Clustering/")
  if (!file.exists(Cluster_dir)){dir.create(Cluster_dir, recursive=T)}
  if(organism == 'human') {
    if(marker.list == 'none'){
      marker.list <- list(Immune=c("PTPRC"),
                          HSC=c("CD34"),
                          Adult_sc=c("PROM1"),
                          T_cell=c("CD3D", "CD3E", "TNFRSF4"),
                          CD4=c("CD4"), T_ProB=c("IL7R"),
                          T_B_CD.subset=c("CCR7"),
                          Treg=c("FOXP3"), CD8=c("CD8A"),
                          ProB=c("MS4A1"), B_cell=c("CD19", "CD22"),
                          Plasma.Cell=c("IGHD", "IGHM"),
                          Immature.B=c("CD79A", "CD38"),
                          NK=c("GNLY", "NKG7"),
                          Monocyte=c("CD14", "LYZ", "CST3", "FCGR3A", "MS4A7"),
                          DC=c("CD1E", "FCER1A"),
                          Platelet=c("PPBP", "ITGA2B"))
    } else {
      signature_dir <- getwd()
      markers <-read.table(paste0(signature_dir,"/",marker.list),sep="\t",header=T)
      marker.list <- as.list(markers)
    }
  }else if(organism == 'mouse'){
    if(marker.list == 'none'){
      marker.list <- list(Immune=c("Ptprc"),
                          HSC=c("Cd34"),
                          Adult_sc=c("Prom1"),
                          T_cell=c("Cd3d", "Cd3e", "Tnfrsf4"),
                          CD4=c("Cd4"), T_ProB=c("Il7r"),
                          T_B_CD.subset=c("Ccr7"),
                          Treg=c("Foxp3"), CD8=c("Cd8a"),
                          ProB=c("Ms4a1"), B_cell=c("Cd19", "Cd22"),
                          Plasma.Cell=c("Ighd", "Ighm"),
                          Immature.B=c("Cd79a", "Cd38"),
                          NK=c("Gnly", "Nkg7"),
                          Monocyte=c("Cd14", "Lyz", "Cst3", "Fcgr3a", "Ms4a7", "Vcan", "Ly6c1", "Ly6c2"),
                          Macrophage=c("Cd80", "Cd68", "C1qa", "C1qb", "Adgre1"),
                          Neutrophil=c("Ly6g", "Cd11b", "S100a8", "S100a9", "Csf3r"),
                          granulocyte=c("Cd66b"),
                          DC=c("Cd1e", "Fcer1a", "Cd208", "Cd265", "Xcr1", "Batf3", "Fscn1"),
                          pDC=c("Siglech", "Clec4b1", "Nrp1", "Ido1"),
                          Basophil=c("Gata2", "Cpa3", "Ms4a2"),
                          Platelet=c("Ppbp", "Itga2b"),
                          Erithrocyte=c("Cd235a", "Gypa", "Hba-a1", "Hba-a2"),
                          Endothelial=c("Mcam", "Vcam1", "Vwf", "Pecam1", "Sele", "Cd93", "Nectin3", "Tek"),
                          Epithelial=c("Cdh1", "Cd326", "Epcam"),
                          Fibroblast=c("Vim", "Pdgfra", "Pdgfrb"))
    } else {
      signature_dir <- getwd()
      markers <-read.table(paste0(signature_dir,"/",marker.list),sep="\t",header=T)
      marker.list <- as.list(markers)
    }
  }
  ### Performing scaling and clustering based on provided inputs
  UMI <- RunUMAP(UMI, dims = 1:dim_pca, verbose=FALSE)
  UMI <- RunTSNE(UMI, dims = 1:dim_pca, verbose=FALSE)   
  UMI <- FindNeighbors(UMI, dims = 1:dim_pca)
  UMI <- FindClusters(UMI, resolution = cluster_res)

  ### clustree
  if (length(cluster_res)>1) { 
  metadata<-UMI@meta.data
  for (res in cluster_res) {
    res_col<-paste0("RNA_snn_res.",res)
    metadata[,res_col]<-as.numeric(as.character(metadata[,res_col]))
  }
  pdf(file.path(Cluster_dir,paste0("04a_clustree_analysis.pdf")), width=14, height=10)
  print(clustree(metadata, prefix = "RNA_snn_res."))
  dev.off()
  }
  
  ### perform also PCA plot
  if(PCA == TRUE){
    algos <- c("umap","tsne","pca")
  }else{
    algos <- c("umap","tsne")
  }
  for(algo in algos){
   algo_plot(Cluster_dir, UMI, algo, "04a_clusters" , cluster_res)
  }
  if ("doublets" %in% colnames(UMI@meta.data)){
  pdf(file.path(Cluster_dir, "04a_umap_doublets.pdf"), width=18, height=10, useDingbats=FALSE)
  p1 <- DimPlot(UMI, reduction="umap", group.by = "doublets", pt.size=1, cols=c("lightgrey","black"))+
  xlab("UMAP 1") + ylab("UMAP 2")
  p2 <- FeaturePlot(UMI, reduction="umap", features="doublets_score", pt.size=1) +scale_colour_gradientn(colours=c("lightgrey", "red", "darkred", "black"))+
  xlab("UMAP 1") + ylab("UMAP 2")
  print(patchwork::wrap_plots(p1 | p2 + plot_layout(guides = 'collect')))
  dev.off()
  pdf(file.path(Cluster_dir, "04a_tsne_doublets.pdf"), width=18, height=10, useDingbats=FALSE)
  p1 <- DimPlot(UMI, reduction="tsne", group.by = "doublets", pt.size=1, cols=c("lightgrey","black"))+
    xlab("tSNE 1") + ylab("tSNE 2")
  p2 <- FeaturePlot(UMI, reduction="tsne", features="doublets_score", pt.size=1) +scale_colour_gradientn(colours=c("lightgrey", "red", "darkred", "black"))+
    xlab("tSNE 1") + ylab("tSNE 2")
  print(patchwork::wrap_plots(p1 | p2 + plot_layout(guides = 'collect')))
  dev.off()  
  }
  
  ### plot con MALAT1?
  
  
  ### Finding differentially expressed features (cluster biomarkers)
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  # Finding cluster markers
  
  ### if only one resolution has been specified
  if (length(cluster_res)==1) {
    umi.markers <- FindAllMarkers(UMI, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
    FC_col <- colnames(umi.markers)[grep("avg_log", colnames(umi.markers))]
    umi.markers <- umi.markers[,c("cluster","gene","pct.1","pct.2",FC_col,"p_val","p_val_adj")]
    umi.markers <- umi.markers[umi.markers$p_val_adj<=0.05,]
    top.markers.2 <- umi.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=2, wt=get(FC_col))
    ftp_h = 5+(0.5*length(levels((umi.markers$cluster))))
   
    ### visualize markers by violin plots
    cat(bold(green("Plotting Top Markers graphs for each cluster \n")))
    
    pdf(paste0(Cluster_dir, "/04b_Violin_top.markers.pdf"), width=12, height=9, useDingbats=FALSE)
    for(cluster in unique(top.markers.2$cluster)) {
      genes <- top.markers.2$gene[top.markers.2$cluster == cluster]
      vg1 <- popsicler:::VLN(UMI, feats=genes[1], point =FALSE) + ggplot2::geom_boxplot(width=0.1, outlier.shape=NA) + NoLegend() + ggplot2::xlab("")
      if (length(genes)==2) {  
      vg2 <- popsicler:::VLN(UMI, feats=genes[2], point =FALSE) +
        ggplot2::geom_boxplot(width=0.1, outlier.shape=NA) +
        NoLegend() +
        ggplot2::xlab(paste("Cluster",cluster, "@ res", cluster_res)) } else {
          vg1<-vg1+xlab(paste("Cluster",cluster, "@ res", cluster_res))
          vg2<-ggplot()
        }
      print(patchwork::wrap_plots(vg1, vg2, ncol=1))
    }
    invisible(dev.off())
    
    ### Visualize TOP markers on a dimensional reduction plot
    pdf(paste0(Cluster_dir, "/04c_UMAP_Top_Markers.pdf"), width=14, height=7, useDingbats=FALSE)
    for(cluster in unique(top.markers.2$cluster)) {
      genes <- top.markers.2$gene[top.markers.2$cluster == cluster]
      vg1 <- FeaturePlot(UMI, features=genes[1], reduction="umap", ncol=1) +
        NoLegend() +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab(paste0("Cluster ",cluster," @ res ", cluster_res,"\nUMAP 2"))
      if (length(genes)==2) {  
      vg2 <- FeaturePlot(UMI, features=genes[2], reduction="umap", ncol=1) +
        NoLegend() +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab("UMAP 2") } else {vg2<-ggplot()}
      print(patchwork::wrap_plots(vg1, vg2, ncol=2))
    }
    invisible(dev.off())
    
    ### TSNE Top Markers
    pdf(paste0(Cluster_dir, "/04c_tSNE_Top_Markers.pdf"), width=14, height=7, useDingbats=FALSE)
    for(cluster in unique(top.markers.2$cluster)) {
      genes <- top.markers.2$gene[top.markers.2$cluster == cluster]
      vg1 <- FeaturePlot(UMI, features=genes[1], reduction="tsne", ncol=1) +
        NoLegend() +
        ggplot2::xlab("tSNE 1") +
        ggplot2::ylab(paste0("Cluster ",cluster," @ res ", cluster_res,"\ntSNE 2"))
      if (length(genes)==2) {  
      vg2 <- FeaturePlot(UMI, features=genes[2], reduction="tsne", ncol=1) +
        NoLegend() +
        ggplot2::xlab("tSNE 1") +
        ggplot2::ylab("tSNE 2")} else {vg2<-ggplot()}
      print(patchwork::wrap_plots(vg1, vg2, ncol=2))
    }
    invisible(dev.off())
    
    ### Plotting cluster biomarkers
    ### expression heatmap
    top.markers.10 <- umi.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = get(FC_col))
    pdf(paste0(Cluster_dir, "/04d_heatmap_top.markers.pdf"), width=18, height=5+(0.5*length(levels((umi.markers$cluster)))), useDingbats=FALSE)
    print(DoHeatmap(UMI, features=top.markers.10$gene, slot="data") + NoLegend())
    invisible(dev.off())
    cat(paste0(silver("Plots saved in: ")),bold(silver("\\04.Clustering\\")), silver("folder \n"))
  }

  ### Visualize common immune markers on a dimensional reduction plot
  ### filtering default immune markers list based on genes in the dataset
  ### old intersect for immune markers list with no cell line classification
  ### new intersection and subset for immune markers list with cell line classification [EDIT]
  marker.list <- lapply(marker.list, function(x) {x[x %in% row.names(UMI)]})
  ### remove empty lists
  marker.list <- marker.list[lapply(marker.list,length)>0]
  #Here we will investigate the immune markers in the analyzed cells
  #A set of default immune markers is available. These markers are:
  
  #Plotting UMAP and TSNE for immune markers
  ftp_h = 8+(0.8*length(marker.list))
  cat(bold(green("Plotting Markers graphs for each cluster \n")))
  popsicler:::FTP(UMI, Cluster_dir, "04e_", "umap", ftp_h, unlist(marker.list), "_markers")
  popsicler:::FTP(UMI, Cluster_dir, "04e_", "tsne", ftp_h, unlist(marker.list), "_markers")
  
  # larghezza va proporzionata al num di cluster. altezza al num di marker
  vln_w<-length(levels(UMI$RNA_snn_res.0.1))*2
  vln_h<-round((length(unlist(marker.list))/4)*3)
  ### visualize common markers by violin plots
  pdf(paste0(Cluster_dir, "/04f_Violin_markers.pdf"), height = vln_h, width= vln_w, useDingbats=FALSE)
  for (res.i in cluster_res){
  print(VlnPlot(UMI, features=unlist(marker.list), ncol=4,pt.size=F, group.by=paste0("RNA_snn_res.",res.i)))
  }         
  dev.off()
 # manca  valore risoluzione nel plot
 
  markers.unlisted <-as.character(unlist(marker.list))
  width_calc<-floor(length(markers.unlisted)/2.5)
  width_sig<-ifelse(width_calc<10, 10, width_calc)
  width_sig<-ifelse(width_sig>16, 16, width_sig)
  
  ### dotplot of immune markers for each cluster
  pdf(paste0(Cluster_dir, "04g_Dotplot_markers.per.cluster.pdf"), width=width_sig, height=10, useDingbats=FALSE)
  for(res.i in cluster_res) {
    print(DTP(UMI, marker.list, paste0("RNA_snn_res.",res.i)))
  }
  invisible(dev.off())
  
  ##### (TEST) assign clusters (before annotation)
  #require(Seurat)

  pdf(paste0(Cluster_dir, "/04h_Clusters_tree.pdf"), width=10, height=10, useDingbats=FALSE)
  for(res.i in cluster_res) {
    Idents(UMI)<-paste0("RNA_snn_res.",res.i)
    clus_time_tree <- BuildClusterTree(UMI)        
    data.tree <- Tool(object = clus_time_tree, slot = "BuildClusterTree")
    plot.phylo(x = data.tree, direction = "rightwards", main=paste0("phylogenetic tree of clusters @ res", res.i))
  }
  invisible(dev.off())
  
  
  cat(paste0(silver("Plots saved in: ")),bold(silver("\\04.Clustering\\")), silver("folder, files are named with the prefixes 4e, 4f, 4g, 4h \n"))
  cat(paste0(cyan("Next suggested step is annotation, run "),bold(cyan("MakeAnnotation \n"))))
  return(UMI)
}





MakeAnnotation <- function(UMI, organism, marker.list='none', thresh=20, cluster_res=NULL, out_folder=getwd()){
  Annot_dir <- paste0(out_folder,"/05.Annotation/")
  if(!is.null(cluster_res)&!all(paste0("RNA_snn_res.",cluster_res)%in%colnames(UMI@meta.data))) { 
    #cat(paste0(red("Clustering at selected resolution has not been calculated yet. Please change resolution or run "),bold(red("CalculateCluster \n"))))
    stop("Clustering at selected resolution has not been calculated yet. Please change resolution or run CalculateCluster")
    }

  if (!file.exists(Annot_dir)){dir.create(Annot_dir, recursive=T)}
  if(organism == 'human'){
    if(marker.list == 'none'){
      marker.list <- list(Immune=c("PTPRC"),
                          HSC=c("CD34"),
                          Adult_sc=c("PROM1"),
                          T_cell=c("CD3D", "CD3E", "TNFRSF4"),
                          CD4=c("CD4"), T_ProB=c("IL7R"),
                          T_B_CD.subset=c("CCR7"),
                          Treg=c("FOXP3"), CD8=c("CD8A"),
                          ProB=c("MS4A1"), B_cell=c("CD19", "CD22"),
                          Plasma.Cell=c("IGHD", "IGHM"),
                          Immature.B=c("CD79A", "CD38"),
                          NK=c("GNLY", "NKG7"),
                          Monocyte=c("CD14", "LYZ", "CST3", "FCGR3A", "MS4A7"),
                          DC=c("CD1E", "FCER1A"),
                          Platelet=c("PPBP", "ITGA2B"))
    } else {
      signature_dir <- getwd()
      markers <-read.table(paste0(signature_dir,"/",marker.list),sep="\t",header=T)
      marker.list <- as.list(markers)
    }
    ### retrieve reference dataset
    # returns a SummarizedExperiment object (scater package) containing matrix of log-expression values with sample-level labels
    # other reference datasets are available, check the vignette
    ##########################################################
    ### Visualize common immune markers on a dimensional reduction plot
    ### filtering default immune markers list based on genes in the dataset
    ### old intersect for immune markers list with no cell line classification
    ### new intersection and subset for immune markers list with cell line classification [EDIT]
    marker.list <- lapply(marker.list, function(x) {x[x %in% row.names(UMI)]})
    ### remove empty lists
    marker.list <- marker.list[lapply(marker.list,length)>0]
    markers.unlisted <-as.character(unlist(marker.list))
    width_calc<-floor(length(markers.unlisted)/2.5)
    width_sig<-ifelse(width_calc<10, 10, width_calc)
    width_sig<-ifelse(width_sig>16, 16, width_sig)

    hpca.se <- suppressMessages(celldex::HumanPrimaryCellAtlasData())
    BpEn.se <- suppressMessages(celldex::BlueprintEncodeData())
    cat(bold(green("Plotting single cell and cluster annotations \n")))
    UMI <- popsicler:::SR_plots("hpca", hpca.se, UMI, Annot_dir, cluster_res)
    UMI <- popsicler:::SR_plots("BpEn", BpEn.se, UMI, Annot_dir, cluster_res)
    cat(paste0(silver("Plots saved in: ")),bold(silver("\\05.Annotation\\")), silver("folder, with the prefixes 5a and 5b \n"))
    annotations <- c("hpca.sc.main.labels","BpEn.sc.main.labels")
    cat(bold(green("Plotting dimensional reduction graphs for each population found in the sample \n")))
    for(single_annot in annotations){
      popsicler:::annotation_plot(Annot_dir, "05c_", UMI, "UMAP", single_annot, "CellPopulations")
      popsicler:::annotation_plot(Annot_dir, "05c_",  UMI, "TSNE", single_annot, "CellPopulations")
      pdf(paste0(Annot_dir, paste0("/05d_DotPlot_Markers_", single_annot, ".pdf")), width=width_sig, height=10, useDingbats=FALSE)
      popsicler:::DTP(UMI, marker.list, single_annot)
      invisible(dev.off())
    }
    cat(paste0(silver("Plots saved in: ")),bold(silver("\\05.Annotation\\")), silver("folder, with the prefixes 05c and 05d \n"))
  } else if(organism == 'mouse') {
    require("scMCA")
    if(marker.list == 'none'){
      marker.list <- list(Immune=c("Ptprc"),
                          HSC=c("Cd34"),
                          Adult_sc=c("Prom1"),
                          T_cell=c("Cd3d", "Cd3e", "Tnfrsf4"),
                          CD4=c("Cd4"), T_ProB=c("Il7r"),
                          T_B_CD.subset=c("Ccr7"),
                          Treg=c("Foxp3"), CD8=c("Cd8a"),
                          ProB=c("Ms4a1"), B_cell=c("Cd19", "Cd22"),
                          Plasma.Cell=c("Ighd", "Ighm"),
                          Immature.B=c("Cd79a", "Cd38"),
                          NK=c("Gnly", "Nkg7"),
                          Monocyte=c("Cd14", "Lyz", "Cst3", "Fcgr3a", "Ms4a7", "Vcan", "Ly6c1", "Ly6c2"),
                          Macrophage=c("Cd80", "Cd68", "C1qa", "C1qb", "Adgre1"),
                          Neutrophil=c("Ly6g", "Cd11b", "S100a8", "S100a9", "Csf3r"),
                          granulocyte=c("Cd66b"),
                          DC=c("Cd1e", "Fcer1a", "Cd208", "Cd265", "Xcr1", "Batf3", "Fscn1"),
                          pDC=c("Siglech", "Clec4b1", "Nrp1", "Ido1"),
                          Basophil=c("Gata2", "Cpa3", "Ms4a2"),
                          Platelet=c("Ppbp", "Itga2b"),
                          Erithrocyte=c("Cd235a", "Gypa", "Hba-a1", "Hba-a2"),
                          Endothelial=c("Mcam", "Vcam1", "Vwf", "Pecam1", "Sele", "Cd93", "Nectin3", "Tek"),
                          Epithelial=c("Cdh1", "Cd326", "Epcam"),
                          Fibroblast=c("Vim", "Pdgfra", "Pdgfrb"))
    }else {
      signature_dir <- getwd()
      markers <-read.table(paste0(signature_dir,"/",marker.list),sep="\t",header=T)
      marker.list <- as.list(markers)
    }
    #################################################################
    ### Visualize common immune markers on a dimensional reduction plot
    ### filtering default immune markers list based on genes in the dataset
    ### old intersect for immune markers list with no cell line classification
    ### new intersection and subset for immune markers list with cell line classification [EDIT]
    marker.list <- lapply(marker.list, function(x) {x[x %in% row.names(UMI)]})
    ### remove empty lists
    marker.list <- marker.list[lapply(marker.list,length)>0]
    markers.unlisted <-as.character(unlist(marker.list))
    width_calc<-floor(length(markers.unlisted)/2.5)
    width_sig<-ifelse(width_calc<10, 10, width_calc)
    width_sig<-ifelse(width_sig>16, 16, width_sig)
    matrice_norm <- as.matrix(GetAssayData(UMI))
    mca_result <- scMCA(scdata = matrice_norm, numbers_plot = 3)
    scMCA_assignment <- mca_result$scMCA
    Ig.se <- suppressMessages(SingleR::ImmGenData())
    ### Loaded reference datasets (ImmGen)
    ### Performing data simplification
    ### simplification of labels
    simple_assignment <- do.call(rbind,strsplit(scMCA_assignment,"(", fixed=T))
    simple_assignment <- unlist(simple_assignment[,1])
    simple_assignment <- toupper(simple_assignment)
    simple_assignment <- gsub("CELLS","CELL", simple_assignment)
    test <- strsplit(simple_assignment, "_")
    test <- sapply(test, function(x,m) c(x, rep(NA, m-length(x))), max(rapply(test,length)))
    simple_assignment <- unlist(test[1,])
    simple_assignment <- gsub("T-CELL","T CELL",simple_assignment)
    simple_assignment <- gsub("MARCROPHAGE","MACROPHAGE",simple_assignment)
    simple_assignment <- gsub("NK CELL","NK",simple_assignment)
    UMI$scMCA <- scMCA_assignment
    UMI$scMCA_simple <- simple_assignment
    ### Plot TSNE scMCA
    cat(bold(green("Plotting tSNE scMCA annotation \n")))
    pdf(paste0(Annot_dir, "/05a_tsne_scMCA_simple.pdf"), width=12, height=10, useDingbats=FALSE)
    print(DimPlot(UMI, reduction="tsne", group.by="scMCA_simple", label=T, pt.size=1) +
            ggplot2::ggtitle(paste0("scMCA simple - "," single cell annotation")) +
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
            ggplot2::guides(col=ggplot2::guide_legend(ncol=1)))
    invisible(dev.off())
    cat(paste0(silver("Plots saved in: ")),bold(silver("\\05.Annotation\\")), silver("folder, with the prefix 05a \n"))
    ### plot each cell type
    clean_labels <- UMI$scMCA_simple
    ultra_rare <- names(table(clean_labels)[table(clean_labels) < thresh])
    clean_labels[clean_labels%in%ultra_rare] <- "other"
    UMI$clean_labels <- clean_labels
    
    ### Plotting annotated populations localization in TSNE, UMAP and PCA. [scMCA]
    cat(bold(green("Plotting scMCA dimensional reduction graphs for each population \n")))
    #popsicler:::annotation_plot(Annot_dir, "05c_", UMI, "PCA", "clean_labels", "scMCA_CellPopulations")
    popsicler:::annotation_plot(Annot_dir, "05c_", UMI, "UMAP", "clean_labels", "scMCA_CellPopulations")
    popsicler:::annotation_plot(Annot_dir, "05c_", UMI, "TSNE", "clean_labels", "scMCA_CellPopulations")
    cat(paste0(silver("Plots saved in: ")),bold(silver("\\05.Annotation\\"), silver("folder, with the prefix 05c \n")))

    ## Custom DotPlot for immune markers
    pdf(paste0(Annot_dir, "/05d_DotPlot_Immune_Markers_scMCA.pdf"), width=width_sig, height=10, useDingbats=FALSE)
    popsicler:::DTP(UMI, marker.list, 'scMCA_simple')
    invisible(dev.off())

    UMI <- popsicler:::SR_plots("ImmGen", Ig.se, UMI, Annot_dir, cluster_res)
     
    # Plotting annotated populations localization in TSNE, UMAP and PCA. [ImmGen
    annotations <- c("ImmGen.sc.main.labels")
    cat(bold(green("Plotting ImmGen dimensional reduction graphs for each population \n")))
    for(single_annot in annotations){
      #popsicler:::annotation_plot(Annot_dir, "05e_", UMI, "PCA", single_annot, "CellPopulations")
      popsicler:::annotation_plot(Annot_dir, "05e_", UMI, "UMAP", single_annot, "CellPopulations")
      popsicler:::annotation_plot(Annot_dir, "05e_", UMI, "TSNE", single_annot, "CellPopulations")
    }
    cat(paste0(silver("Plots saved in: ")),bold(silver("\\05.Annotation\\")), silver("folder, with the prefix 05e \n"))
  }
  return(UMI)
}
