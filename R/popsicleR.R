###########################################################################################
### popsicleR
###########################################################################################

plotGene <- function(genelist, umi, dir){
  ### Density plot
  suppressWarnings({pdf(file.path(dir, paste0("01d_QC_Hist_Check.pdf")), useDingbats=FALSE)
    cat(bold(green("Plotting QC per gene Histograms \n")))
    for(gene in genelist)
    {
      if(gene %in% row.names(umi@assays$RNA@counts)) {
        expr_gene <- paste0(gene, "_expressed")
        umi@meta.data[, expr_gene] <- ifelse(GetAssayData(object=umi, slot="counts")[gene,]>0, "TRUE", "FALSE")
        plot.title <- "Density total genes"
        p1 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=get(expr_gene), fill = get(expr_gene))) + ggplot2::geom_density(size=0.5, alpha=0.2) + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
        p1 <- p1 + guides(color=guide_legend("Expressed:"), fill =guide_legend("Expressed:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)
        plot.title <- "Density total genes zoom"
        p2 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=get(expr_gene), fill = get(expr_gene))) + ggplot2::geom_density(size=0.5, alpha=0.2) + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + ggplot2::xlim(0,2000) + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
        p2 <- p2 + guides(color=guide_legend("Expressed:"), fill =guide_legend("Expressed:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)
        plot.title <- "Density total UMI"
        p3 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=get(expr_gene), fill = get(expr_gene))) + ggplot2::geom_density(size=0.5, alpha=0.2) + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
        p3 <- p3 + guides(color=guide_legend("Expressed:"), fill =guide_legend("Expressed:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)
        plot.title <- "Density total UMI zoom"
        p4 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=get(expr_gene), fill = get(expr_gene))) + ggplot2::geom_density(size=0.5, alpha=0.2) + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + ggplot2::xlim(0,5000) + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
        p4 <- p4 + guides(color=guide_legend("Expressed:"), fill =guide_legend("Expressed:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5)
        print(patchwork::wrap_plots((p1+p2)/(p3+p4)) + plot_layout(guides = 'collect')
              + plot_annotation(title = gene, theme = theme(plot.title = element_text(size = 16, hjust = 0.4, face="bold"))))
      }
    }
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("01.QC_Plots\\01d_QC_Hist_Check.pdf \n")))
  ### Scatter Plot
  suppressWarnings({pdf(file.path(dir, paste0("01e_QC_Scatter_Check.pdf")), useDingbats=FALSE)
    cat(bold(green("Plotting QC per gene Scatter plots \n")))
    for(gene in genelist)
    {
      if(gene %in% row.names(umi@assays$RNA@counts)) {

        expr_gene <- paste0(gene, "_expressed")
        umi@meta.data[, expr_gene] <- ifelse(GetAssayData(object=umi, slot="counts")[gene,]>0, "TRUE", "FALSE")

        x.zoom.genes<-2000
        x.zoom.umi<-5000
        #genes vs mt%
        vars1 <- table(umi@meta.data[, expr_gene]=="TRUE")["TRUE"][[1]]
        gp3<-FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_mt", pt.size=0.2, cells=colnames(umi)[umi@meta.data[, expr_gene]=="TRUE"], col="#00BFC4") +
          #ggplot2::ylim(0,100) +
          ggplot2::ggtitle(" ")+
          theme(text=element_text(size=10),axis.title=element_text(size=10)) + NoLegend()
        gp4<-FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_mt", pt.size=0.2, cells=colnames(umi)[umi@meta.data[, expr_gene]=="TRUE"], col="#00BFC4") +
          ggplot2::xlim(0,x.zoom.genes) +
          #ggplot2::ylim(0,100) +
          ggplot2::ggtitle(" ")+
          theme(text=element_text(size=10),axis.title=element_text(size=10))+ NoLegend()

        up1<-FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_mt", pt.size=0.2,  cells=colnames(umi)[umi@meta.data[, expr_gene]=="TRUE"], col="#00BFC4") +
          ggplot2::ggtitle(" ") +
          guides(colour = guide_legend(paste0(as.character(gene),"+")))+
          theme(text=element_text(size=10),axis.title=element_text(size=10))+ NoLegend()
        up2<-FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_mt", pt.size=0.2,  cells=colnames(umi)[umi@meta.data[, expr_gene]=="TRUE"], col="#00BFC4") +
          ggplot2::xlim(0,x.zoom.umi) +
          ggplot2::ggtitle(" ") +
          guides(colour = guide_legend(paste0(as.character(gene),"+")))+
          theme(text=element_text(size=10),axis.title=element_text(size=10))+ NoLegend()
        print(patchwork::wrap_plots(gp3+gp4+up1+up2+plot_layout(guides = 'collect'))
              + plot_annotation(title = paste0(as.character(gene),"\n(expressed in ", vars1, " cells)"),
                                theme = theme(plot.title = element_text(hjust = 0.5, face="bold", size =16))))

        gu1<-FeatureScatter(umi, feature1="nFeature_RNA", feature2="nCount_RNA", pt.size=0.2, cells=colnames(umi)[umi@meta.data[, expr_gene]=="TRUE"], col="#00BFC4") +
          ggplot2::ggtitle(" ") +
          theme(text=element_text(size=10),axis.title=element_text(size=10))+ NoLegend()

        # number of genes vs single marker expression

        mp1<-FeatureScatter(umi, feature1="nFeature_RNA", feature2=gene, pt.size=0.2, group.by=expr_gene) +
          ggplot2::ylab(paste0(as.character(gene)," counts")) +
          ggplot2::ggtitle(" ") +
          guides(colour = guide_legend(paste0(as.character(gene),"+")))+
          theme(text=element_text(size=10),axis.title=element_text(size=10))
        mp2<-FeatureScatter(umi, feature1="nFeature_RNA", feature2=gene, pt.size=0.2, group.by=expr_gene) +
          ggplot2::ylab(paste0(as.character(gene)," counts")) +
          ggplot2::xlim(0,x.zoom.genes) +
          ggplot2::ggtitle(" ") +
          guides(colour = guide_legend(paste0(as.character(gene),"+")))+
          theme(text=element_text(size=10),axis.title=element_text(size=10))
        print(patchwork::wrap_plots(mp1+mp2+gu1+ guide_area()+ plot_layout(guides = 'collect'))
              + plot_annotation(title = paste0(as.character(gene),"\n(expressed in ", vars1, " cells)"),
                                theme = theme(plot.title = element_text(hjust = 0.5, face="bold", size =16))))
      }
    }
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("01.QC_Plots\\01e_QC_Scatter_Check.pdf \n")))
  cat(paste0(cyan("\nNow check the graphs, choose your thresholds and then run")),bold(cyan("FilterPlots \n")))

}


algo_plot_features <- function(dir, data, type, plot_name, organism){

  if (type=="pca"){prefix <- "pc"}else{prefix <- type}

  pdf(file.path(paste0(dir, plot_name, "_", type, "_cell_features.pdf")), width=14, height=12, useDingbats=FALSE)
  four_plots(data, type)


  PreA <- FeaturePlot(data, pt.size=1, features="percent_ribo", reduction=type) + ggplot2::xlab(paste0(toupper(prefix), " 1")) + ggplot2::ylab(paste0(toupper(prefix), " 2"))
  PreB <- FeaturePlot(data, pt.size=1, features="percent_disso", reduction=type) + ggplot2::xlab(paste0(toupper(prefix), " 1")) + ggplot2::ylab(paste0(toupper(prefix), ' 2'))
  if (organism == "human") {
    if("MALAT1"%in%rownames(data))
      PreC <- FeaturePlot(data, pt.size=1, features="MALAT1", reduction=type, min.cutoff="q75") + ggplot2::xlab(paste0(toupper(prefix), " 1")) + ggplot2::ylab(paste0(toupper(prefix), " 2"))+
        ggplot2::ggtitle("MALAT1 top25%") +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) else PreC<-ggplot()
        if("GAPDH"%in%rownames(data))
          PreD <- FeaturePlot(data, pt.size=1, features="GAPDH", reduction=type) + ggplot2::xlab(paste0(toupper(prefix), " 1")) + ggplot2::ylab(paste0(toupper(prefix), " 2")) else PreD<-ggplot()
  } else if (organism == "mouse") {
    if("Malat1"%in%rownames(data))
      PreC <- FeaturePlot(data, pt.size=1, features="Malat1", reduction=type, min.cutoff="q75") + ggplot2::xlab(paste0(toupper(prefix), " 1")) + ggplot2::ylab(paste0(toupper(prefix), " 2")) +
        ggplot2::ggtitle("Malat1 top25%") +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) else PreC<-ggplot()
        if("Gapdh"%in%rownames(data))
          PreD <- FeaturePlot(data, pt.size=1, features="Gapdh", reduction=type) + ggplot2::xlab(paste0(toupper(prefix), " 1")) + ggplot2::ylab(paste0(toupper(prefix), " 2")) else PreD<-ggplot()
  }
  print(patchwork::wrap_plots(PreA, PreB, PreC, PreD, ncol=2))

  if ("doublets" %in% colnames(data@meta.data)){
    p1 <- DimPlot(data, reduction=type, group.by = "doublets", pt.size=1, cols=c("lightgrey","black"))+
      xlab(paste0(toupper(prefix), " 1")) + ylab(paste0(toupper(prefix), " 2"))
    p2 <- FeaturePlot(data, reduction=type, features="doublets_score", pt.size=1) +scale_colour_gradientn(colours=c("lightgrey", "red", "darkred", "black"))+
      xlab(paste0(toupper(prefix), " 1")) + ylab(paste0(toupper(prefix), " 2"))
    print(patchwork::wrap_plots(p1 | p2 + plot_layout(guides = 'collect'),nrow=2))
  }

  invisible(dev.off())
}


algo_plot_clusters<- function(dir, data, type, plot_name, res){
  pdf(file.path(paste0(dir, plot_name, "_", type,"_clusters.pdf")), width=10, height=10, useDingbats=FALSE)
  if (type=="pca"){prefix <- "pc"}else{prefix <- type}
  for (res.i in res){
    print(DimPlot(data, reduction=type, group.by=paste0("RNA_snn_res.",res.i), label=T, pt.size=1, repel=TRUE) +
            ggplot2::xlab(paste0(toupper(prefix), " 1")) +
            ggplot2::ylab(paste0(toupper(prefix), " 2")) +
            ggplot2::ggtitle(paste0("Clusters @ res ", res.i)) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  }
  invisible(dev.off())
}



SR_plots <- function(db_name, annot_db, data, directory, cluster_res) { ### aggiungere BBpar per settare diversi core
  sc_anal <- paste0(db_name,".sc.main.labels")
  var.sc <- SingleR(test=data@assays$RNA@data, ref=annot_db, labels=annot_db$label.main, method="single")
  data[[sc_anal]] <- var.sc$labels

  pdf(paste0(directory, "/04a_UMAP_", sc_anal, ".pdf"), width=12, height=10, useDingbats=FALSE)
  print(UMAPPlot(data, label=T, pt.size=1, group.by=sc_anal, repel=TRUE) +
          ggplot2::ggtitle(paste0(db_name," single cell annotation")) +
          ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
          ggplot2::guides(col=ggplot2::guide_legend(ncol=1))+ xlab("UMAP 1") + ylab("UMAP 2"))
  invisible(dev.off())

  pdf(paste0(directory, "/04a_TSNE_", sc_anal, ".pdf"), width=12, height=10, useDingbats=FALSE)
  print(DimPlot(data, reduction="tsne", group.by=sc_anal, label=T, pt.size=1, repel=T) +
          ggplot2::ggtitle(paste0(db_name," single cell annotation")) +
          ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
          ggplot2::guides(col=ggplot2::guide_legend(ncol=1))+xlab("TSNE 1") + ylab("TSNE 2"))
  invisible(dev.off())


  if(is.null(cluster_res)) {
    sel_cluster<-"seurat_clusters"
    cl_anal <- paste0(db_name,".cl.main.labels")
    var.cl <- SingleR(test=data@assays$RNA@data, ref=annot_db, labels=annot_db$label.main, method="cluster", clusters=data@meta.data[[sel_cluster]])
    data[[cl_anal]] <- paste0(data[[]][[sel_cluster]], ":", var.cl$labels[match(data[[]][[sel_cluster]], rownames(var.cl))])
    data@meta.data[[cl_anal]]<-factor(data@meta.data[[cl_anal]], levels=mixedsort(unique(data@meta.data[[cl_anal]])))

    pdf(paste0(directory, "/04b_UMAP_", paste0(db_name,".cl.main.labels"), ".pdf"), width=12, height=10, useDingbats=FALSE)

    print(UMAPPlot(data, label=T, pt.size=1, group.by=cl_anal, repel=T) +
            NoLegend() +
            ggplot2::ggtitle(paste0(db_name," cluster annotation")) +
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+ xlab("UMAP 1") + ylab("UMAP 2"))

    invisible(dev.off())

    pdf(paste0(directory, "/04b_TSNE_", paste0(db_name,".cl.main.labels"), ".pdf"), width=12, height=10, useDingbats=FALSE)
    print(DimPlot(data, reduction="tsne", group.by=cl_anal, label=T, pt.size=1, repel=T) +
            NoLegend() +
            ggplot2::ggtitle(paste0(db_name," cluster annotation")) +
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+xlab("TSNE 1") + ylab("TSNE 2"))
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

    pdf(paste0(directory, "/04b_UMAP_", paste0(db_name,".cl.main.labels"), ".pdf"), width=12, height=10, useDingbats=FALSE)
    for (cl_anal in cl_analyses) {
      print(UMAPPlot(data, label=T, pt.size=1, group.by=cl_anal, repel=T) +
              NoLegend() +
              ggplot2::ggtitle(paste0(db_name," cluster ",gsub(paste0(db_name,".cl.main.labels."),"",cl_anal)," annotation")) +
              ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+ xlab("UMAP 1") + ylab("UMAP 2"))
    }
    invisible(dev.off())

    pdf(paste0(directory, "/04b_TSNE_", paste0(db_name,".cl.main.labels"), ".pdf"), width=12, height=10, useDingbats=FALSE)
    for (cl_anal in cl_analyses) {
      print(DimPlot(data, reduction="tsne", group.by=cl_anal, label=T, pt.size=1, repel=T) +
              NoLegend() +
              ggplot2::ggtitle(paste0(db_name," cluster ",gsub(paste0(db_name,".cl.main.labels."),"",cl_anal)," annotation")) +
              ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))+xlab("TSNE 1") + ylab("TSNE 2"))
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
  list_to_plot <- lapply(to_be_plotted, FUN=function(x){
  plot <- FeaturePlot(data, features = x, reduction=dimensional_redux) +
          ggplot2::xlab(paste0(toupper(dimensional_redux), " 1")) +
          ggplot2::ylab(paste0(toupper(dimensional_redux), " 2"))
  return(plot)
  })
  print(wrap_plots(list_to_plot, ncol=4))
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
  names(markers.to.plot) <- NULL
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
           ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.5,3.5))) + # ggplot2 3.3
           ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.5,2))) + # ggplot2 3.3
           RotatedAxis() + ggplot2::xlab("marker") + ggplot2::ylab(annotation) +
           ggplot2::annotate("text", x = mrk_x, y = mrk_y, label=names(markers), angle=head.angle, hjust = 0))
}

four_plots <- function(data, redux) {
  if (redux=="pca"){prefix <- "pc"}else{prefix <- redux}
  PreA <- FeaturePlot(data, pt.size=1, features="nCount_RNA", reduction=redux) + ggplot2::xlab(paste0(toupper(prefix), " 1")) + ggplot2::ylab(paste0(toupper(prefix), " 2"))
  PreB <- FeaturePlot(data, pt.size=1, features="nFeature_RNA", reduction=redux) + ggplot2::xlab(paste0(toupper(prefix), " 1")) + ggplot2::ylab(paste0(toupper(prefix), ' 2'))
  PreC <- FeaturePlot(data, pt.size=1, features="percent_mt", reduction=redux) + ggplot2::xlab(paste0(toupper(prefix), " 1")) + ggplot2::ylab(paste0(toupper(prefix), " 2"))
  PreD <- DimPlot(data, pt.size=1, group.by="Phase", reduction=redux) + ggplot2::xlab(paste0(toupper(prefix), " 1")) + ggplot2::ylab(paste0(toupper(prefix), " 2")) +
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
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) + NoLegend() +
            ggplot2::xlab(paste0(redux, " 1")) + ggplot2::ylab(paste0(redux, " 2")))
  }
  invisible(dev.off())
}

#' PrePlots
#'
#' @description
#'
#' Perfoms preprossing analyses and quality controls on the provided dataset.
#' Reports information on the dataset and outputs plots describing quality of the data to be further analyzed.
#'
#' @usage
#'
#' PrePlots(sample_name, input_data, genelist=NULL, percentage=0.1, gene_filter=200, cellranger=TRUE, organism=c("human","mouse"), out_folder=getwd())
#'
#' @param sample_name Name of the sample investigated
#' @param input_data Input data for popsicleR. When data are generated via Cell Ranger input_data must be a path to a folder containing the approriate output. Otherwise, input_data must be a path to a .txt matrix structured with gene names in the first column and cell names in the first row.
#' @param genelist List of gene for singular gene plots. By default is NULL and a proprietary genelist is used to generate plots. A personal genelist can be fed as characters vector
#' @param percentage Minimum percentage of cells in which a gene must be expressed to be retained for the subsequent analysis. Default is 0.1
#' @param gene_filter Minimum number of genes detected to retained a cell for the subsequent analysis. Default is 200
#' @param cellranger Specify if data were generated via Cell Ranger. Default is TRUE
#' @param organism The organism on which perfom the analysis. Can be human or mouse.
#' @param out_folder Output folder. Default is working directory
#'
#' @details
#' PrePlots returns several graphs in the "01.QC_Plots" folder:
#'
#' violin (1a), density (1b), and scatter (1c) plots on number of genes (nFeature_RNA), number of UMI (nCount_RNA), the percentage of reads that map to the mitochondrial genome (percent_mt), to ribosomal (percent_ribo) and dissociation genes (percent_disso).
#'
#' PrePlots also returns the distributions (1d) and the cell abundance (1e) of specific cell populations identified by a list of user-defined marker genes.
#'
#' @return Returns a Seurat Object generated from input data.
#'
#' Genes expressed in less than percentage of cells the value specified by the "percentage" parameter and cells that express less than number of genes specified by "gene_filter" parameter are filtered-out from input data.
#'
#' @examples
#'
#' PrePlots('breast_single_cell', input_data="/path/to/cellranger_output_directory/", genelist=c('TP53','PTEN'), percentage=0.1, gene_filter=200, cellranger = TRUE, organism="human")
#'
#' PrePlots("mouse_sample", input_data="/path/to/matrix.txt", genelist=c('Tp53','Pten'), cellranger=FALSE, organism= "mouse")
#'
#' @author Jimmy Caroli, Francesco Grandi
#'
#'
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom magrittr "%>%"
#' @export
#'
PrePlots <- function(sample_name, input_data, genelist=NULL, percentage=0.1, gene_filter=200, cellranger=TRUE, organism=c("human","mouse"), out_folder=getwd()){
  ### Root directory
  root <- out_folder
  if (!file.exists(root)){dir.create(root, recursive=T)}
  ### QC Plots directory
  QC_dir <- paste0(root,"/01.QC_Plots/")
  if (!file.exists(QC_dir)){dir.create(QC_dir, recursive=T)}
  ### PreProcessing directory
  ###########################
  ###   dataset loading   ###
  ###########################
  if (cellranger == TRUE){
    umi_data <- Read10X(data.dir=input_data) #check if CellRanger output
  } else {
    umi_data <- read.table(file=input_data, header=T, sep="\t", row.names=1)
  } #else stop(red("Please check your input file. Must be a CellRanger output or a .txt matrix with gene names in the first column and cell names in the first row \n"))
  #### creation of a Seurat object
  # Keep all genes expressed in >= ~0.1% of cells
  # Keep all cells with at least 200 detected genes
  n_cells <- round((ncol(umi_data)*percentage)/100)
  n_genes <- gene_filter
  # min.cells: Include features detected in at least this many cells
  # min.features: Include cells where at least this many features are detected
  umi <- CreateSeuratObject(counts=umi_data, min.cells=n_cells, min.features=n_genes, project=sample_name)
  Starting_Cells <- ncol(umi)
  #####################################
  ###	  calculation of QC metrics   ###
  #####################################
  ### calculate the percentage of mitochondria, ribosomal and dissociation genes
  dissociation_genes <- c("ACTG1","ANKRD1","ARID5A","ATF3","ATF4","BAG3","BHLHE40","BRD2","BTG1","BTG2","CCNL1","CCRN4L","CEBPB","CEBPD","CEBPG","CSRNP1","CXCL1","CYR61","DCN","DDX3X","DDX5","DES","DNAJA1","DNAJB1","DNAJB4","DUSP1","DUSP8","EGR1","EGR2","EIF1","EIF5","ERF","ERRFI1","FAM132B","FOS","FOSB","FOSL2","GADD45A","GADD45G","GCC1","GEM","H3F3B","HIPK3","HSP90AA1","HSP90AB1","HSPA1A","HSPA1B","HSPA5","HSPA8","HSPB1","HSPE1","HSPH1","ID3","IDI1","IER2","IER3","IER5","IFRD1","IL6","IRF1","IRF8","ITPKC","JUN","JUNB","JUND","KCNE4","KLF2","KLF4","KLF6","KLF9","LITAF","LMNA","MAFF","MAFK","MCL1","MIDN","MIR22HG","MT1","MT2","MYADM","MYC","MYD88","NCKAP5L","NCOA7","NFKBIA","NFKBIZ","NOP58","NPPC","NR4A1","ODC1","OSGIN1","OXNAD1","PCF11","PDE4B","PER1","PHLDA1","PNP","PNRC1","PPP1CC","PPP1R15A","PXDC1","RAP1B","RASSF1","RHOB","RHOH","RIPK1","SAT1","SBNO2","SDC4","SERPINE1","SKIL","SLC10A6","SLC38A2","SLC41A1","SOCS3","SQSTM1","SRF","SRSF5","SRSF7","STAT3","TAGLN2","TIPARP","TNFAIP3","TNFAIP6","TPM3","TPPP3","TRA2A","TRA2B","TRIB1","TUBB4B","TUBB6","UBC","USP2","WAC","ZC3H12A","ZFAND5","ZFP36","ZFP36L1","ZFP36L2","ZYX")
  if(organism == 'human'){
    dissociation_genes <- paste('^', paste(dissociation_genes, collapse='$|^'), '$', sep = '')
    if (is.null(genelist)) genelist<-c("EPCAM","VIM", "COL1A1", "PECAM1", "PTPRC", "CD3D", "CD14")
    umi[["percent_mt"]] <- PercentageFeatureSet(umi, pattern = "^MT-")
    umi[["percent_ribo"]] <- PercentageFeatureSet(umi, pattern = '^RPL|^RPS|^MRPL|^MRPS')
    umi[["percent_disso"]] <- PercentageFeatureSet(umi, pattern = dissociation_genes)
  } else if(organism == 'mouse') {
    dissociation_genes <- tools::toTitleCase(tolower(dissociation_genes))
    dissociation_genes <- paste('^', paste(dissociation_genes, collapse='$|^'), '$', sep = '')
    if (is.null(genelist)) genelist<-c("Epcam","Vim", "Col1a1", "Pecam1", "Ptprc", "Cd3d", "Cd14")
    umi[["percent_mt"]] <- PercentageFeatureSet(umi, pattern = "^mt-")
    umi[["percent_ribo"]] <- PercentageFeatureSet(umi, pattern ='^Rpl|^Rps|^Mrpl|^Mrps')
    umi[["percent_disso"]] <- PercentageFeatureSet(umi, pattern = dissociation_genes)
  } else {stop("organism must be human or mouse")}
  ### Violin Plot on number of genes, number of UMI and fraction of mitochondrial genes
  cat(bold(green("\nPlotting QC Violin plots \n")))
  suppressWarnings({pdf(paste0(QC_dir, "01a_QC_violin_plots.pdf"), width=24, useDingbats=FALSE)
    vln1 <- popsicleR:::VLN(umi, 10, feats="nFeature_RNA", colours= "tomato", 0.01)+ NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggplot2::xlab("")
    vln2 <- popsicleR:::VLN(umi, 10, feats="nCount_RNA", colours= "tomato", 0.01)+ NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggplot2::xlab("")
    vln3 <- popsicleR:::VLN(umi, 10, feats="percent_mt", colours= "dodgerblue", 0.01)+ NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggplot2::xlab("")
    vln4 <- popsicleR:::VLN(umi, 10, feats="percent_ribo", colours= "yellow", 0.01)+ NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggplot2::xlab("")
    vln5 <- popsicleR:::VLN(umi, 10, feats="percent_disso", colours= "forestgreen", 0.01)+ NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5)) + ggplot2::xlab("")
    print(patchwork::wrap_plots(vln1 | vln2 | vln3 | vln4 | vln5 + plot_layout(guides = 'collect') + NoLegend()) +
            plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))))
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("01.QC_Plots\\01a_QC_violin_plots.pdf \n")))
  ### Density plot
  cat(bold(green("Plotting QC Density plots \n")))
  suppressWarnings({pdf(file.path(QC_dir, "01b_QC_Hist_nGene_nUMI_MTf_Ribo.pdf"), useDingbats=FALSE)
    plot.title <- "Density total genes"
    nGene <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    plot.title <- "Density total genes zoom"
    nGene_zoom <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title) + ggplot2::xlim(0,2000) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    plot.title <- "Density total UMI"
    nUMI <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    plot.title <- "Density total UMI zoom"
    nUMI_zoom <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title) + ggplot2::xlim(0,5000) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    print(patchwork::wrap_plots(nGene + nGene_zoom + nUMI + nUMI_zoom + plot_layout(guides = 'collect')+ NoLegend()) +
            plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))) + NoLegend())
    plot.title <- "Density mitochondrial fraction"
    mt_fraction <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_mt, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="blue", fill="dodgerblue") + ggplot2::ggtitle(plot.title) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    plot.title <- "Density ribosomal fraction"
    ribosomal_fraction <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_ribo, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="darkgoldenrod3", fill="yellow") + ggplot2::ggtitle(plot.title) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    plot.title <- "Density dissociation fraction"
    dissociation_fraction <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_disso, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="darkgreen", fill="forestgreen") + ggplot2::ggtitle(plot.title) + guides(colour = guide_legend("Sample ID:"), fill = guide_legend("Sample ID:")) + theme(legend.title = element_text(face = "bold"),legend.title.align = 0.5) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    print(patchwork::wrap_plots(mt_fraction / (ribosomal_fraction + dissociation_fraction) + plot_layout(guides = 'collect') + NoLegend()) +
            plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))) + NoLegend())
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("01.QC_Plots\\01b_QC_Hist_nGene_nUMI_MTf_Ribo.pdf \n")))
  ### Scatter Plot
  cat(bold(green("Plotting QC Scatter plots \n")))
  suppressWarnings({pdf(file.path(QC_dir, "01c_QC_Scatter_nGene_nUMI_MTf.pdf"), width=18, height=12, useDingbats=FALSE)
    plotA <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="nCount_RNA", pt.size=0.3, cols="tomato") +  theme(plot.title=element_blank()) + NoLegend()
    plotB <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_mt", pt.size=0.3, cols="dodgerblue") +  theme(plot.title=element_blank()) + NoLegend()
    plotC <- FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_mt", pt.size=0.3, cols="dodgerblue") +  theme(plot.title=element_blank()) + NoLegend()
    plotD <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_ribo", pt.size=0.3, cols="yellow") +  theme(plot.title=element_blank())+ NoLegend()
    plotE <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_disso", pt.size=0.3, cols="forestgreen") +  theme(plot.title=element_blank())+ NoLegend()
    print(patchwork::wrap_plots(plotA + ggExtra::ggMarginal(plotB, type="density", color="blue", fill="dodgerblue") + ggExtra::ggMarginal(plotC, type="density", color="blue", fill="dodgerblue") + plot_spacer() + ggExtra::ggMarginal(plotD, type="density", color="darkgoldenrod3", fill="yellow") + ggExtra::ggMarginal(plotE, type="density", color="darkgreen", fill="forestgreen") + plot_layout(guides = 'collect')+ NoLegend()))
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("01.QC_Plots\\01c_QC_Scatter_nGene_nUMI_MTf.pdf \n")))
  #
  ###########################################################
  ###   check for presence and number of selected genes   ###
  ###########################################################
  #
  #umi@meta.data$[, paste0(opt$gene, '_expressed')]
  if(length(genelist)!=0){
    #for (gene in genelist)
    #{
    plotGene(genelist, umi, QC_dir)
    # }
  }
  ### write log
  write.table(t(c(as.character(as.POSIXct(Sys.time())),"PrePlots:","percentage",percentage,"gene_filter", gene_filter)),file=file.path(out_folder,"popsicleR.log"), sep="\t", row.names=F, col.names=F, quote=F, append=T)
  return(umi)
}

#' FilterPlots
#'
#' @description
#'
#' Apply filters to input data. Creates plot with provided thresholds as input.
#'
#' @usage
#'
#' FilterPlots(UMI, G_RNA_low = 0, G_RNA_hi = Inf, U_RNA_low = 0, U_RNA_hi = Inf, percent_mt_hi = 100, percent_ribo_hi = 100, percent_disso_hi = 100, out_folder=getwd())
#'
#' @param UMI Input UMI object generated by PrePlots function
#' @param G_RNA_low Minimum number of genes detected in each cell. Default is 0
#' @param G_RNA_hi Maximum number of genes detected in each cell. Default is Infinite
#' @param U_RNA_low Minimum number of molecules detected within a cell. Default is 0
#' @param U_RNA_hi Maximum number of molecules detected within a cell. Default is Infinite
#' @param percent_mt_hi Maximum percentage of mitochondrial genes. Default is 100
#' @param percent_ribo_hi Maximum percentage of ribosomal genes. Default is 100
#' @param percent_disso_hi Maximum percentage of dissociation genes. Default is 100
#' @param out_folder Output folder. Default is working directory
#'
#' @details
#'
#' FilterPlots returns two graphs in the "01.QC_Plots" folder:
#'
#' the first plot (01f) displays data distributions after filtering; filtering thresholds are displayed as gray lines.
#'
#' the second graph (01g) shows data correlations, that highlight with colors, cells passing the filtering procedure; filtering thresholds are displayed as gray lines.
#'
#' @return Returns a Seurat Object filtered according to the user-defined set of thresholds.
#'
#' @examples
#'
#' FilterPlots(umi_object)
#'
#' FilterPlots(umi_object, G_RNA_low = 500, G_RNA_hi = 10000, U_RNA_low = 1000, U_RNA_hi = 70000, percent_mt_hi = 20, percent_ribo_hi = 50, percent_disso_hi = 10)
#'
#' @author Jimmy Caroli, Francesco Grandi
#'
#' @export
#'

FilterPlots <- function(UMI, G_RNA_low = 0, G_RNA_hi = Inf, U_RNA_low = 0, U_RNA_hi = Inf, percent_mt_hi = 100, percent_ribo_hi = 100, percent_disso_hi = 100, out_folder=getwd()){
  #######################################
  ###   define and apply thresholds   ###
  #######################################
  QC_dir <- paste0(out_folder,"/01.QC_Plots/")
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
  suppressWarnings({pdf(paste0(QC_dir, "/01f_final_Hist_plots.pdf"), useDingbats=FALSE)
    plot.title <- "Density total genes"
    final <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) +
      genes.dn.lim + genes.up.lim + NoLegend() +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    plot.title <- "Density total genes zoom"
    final1 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title)  + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + ggplot2::xlim(0,2000) + theme(axis.title.y = element_blank()) +
      genes.dn.lim + genes.up.lim + NoLegend() +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    plot.title <- "Density total UMI"
    final2 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title)  + theme(plot.title = element_text(hjust = 0.5, face ="bold")) +
      umi.dn.lim + umi.up.lim + NoLegend() +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    plot.title <- "Density total UMI zoom"
    final3 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=nCount_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="red", fill="tomato") + ggplot2::ggtitle(plot.title)  + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + ggplot2::xlim(0,5000) + theme(axis.title.y = element_blank()) +
      umi.dn.lim + umi.up.lim + NoLegend() +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
    print(patchwork::wrap_plots(final + final1 + final2 + final3 + plot_layout(guides = 'collect')))
    plot.title <- "Density mitochondrial fraction"
    final4 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_mt, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="blue", fill="dodgerblue") + ggplot2::ggtitle(plot.title) +
      theme(plot.title = element_text(hjust = 0.5, face ="bold")) + NoLegend() + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + mt.dn.lim
    plot.title <- "Density ribosomal fraction"
    final5 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_ribo, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="darkgoldenrod3", fill="yellow") + ggplot2::ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + NoLegend() +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + ribo.dn.lim
    plot.title <- "Density dissociation fraction"
    final6 <- ggplot2::ggplot(umi@meta.data, ggplot2::aes(x=percent_disso, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="darkgreen", fill="forestgreen") + ggplot2::ggtitle(plot.title) +
      theme(plot.title = element_text(hjust = 0.5, face ="bold")) + NoLegend() + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank()) + disso.dn.lim
    print(patchwork::wrap_plots(final4/(final5 + final6) + plot_layout(guides = 'collect')))})
  invisible(dev.off())
  cat(paste0(silver("\nPlots saved in: ")),bold(silver("01.QC_Plots\\01f_final_Hist_plots.pdf \n")))

  cat(bold(green("Plotting QC final scatter plots \n")))
  suppressWarnings({pdf(paste0(QC_dir, "/01g_final_Scatter_plots.pdf"),width=18, height=12, useDingbats=FALSE)
    plot1 <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="nCount_RNA", pt.size = 0.3, group.by = "filtered")+ scale_color_manual(values = c("tomato", "#666666")) +  theme(plot.title=element_blank()) +
      genes.dn.lim + genes.up.lim + umi.dn.lim.y + umi.up.lim.y + NoLegend()
    plot2 <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_mt", pt.size = 0.3, group.by = "filtered")+ scale_color_manual(values = c("dodgerblue", "#666666")) + theme(plot.title=element_blank()) +
      genes.dn.lim + genes.up.lim + mt.dn.lim.y + NoLegend()
    plot3 <- FeatureScatter(umi, feature1="nCount_RNA", feature2="percent_mt", pt.size = 0.3, group.by = "filtered")+ scale_color_manual(values = c("dodgerblue", "#666666")) + theme(plot.title=element_blank()) +
      umi.dn.lim + umi.up.lim + mt.dn.lim.y + NoLegend()
    plot4 <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_ribo", pt.size = 0.3, group.by = "filtered")+ scale_color_manual(values = c("yellow", "#666666")) + theme(plot.title=element_blank()) +
      genes.dn.lim + genes.up.lim + ribo.dn.lim.y  + NoLegend()
    plot5 <- FeatureScatter(umi, feature1="nFeature_RNA", feature2="percent_disso", pt.size = 0.3, group.by = "filtered")+ scale_color_manual(values = c("forestgreen", "#666666")) + theme(plot.title=element_blank()) +
      genes.dn.lim + genes.up.lim + disso.dn.lim.y  + NoLegend()
    print(patchwork::wrap_plots(plot1 + plot2 + plot3 + plot_spacer() + plot4 + plot5+ plot_layout(guides = 'collect') + NoLegend()))
    invisible(dev.off())})
  cat(paste0(silver("Plots saved in: ")),bold(silver("01.QC_Plots\\01g_final_Scatter_plots.pdf \n")))
  umi <- subset(umi, subset = nFeature_RNA > G_RNA_low &
                  nFeature_RNA < G_RNA_hi &
                  nCount_RNA   > U_RNA_low &
                  nCount_RNA   < U_RNA_hi &
                  percent_mt   < percent_mt_hi &
                  percent_ribo < percent_ribo_hi &
                  percent_disso < percent_disso_hi)

  cat(paste0(cyan("\nNext suggested step is Doublets Calculation, run")),bold(cyan("CalculateDoublets \n")))
  cat(paste0(bold(cyan("\nWARNING: \n")),cyan("It is recommended to first run CalculateDoublets step setting "), bold(cyan("dbs_thr ='none' ")),cyan("or "),bold(cyan("dbs_rate =NULL ")),cyan("and "),  bold(cyan("dbs_remove= FALSE. \n")), cyan("Once checked the graphs it is possible to re-run this step specifying a custom threshold \nthrough the 'dbs_thr' or 'dbs_rate' parameter and removing all the cells identified as doublets setting 'dbs_remove' parameter as TRUE. \n")))

  ### write log
  write.table(t(c(as.character(as.POSIXct(Sys.time())),"FilterPlots:","G_RNA_low",G_RNA_low,"G_RNA_hi", G_RNA_hi, "U_RNA_low",U_RNA_low,"U_RNA_hi",U_RNA_hi,"percent_mt_hi",percent_mt_hi,"percent_ribo_hi",percent_ribo_hi,"percent_disso_hi",percent_disso_hi)),file=file.path(out_folder,"popsicleR.log"), sep="\t", row.names=F, col.names=F, quote=F, append=T)

  return(umi)
}

#' CalculateDoublets
#'
#' @description
#' Calculates, highlights and optionally removes doublets in the provided dataset.
#' This function can be run twice sequentially, first time to explore doublets graphs and subsequently to set the appropriate threshold for doublets identification.
#'
#' @usage
#'
#' CalculateDoublets(UMI, method=c("scrublet","scDblFinder"), dbs_thr='none', dbs_rate=NULL, dbs_remove=TRUE, out_folder=getwd())
#'
#' @param UMI Input UMI object generated either via previous popsicleR funtions or CalculateDoublets function
#' @param method Method to use for doublets calculation. It can be "scrublet" or "scDblFinder". Default is "scrublet"
#' @param dbs_thr Doublets threshold. Default is auto calculated. Used only if method = "scrublet"
#' @param dbs_rate Expected doublets rate. Default is auto calculated. Used only if method = "scDblFinder"
#' @param dbs_remove Remove doublets if TRUE. Default is TRUE
#' @param out_folder Output folder. Default is the working directory
#'
#' @details
#'
#' CalculateDoublets returns two graphs in the "01.QC_Plots" folder:
#'
#' the first (01h) allows the visualization of the cell doublets and of the doublet scores on a UMAP projection.
#'
#' the second (01i), computed only if method = "scrublet", returns histograms of observed and simulated doublet score distributions and the value of the automatically detected threshold for the simulated doublet score.
#'
#' @return Returns a Seurat Object after doublets calculation.
#'
#' If "dbs_remove = TRUE", cells labeled as doublets are removed from the object.
#'
#' @references
#'
#' The main function used to identify doublets with Scrublet has been modified from the R code provided at https://rdrr.io/github/ChengxiangQiu/rscrublet/ by Chengxiang Qiu \email{<cxqiu@uw.edu>} based on the method implemented in the python module Scrublet
#'
#' Samuel L. Wolock, Romain Lopez, Allon M. Klein, "Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data" (2019). \emph{Cell System}, Volume \bold{8},--281-291--.e9,ISSN 2405-4712,\url{https://doi.org/10.1016/j.cels.2018.11.005.}
#'
#' Germain PL, Lun A, Macnair W and Robinson MD. Doublet identification in single-cell sequencing data using scDblFinder. \emph{F1000Research} (2021), 10:979 \url{https://doi.org/10.12688/f1000research.73600.1}
#'
#' @examples
#'
#' ## first run using scrublet:
#' umi_object <- CalculateDoublets(umi_object, method = "scrublet", dbs_thr='none', dbs_remove=FALSE)
#'
#' ## second run using scrublet:
#' CalculateDoublets(umi_object, method = "scrublet", dbs_thr=0.22, dbs_remove=TRUE)
#'
#' ## first run using scDblFinder:
#' umi_object <- CalculateDoublets(umi_object, method = "scDblFinder", dbs_rate=NULL, dbs_remove=FALSE)
#'
#' ## second run using scDblFinder:
#' CalculateDoublets(umi_object, method = "scDblFinder", dbs_rate=0.22, dbs_remove=TRUE)
#'
#' @author Jimmy Caroli, Francesco Grandi, Chengxiang Qiu
#'
#' @export
#'

CalculateDoublets <- function(UMI, method=c("scrublet","scDblFinder"), dbs_thr='none', dbs_rate=NULL, dbs_remove=TRUE, out_folder=getwd()){
  QC_dir <- paste0(out_folder,"/01.QC_Plots/")
  suppressWarnings(if (!file.exists(QC_dir)){dir.create(QC_dir, recursive=T)})
  method<- match.arg(method)

  if(method=="scrublet"){
  ### calculate doublets with scrublet

  if(dbs_thr == 'none'| !"scrublet_score" %in% colnames(UMI@meta.data)){
    doublets <- scrubDoublets(as.matrix(UMI@assays$RNA@counts), directory=QC_dir, expected_doublet_rate=0.1)
    names(doublets) <- c("predicted", "score_predicted", "score_simulated")
    ### if doublets exist, calculate, set to 0 otherwise
    ifelse(TRUE %in% doublets$predicted, dbs_found <- table(doublets$predicted)[[2]], dbs_found <- 0)
    thr <- round(mean(c(min(doublets$score_predicted[doublets$predicted]),max(doublets$score_predicted[!doublets$predicted]))), 3)
    ### Threshold definition
    doublets$predicted <- c(doublets$score_predicted > thr)
    ### add features to umi object
    UMI$scrublet_score <- doublets$score_predicted
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
    pdf(paste0(QC_dir,"/01h_doublets_umap.pdf"),15,8, useDingbats=FALSE)
    p1 <- DimPlot(UMI2, reduction="umap", group.by = "doublets", pt.size=0.5, cols=c("lightgrey"), cells.highlight = highlight_labels, cols.highlight = "black")+ ggplot2::xlab("UMAP 1") +ggplot2::ylab("UMAP 2")
    p2 <- FeaturePlot(UMI2, reduction="umap", features="doublets_score", pt.size=0.5) +scale_colour_gradientn(colours=c("lightgrey", "red", "darkred", "black"))+ ggplot2::xlab("UMAP 1") +ggplot2::ylab("UMAP 2")
    print(patchwork::wrap_plots(p1 | p2 + plot_layout(guides = 'collect')))
    dev.off()
    cat(paste0(silver("Plots saved in: ")),bold(silver("01.QC_Plots\\01h_doublets_umap.pdf \n")))
    cat(paste0(cyan("Once checked the graphs it is possible to re-run this step specifying a custom threshold through the 'dbs_thr' parameter and removing all the cells identified as doublets setting 'dbs_remove' parameter as TRUE. \n")))
    if(dbs_remove == TRUE){
      cat(paste0(green("Removing", sum(UMI$doublets) ,"Doublets \n")))
      UMI <- UMI[,!UMI$doublets]
    }
  }else{
    UMI$doublets <- ifelse(UMI$scrublet_score > dbs_thr, TRUE, FALSE)
    if(dbs_remove == TRUE){
      cat(paste0(green("Removing", sum(UMI$doublets) ,"Doublets \n")))
      UMI <- UMI[,!UMI$doublets]
      cat(paste0(cyan("Next suggested step is data normalization, run")),bold(cyan("Normalize \n")))
    }
  }
  } else {
    ### calculate doublets with scDblFinder
    sample.sce <- as.SingleCellExperiment(UMI)
    set.seed(123)
    sample.sce <- scDblFinder(sample.sce, dbr=dbs_rate)
    UMI$scDblFinder_score <- sample.sce$scDblFinder.score
    UMI$doublets_score <- sample.sce$scDblFinder.score
    UMI$doublets <- sample.sce$scDblFinder.class
    UMI$doublets <- UMI$doublets=="doublet"

    UMI2 <- NormalizeData(UMI, verbose=FALSE)
    UMI2 <- FindVariableFeatures(UMI2, selection.method="vst", nfeatures=2000, verbose=FALSE)
    UMI2 <- ScaleData(UMI2, vars.to.regress=NULL, verbose=FALSE)
    UMI2 <- RunPCA(UMI2, n_pcs=30, verbose=FALSE)
    UMI2 <- suppressWarnings(RunUMAP(UMI2, dims=1:10, verbose=FALSE))
    Idents(UMI2) <- "doublets"
    highlight_labels <- list("doublet"= WhichCells(UMI2, idents = TRUE), "singlet"= WhichCells(UMI2, idents = FALSE))
    cat(bold(green("Plotting doublets UMAP \n")))
    pdf(paste0(QC_dir,"/01h_doublets_umap.pdf"),15,8, useDingbats=FALSE)
    p1 <- DimPlot(UMI2, reduction="umap", group.by = "doublets", pt.size=0.5, cols=c("lightgrey"), cells.highlight = highlight_labels, cols.highlight = "black")+ ggplot2::xlab("UMAP 1") +ggplot2::ylab("UMAP 2")
    p2 <- FeaturePlot(UMI2, reduction="umap", features="doublets_score", pt.size=0.5) +scale_colour_gradientn(colours=c("lightgrey", "red", "darkred", "black"))+ ggplot2::xlab("UMAP 1") +ggplot2::ylab("UMAP 2")
    print(patchwork::wrap_plots(p1 | p2 + plot_layout(guides = 'collect')))
    dev.off()
    cat(paste0(silver("Plots saved in: ")),bold(silver("01.QC_Plots\\01h_doublets_umap.pdf \n")))
    cat(paste0(cyan("Once checked the graphs it is possible to re-run this step specifying a custom doublet rate through the 'dbs_rate' parameter and removing all the cells identified as doublets setting 'dbs_remove' parameter as TRUE. \n")))
    if(dbs_remove == TRUE){
    cat(paste0(green("Removing", sum(UMI$doublets) ,"Doublets \n")))
    UMI <- UMI[,!UMI$doublets]
    }
   cat(paste0(cyan("Next suggested step is data normalization, run")),bold(cyan("Normalize \n")))
   }
  ### write log
  dbs_rate<-ifelse(is.null(dbs_rate),"NULL",dbs_rate)
  write.table(t(c(as.character(as.POSIXct(Sys.time())),"CalculateDoublets:","method",method,"dbs_thr", dbs_thr, "dbs_rate",dbs_rate,"dbs_remove",dbs_remove)),file=file.path(out_folder,"popsicleR.log"), sep="\t", row.names=F, col.names=F, quote=F, append=T)

  return(UMI)
}

#' Normalize
#'
#' @description
#'
#' Perform data normalization on the provided input. User must define the number of variable genes to be investigated
#' after normalization.
#'
#' @usage
#'
#' Normalize(UMI, variable_genes=2000, out_folder=getwd())
#'
#' @param UMI Input UMI object generated via previous popsicleR functions
#' @param variable_genes Number of highly variable genes investigated after normalization. Default is 2000
#' @param out_folder Output folder. Default is the working directory
#'
#' @details
#' Normalize returns two graphs in the "02.PreProcessing" folder:
#'
#' the first graph (02a) shows expression before and after the normalization step,
#'
#' while the second (02b) highlights the highly variable genes.
#'
#' @return Returns a Seurat Object after normalization process.
#'
#' @examples
#'
#' Normalize(seurat_object, variable_genes=2000)
#'
#' @author Jimmy Caroli, Francesco Grandi
#'
#' @export
#'

Normalize <- function(UMI, variable_genes=2000, out_folder=getwd()){
  PP_dir <- paste0(out_folder,"/02.PreProcessing/")
  if (!file.exists(PP_dir)){dir.create(PP_dir, recursive=T)}
  suppressWarnings({umi <- NormalizeData(object=UMI, normalization.method="LogNormalize", scale.factor=1e4)
  cat(bold(green("\nPlotting Normalization graphs \n")))
  pdf(paste0(PP_dir, "/02a_total_expression_after_before_norm.pdf"), useDingbats=FALSE)
  par(mfrow = c(2,1))
  hist(colSums(as.matrix(umi@assays$RNA@counts)), breaks=100, main="Total expression before normalization", xlab="Sum of expression")
  hist(colSums(as.matrix(umi@assays$RNA@data)), breaks=100, main="Total expression after normalization", xlab="Sum of expression")
  invisible(dev.off())
  cat(paste0(silver("Plots saved in: ")),bold(silver("02.PreProcessing\\02a_total_expression_after_before_norm.pdf \n")))
  umi <- FindVariableFeatures(umi, selection.method = "vst", nfeatures = variable_genes)
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(umi), 10)
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(umi)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # When using repel, set xnudge and ynudge to 0 for optimal results
  cat(bold(green("Plotting High Variables Genes \n")))
  pdf(paste0(PP_dir, "/02b_plot_FindVariableGenes.pdf"), useDingbats=FALSE)
  print(patchwork::wrap_plots(plot1 / plot2 + plot_layout(guides = 'collect')), ncol=1, nrow=2)
  cat(paste0(silver("Plots saved in: ")),bold(silver("02.PreProcessing\\02b_plot_FindVariableGenes.pdf \n")))
  invisible(dev.off())})
  cat(paste0(cyan("Next suggested step is regression, run ")),bold(cyan("ApplyRegression")),cyan("\nIt is suggested to apply regression without providing any variables to regress and without exploring PCs (explore_PC=FALSE) to save time and computational effort.\nRegression variables can be chosen through visual inspection of this function outputs \nOnce selected, set explore_PC as TRUE and explore graphs to identify the PCs number to use in your analysis \n"))
  ### write log
  write.table(t(c(as.character(as.POSIXct(Sys.time())),"Normalize:","variable_genes",variable_genes)),file=file.path(out_folder,"popsicleR.log"), sep="\t", row.names=F, col.names=F, quote=F, append=T)

  return(umi)
}

#' ApplyRegression
#'
#' @description
#' Perform data scaling and optionally regression on normalized data. User can specify variable/s on which calculate the regression.
#' Before scaling cell cycle scores are calculated. When regression is required, this function should be run twice sequentially, first time to explore variables distribution and subsequently to regress on appropriate variables.
#'
#' @usage
#'
#' ApplyRegression(UMI, organism=c("human","mouse"), variables='none', explore_PC=FALSE, out_folder=getwd())
#'
#' @param UMI Input UMI object, generated via Normalize function
#' @param organism The organism on which perfom the analysis. Can be human or mouse.
#' @param variables String or character vector specifying variable/s on which perform regression. Commonly, variables used for regression are: nFeature_RNA, nCount_RNA, percent_mt, S.Score, G2M.Score. Default is 'none' as to not apply any regression
#' @param explore_PC When TRUE performs additional investigations and generates further exploration graphs. Default is FALSE.
#' @param out_folder Output folder. Default is the working directory
#'
#' @details
#'
#' ApplyRegression returns several graphs in dedicated subfolders inside the "02.PreProcessing" directory. When no regression variables are specified, a "No_Regression" subfolder is generated. Otherwise, a dedicated subfolder is generated for each combination of regression variables.
#'
#' graph 02c, returns projections on reduced spaces (PCA, t-SNE, and UMAP using 20 principal components) to observe variables distribution.
#'
#' Furthermore, when explore_PC=TRUE graphs from 02d to 02g are generated in order to choose the optimal number of principal components to consider for further analysis.
#'
#' graph 02d displays top genes associated with the first two principal components.
#'
#' graph 02e shows a heatmap focusing on a principal component. Both cells and genes are sorted by their loadings for each principal component.
#'
#' graph 02f returns the results of the JackStraw analysis for PCA significance.
#'
#' graph 02g reports the standard deviation observed for each principal component.
#'
#' @return Returns a Seurat Object after regression.
#'
#' @examples
#'
#' ApplyRegression(UMI= umi_object, organism= "human", variables= 'none', explore_PC=FALSE)
#'
#' ApplyRegression(UMI= umi_object, organism= "human", variables= c("nFeature_RNA","nCount_RNA", "percent_mt", "S.Score", "G2M.Score"), explore_PC=TRUE)
#'
#' @author Jimmy Caroli, Francesco Grandi
#'
#' @export
#'

ApplyRegression <- function(UMI, organism=c("human","mouse"), variables='none', explore_PC=FALSE,  out_folder=getwd()){

  if(all(variables!='none') & !all(variables %in% colnames(UMI@meta.data))) {stop(red("ERROR: Specified variables are not in metadata field"))}
  PP_dir <- file.path(out_folder,"02.PreProcessing")
  suppressWarnings(if (!file.exists(PP_dir)){dir.create(PP_dir, recursive=T)})

  if(!all(c("S.Score", "G2M.Score")%in%colnames(UMI@meta.data))){
    cat(bold(green("Calculating Cell Cycle Score \n")))
    if(organism == 'human') {
      cc.genes$s.genes <- intersect(cc.genes$s.genes, row.names(UMI@assays$RNA@counts))
      cc.genes$g2m.genes <- intersect(cc.genes$g2m.genes, row.names(UMI@assays$RNA@counts))
      UMI <- CellCycleScoring(UMI, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=T)
    } else if(organism == 'mouse'){
      m.s.genes <- c("Mcm4", "Exo1", "Slbp", "Gmnn", "Cdc45", "Msh2", "Mcm6", "Rrm2", "Pold3", "Blm", "Ubr7", "Mcm5", "Clspn", "Hells", "Nasp", "Rpa2", "Rad51ap1", "Tyms", "Rrm1", "Rfc2", "Prim1", "Brip1", "Usp1", "Ung", "Pola1", "Mcm2", "Fen1", "Tipin", "Pcna", "Cdca7", "Uhrf1", "Casp8ap2", "Cdc6", "Dscc1", "Wdr76", "E2f8", "Dtl", "Ccne2", "Atad2", "Gins2", "Chaf1b", "Pcna-ps2")
      m.g2m.genes <- c("Nuf2", "Psrc1", "Ncapd2", "Ccnb2", "Smc4", "Lbr", "Tacc3", "Cenpa", "Kif23", "Cdca2", "Anp32e", "G2e3", "Cdca3", "Anln", "Cenpe", "Gas2l3", "Tubb4b", "Cenpf", "Dlgap5", "Hjurp", "Cks1brt", "Gtse1", "Bub1", "Birc5", "Ube2c", "Rangap1", "Hmmr", "Ect2", "Tpx2", "Ckap5", "Cbx5", "Nek2", "Ttk", "Cdca8", "Nusap1", "Ctcf", "Cdc20", "Cks2", "Mki67", "Tmpo", "Ckap2l", "Aurkb", "Kif2c", "Cdk1", "Kif20b", "Top2a", "Aurka", "Ckap2", "Hmgb2", "Cdc25c", "Ndc80", "Kif11")
      m.s.genes <- intersect(m.s.genes, row.names(UMI@assays$RNA@counts))
      m.g2m.genes <- intersect(m.g2m.genes, row.names(UMI@assays$RNA@counts))
      UMI <- CellCycleScoring(UMI, s.features=m.s.genes, g2m.features=m.g2m.genes, set.ident=T)
    } else {stop("organism must be human or mouse")}
  }#end if

  ### scaling and regression
  all.genes <- rownames(UMI)
  if (all(variables == 'none')){
    cycle.dir <- file.path(PP_dir,paste0("No_Regression"))
    suppressWarnings(if (!file.exists(cycle.dir)){dir.create(cycle.dir)})
    UMI <- ScaleData(object=UMI, features=all.genes)
    ### perform PCA on the scaled data.
    UMI <- suppressMessages(RunPCA(UMI, features=UMI[["RNA"]]@var.features, npcs=30, do.print=F, verbose=FALSE))
    UMI2 <- suppressMessages(RunTSNE(UMI, dims = 1:20))
    UMI2 <- suppressWarnings(RunUMAP(UMI2, dims = 1:20, verbose=FALSE))
    ### Single pdf with 3 pages, 4 plots per page
    cat(bold(green("Plotting dimensional reduction graphs with no regression \n")))
    pdf(file.path(cycle.dir, "/02c_DimReduction_NoRegression.pdf"), width=14, height=12, useDingbats=FALSE)
    popsicleR:::four_plots(UMI2, "pca")
    popsicleR:::four_plots(UMI2, "tsne")
    popsicleR:::four_plots(UMI2, "umap")
    invisible(dev.off())
    cat(paste0(silver("Plots saved in: ")),bold(silver("02.PreProcessing dedicated subfolder \n")))
  } else {
    cycle.dir <- file.path(PP_dir,paste0("Regression_on_",paste(unlist(variables), collapse='_')))
    suppressWarnings(if (!file.exists(cycle.dir)){dir.create(cycle.dir)})
    UMI <- ScaleData(object=UMI, vars.to.regress=variables, features=all.genes)
    ### perform PCA on the scaled data.
    UMI <- RunPCA(UMI, features=UMI[["RNA"]]@var.features, npcs=30, do.print=F, verbose=FALSE)
    UMI2 <- suppressMessages(RunTSNE(UMI, dims = 1:20))
    UMI2 <- suppressWarnings(RunUMAP(UMI2, dims = 1:20, verbose=FALSE))
    ### PCA plots
    cat(bold(green("Plotting dimensional reduction graphs after regression \n")))
    pdf(file.path(cycle.dir, "/02c_DimReduction_PostRegression.pdf"), width=14, height=12, useDingbats=FALSE)
    popsicleR:::four_plots(UMI2, "pca")
    popsicleR:::four_plots(UMI2, "tsne")
    popsicleR:::four_plots(UMI2, "umap")
    invisible(dev.off())
    cat(paste0(silver("Plots saved in: ")),bold(silver("02.PreProcessing dedicated subfolder \n")))
  }

  ### PC exploration
  if (explore_PC == TRUE){
    cat(bold(green("Plotting graphs to explore PCs \n")))
    pdf(paste0(cycle.dir, "/02d_VizPCA_HVG.pdf"), useDingbats=FALSE)
    print(VizDimLoadings(UMI, dims = 1:2, reduction = "pca"))
    invisible(dev.off())
    ### Heatmap
    pdf(paste0(cycle.dir, "/02e_PCHeatmap.pdf"), width=9, height=12, useDingbats=FALSE)
    DimHeatmap(UMI, dims=1:9, cells = 500, balanced = TRUE)
    DimHeatmap(UMI, dims=10:18, cells = 500, balanced = TRUE)
    DimHeatmap(UMI, dims=19:27, cells = 500, balanced = TRUE)
    invisible(dev.off())
    ### Jackstraw
    UMI2 <- suppressWarnings(JackStraw(object=UMI, dims=30, num.replicate=100))
    UMI2 <- suppressWarnings(ScoreJackStraw(UMI2, dims = 1:30))
    pdf(paste0(cycle.dir, "/02f_JackStrawPlot.pdf"), width=15, height=15, useDingbats=FALSE)
    print(JackStrawPlot(object=UMI2, dims=1:30))
    invisible(dev.off())
    ### Elbow Plot
    pdf(paste0(cycle.dir, "/02g_PCElbowPlot.pdf"), useDingbats=FALSE)
    print(ElbowPlot(object=UMI, ndims=30))
    invisible(dev.off())
    cat(paste0(silver("Plots saved in: ")),bold(silver("02.PreProcessing dedicated subfolder \n")))
    cat(paste0(cyan("\nOnce identified the PCs number to use in your analysis, perform clustering running "),bold(cyan("CalculateCluster \n"))))
  }
  ### write log
  write.table(t(c(as.character(as.POSIXct(Sys.time())),"ApplyRegression:","variables",variables,"explore_PC",explore_PC)),file=file.path(out_folder,"popsicleR.log"), sep="\t", row.names=F, col.names=F, quote=F, append=T)

  return(UMI)
}

#' CalculateCluster
#'
#' @description
#' Performs a Louvain clustering exploiting Seurat functions and generate the reduced dimensional embeddings (tSNE and UMAP). Clusters and embeddings are calculated according to the user-provided number of principal components.
#' first time to explore variables distribution and subsequently to regress on appropriate variables.
#'
#' @usage
#'
#' CalculateCluster(UMI, dim_pca, organism=c("human","mouse"), marker.list='none', PCA=TRUE, cluster_res=0.8, out_folder=getwd())
#'
#' @param UMI Input UMI object generated via ApplyRegression function
#' @param dim_pca Number of principal component to use for clustering and embeddings calculation.
#' @param organism Input organism to define the marker list to use. Can be human or mouse
#' @param marker.list List of markers. Default is 'none' and a proprietary marker list is applied for the selected organism. User can fed a custom genelist specifying the file_name of a tab delimited file stored in the working directory; first row must be a header.
#' @param PCA If TRUE plots also the PCA embeddings, else only tSNE and UMAP are shown. Default is TRUE
#' @param cluster_res Cluster resolution. Default is 0.8. User can specify, through a numeric vector, multiple clustering resolution to explore. When more than one resolution are specified, time consuming steps (e.g. Seurat::FindAllMarkers) are not performed.
#' @param out_folder Output folder. Default is the working directory.
#'
#' @details
#'
#' CalculateCluster returns several graphs in the "03.Clustering" dedicated subfolder:
#
#' graphs 03a, show for each user-provided clustering resolution projections on reduced spaces (PCA, t-SNE, and UMAP) colored according to the clustering.
#'
#' graph 03b, computed only when more clustering resolution are provided, displays a clustering tree that shows how cells are assigned to the clusters at the various clustering resolutions
#'
#' graph 03c, reports a phylogenetic tree showing correlation between the identified clusters for each clustering resolution.
#'
#' graphs 03d, explore data distribution based on QC features, MALAT1 and GAPDH expression and doublets estimation through projections on reduced spaces (PCA, t-SNE, and UMAP) final embeddings.
#'
#' graph 03e, returns violin plots reporting the expression of the 2 top marker genes for each cluster.
#'
#' graphs 03f, report the expression of the 2 top marker genes for each cluster on reduced spaces projections (t-SNE and UMAP).
#'
#` graph 03g, plots a heatmap that resume for each cluster the top10-markers expression values.
#'
#' graphs 03h, show marker genes expression on reduced spaces projections (t-SNE and UMAP).
#'
#' graph 03i, returns violin plots of marker genes expression for each clustering resolution provided.
#'
#' graph 03j, displays dot plots reporting marker genes expression in the identified clusters for each user-provided clustering resolution.
#'
#' Graphs from 03e to 03g are computed only when a single clustering resolution is provided by the user. popsicleR also provides a .txt file containing a marker list for each cluster. This file is computed only when a single clustering resolution is specified.
#'
#' @return Returns a Seurat Object after clustering.
#'
#' @examples
#'
#' CalculateCluster(umi_object, dim_pca= 12, organism= "human", cluster_res= c(0.4, 0.6, 0.8))
#'
#' CalculateCluster(umi_object, dim_pca= 12, organism= "human", cluster_res= 0.8)
#'
#' @author Jimmy Caroli, Francesco Grandi
#'
#' @export
#'

CalculateCluster <- function(UMI, dim_pca, organism=c("human","mouse"), marker.list='none', PCA=TRUE, cluster_res=0.8, out_folder=getwd()){
  #dimpca --> number of principal components
  #cluster_res --> cluster resolution (default 0.8)


  Cluster_dir <- paste0(out_folder,"/03.Clustering/")
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
  }else {stop("organism must be human or mouse")}
  ### Performing scaling and clustering based on provided inputs
  UMI <- RunUMAP(UMI, dims = 1:dim_pca, verbose=FALSE)
  UMI <- RunTSNE(UMI, dims = 1:dim_pca, verbose=FALSE)
  UMI <- FindNeighbors(UMI, dims = 1:dim_pca)
  UMI <- FindClusters(UMI, resolution = cluster_res)

  ### plot clusters
  if(PCA == TRUE){
    algos <- c("umap","tsne","pca")
  }else{
    algos <- c("umap","tsne")
  }

  for(algo in algos){
    algo_plot_clusters(Cluster_dir, UMI, algo, "03a", cluster_res)
  }

  ### clustree analysis
  if (length(cluster_res)>1) {
    metadata<-UMI@meta.data
    for (res in cluster_res) {
      res_col<-paste0("RNA_snn_res.",res)
      metadata[,res_col]<-as.numeric(as.character(metadata[,res_col]))
    }
    pdf(file.path(Cluster_dir,paste0("03b_clustree_resolution_analysis.pdf")), width=14, height=10, useDingbats=FALSE)
    print(clustree(metadata, prefix = "RNA_snn_res."))
    dev.off()
  }

  ### phylo tree
  pdf(paste0(Cluster_dir, "/03c_clusters_phylogenetic_tree.pdf"), width=10, height=10, useDingbats=FALSE)
  for(res.i in cluster_res) {
    Idents(UMI)<-paste0("RNA_snn_res.",res.i)
    clus_time_tree <- BuildClusterTree(UMI)
    data.tree <- Tool(object = clus_time_tree, slot = "BuildClusterTree")
    plot.phylo(x = data.tree, direction = "rightwards", main=paste0("phylogenetic tree of clusters @ res", res.i))
  }
  invisible(dev.off())

  ### plot cell features
  for(algo in algos){
    algo_plot_features(Cluster_dir, UMI, algo, "03d", organism)
  }

  ### Finding differentially expressed features (cluster biomarkers)
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  # Finding cluster markers

  ### if only one resolution has been specified
  if (length(cluster_res)==1) {
    umi.markers <- FindAllMarkers(UMI, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
    FC_col <- colnames(umi.markers)[grep("avg_log", colnames(umi.markers))]
    umi.markers <- umi.markers[,c("cluster","gene","pct.1","pct.2",FC_col,"p_val","p_val_adj")]
    umi.markers <- umi.markers[umi.markers$p_val_adj<=0.05,]

    write.table(umi.markers, file=paste0(Cluster_dir,"03_markers@res",cluster_res,".txt"), sep="\t", col.names=T, row.names=F, quote=F)

    top.markers.2 <- umi.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n=2, wt=get(FC_col))
    ftp_h = 5+(0.5*length(levels((umi.markers$cluster))))

    ### visualize markers by violin plots
    cat(bold(green("Plotting Top Markers graphs for each cluster \n")))

    pdf(paste0(Cluster_dir, "/03e_Violin_Top_marker_genes.pdf"), width=12, height=9, useDingbats=FALSE)
    for(cluster in unique(top.markers.2$cluster)) {
      genes <- top.markers.2$gene[top.markers.2$cluster == cluster]
      vg1 <- popsicleR:::VLN(UMI, feats=genes[1], point =FALSE) + ggplot2::geom_boxplot(width=0.1, outlier.shape=NA) + NoLegend() + ggplot2::xlab("")
      if (length(genes)==2) {
        vg2 <- popsicleR:::VLN(UMI, feats=genes[2], point =FALSE) +
          ggplot2::geom_boxplot(width=0.1, outlier.shape=NA) +
          NoLegend() +
          ggplot2::xlab(paste("Top markers for cluster",cluster, "@ res", cluster_res)) } else {
            vg1<-vg1+xlab(paste("Top markers for cluster",cluster, "@ res", cluster_res))
            vg2<-ggplot()
          }
      print(patchwork::wrap_plots(vg1, vg2, ncol=1))
    }
    invisible(dev.off())

    ### Visualize TOP markers on a dimensional reduction plot
    pdf(paste0(Cluster_dir, "/03f_umap_Top_marker_genes.pdf"), width=14, height=7, useDingbats=FALSE)
    for(cluster in unique(top.markers.2$cluster)) {
      genes <- top.markers.2$gene[top.markers.2$cluster == cluster]
      vg1 <- FeaturePlot(UMI, features=genes[1], reduction="umap", ncol=1) +
        NoLegend() +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab(paste0("Top markers for cluster ",cluster," @ res ", cluster_res,"\nUMAP 2"))
      if (length(genes)==2) {
        vg2 <- FeaturePlot(UMI, features=genes[2], reduction="umap", ncol=1) +
          NoLegend() +
          ggplot2::xlab("UMAP 1") +
          ggplot2::ylab("UMAP 2") } else {vg2<-ggplot()}
      print(patchwork::wrap_plots(vg1, vg2, ncol=2))
    }
    invisible(dev.off())

    ### TSNE Top Markers
    pdf(paste0(Cluster_dir, "/03f_tsne_Top_marker_genes.pdf"), width=14, height=7, useDingbats=FALSE)
    for(cluster in unique(top.markers.2$cluster)) {
      genes <- top.markers.2$gene[top.markers.2$cluster == cluster]
      vg1 <- FeaturePlot(UMI, features=genes[1], reduction="tsne", ncol=1) +
        NoLegend() +
        ggplot2::xlab("tSNE 1") +
        ggplot2::ylab(paste0("Top markers for cluster ",cluster," @ res ", cluster_res,"\ntSNE 2"))
      if (length(genes)==2) {
        vg2 <- FeaturePlot(UMI, features=genes[2], reduction="tsne", ncol=1) +
          NoLegend() +
          ggplot2::xlab("TSNE 1") +
          ggplot2::ylab("TSNE 2")} else {vg2<-ggplot()}
      print(patchwork::wrap_plots(vg1, vg2, ncol=2))
    }
    invisible(dev.off())

    ### Plotting cluster biomarkers
    ### expression heatmap
    top.markers.10 <- umi.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = get(FC_col))
    pdf(paste0(Cluster_dir, "/03g_heatmap_top.markers.pdf"), width=18, height=5+(0.5*length(levels((umi.markers$cluster)))), useDingbats=FALSE)
    print(DoHeatmap(UMI, features=top.markers.10$gene, slot="scale.data") + NoLegend())
    invisible(dev.off())
  }

  ### Visualize given marker list on a dimensional reduction plot

  marker.list <- lapply(marker.list, function(x) {x[x %in% row.names(UMI)]})
  ### remove empty lists
  marker.list <- marker.list[lapply(marker.list,length)>0]

  #Plotting UMAP and TSNE for immune markers
  ftp_h = 4*(ceiling(length(unlist(marker.list))/4))
  cat(bold(green("Plotting Markers graphs for each cluster \n")))
  popsicleR:::FTP(UMI, Cluster_dir, "03h_", "umap", ftp_h, unlist(marker.list), "_marker_list")
  popsicleR:::FTP(UMI, Cluster_dir, "03h_", "tsne", ftp_h, unlist(marker.list), "_marker_list")

  # width proportional to the number of clusters, height proportional to the markers number
  vln_w<-length(levels(UMI$seurat_clusters))*2
  vln_h<-round((length(unlist(marker.list))/4)*3)
  ### visualize common markers by violin plots
  pdf(paste0(Cluster_dir, "/03i_Violin_marker_list.pdf"), height = vln_h + 1, width= vln_w, useDingbats=FALSE)
  for (res.i in cluster_res){
    vln <- VlnPlot(UMI, features=unlist(marker.list), ncol=4,pt.size=F, group.by=paste0("RNA_snn_res.",res.i))
    print(vln + plot_annotation(title = paste0("Markers expression for each cluster at res. ", res.i),
                      theme = theme(plot.title = element_text(hjust = 0.5, face="bold", size =24))))
  }
  dev.off()
  # manca  valore risoluzione nel plot


  marker_list <- list()
    for (m in marker.list) {marker_list <- c(marker_list,list(m[m%in% rownames(UMI)]))}
  names(marker_list) <- names(marker.list)
  markers_sub <- marker_list[lapply(marker_list,length)>0]
  markers.to.plot <- as.character(unlist(markers_sub))
  names(markers.to.plot) <- NULL
  ann_unique <- unique(as.character(UMI@meta.data[,paste0("RNA_snn_res.",max(cluster_res))]))
  len_lab <- length(ann_unique)
  max_markers <- max(nchar(markers.to.plot))
  len_markers <- length(markers.to.plot)
  k_height <- ifelse(len_markers>10, 2, 1)
  width_sig <- len_markers*0.35 + 3
  height_sig <- len_lab*0.35 + max_markers*0.3 + k_height
  ### dotplot of immune markers for each cluster
  pdf(paste0(Cluster_dir, "03j_Dotplot_marker_list.pdf"), width=width_sig, height=height_sig, useDingbats=FALSE)
  for(res.i in cluster_res) {
    print(DTP(UMI, markers_sub, paste0("RNA_snn_res.",res.i)))
  }
  invisible(dev.off())

  cat(paste0(silver("Plots saved in: ")),bold(silver("03.Clustering")), silver("folder \n"))
  cat(paste0(cyan("Next suggested step is annotation, run "),bold(cyan("MakeAnnotation \n"))))
  ### write log
  write.table(t(c(as.character(as.POSIXct(Sys.time())),"CalculateCluster:","dim_pca",dim_pca,"cluster_res",cluster_res)),file=file.path(out_folder,"popsicleR.log"), sep="\t", row.names=F, col.names=F, quote=F, append=T)

  return(UMI)
}

#' MakeAnnotation
#'
#' @description
#'
#' MakeAnnotation annotates cells using SingleR and scMCA built-in references.
#'
#' For human samples, annotation is performed basing on the Human Primary Cell Atlas (HPCA) and the Blueprint Encode (BpEn) references from celldex and SingleR packages.
#'
#' For mouse samples, annotation is performed according to mouse bulk RNA-seq datasets (Mouse RNA-seq) and Immunological Genome Project (ImmGen) references from SingleR, and the single-cell Mouse Cell Atlas (scMCA) reference.
#'
#' Annotations are displayed on PCA, UMAP and TSNE embeddings, and evaluated through further graphs.
#'
#' @usage
#'
#' MakeAnnotation(UMI, organism=c("human","mouse"), marker.list='none', thresh=20, cluster_res=NULL, out_folder=getwd())
#'
#' @param UMI The input umi object generated via CalculateCluster function
#' @param organism The organism on which perfom the analysis. Can be human or mouse.
#' @param marker.list Markers list used for the plots. Must be the same used in CalculateCluster. Default uses proprietary lists.
#' @param thresh Rare population number threshold. Default is 20
#' @param cluster_res Defines on which cluster resolution perform the annotation. Must be a clustering resolution already calculated by CalculateCluster in the previous step. Default is NULL.
#' @param out_folder Output folder. Default is the working directory
#'
#' @details
#' MakeAnnotation returns several graphs in the "04.Annotation" dedicated subfolder:
#'
#' graphs 04a, show single-cell annotations for each reference on reduced spaces projections (t-SNE and UMAP). Cells are colored according to the assigned cell type.
#'
#' graphs 04b, show cluster annotations for SingleR references on reduced spaces projections (t-SNE and UMAP). Cells are colored according to the clustering and labeled basing on the assigned cell type.
#'
#' graphs 04c, report single-cell annotations for each reference on reduced spaces projections (t-SNE and UMAP). Cells are colored in red according to the assigned cell type. Each cell type is displayed individually.
#'
#' graphs 04d, return dot plots reporting, in the cell types assigned using the different references, the markers genelist expression.
#'
#' graphs 04e, are pie-chart plots displaying, for each cluster, the proportion of cells annotated with a specific label based on the different references.
#'
#' @return Returns a Seurat Object after annotation.
#'
#' @examples
#'
#' MakeAnnotation(umi_object, organism="human", marker.list= "markers.txt", cluster_res=0.8)
#'
#' @author Jimmy Caroli, Francesco Grandi
#'
#' @export
#'

MakeAnnotation <- function(UMI, organism=c("human","mouse"), marker.list='none', thresh=20, cluster_res=NULL, out_folder=getwd()){
  Annot_dir <- paste0(out_folder,"/04.Annotation/")
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
    UMI <- SR_plots("hpca", hpca.se, UMI, Annot_dir, cluster_res)
    UMI <- SR_plots("BpEn", BpEn.se, UMI, Annot_dir, cluster_res)

    annotations <- c("hpca.sc.main.labels","BpEn.sc.main.labels")
    cat(bold(green("Plotting dimensional reduction graphs for each population found in the sample \n")))
    for(single_annot in annotations){
      annotation_plot(Annot_dir, "04c_", UMI, "UMAP", single_annot, "CellPopulations")
      annotation_plot(Annot_dir, "04c_",  UMI, "TSNE", single_annot, "CellPopulations")
      marker_list <- list()
        for (m in marker.list) {marker_list <- c(marker_list,list(m[m%in% rownames(UMI)]))}
      names(marker_list) <- names(marker.list)
      markers_sub <- marker_list[lapply(marker_list,length)>0]
      markers.to.plot <- as.character(unlist(markers_sub))
      names(markers.to.plot) <- NULL
      ann_unique <- unique(as.character(UMI@meta.data[,single_annot]))
      len_lab <- length(ann_unique)
      max_lab <- max(nchar(ann_unique))
      max_markers <- max(nchar(markers.to.plot))
      len_markers <- length(markers.to.plot)
      k_height <- ifelse(len_markers>10, 2, 1)
      width_sig <- len_markers*0.35 + max_lab*0.3 + 3
      height_sig <- len_lab*0.35 + max_markers*0.3 + k_height
      pdf(paste0(Annot_dir, paste0("04d_DotPlot_Markers_", single_annot, ".pdf")), width=width_sig, height=height_sig, useDingbats=FALSE)
      print(popsicleR:::DTP(UMI, markers_sub, single_annot))
      invisible(dev.off())
    }
    ### corrplot
    for(single_annot in annotations){

      pdf(paste0(Annot_dir, paste0("04e_Corrplot_clusters_vs_", single_annot, ".pdf")), width=8, height=7, useDingbats=FALSE)
      for (res in cluster_res ) {
        res.col<-paste0("RNA_snn_res.",res)
        clusters<-UMI@meta.data[,res.col]
        labels<-UMI@meta.data[,single_annot]
        a <- table(clusters,labels)
        a_freq<-a/rowSums(a)
        corrplot(t(a_freq),is.corr=F,tl.col='black',title = paste('clusters @',res, "vs", single_annot),  mar=c(0,0,2,0), method="pie", tl.srt=45)
      }
      invisible(dev.off())
    }


    cat(paste0(silver("Plots saved in: ")),bold(silver("\\04.Annotation\\")), silver("folder \n"))
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

    ### singleR with ImmGen and MouseRNAseq
    Ig.se <- suppressMessages(celldex::ImmGenData())
    mouseRNA.se <- suppressMessages(celldex::MouseRNAseqData())
    cat(bold(green("Plotting single cell and cluster annotations \n")))
    UMI <- popsicleR:::SR_plots("ImmGen", Ig.se, UMI, Annot_dir, cluster_res)
    UMI <- popsicleR:::SR_plots("MouseRNAseq", mouseRNA.se, UMI, Annot_dir, cluster_res)

    ### run scMCA
    matrice_norm <- as.matrix(GetAssayData(UMI))
    mca_result <- scMCA(scdata = matrice_norm, numbers_plot = 3)
    scMCA_assignment <- mca_result$scMCA

    ### Performing data simplification
    ### simplification of labels
    simple_assignment <- do.call(rbind,strsplit(scMCA_assignment,"(", fixed=T))
    simple_assignment <- unlist(simple_assignment[,1])
    #simple_assignment <- toupper(simple_assignment)
    simple_assignment <- gsub("T-","T",simple_assignment)
    simple_assignment <- gsub(" *(C|c)ells*","", simple_assignment)
    test <- strsplit(simple_assignment, "_")
    test <- sapply(test, function(x,m) c(x, rep(NA, m-length(x))), max(rapply(test,length)))
    simple_assignment <- unlist(test[1,])
    simple_assignment <- gsub("arcrophage","acrophage",simple_assignment)
    simple_assignment <- gsub("^T$","T cells", simple_assignment)
    simple_assignment <- gsub("^NK$","NK cells", simple_assignment)
    simple_assignment <- gsub("^B$","B cells", simple_assignment)
    simple_assignment <- gsub("Neutrophil granulocyte","Neutrophil", simple_assignment)
    simple_assignment <- gsub("Basophil granulocyte","Basophil", simple_assignment)
    simple_assignment <- gsub("Eosinophil granulocyte","Eosinophil", simple_assignment)
    UMI$scMCA <- scMCA_assignment
    UMI$scMCA_simple <- simple_assignment
    ### Plot TSNE and UMAP with scMCA annotation
    pdf(paste0(Annot_dir, "/04a_TSNE_scMCA_simple.pdf"), width=20, height=10, useDingbats=FALSE)
    dimplot <- DimPlot(UMI, reduction="tsne", group.by="scMCA_simple", label=T, pt.size=1, repel=T) +
            ggplot2::ggtitle(paste0("scMCA simple - "," single cell annotation")) +
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
            ggplot2::guides(col=ggplot2::guide_legend(ncol=1)) + xlab("TSNE 1") + ylab("TSNE 2")
    print(patchwork::wrap_plots(dimplot + guide_area() + plot_layout(guides = 'collect')))
    invisible(dev.off())
    pdf(paste0(Annot_dir, "/04a_UMAP_scMCA_simple.pdf"), width=20, height=10, useDingbats=FALSE)
    dimplot <- DimPlot(UMI, reduction="umap", group.by="scMCA_simple", label=T, pt.size=1, repel=T) +
            ggplot2::ggtitle(paste0("scMCA simple - "," single cell annotation")) +
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5)) +
            ggplot2::guides(col=ggplot2::guide_legend(ncol=1)) + xlab("UMAP 1") + ylab("UMAP 2")
    print(patchwork::wrap_plots(dimplot + guide_area() + plot_layout(guides = 'collect')))
    invisible(dev.off())

    ### plot each cell type
    clean_labels <- UMI$scMCA_simple
    ultra_rare <- names(table(clean_labels)[table(clean_labels) < thresh])
    clean_labels[clean_labels%in%ultra_rare] <- "Other"
    UMI$clean_labels <- clean_labels

    ### Plotting annotated populations localization in TSNE, UMAP and PCA. [scMCA]
    cat(bold(green("Plotting dimensional reduction graphs for each population \n")))
    popsicleR:::annotation_plot(Annot_dir, "04c_", UMI, "UMAP", "clean_labels", "CellPopulations_scMCA")
    popsicleR:::annotation_plot(Annot_dir, "04c_", UMI, "TSNE", "clean_labels", "CellPopulations_scMCA")

    # Plotting annotated populations localization in TSNE, UMAP and PCA. [ImmGen]
    annotations <- c("ImmGen.sc.main.labels", "MouseRNAseq.sc.main.labels")
    for(single_annot in annotations){
      popsicleR:::annotation_plot(Annot_dir, "04c_", UMI, "UMAP", single_annot, "CellPopulations")
      popsicleR:::annotation_plot(Annot_dir, "04c_", UMI, "TSNE", single_annot, "CellPopulations")
      marker_list <- list()
        for (m in marker.list) {marker_list <- c(marker_list,list(m[m%in% rownames(UMI)]))}
      names(marker_list) <- names(marker.list)
      markers_sub <- marker_list[lapply(marker_list,length)>0]
      markers.to.plot <- as.character(unlist(markers_sub))
      names(markers.to.plot) <- NULL
      ann_unique <- unique(as.character(UMI@meta.data[,single_annot]))
      len_lab <- length(ann_unique)
      max_lab <- max(nchar(ann_unique))
      max_markers <- max(nchar(markers.to.plot))
      len_markers <- length(markers.to.plot)
      k_height <- ifelse(len_markers>10, 2, 1)
      width_sig <- len_markers*0.35 + max_lab*0.3 + 3
      height_sig <- len_lab*0.35 + max_markers*0.3 + k_height
      pdf(paste0(Annot_dir, paste0("04d_DotPlot_Markers_", single_annot, ".pdf")), width=width_sig, height=height_sig, useDingbats=FALSE)
      print(popsicleR:::DTP(UMI, markers_sub, single_annot))
      invisible(dev.off())
    }

    ## Custom DotPlot for markers
    single_annot <- "scMCA_simple"
    marker_list <- list()
      for (m in marker.list) {marker_list <- c(marker_list,list(m[m%in% rownames(UMI)]))}
    names(marker_list) <- names(marker.list)
    markers_sub <- marker_list[lapply(marker_list,length)>0]
    markers.to.plot <- as.character(unlist(markers_sub))
    names(markers.to.plot) <- NULL
    ann_unique <- unique(as.character(UMI@meta.data[,single_annot]))
    len_lab <- length(ann_unique)
    max_lab <- max(nchar(ann_unique))
    max_markers <- max(nchar(markers.to.plot))
    len_markers <- length(markers.to.plot)
    k_height <- ifelse(len_markers>10, 2, 1)
    width_sig <- len_markers*0.35 + max_lab*0.3 + 3
    height_sig <- len_lab*0.35 + max_markers*0.3 + k_height
    pdf(paste0(Annot_dir, "/04d_DotPlot_Markers_scMCA.pdf"), width=width_sig, height=height_sig, useDingbats=FALSE)
    print(popsicleR:::DTP(UMI, marker.list, "scMCA_simple"))
    invisible(dev.off())

    cat(paste0(silver("Plots saved in: "),bold(silver("\\04.Annotation\\")), silver(" folder \n")))


    ### corrplot
    for(single_annot in c(annotations,"scMCA_simple")){

      pdf(paste0(Annot_dir, paste0("04e_Corrplot_clusters_vs_", single_annot, ".pdf")), width=8, height=7, useDingbats=FALSE)
      for (res in cluster_res ) {
        res.col<-paste0("RNA_snn_res.",res)
        clusters<-UMI@meta.data[,res.col]
        labels<-UMI@meta.data[,single_annot]
        a <- table(clusters,labels)
        a_freq<-a/rowSums(a)
        corrplot(t(a_freq),is.corr=F,tl.col='black',title = paste('clusters @',res, "vs", single_annot),  mar=c(0,0,2,0), method="pie", tl.srt=45)
      }
      invisible(dev.off())
    }

  }else {stop("organism must be human or mouse")}
  ### write log
  write.table(t(c(as.character(as.POSIXct(Sys.time())),"MakeAnnotation:","thresh",thresh, "cluster_res",cluster_res)),file=file.path(out_folder,"popsicleR.log"), sep="\t", row.names=F, col.names=F, quote=F, append=T)

  return(UMI)
}


###################################################
### The main function used to identify doublets ###
### modified from the R code provided at https://rdrr.io/github/ChengxiangQiu/rscrublet/
### by Chengxiang Qiu <cxqiu@uw.edu> based on the method
### implemented in the python module Scrublet
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
  p_obs <- p_obs + labs(x="Doublet scores", y="Counts", title="Observed Cells") +
           theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
  p_sim <- ggplot2::ggplot(dat_sim, aes(x = doublet_scores))
  p_sim <- p_sim + geom_histogram(aes(fill = clust), binwidth = 0.02, color = "grey50")
  p_sim <- p_sim + geom_vline(xintercept = predicted_threshold, color = "blue")
  p_sim <- p_sim + labs(x="Doublet scores", y="Counts", title="Simulated Doublets") + theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))

  p_obs2 <- ggplot2::ggplot(dat_obs, aes(x = doublet_scores)) + stat_density(geom="line", color="red") + geom_vline(xintercept = predicted_threshold, color = "blue")
  p_obs2 <- p_obs2 + labs(x="Doublet scores", y="Density", title="") + theme_classic(base_size = 10) +
            theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
  p_sim2 <- ggplot2::ggplot(dat_sim, aes(x = doublet_scores)) + stat_density(geom="line", color="red") + geom_vline(xintercept = predicted_threshold, color = "blue")
  p_sim2 <- p_sim2 + labs(x="Doublet scores", y="Density", title="") + theme_classic(base_size = 10) + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
  #
  pdf(paste0(directory,"/01i_histogram of doublet scores.pdf"),8,8, useDingbats=FALSE)
  gridExtra::grid.arrange(p_obs, p_sim, p_obs2, p_sim2, nrow = 2, ncol = 2)
  dev.off()
  cat(paste0(silver("Plots saved in: ")),bold(silver("01.QC_Plots\\01i_histogram of doublet scores.pdf \n")))

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
