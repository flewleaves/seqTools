
#' @keywords internal
#' @noRd
sc_filter <- function(scRNA, s = 0.02) {
  
  umi_filter <- function(x, s){
    xl = log10(x[x > 0])
    q = quantile(xl, seq(.01, .99, s), na.rm = T)
    d = as.numeric(na.omit(stats::filter(diff(q), rep(1/3, 3))))
    d2 = diff(d)
    m = length(d2) %/% 2
    b = median(abs(d2)) + 1e-8
    
    cl = max(0, d2[1:m]) / b
    ch = max(0, -d2[(m+1):length(d2)]) / b
    
    md = median(xl)
    w = mad(xl, constant = 1)
    
    return(c(lower = 10^(md - max(2.5, 5 - max(0, cl - 2)) * w),
      upper = 10^(md + max(2.5, 4 - max(0, ch - 2)) * w)))
  }

  mt_filter <- function(x, ks = c(10, 15, 20, 25)) {
    frac <- max(x, na.rm = TRUE) <= 1
    if (frac) x <- x * 100
    
    kp <- sapply(ks, function(k) mean(x <= k))
    dkp <- diff(kp)
    idx <- which(dkp < 0.05)[1]
    k <- if (is.na(idx)) ks[length(ks)] else ks[idx]
    
    # 返回原始尺度的阈值
    if (frac) k / 100 else k
  }

  umi <- scRNA$nCount_RNA
  nf <- scRNA$nFeature_RNA
  mt <- scRNA$percent.mt
  filter <- c(umi_filter(umi, s), umi_filter(nf, s), mt_filter(mt))
  print(paste("Using filter parameters as below: nCount_RNA min", filter[1], "max", filter[2], "nFeature_RNA min", filter[3], "max", filter[4], "mt.percent less than", filter[5], sep = " "))
  scRNA <- subset(scRNA, subset = nCount_RNA > filter[1] & nCount_RNA < filter[2] & nFeature_RNA > filter[3] & nFeature_RNA < filter[4] & percent.mt < filter[5])
  return(scRNA)
}

#'
#' @param scRNA
#' @keywords internal
#' @noRd
find_doublet <- function(scRNA, samples){

  scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
  scRNA<- FindVariableFeatures(object = scRNA, mean.function = ExpMean, dispersion.function = LogVMR)
  scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))

  temp <- as.SingleCellExperiment(scRNA)
  samples <- scRNA[[samples]][,1]
  temp <- scDblFinder::scDblFinder(temp, BPPARAM=SnowParam(3),samples = samples)
  scRNA <- CreateSeuratObject(
    counts = assay(temp, "counts"),  # 提取基因表达矩阵
    assay = "RNA",                  # 命名数据层为RNA
    meta.data = as.data.frame(colData(temp))  # 携带质检结果
  )

  scRNA <- subset(scRNA, scDblFinder.class == "singlet")

  return(scRNA)

}


#' Preprocessing of single-cell data
#'
#' @param dataList Must be a list containing Seurat objects.
#' @param species "Hs" for human and "Mm" for mouse, used for estimating mitochondria DNA percent.
#' @param DoubletFind Whether to remove doublets.
#' @param filter A numeric vector containing filtering criteria for min features, max features, min UMI count, max UMI count, and mitochondria DNA percent.
#' @param progress_saving Whether to save the object after running the function.
#' @param samples A string of the name in meta.data for grouping cells, only works when DoubletFind = T.
#'
#' @returns A dataList with Seurat objects preprocessed.
#' @export
#'
#' @examples
#' dataList <- scRNA_preprocessing(dataList = dataList, species = "Hs")
scRNA_preprocessing <- function(dataList, species, DoubletFind = T, filter = NULL, progress_saving = F, samples = NULL){

    cal = is.null(filter) || length(filter) != 5
    for(i in 1:length(dataList)){
      dataList[[i]][["percent.mt"]] <- PercentageFeatureSet( dataList[[i]], pattern = ifelse(species == "Mm", "^mt-", "^MT-"))
      if(cal){
        dataList[[i]] <- sc_filter(dataList[[i]])
      }else{
        dataList[[i]] <- subset(dataList[[i]], subset = nCount_RNA > filter[1] & nCount_RNA < filter[2] & nFeature_RNA > filter[3] & nFeature_RNA < filter[4] & percent.mt < filter[5])
      }
      print(Seurat::VlnPlot(dataList[[i]], features = c("nFeature_RNA","nCount_RNA","percent.mt"),layer = "counts", pt.size = 0.01, ncol = 3))
    }


  if(DoubletFind){

    print("Removing Doublets...")
    dataList <- lapply(dataList, function(x) {
      find_doublet(scRNA = x, samples = samples)
    })
    cat("Done\n")

  }

  if(progress_saving) saveRDS(dataList, "preprocessed_list.rds")

  return(dataList)
}

#' Perform standard normalization and dimension reduction for Seurat object
#'
#' @param dataList The result returned by scRNA_preprocessing()
#' @param species "Hs" for human and "Mm" for mouse, used for Cell Cycle Scoring.
#' @param CellCycleScoring Whether to regress the factor of cell cycle gene expression.
#' @param vars.to.regress Factors expected to be regressed. If CellCycleScoring = T,"S.Score", "G2M.Score" will be automatically added.
#' @param progress_saving Whether to save the object after running the function.
#'
#' @returns A merged Seurat Object which has performed normalization, finding variable features, scale and PCA.
#' @export
#'
#' @examples
#' scRNA <- scRNA_Normalization_Reduction(dataList = dataList, species = "Hs")
scRNA_Normalization_Reduction <- function(dataList, species, CellCycleScoring = F, vars.to.regress = c('percent.mt'), nfeatures = 2000 ,progress_saving = F){

  print("Merging objects...")
  if(length(dataList) > 1) {
    Merge.Seurat <- merge(x = dataList[[1]], y = dataList[-1])
  } else {
    Merge.Seurat <- dataList[[1]]
  }
  cat("Done\n")

  print("Normalizing and scaling...")
  Merge.Seurat <- NormalizeData(Merge.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  Merge.Seurat<- FindVariableFeatures(object = Merge.Seurat, nfeatures = nfeatures)
  #VariableFeaturePlot(Merge.Seurat)

  if(CellCycleScoring){

    print("Regressing cell cycle scoring")
    if(species == "Hs") {
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
    } else {
      s.genes <- stringr::str_to_title(cc.genes$s.genes)
      g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)
    }
    Merge.Seurat[['RNA']] <- JoinLayers(Merge.Seurat[['RNA']])
    gc()
    Merge.Seurat <- CellCycleScoring(Merge.Seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    # Visualize the distribution of cell cycle markers across
    #RidgePlot(Merge.Seurat, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
    Merge.Seurat <- ScaleData(Merge.Seurat, vars.to.regress = unique(c("S.Score", "G2M.Score", vars.to.regress)), features = VariableFeatures(Merge.Seurat))

  }else{

    Merge.Seurat <- ScaleData(Merge.Seurat, vars.to.regress = vars.to.regress, features = VariableFeatures(Merge.Seurat))

  }

  Merge.Seurat <- RunPCA(Merge.Seurat,features=VariableFeatures(object = Merge.Seurat))

  if(progress_saving){

    saveRDS(Merge.Seurat, "pca_scRNA.rds")
  }

  return(Merge.Seurat)
}


#' Perform SCTransform for Seurat object
#'
#' @param dataList The result returned by scRNA_preprocessing()
#' @param vars.to.regress Factors expected to be regressed.
#' @param progress_saving Whether to save the object after running the function.
#'
#' @returns A merged Seurat Object which has performed SCTransform. Not compatible with scRNA_Normalization_Reduction()
#' @export
#'
#' @examples
#' scRNA <- scRNA_SCTransform(dataList = dataList)
scRNA_SCTransform <- function(dataList, vars.to.regress = c('percent.mt'), nfeatures = 3000 ,progress_saving = F){

  print("Running SCTransform...")
  if(length(dataList) > 1) {
    dataList <- lapply(dataList, function(s){SCTransform(s, vars.to.regress = vars.to.regress, variable.features.n = nfeatures,verbose = TRUE)})
    Merge.Seurat <- merge(x = dataList[[1]], y = dataList[-1])
  } else {
    temp <- dataList[[1]]
    Merge.Seurat <- SCTransform(temp, vars.to.regress = vars.to.regress, verbose = TRUE)
  }
  cat("Done\n")
  DefaultAssay(Merge.Seurat) <- "SCT"
  Merge.Seurat <- RunPCA(Merge.Seurat)
  if(progress_saving){

    saveRDS(Merge.Seurat, "sct_scRNA.rds")
  }
  return(Merge.Seurat)
}


#' Integrate single-cell data
#'
#' @param Merge.Seurat The result returned by scRNA_Normalization_Reduction() or scRNA_SCTransform(). Must have multiple layers.
#' @param SCTransform Whether SCTransform is performed.
#' @param progress_saving Whether to save the object after running the function.
#'
#' @returns A Seurat Object which has cells mapped together for downstream visualization.
#' @export
#'
#' @examples
#' scRNA <- scRNA_Integration(scRNA)
scRNA_Integration <- function(Merge.Seurat, SCTransform = F, theta = NULL,progress_saving = F){

  scRNA <- IntegrateLayers(object = Merge.Seurat, method = HarmonyIntegration,
                         orig.reduction = "pca", new.reduction = 'harmony',
                         assay = ifelse(SCTransform, "SCT", "RNA"), theta = theta,verbose = TRUE)
  gc()
  scRNA[['RNA']] <- JoinLayers(scRNA[['RNA']])
  if(progress_saving){

    saveRDS(scRNA, "scRNA_harmony.rds")
  }
  return(scRNA)
}


#' Clustering cells
#'
#' @param scRNA Result of scRNA_Integration or single Seurat object after normalization.
#' @param PCs Dims used for clustering. If NULL, automatically choose the dim that can explain 90 percent variance.
#' @param resolution Resolution for clustering. Higher resolution means more communities. It can be a numeric vector, and if so, the function will ask the optimal choice.
#' @param reduction "harmony" if Integration is performed, otherwise "pca" can be used.
#' @param SCTransform Whether SCTransform is performed.
#' @param progress_saving Whether to save the object after running the function.
#'
#' @return A Seurat object with cells clustered.
#' @export
#'
#' @examples
#' scRNA <- scRNA_clustering(scRNA, PCs = 25, resolution = seq(0.1, 1.5, 0.1))
scRNA_clustering <- function(scRNA, PCs = NULL, resolution = seq(0.2,1,0.2), reduction = "harmony", SCTransform = F, k.param = 20,n.neighbors = 30L , min.dist = 0.3, method = "UMAP", progress_saving = F){

  if(!method %in% c("UMAP", "TSNE")){
    stop("Please check your dimension reduction method")
  }

  if(is.null(PCs)){
    stdev <- scRNA@reductions$pca@stdev
    var <- stdev^2
    EndVar = 0
    total <- sum(var)
    for(i in 1:length(var)){
      numerator <- sum(var[1:i])
      expvar <- numerator/total
      if(EndVar == 0){
        if(expvar > 0.9){
          EndVar <- EndVar + 1
          PCNum <- i
        }
      }
    }
    #Confirm #PC's determined explain > 90% of variance
    print(paste("Choose dims 1:",PCNum,"PC explains", sum(var[1:PCNum])/ sum(var), "variance"))

  }else{
    PCNum <- PCs[length(PCs)]
  }

  scRNA <- FindNeighbors(scRNA, dims = 1:PCNum, reduction = reduction, k.param = k.param)
  scRNA <- FindClusters(scRNA, verbose = TRUE, resolution = resolution)
  if(method == "UMAP"){
    scRNA <- RunUMAP(scRNA, dims = 1:PCNum, reduction = reduction, n.neighbors = n.neighbors, min.dist = min.dist)
  }else if(method == "TSNE"){
    scRNA <- RunTSNE(scRNA, dims = 1:PCNum, reduction = reduction)
  }

  if(length(resolution) > 1){
    print(clustree::clustree(scRNA@meta.data, prefix = ifelse(SCTransform,"SCT_snn_res.","RNA_snn_res.")))
  }

  print(DimPlot(scRNA))
  if(progress_saving){
    saveRDS(scRNA, "Clustered_scRNA.rds")
  }
  return(scRNA)
}


#' A simple pipeline for quick single cell analysis
#'
#' @param dataList Refer to corresponding functions
#' @param species 
#' @param SCTransform
#' @param DoubletFind
#' @param filter
#' @param progress_saving
#' @param CellCycleScoring
#' @param vars.to.regress
#' @param PCs
#' @param resolution 
#'
#' @returns
#' @export
#'
#' @examples
#' scRNA <- quick_single_cell(dataList = dataList, species = "Hs", SCTransform = T)
quick_single_cell <- function(dataList, species, SCTransform = F, DoubletFind = T, filter = NULL,
                              progress_saving = F, CellCycleScoring = F, vars.to.regress = c('percent.mt'),
                              PCs = NULL, resolution = seq(0.2,1,0.2), samples = NULL){

  dataList <- scRNA_preprocessing(dataList = dataList, species = species,  DoubletFind = DoubletFind, filter = filter, progress_saving = progress_saving, samples = samples)
  if(SCTransform){

    scRNA <- scRNA_SCTransform(dataList = dataList, vars.to.regress = vars.to.regress, progress_saving = progress_saving)
  }else{
    scRNA <- scRNA_Normalization_Reduction(dataList = dataList, species = species,vars.to.regress = vars.to.regress, CellCycleScoring = CellCycleScoring, progress_saving = progress_saving)
  }

  if(length(dataList) > 1){
    scRNA <- scRNA_Integration(Merge.Seurat = scRNA, SCTransform = SCTransform, progress_saving = progress_saving)
  }
    scRNA <- scRNA_clustering(scRNA = scRNA, PCs = PCs, resolution = resolution, reduction = ifelse(length(dataList) > 1, "harmony", "pca"), SCTransform = SCTransform, progress_saving = progress_saving)

  return(scRNA)

}



#' Sample Seurat Object
#'
#' @param scRNA A Seurat Object of a list of Seurat Objects
#' @param proportion The proportion of cells in the sample
#'
#' @returns
#' @export
#'
#' @examples
sample_object <- function(scRNA, proportion){

  if(class(scRNA) == "list"){

    for(i in 1:length(scRNA)){

      sce.tmp <- scRNA[[i]]
      sce.tmp <- sce.tmp[, sample(1:ncol(sce.tmp),round(ncol(sce.tmp)*proportion)) ]
      scRNA[[i]] <- sce.tmp
      rm(sce.tmp)

    }

  }else{

    scRNA <- scRNA[, sample(1:ncol(scRNA),round(ncol(scRNA)*proportion)) ]
  }

  return(scRNA)

}


#' Extract expression matrix from Seurat Object
#'
#' @param scRNA
#' @param Normalized
#'
#' @returns
#' @export
#'
#' @examples
Extract_Expression_Matrix <- function(scRNA, Normalized = F){

  lifecycle::deprecate_warn(what = "Extract_Expression_Matrix()", "GetAssayData(scRNA, assay = 'RNA', layer = 'counts')")

  require(Seurat)
  require(dplyr)

  version <- strsplit(as.character(scRNA@version), split = "[.]")[[1]]

  if( version[1] == "5"){

    scRNA[['RNA']] <- JoinLayers(scRNA[['RNA']])
  }

  if(Normalized){

    raw.data <- GetAssayData(scRNA, assay = "RNA", layer = "data")

  }else{

    raw.data <- GetAssayData(scRNA, assay = "RNA", layer = "counts")

  }


  Features(scRNA) %>% rownames() -> a
  row.names(raw.data) <- a
  rm(a)

  Cells(scRNA) %>% rownames() -> b
  colnames(raw.data) <- b
  rm(b)

  if(TRUE %in% lapply(dimnames(raw.data), is.null)){

    dimnames(raw.data) <- list(Features(scRNA), Cells(scRNA))
  }


  return(raw.data)

}


#' Score the specified genesets in Single cell object
#'
#' @param scRNA A Seurat object.
#' @param geneSets A list of genesets or a whole geneset dataframe categorized by column "term"
#' @param method A character either "AUCell" or "ModuleScore". "AUCell" is an extra R package and it is better to set seeds because it use random numbers.
#' "ModuleScore" use internal Seurat methods.
#' @import AUCell
#' @importFrom dplyr %>%
#'
#' @returns A Seurat object with activity in meta.data.
#' @export
#'
#' @examples
score_geneset_activity <- function(scRNA, geneSets, method = "AUCell"){

  if(class(geneSets) == "data.frame"){

    geneSets <- split(geneSets, geneSets$term)

  }

  if(!method %in% c("AUCell", "ModuleScore")){
    stop("Such method is not supported.")
  }

  geneSets_name <- lapply(geneSets, function(x){unique(x$term)})
  geneSets_name <- as.character(unlist(geneSets_name))
  geneSets_gene <- lapply(geneSets, function(x){x$gene})
  geneSets <- geneSets_gene
  names(geneSets) <- geneSets_name


  if(method == "AUCell"){

  cells_rankings <- AUCell_buildRankings(Extract_Expression_Matrix(scRNA, Normalized = T), plotStats=TRUE)

  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,aucMaxRank = nrow(cells_rankings)*0.05)


  set.seed(123)
  #cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)

  for(i in 1:length(geneSets_name)){

    AUCell_Score <- as.numeric(getAUC(cells_AUC)[geneSets_name[i], ])
    scRNA@meta.data[[geneSets_name[i]]] <- AUCell_Score
    print(geneSets_name[i])
  }

  return(scRNA)

  }

  if(method == "ModuleScore"){

    p <- ncol(scRNA@meta.data)

    require(Seurat)
    scRNA <- AddModuleScore(scRNA,
                            features = geneSets,
                            ctrl = 100,
                            name = "GeneSet")

    colnames(scRNA@meta.data)[(p+1):(p+length(geneSets))] <- geneSets_name



  }


  return(scRNA)

}


#' Running pseudobulk to find DEGs in scRNA-seq
#'
#' @param scRNA A Seurat object.
#' @param idents_to_check The target idents to be checked. Name is decided by group.
#' For example, idents are grouped by "feature1" (including "A" and "B") and "feature2" (including "a" and "b"),
#' then the idents are "A_a", "A_b", "B_a" and "B_b".
#' @param group A character vector specifying how to group cells. Elements must exist in metadata.
#' @param replicates A character vector specifying for each group, what are the replicates. Note that replicates should not be specified in "idents_to_check"
#' @param aggregate_expr Whether to aggregate the expression matrix as a pseudobulk matrix. TRUE is recommended, especially when replicates exist.
#' @param p.value The p.value cutoff for result. Default is 1.
#' @param min.pct The minimum percent of cells that express one gene. Default is 0.
#' @param logfc The logFC cutoff for result. Default is 0.
#' @param test.use Refer to FindMarkers().
#'
#' @returns
#' @export
#'
#' @examples
pseudoBulk <- function(scRNA, idents_to_check, group = NULL, replicates = NULL,aggregate_expr = T, p.value = 1, min.pct = 0,logfc = 0, test.use = "DESeq2"){

  if(!aggregate_expr){

    p <- do.call(paste, c(scRNA@meta.data[, group, drop = FALSE], sep = "_")) %>% as.character()
    scRNA$celltype.stim <- p
    Idents(scRNA) <- scRNA$celltype.stim
    bulk.mono.de <- FindMarkers(scRNA, ident.1 = idents_to_check[1], ident.2 = idents_to_check[2], min.pct = min.pct, logfc.threshold = logfc)

  }else{

    pseudo_ifnb <- AggregateExpression(scRNA,
                                       assays = "RNA",
                                       return.seurat = T,
                                       group.by = c(group, replicates))

    p <- do.call(paste, c(pseudo_ifnb@meta.data[, group, drop = FALSE], sep = "_")) %>% as.character()
    pseudo_ifnb$celltype.stim <- p
    Idents(pseudo_ifnb) <- "celltype.stim"
    bulk.mono.de <- FindMarkers(pseudo_ifnb,
                                ident.1 = idents_to_check[1],
                                ident.2 = idents_to_check[2],
                                test.use = test.use, min.pct = min.pct, logfc.threshold = logfc)
    bulk.mono.de <- na.omit(bulk.mono.de)
  }

  bulk.mono.de_sig <- bulk.mono.de[bulk.mono.de$p_val_adj <= p.value,]
  colnames(bulk.mono.de_sig)[2] <- "log2FoldChange"
  colnames(bulk.mono.de_sig)[5] <- "padj"
  return(bulk.mono.de_sig)
}


#' Quick annotate cell clusters.
#'
#' @param scRNA A Seurat object.
#' @param celltype The specific celltype.
#' @param pos The number of target seurat cluster (start from 0).
#'
#' @returns A Seurat object with metadata having a column "cell_type".
#' @export
#'
#' @examples
quick_manual_annotation <- function(scRNA, celltype, pos){


  if(is.null(scRNA@meta.data$cell_type)){

    scRNA@meta.data$cell_type <- rep("unclassified", nrow(scRNA@meta.data))

  }

  levels(scRNA@meta.data$cell_type) <- unique(c(levels(scRNA@meta.data$cell_type), celltype))

  scRNA@meta.data$cell_type[scRNA@meta.data$seurat_clusters %in% pos] <- celltype
  return(scRNA)
}


create_SingleR_ref <- function(scRNA, group.by){

  av <- AggregateExpression(scRNA , assays = "RNA", group.by = group.by)
  Ref = av[[1]]

  require(SingleCellExperiment)
  ref_sce = SingleCellExperiment::SingleCellExperiment(assays=list(counts=Ref))
  require(scater)
  ref_sce=scater::logNormCounts(ref_sce)
  colData(ref_sce)$Type=colnames(Ref)
  return(ref_sce)

}

#' Quick annotate cell clusters using SingleR
#'
#' @param scRNA A Seurat object.
#' @param ref Reference data.
#' @param labels Choose the labels in the reference data.
#' @param return_labeled If TRUE, return a Seurat object that is already annotated. If FALSE, return a table containing the annotation.
#'
#' @importFrom dplyr %>%
#' @returns
#' @export
#'
#' @examples
predict_SingleR <- function(scRNA, ref, labels, return_labeled = F){

  assay_for_SingleR <- GetAssayData(scRNA, assay = "RNA", layer = "data")
  pred <- SingleR::SingleR(test = assay_for_SingleR, ref = ref, labels = labels, assay.type.test = 1)

  cellType <- data.frame(seurat = scRNA@meta.data$seurat_clusters, pred = pred$labels)
  sort(table(cellType[,1]))

  temp <- as.data.frame(table(cellType[,1:2]))
  temp <- dplyr::group_by(temp,seurat)
  temp <- dplyr::top_n(x = temp, n = 1, wt = Freq)

  finalmap <- temp[order(temp$seurat),]$pred

  if(return_labeled){

    names(finalmap) <- levels(scRNA)
    scRNA <- RenameIdents(scRNA, finalmap)
    scRNA$cell_type <- Idents(scRNA)
    return(scRNA)
  }

  return(finalmap)

}

#' Draw ratio of different cell populations.
#'
#' @param Anno.by A vector from meta.data showing the annotation of cells.
#' @param group.by A vector from meta.data showing the condition of cells.
#' @param y.axis The title on the y.axis.
#' @param width width of columns in the figure.
#' @param return_table If TRUE, return a table of cell proportions. If FALSE, return a figure showing the ratio.
#'
#' @import ggplot2
#' @returns
#' @export
#'
#' @examples
#' p <- draw_ratio(scRNA$cell_type, scRNA$treatment)
draw_ratio <- function(anno.by, group.by, y.axis = "Fraction of clusters in each group", width = 0.5, return_table = F){

  Cellratio <- prop.table(table(anno.by, group.by), margin = 2)
  Cellratio <- as.data.frame(Cellratio)
  colnames(Cellratio)[1:2] <- c("cluster", "group")
  if(return_table) return(Cellratio)
  colourCount = length(unique(Cellratio$cluster))
  require(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  return(ggplot(data = Cellratio, aes(x =group, y = Freq, fill =  cluster)) +
           geom_bar(stat = "identity", width=width,position="fill")+
           scale_fill_manual(values=col_vector[1:20]) +
           theme_bw()+
           theme(panel.grid =element_blank()) +
           labs(x="",y= y.axis)+
           ####用来将y轴移动位置
           theme(axis.text.y = element_text(size=12, colour = "black"))+
           theme(axis.text.x = element_text(size=12, colour = "black"))+
           theme(
             axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
           ) )

}

#' Find markers for each cell cluster
#'
#' @param scRNA A Seurat object.
#' @param n Top n markers.
#' @param draw If TRUE, return a figure showing the markers. If FALSE, return the marker genes.
#' @param high_expr Whether the markers genes must have a high log fold change compared to other clusters.
#' @param SCTransform Whether the object is SCTransformed.
#' @param ... Refer to FindAllMarkers()
#'
#' @returns
#' @export
#'
#' @examples
Find_topn_markers <- function(scRNA, n = 10, draw = T, high_expr = T, SCTransform = F,...){

  if(SCTransform){
    scRNA <- PrepSCTFindMarkers(scRNA)
  }
  markers <- FindAllMarkers(scRNA, ...)
  require(dplyr)
  if(high_expr){

     markers <- markers[(markers$avg_log2FC > 0) & (markers$p_val_adj < 0.05),]

  }else{

    markers <- markers[markers$p_val_adj < 0.05,]
  }
  pct <- markers$pct.1
  logFC <- markers$avg_log2FC
  p_val_number <- -log10(markers$p_val_adj+1e-10)
  order <- pct^2 * logFC^2 * p_val_number
  markers$order <- order
  markers <- dplyr::group_by(markers, cluster)
  markers <- dplyr::top_n(x = markers, n = n, wt = order)

  markers <- markers[!duplicated(markers$gene), ]

  if(draw){

    return(DotPlot(scRNA, features = markers$gene)+theme(
      axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
    ))

  }else{

    return(markers$gene)
  }
}

SAVER_imputation <- function(scRNA, genes_to_check = NULL, estimate_size = deprecated()){

  require(SAVER)
  require(dplyr)

  if (lifecycle::is_present(estimate_size)) {

    lifecycle::deprecate_warn(what = "estimate_size", details = "Now this function uses 2000 high variable genes instead")

  }

  #expr <- Extract_Expression_Matrix(scRNA)
  expr <- GetAssayData(scRNA, assay = "RNA", layer = "counts")
  if(is.null(VariableFeatures(scRNA))){

    cat("No variable features detected. Finding variable features...")
    scRNA <- NormalizeData(scRNA,assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
    scRNA <- FindVariableFeatures(object = scRNA,assay = "RNA", mean.function = ExpMean, dispersion.function = LogVMR)

  }

  hvg <- VariableFeatures(scRNA)
  if(!is.null(genes_to_check)){

    gene_pos <- which(row.names(expr) %in% unique(c(genes_to_check, hvg)))
    imputed_data <- saver(expr, pred.genes = gene_pos, estimates.only = T)
  }else{

    imputed_data <- saver(expr, estimates.only = T)
  }

  gc()
  return(imputed_data)
}

#' Single Cell Draw Boxplot
#'
#' @param scRNA A Seurat object.
#' @param n Genes to check.
#' @param group.by A character denoting the group.
#'
#' @returns
#' @export
#'
#' @examples
sc_BoxPlot <- function(scRNA, genes_to_check, group.by){

  require(ggplot2)
  require(dplyr)
  #expr <- Extract_Expression_Matrix(scRNA, Normalized = T) %>% .[genes_to_check,] %>% as.data.frame() %>% t()
  expr <- GetAssayData(scRNA, assay = "RNA", layer = "data") %>% .[genes_to_check,] %>% as.data.frame() %>% t()
  group <- scRNA[,group.by,drop = TRUE]
  return(tinyarray::draw_boxplot(expr, group)&labs(x='')&
           theme(axis.title.x = element_blank(),
                 axis.text.x = element_text(color = 'black',face = "bold", size = 12),
                 axis.text.y = element_text(color = 'black', face = "bold"),
                 axis.title.y = element_text(color = 'black', face = "bold", size = 15),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(color="black",linewidth = 1.2, linetype="solid"),
                 panel.spacing = unit(0.12, "cm"),
                 plot.title = element_text(hjust = 0.5, face = "bold.italic")))
}


#' Quick Reader
#'
#' @param dir_path A file path denoting the position of single cell raw data, including files like "matrix", "feature", "barcode".
#' @param sep The symbol separating the files, e.g., in "GSM123456_Sample1_barcode", the sep is "_"
#' @param pattern_index The index of character that denotes the sample name after separation. e.g., in "GSM123456_Sample1_barcode", index is 2 
#'
#' @returns
#' @export
#'
#' @examples
quick_data_reader <- function(dir_path, sep, pattern_index, min.cells = 3, min.features = 200) {
  library(dplyr)
  library(fs)

  files <- list.files(dir_path)
  groups <- strsplit(files, split = sep) |>
    sapply(function(x) x[[pattern_index]]) |>
    unique()


  dir_create(path = file.path(dir_path, groups))
  dataList <- list()
  for (group in groups) {

    print(paste0("processing object ", group))
    files_to_move <- dir_ls(dir_path, regexp = group, type = "file")

    if(length(files_to_move) != 3) {  
      stop(paste0(group, " has ", length(files_to_move), " files (expected 3)"))
    }

    new_names <- sapply(files_to_move, function(f) {
      base <- basename(f)
      if (grepl("barcode", base)) return("barcodes.tsv.gz")
      if (grepl("matrix", base)|grepl("mtx", base))  return("matrix.mtx.gz")
      return("features.tsv.gz")
    })

    new_paths <- file.path(dir_path, group, new_names)
    file_move(files_to_move, new_paths)
    counts <- Seurat::Read10X(data.dir = file.path(dir_path, group))
    sce <- Seurat::CreateSeuratObject(
      counts = counts,
      min.cells = min.cells,
      min.features = min.features,
      project = group
    )
    dataList[[group]] <- sce
  }
  return(dataList)
}


#' Get correlation in single cell data
#'
#' @returns
#' @export
#'
#' @examples

sc_cor <- function(scRNA, genes.check, group.by, ident.group, reduction = NULL, target.group,k = 25, max_shared = 10, min.cells = 100, name = "test", method = "pearson"){

    if(is.null(reduction)) reduction <- if ("harmony" %in% names(scRNA@reductions)) "harmony" else "pca"
    if(length(genes.check) < 2){
      stop("Please check the gene.check input")
    }

    create_meta_cells <- function(scRNA, group.by, ident.group, k, max_shared, min.cells){
      wg <- hdWGCNA::SetupForWGCNA(
        scRNA,
        gene_select = "fraction", 
        fraction = 0.05, 
        wgcna_name = name 
      )

      wg <- hdWGCNA::MetacellsByGroups(
        seurat_obj = wg,
        group.by = group.by,
        reduction = reduction, 
        k = k, 
        max_shared = max_shared, 
        ident.group = ident.group, 
        min_cells = min.cells 
      )

      wg <- hdWGCNA::NormalizeMetacells(wg)
      tp <- hdWGCNA::GetMetacellObject(wg)
      saveRDS(tp, paste0(name, "_Wgcna.rds"))
      return(tp)
    }

  if(file.exists(paste0(name, "_Wgcna.rds"))){
    tp <- readRDS(paste0(name, "_Wgcna.rds"))
  }else{
    tp <- create_meta_cells(scRNA, group.by, ident.group, k = k, max_shared = max_shared, min.cells = min.cells)
  }
  
  cells<- WhichCells(tp, expression = !!rlang::sym(ident.group) %in% target.group)
  inter = intersect(Features(tp), genes.check)
  mat<- LayerData(tp, assay = "RNA", layer = "data")[inter, cells] %>% as.matrix() %>% t()
  if(length(genes.check) == 2){
    print(ggplot2::ggplot(data = mat, ggplot2::aes(x = .data[[genes.check[1]]], y = .data[[genes.check[2]]])) +
    ggplot2::geom_smooth(method = "lm", color="#4D79A6", formula = y~x, fill = "#A1CAE6") +
    ggplot2::theme_bw()+
    ggplot2::geom_point(colour = '#669BD2', size = 2)+
    ggpubr::stat_cor(method = method))
    return(mat)
  }else{
    cor_res<- Hmisc::rcorr(mat, type = method) 
    cor_res$r[is.na(cor_res$r)]<- 0
    corrplot::corrplot(cor_res$r, diag = FALSE, col = rev(corrplot::COL2('RdBu', 200)), method = "number",
         p.mat = cor_res$p, 
         insig = "pch",
         pch = 4)
    return(cor_res)
  }
  
}