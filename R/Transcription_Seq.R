
#' Analyze Differentially Expressed Genes From RNA-seq Data
#'
#' @description `r lifecycle::badge("deprecated")`
#' Please use [new_function()] instead.
#'
#' @param countData A standard expression matrix with rows the genes and columns the samples. The matrix can be raw count or normalized according to the method chose.
#' @param group A vector of groups of the sample. Note that only 2 groups are supported.
#' @param contrast A string with the format "A-B", in which "A" and "B" are to groups in the group parameter.
#'  The result demonstrates the gene expression in group A vs. that in group B if input is "A-B".
#'  If NULL, the result will denote unique(group)[1] vs. unique(group)[2].
#' @param p.value A number denotes the largest adjusted p value in the result table.
#' @param logFC A number denotes the minimum absolute log2 fold change in the result table.
#' @param method A string tells the method used for DEG analysis, can be "DESeq2", "edgeR", "limma", "wilcox" or "t".
#'  DESeq2 uses a negative binomial model. It accepts raw count data and is the default method.
#'  edgeR uses a negative binomial model as well but uses a different normalization method. It accepts raw count data.
#'  limma uses a linear model and uses Bayesian methods to reduce false positive results. It accepts log2-transformed raw count or TPM data.
#'  wilcox performs Wilcoxn test to every gene. Suitable for data with a large sample size. No require on data distribution.
#'  t performs t test to every gene. Suitable for data following a normal distribution and variances are homogeneous. Not recommended for most cases.
#' @param log_transformed A logical value denotes whether the input is log-transformed. Default is NULL and the function detects it automatically. Only when the detection fails should it be specified.
#'
#' @param FilterGene A logical value denotes whether filter genes with low expression to reduce false positive rate. Only works for limma, wilcox and t methods.
#' @returns A result table containing log2 fold change, p value and p value adjusted by FDR method.
#' @export
#'
#' @examples
#' res <- DEG_analysis(test, group = c("A", "A", "A", "B", "B", "B"), contrast = "B-A", method = "wilcox", Normalized = T, FilterGene = F)
DEG_analysis <- function(
  countData,
  group,
  contrast = NULL,
  p.value = 1,
  logFC = 0,
  method = "DESeq2",
  FilterGene = T,
  log_transformed = NULL
) {
  lifecycle::deprecate_warn(what = "DEG_analysis()", "DEG_analysis_v2()")

  degSig <- NULL

  if (is.null(contrast)) {
    cat("There is no contrast. Using")
    cat(paste(unique(group)[1], "vs.", unique(group)[2], sep = " "))
    contrast = paste(unique(group)[1], "-", unique(group)[2], sep = "")
  }

  if (method == "DESeq2") {
    print("Using DESeq2 for analysis...")
    if (any(countData %% 1 != 0)) {
      print(
        "Inputs for DESeq2 should be raw count. Please check your input and be careful with the result."
      )
    }

    contrast = strsplit(contrast, split = "-")
    contrast <- contrast[[1]]
    countData <- countData[rowMeans(countData) > 1, ]
    condition <- factor(group)
    colData <- data.frame(row.names = colnames(countData), condition)

    countData <- round(countData)
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = countData,
      colData = colData,
      design = ~condition
    )
    dds1 <- DESeq2::DESeq(
      dds,
      fitType = 'mean',
      minReplicatesForReplace = 7,
      parallel = FALSE
    )
    res <- DESeq2::results(
      dds1,
      contrast = c('condition', contrast[1], contrast[2])
    )
    res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

    res1 <- na.omit(res1)
    degSig <- res1[abs(res1$log2FoldChange) >= logFC & res1$padj <= p.value, ]
    degSig <- degSig[order(degSig$padj), ]

    #关于结果中的NA:如果在一行中，所有样本的计数都为零，则基础平均值（baseMean）列将为零，log2 FC、p值和调整后的p值都将被设置为NA
    #详见https://cloud.tencent.com/developer/article/2327198
  }

  if (method == "limma") {
    print("Using limma for analysis...")
    if (is.null(log_transformed)) {
      log_transformed <- median(countData[, 1, drop = TRUE]) < 5 &&
        max(countData[, 1, drop = TRUE]) < 50
    }

    if (!log_transformed) {
      print(
        "Inputs for limma should be log2-transformed. Please check your input and be careful with the result."
      )
    }

    design_df <- data.frame(colnames(countData), group)
    TS <- factor(design_df$group, unique(design_df$group))
    design <- stats::model.matrix(~ 0 + TS)
    colnames(design) <- c(unique(group)[1], unique(group)[2])

    if (FilterGene) {
      # keep = rowSums( countData >= 10 ) >= min(table(gs))
      # countData <- countData[keep,]
      countData <- countData[rowSums(countData) > 1, ]
    }

    fit <- limma::lmFit(countData, design)

    cont.matrix <- limma::makeContrasts(contrasts = contrast, levels = design)
    fit2 <- limma::contrasts.fit(fit, cont.matrix)
    fit2 <- limma::eBayes(fit2)
    degSig <- limma::topTable(
      fit2,
      number = Inf,
      adjust.method = "BH",
      p.value = p.value
    )

    degSig <- degSig[abs(degSig$logFC) >= logFC, ]
    colnames(degSig)[colnames(degSig) == "logFC"] <- "log2FoldChange"
    colnames(degSig)[colnames(degSig) == "P.Value"] <- "padj"
  }

  if (method == "edgeR") {
    print("Using edgeR for analysis...")
    if (any(countData %% 1 != 0)) {
      print(
        "Inputs for edgeR should be raw count. Please check your input and be careful with the result."
      )
    }

    contrast = strsplit(contrast, split = "-")
    contrast <- contrast[[1]]
    group <- factor(group, levels = rev(contrast))

    dge <- edgeR::DGEList(counts = countData, group = group)
    keep <- edgeR::filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge$samples$lib.size <- base::colSums(dge$counts)
    dge <- edgeR::calcNormFactors(dge)

    design <- stats::model.matrix(~ 0 + group)
    rownames(design) <- colnames(dge)
    colnames(design) <- rev(contrast)

    dge <- edgeR::estimateDisp(dge, design)
    fit <- edgeR::glmQLFit(dge, design)
    lrt <- edgeR::glmQLFTest(fit, contrast = c(-1, 1))
    nrDEG <- edgeR::topTags(lrt, n = Inf)
    deg <- as.data.frame(nrDEG)
    degSig <- deg[abs(deg$logFC) >= logFC & deg$FDR <= p.value, ]

    colnames(degSig)[1] <- "log2FoldChange"
    colnames(degSig)[5] <- "padj"
  }

  if (method == "wilcox" | method == "t") {
    print("Inputs for wilcox or t should be normalized.")
    if (is.null(log_transformed)) {
      log_transformed <- median(countData[, 1, drop = TRUE]) < 5 &&
        max(countData[, 1, drop = TRUE]) < 50
      if (log_transformed) cat("Inputs seem to be log2-transformed.")
    }

    #Independent filtering
    countData <- edgeR::DGEList(counts = countData, group = group)

    if (FilterGene) {
      keep <- edgeR::filterByExpr(countData)
      countData <- countData[keep, , keep.lib.sizes = FALSE]
    }

    countData <- as.matrix(countData)
    require(future.apply)
    contrast = strsplit(contrast, split = "-")
    contrast <- contrast[[1]]

    groups <- unique(group)
    if (length(groups) != 2) {
      stop("Two groups are required")
    }
    group_order <- if (groups[1] == contrast[1]) groups else rev(groups)

    group1_idx <- which(group == group_order[1])
    group2_idx <- which(group == group_order[2])
    exact <- min(length(group1_idx), length(group2_idx)) <= 50

    if (!log_transformed) {
      countData <- log2(countData + 1)
    }
    logFC_vec <- rowMeans(countData[, group1_idx]) -
      rowMeans(countData[, group2_idx])

    if (method == "wilcox") {
      print("Using wilcox for analysis...")
      suppressWarnings({
        pvals <- vapply(
          1:nrow(countData),
          function(i) {
            stats::wilcox.test(
              countData[i, group1_idx],
              countData[i, group2_idx],
              exact = exact
            )$p.value
          },
          numeric(1)
        )
      })
    } else {
      print("Using student's t for analysis...")
      suppressWarnings({
        pvals <- vapply(
          1:nrow(countData),
          function(i) {
            stats::t.test(
              countData[i, group1_idx],
              countData[i, group2_idx]
            )$p.value
          },
          numeric(1)
        )
      })
    }

    results <- data.frame(p.value = pvals, log2FoldChange = logFC_vec)
    results$padj <- stats::p.adjust(results$p.value, method = "fdr")
    results <- na.omit(results)
    degSig <- results[
      abs(results$log2FoldChange) >= logFC & results$padj <= p.value,
    ]
    degSig <- degSig[order(degSig$padj), ]
  }

  return(degSig)
}


#' Perform quick GSEA analysis
#'
#' @param geneList A named geneList including the ranking metrics of all genes with names the gene name
#' @param genesets The geneset for GSEA analysis. Can be common genesets including "GO", "KEGG", "Reactome", or custom gmtfile.
#' @param species "Hs" for human and "Mm" for mouse.
#' @param from The gene format of the input data.
#' @param p.value The cut off of GSEA analysis.
#' @param minGSSize The minimum geneset size for analysis.
#' @param maxGSSize The maximum geneset size for analysis.
#'
#' @returns
#' @export
#'
#' @examples
quick_GSEA <- function(
  geneList,
  genesets = "GO",
  species,
  from,
  p.value = 0.05,
  minGSSize = 10,
  maxGSSize = 500
) {
  if (!species %in% c("Mm", "Hs")) {
    stop("Such species is not supported")
  }

  if (class(genesets) != "character") {
    if (from != "SYMBOL") {
      names(geneList) <- Quick_ID_conversion(
        names(geneList),
        species = species,
        from = from,
        to = "SYMBOL",
        matrix = F
      )
    }
    result <- clusterProfiler::GSEA(
      geneList,
      TERM2GENE = genesets,
      pvalueCutoff = p.value
    )
  } else {
    if (from != "ENTREZID") {
      names(geneList) <- Quick_ID_conversion(
        names(geneList),
        species = species,
        from = from,
        to = "ENTREZID",
        matrix = F
      )
    }
    if (genesets == "GO") {
      result <- clusterProfiler::gseGO(
        geneList,
        ifelse(
          species == "Hs",
          'org.Hs.eg.db::org.Hs.eg.db',
          'org.Mm.eg.db::org.Mm.eg.db'
        ),
        keyType = "ENTREZID",
        ont = "ALL",
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        pvalueCutoff = p.value
      )
      #使用GSEA进行GO富集分析（'org.Hs.eg.db::org.Hs.eg.db'：对应物种的数据库；ont：选择输出条目，可选“BP,MF,CC或者ALL”，pvalueCutoff：设置P的阈值）
    } else if (genesets == "KEGG") {
      result <- clusterProfiler::gseKEGG(
        geneList,
        organism = ifelse(species == "Hs", "hsa", "mmu"),
        minGSSize = minGSSize, #默认 10
        maxGSSize = maxGSSize,
        pvalueCutoff = p.value
      ) #使用GSEA进行KEGG富集分析
    } else if (genesets == "Reactome") {
      result <- clusterProfiler::gsePathway(
        gene = geneList,
        organism = ifelse(species == "Mm", 'mouse', 'human'), #物种选择 "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
        minGSSize = minGSSize, #默认 10
        maxGSSize = maxGSSize, #默认 500
        pAdjustMethod = "fdr",
        pvalueCutoff = p.value
      )
    } else {
      stop("Please Check the geneset you inputed.")
    }
  }
  return(result)
}

#' Perform quick enrich analysis
#'
#' @param geneList A geneList including names of significant differentially expressed genes
#' @param genesets The geneset for enrich analysis. Can be "GO", "KEGG", and "Reactome"
#' @param species "Hs" for human and "Mm" for mouse.
#' @param from The gene format of the input data.
#' @param universe The background geneset for ORA analysis.
#' @param q.value The cut off adjusted p for enrich analysis.
#' @param minGSSize The minimum geneset size for analysis.
#' @param maxGSSize The maximum geneset size for analysis.
#'
#' @returns
#' @export
#'
#' @examples
quick_enrich <- function(
  geneList,
  genesets = "GO",
  species,
  from,
  universe = NULL,
  p.value = 0.05,
  minGSSize = 10,
  maxGSSize = 500
) {
  if (is.null(intersect(geneList, universe))) {
    stop("genes in the geneList should be part of the universe genes")
  }

  if (from != "ENTREZID") {
    geneList <- Quick_ID_conversion(
      geneList,
      species = species,
      from = from,
      to = "ENTREZID",
      matrix = F
    )
    universe <- Quick_ID_conversion(
      universe,
      species = species,
      from = from,
      to = "ENTREZID",
      matrix = F
    )
  }

  if (genesets == "GO") {
    GO_database <- ifelse(
      species == "Mm",
      'org.Mm.eg.db::org.Mm.eg.db',
      'org.Hs.eg.db::org.Hs.eg.db'
    ) #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb

    return(clusterProfiler::enrichGO(
      geneList, #GO富集分析
      OrgDb = GO_database,
      keyType = "ENTREZID", #设定读取的gene ID类型
      ont = "ALL", #(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
      universe = universe,
      minGSSize = minGSSize, #默认 10
      maxGSSize = maxGSSize, #默认 500
      pAdjustMethod = "fdr",
      pvalueCutoff = p.value, #设定q值阈值
      readable = T
    ))
  } else if (genesets == "KEGG") {
    KEGG_database <- ifelse(species == "Mm", 'mmu', 'hsa') #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
    return(clusterProfiler::enrichKEGG(
      geneList, #KEGG富集分析
      organism = KEGG_database,
      universe = universe,
      minGSSize = minGSSize, #默认 10
      maxGSSize = maxGSSize, #默认 500
      pAdjustMethod = "fdr",
      pvalueCutoff = p.value
    ))
  } else if (genesets == "Reactome") {
    return(clusterProfiler::enrichPathway(
      gene = geneList, #表示前景基因，即待富集的基因列表;[,1]表示对entrezid_all数据集的第1列进行处理
      organism = ifelse(species == "Mm", 'mouse', 'human'), #物种选择 "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
      universe = universe,
      minGSSize = minGSSize, #默认 10
      maxGSSize = maxGSSize, #默认 500
      pAdjustMethod = "fdr", # 指定多重假设检验矫正的方法，选项包含 "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
      pvalueCutoff = p.value, #指定 p 值阈值（指定 1 以输出全部），默认为0.05
      readable = T
    ))
  } else {
    stop("This Geneset is Not supported")
  }
}


#' Perform Normalization to RNAseq data
#'
#' @param countData A standard count matrix with rows the genes and columns the samples.
#' @param method Can be "cpm", "rpkm", "tpm", "tmm" or "vst/rlog".
#' @param species "Hs" for human and "Mm" for mouse. Only works when gene.length = NULL.
#' @param group Used for vst and tmm method.
#' @param from The gene format of the input data. Only works when gene.length = NULL.
#' @param gene.length The length of genes that required to be normalized. If NULL, use internal data to get gene length.
#'
#' @returns
#' @export
#'
#' @examples
normalization <- function(countData, method, group = NULL, gene.length = NULL) {
  if (!method %in% c("cpm", "rpkm", "tpm", "tmm", "vst", "rlog", "vst/rlog")) {
    stop("Such normalization is not supported")
  }

  if (method == "cpm") {
    return(edgeR::cpm(countData))
  } else if (method == "tmm") {
    y <- edgeR::DGEList(counts = countData, group = group)
    y <- edgeR::calcNormFactors(y, method = "TMM")
    tmm <- edgeR::cpm(y, log = F, prior.count = 0)
    return(tmm)
  } else if (method == "vst" | method == "rlog" | method == "vst/rlog") {
    blind = FALSE
    if (is.null(group)) {
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData,
        data.frame(row.names = colnames(countData)),
        design = ~1
      )
      blind <- TRUE
    } else {
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData,
        data.frame(condition = factor(group), row.names = colnames(countData)),
        design = ~condition
      )
    }

    if (ncol(countData) <= 30) {
      print("Small sample size. Performing rlog normalization.")
      rg <- DESeq2::rlog(dds, blind = blind)
      return(SummarizedExperiment::assay(rg))
    } else {
      print("Big sample size. Performing vst normalization.")
      vsd <- DESeq2::vst(dds, blind = blind)
      return(SummarizedExperiment::assay(vsd))
    }
  } else if (method == "rpkm") {
    if (is.null(gene.length)) {
      stop("rpkm requires gene.length input")
    }
    return(edgeR::rpkm(countData, gene.length = gene.length))
  } else if (method == "tpm") {
    if (is.null(gene.length)) {
      stop("tpm requires gene.length input")
    }
    rpkm <- edgeR::rpkm(countData, gene.length = gene.length)
    return(t(t(rpkm) / colSums(rpkm, na.rm = TRUE) * 1e6))
  }
}

#' To get gene.length from hg39 when the RNAseq length info is missing (not recommended)
get_standard_gene_length <- function(countData, species, from) {
  geneid_efflen <- if (species == "Mm") Mm_ref else Hs_ref
  colnames(geneid_efflen)[1] <- "Geneid"
  if (from != "ENSEMBL") {
    geneid_efflen <- Quick_ID_conversion(
      geneid_efflen,
      species = species,
      from = "ENSEMBL",
      to = from,
      matrix = T
    )
  }
  gene.length <- geneid_efflen[
    match(row.names(countData), geneid_efflen$Geneid),
    "efflen"
  ]
  if (any(is.na(gene.length))) {
    warning("There are NAs in the result. Please be careful when using it.")
  }

  return(gene.length)
}


#' Get Signal-to-Noise ranking metrics for GSEA
#'
#' @param tpm A tpm-normalized expression matrix
#' @param group A vector of groups of the sample. Note that only 2 groups are supported.
#' @param contrast A string with the format "A-B", in which "A" and "B" are to groups in the group parameter.
#'  The result demonstrates the gene expression in group A vs. that in group B if input is "A-B".
#'  If NULL, the result will denote unique(group)[1] vs. unique(group)[2].
#' @param abs Whether to get the absolute values.
#'
#' @returns A named, decreasingly ordered numeric vector.
#' @export
#'
#' @examples
s2nCalulator <- function(tpm, group, contrast = NULL, abs = F) {
  single_gene_s2n <- function(expr, group, contrast, abs) {
    group <- factor(group, levels = unique(group))

    if (is.null(contrast)) {
      message(
        "No contrast specified. Using ",
        levels(group)[1],
        " vs. ",
        levels(group)[2]
      )
      contrast <- paste(levels(group)[1], "-", levels(group)[2], sep = "")
    }

    contrast_parts <- strsplit(contrast, "-")[[1]]
    if (length(contrast_parts) != 2) {
      stop("contrast must be in format 'group1-group2'")
    }

    if (!all(contrast_parts %in% levels(group))) {
      stop("All groups in contrast must be present in the group factor")
    }

    g1 <- expr[which(group == contrast_parts[1])]
    g2 <- expr[which(group == contrast_parts[2])]

    diff <- mean(g1) - mean(g2)
    a <- ifelse(abs, abs(diff), diff)
    b <- max(0.2 * max(mean(g1), 1), sd(g1)) +
      max(0.2 * max(mean(g2), 1), sd(g2))

    return(a / b)
  }

  s2n <- apply(
    tpm,
    1,
    single_gene_s2n,
    group = group,
    contrast = contrast,
    abs = abs
  )
  names(s2n) <- row.names(tpm)
  s2n <- s2n[order(s2n, decreasing = T)]
  return(s2n)
}


#' Draw Volcano plot for RNAseq result
#'
#' @param deg Result of DEG_analysis(), or a data.frame with at least two columns "log2FoldChange" and "padj".
#' @param p.value The p.value cutoff for identifying the significant DEGs.
#' @param logFC The logFC cutoff for identifying the significant DEGs.
#' @param color A vector of three colors denoting the down regulated, not significant, and up regulated dots.
#' @param label A vector of genes that require to show their names.
#' @param highlight A vector of genes that require to highlight their dots.
#' @param highlight_color Specify the color of highlight.
#' @param max.overlaps refer to ggrepel::geom_label_repel()
#'
#' @returns
#' @export
#'
#' @examples
draw_Volcano <- function(
  deg,
  p.value = 0.05,
  logFC = 1,
  color = c("#6a82ed", "grey", "#ed6f6f"),
  alpha = 0.3,
  size = 2,
  label = NULL,
  highlight = NULL,
  highlight_color = "#07a818",
  max.overlaps = 10
) {
  deg <- as.data.frame(deg)
  deg$regulate <- with(
    deg,
    ifelse(
      padj < p.value & abs(log2FoldChange) >= logFC,
      ifelse(log2FoldChange > 0, "up-regulated", "down-regulated"),
      "not significant"
    )
  )

  p <- ggplot(deg, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = alpha, size = size, aes(color = regulate)) +
    scale_color_manual(values = color) +
    ylab("-log10(q-value)") +
    geom_vline(
      xintercept = c(-logFC, logFC),
      lty = 4,
      col = "black",
      lwd = 0.8
    ) +
    geom_hline(
      yintercept = -log10(p.value),
      lty = 4,
      col = "black",
      lwd = 0.8
    ) +
    theme_bw()

  if (!is.null(highlight)) {
    highlight_data <- deg[row.names(deg) %in% highlight, ]
    p <- p +
      geom_point(
        data = highlight_data,
        aes(x = log2FoldChange, y = -log10(padj)),
        colour = highlight_color,
        size = size
      )
  }

  if (!is.null(label)) {
    label_data <- deg[row.names(deg) %in% label, ]
    p <- p +
      ggrepel::geom_text_repel(
        min.segment.length = 0,
        data = label_data,
        aes(label = row.names(label_data)),
        box.padding = unit(0.4, "lines"),
        segment.color = "black",
        segment.size = 0.6,
        max.overlaps = max.overlaps
      )
  }

  return(p)
}


#' Draw pca plot
#'
#' @param data A data.frame with row the variables and the column samples.
#' @param group A vector specifies the group of the samples in data.
#' @param ncp A number specifies the principal components in the result. Default is 2.
#' @param label,mean.point,addEllipses,palette refer to factoextra::fviz_pca_ind
#'
#' @returns
#' @export
#'
#' @examples
draw_pca <- function(
  data,
  group = NULL,
  ncp = 2,
  draw = T,
  label = "none",
  addEllipses = T,
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  geom.ind = c("point")
) {
  if (ncp < 2) {
    stop("ncp must be larger than 1")
  }

  if (all(abs(data - round(data)) < 1e-10)) {
    print("Inputs seem to be count matrix. Performing normalization...")
    data <- normalization(data, "vst/rlog", group = group)
  }

  data <- t(data)
  pca <- FactoMineR::PCA(data, ncp = ncp, scale.unit = TRUE, graph = FALSE)
  eig1 <- round(pca$eig[1, 2], 2)
  eig2 <- round(pca$eig[2, 2], 2)
  print(paste0("The first pc explains ", eig1, " variance"))
  print(paste0("The second pc explains ", eig2, " variance"))

  if (!draw) {
    return(pca)
  }

  return(factoextra::fviz_pca_ind(
    pca,
    label = label,
    col.ind = group,
    mean.point = F,
    palette = palette,
    addEllipses = addEllipses,
    geom.ind = geom.ind
  ))
}

#' Draw heatmap
#'
#' A quick tool for drawing heatmap.
#'
#' @param data A data.frame/matrix with row the variables and the column samples.
#' @param group A vector specifies the group of the samples in data.
#' @param scale Whether to perform Z-score scale.
#' @param cutoff The cutoff for color shown. e.g., 0.1 means the lowest 10 percent/highest 10 percent will be shown as the darkest color.
#' @param ... Other parameters for ComplexHeatmap::Heatmap
#'
#' @return A Heatmap object.
#' @export
#'
#' @examples
draw_heatmap <- function(data, group = NULL, scale = TRUE, cutoff = 0.1, ...) {
  if (all(abs(data - round(data)) < 1e-10)) {
    print("Inputs seem to be count matrix. Performing normalization...")
    data <- normalization(data, "vst/rlog", group = group)
  }

  if (scale) {
    data <- t(scale(t(data)))
  }
  data[is.na(data)] <- 0
  max.c <- quantile(data, 1 - cutoff)
  min.c <- quantile(data, cutoff)
  col_fun = circlize::colorRamp2(
    c(min.c, 0, max.c),
    c("#4DBBD5CC", "white", "#E64B35CC")
  )
  hp <- ComplexHeatmap::Heatmap(
    data,
    col = col_fun,
    column_dend_gp = grid::gpar(lwd = 1.5, col = "black"),
    column_title_gp = grid::gpar(fontface = "bold", fontsize = 14),
    row_title_gp = grid::gpar(fontface = "bold", fontsize = 14),
    row_names_gp = grid::gpar(fontface = "bold"),
    ...
  )
  ComplexHeatmap::draw(hp)
  return(hp)
}


#' Analyze Differentially Expressed Genes From RNA-seq Data
#'
#' @description
#' Get DE genes using one of several common methods.
#'
#' @param Input A standard expression matrix with rows the genes and columns the samples. The matrix can be raw count or normalized according to the method chose.
#' @param group A vector of groups of the sample. Note that only 2 groups are supported.
#' @param contrast A string with the format "A-B", in which "A" and "B" are to groups in the group parameter.
#'  The result demonstrates the gene expression in group A vs. that in group B if input is "A-B".
#'  If NULL, the result will denote unique(group)[1] vs. unique(group)[2].
#' @param p.value A number denotes the largest adjusted p value in the result table.
#' @param logFC A number denotes the minimum absolute log2 fold change in the result table.
#' @param method A string tells the method used for DEG analysis, can be "DESeq2", "edgeR", "limma", "wilcox" or "t".
#'  DESeq2 uses a negative binomial model. It accepts raw count data and is the default method.
#'  edgeR uses a negative binomial model as well but uses a different normalization method. It accepts raw count data.
#'  limma uses a linear model and uses Bayesian methods to reduce false positive results. It accepts log2-transformed raw count or TPM data.
#'  wilcox performs Wilcoxn test to every gene. Suitable for data with a large sample size. No require on data distribution.
#'  t performs t test to every gene. Suitable for data following a normal distribution and variances are homogeneous. Not recommended for most cases.
#' @param FilterGene A logical value denotes whether filter genes with low expression to reduce false positive rate. Not work for edgeR.
#' @param log_transformed A logical value denotes whether the input is log-transformed. Default is NULL and the function detects it automatically. Only when the detection fails should it be specified.
#' @returns A result table containing log2 fold change, p value and p value adjusted by FDR method.
#' @export
#'
#' @examples
#' res <- DEG_analysis_v2(test, group = c("A", "A", "A", "B", "B", "B"), contrast = "B-A", method = "wilcox", FilterGene = F)
DEG_analysis_v2 <- function(
  Input,
  group,
  contrast = NULL,
  p.value = 1,
  logFC = 0,
  method = "DESeq2",
  FilterGene = T,
  log_transformed = NULL
) {
  #Check parameters
  if (!is.matrix(Input) && !is.data.frame(Input)) {
    stop("Input must be a matrix or data.frame")
  }

  if (length(group) != ncol(Input)) {
    stop("Length of group must equal number of columns in Input")
  }

  if (!method %in% c("DESeq2", "limma", "edgeR", "wilcox", "t")) {
    stop("method must be one of: 'DESeq2', 'limma', 'edgeR', 'wilcox', 't'")
  }

  group <- factor(group, levels = unique(group))

  if (is.null(contrast)) {
    message(
      "No contrast specified. Using ",
      levels(group)[1],
      " vs. ",
      levels(group)[2]
    )
    contrast <- paste(levels(group)[1], "-", levels(group)[2], sep = "")
  }

  contrast_parts <- strsplit(contrast, "-")[[1]]
  if (length(contrast_parts) != 2) {
    stop("contrast must be in format 'group1-group2'")
  }

  if (!all(contrast_parts %in% levels(group))) {
    stop("All groups in contrast must be present in the group factor")
  }

  Input <- check_input(Input, method, log_transformed)
  raw_results <- switch(
    method,
    DESeq2 = run_deseq2(Input, group, contrast_parts, FilterGene),
    limma = run_limma(Input, group, contrast, FilterGene),
    edgeR = run_edger(Input, group, contrast_parts),
    wilcox = run_wilcox(Input, group, contrast_parts, FilterGene),
    t = run_t(Input, group, contrast_parts, FilterGene)
  )

  return(format_results(raw_results, method, p.value, logFC))
}

check_input <- function(Input, method, log_transformed) {
  is_integer_matrix <- function(mat) {
    diff <- mat - round(mat)
    return(all(abs(diff) < 1e-10, na.rm = TRUE))
  }
  is_raw_count <- is_integer_matrix(Input) && max(Input) > 1000

  if (is.null(log_transformed)) {
    qx <- as.numeric(quantile(
      Input,
      c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
      na.rm = T
    ))
    log_transformed <- !((qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0))
  }
  print(paste(
    "Inputs seem to be",
    ifelse(
      is_raw_count,
      "raw counts",
      ifelse(log_transformed, "log-transformed", "normalized")
    )
  ))

  if (method %in% c("DESeq2", "edgeR") && !is_raw_count) {
    warning(
      "Inputs for ",
      method,
      " should be raw counts. ",
      "Please check your input data."
    )
  }

  if (method == "limma" && !log_transformed) {
    warning(
      "Inputs for limma should be log2-transformed. Perform log transformation."
    )
    Input <- log2(Input + 1)
  }

  if (method == "wilcox") {
    if (is_raw_count) {
      warning(
        "Inputs for ",
        method,
        " are recommended to be normalized or log transformed."
      )
    }
    if (!log_transformed) Input <- log2(Input + 1)
  }

  if (method == "t" && !log_transformed) {
    warning(
      "Inputs for ",
      method,
      " are recommended to be log transformed. Perform log transformation"
    )
    Input <- log2(Input + 1)
  }

  return(Input)
}

format_results <- function(results, method, p.value, logFC) {
  col_map <- list(
    DESeq2 = c(log2FoldChange = "log2FoldChange", padj = "padj"),
    limma = c(log2FoldChange = "logFC", padj = "adj.P.Val"),
    edgeR = c(log2FoldChange = "logFC", padj = "FDR"),
    wilcox = c(log2FoldChange = "log2FoldChange", padj = "padj"),
    t = c(log2FoldChange = "log2FoldChange", padj = "padj")
  )

  if (method %in% names(col_map)) {
    mapping <- col_map[[method]]
    for (new_name in names(mapping)) {
      old_name <- mapping[new_name]
      if (old_name %in% colnames(results)) {
        colnames(results)[colnames(results) == old_name] <- new_name
      }
    }
  }

  if (!is.null(results)) {
    results <- results[
      abs(results$log2FoldChange) >= logFC &
        results$padj <= p.value,
      ,
      drop = FALSE
    ]
    results <- results[order(results$padj, -abs(results$log2FoldChange)), ]
  }

  return(results)
}

run_deseq2 <- function(Input, group, contrast_parts, FilterGene) {
  print("Using DESeq2 for analysis...")
  if (FilterGene) {
    Input <- Input[rowMeans(Input) > 1, ]
  }
  condition <- factor(group)
  colData <- data.frame(row.names = colnames(Input), condition)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = Input,
    colData = colData,
    design = ~condition
  )
  dds1 <- DESeq2::DESeq(
    dds,
    fitType = 'mean',
    minReplicatesForReplace = 7,
    parallel = FALSE
  )
  res <- DESeq2::results(
    dds1,
    contrast = c('condition', contrast_parts[1], contrast_parts[2])
  )
  return(data.frame(res, stringsAsFactors = FALSE, check.names = FALSE))
}

run_limma <- function(Input, group, contrast, FilterGene) {
  print("Using limma for analysis...")
  if (FilterGene) {
    Input <- Input[rowMeans(Input) > 1, ]
  }

  design_df <- data.frame(colnames(Input), group)
  TS <- factor(design_df$group, unique(design_df$group))
  design <- stats::model.matrix(~ 0 + TS)
  colnames(design) <- c(unique(group)[1], unique(group)[2])

  fit <- limma::lmFit(Input, design)
  cont.matrix <- limma::makeContrasts(contrasts = contrast, levels = design)
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2)
  return(limma::topTable(fit2, number = Inf, adjust.method = "BH", p.value = 1))
}

run_edger <- function(Input, group, contrast_parts) {
  print("Using edgeR for analysis...")
  group <- factor(group, levels = rev(contrast_parts))
  dge <- edgeR::DGEList(counts = Input, group = group)
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge$samples$lib.size <- base::colSums(dge$counts)
  dge <- edgeR::calcNormFactors(dge)
  design <- stats::model.matrix(~ 0 + group)
  rownames(design) <- colnames(dge)
  colnames(design) <- rev(contrast_parts)

  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design)
  lrt <- edgeR::glmQLFTest(fit, contrast = c(-1, 1))
  nrDEG <- edgeR::topTags(lrt, n = Inf)
  return(as.data.frame(nrDEG))
}

run_wilcox <- function(Input, group, contrast_parts, FilterGene) {
  if (FilterGene) {
    Input <- Input[rowMeans(Input) > 1, ]
  }

  group1_idx <- which(group == contrast_parts[1])
  group2_idx <- which(group == contrast_parts[2])

  group1_means <- rowMeans(Input[, group1_idx, drop = FALSE])
  group2_means <- rowMeans(Input[, group2_idx, drop = FALSE])
  logFC_vec <- group1_means - group2_means
  exact <- min(length(group1_idx), length(group2_idx)) <= 50

  print("Using wilcox for analysis...")
  suppressWarnings({
    pvals <- vapply(
      1:nrow(Input),
      function(i) {
        stats::wilcox.test(
          Input[i, group1_idx],
          Input[i, group2_idx],
          exact = exact
        )$p.value
      },
      numeric(1)
    )
  })
  res <- data.frame(p.value = pvals, log2FoldChange = logFC_vec)
  res$padj <- stats::p.adjust(res$p.value, method = "fdr")
  return(res)
}

run_t <- function(Input, group, contrast_parts, FilterGene) {
  if (FilterGene) {
    Input <- Input[rowMeans(Input) > 1, ]
  }
  temp <- genefilter::rowttests(Input, factor(group, levels = contrast_parts))
  logFC_vec <- temp$dm
  pvals <- temp$p.value
  res <- data.frame(
    p.value = pvals,
    log2FoldChange = logFC_vec,
    row.names = row.names(temp)
  )
  res$padj <- stats::p.adjust(res$p.value, method = "fdr")
  return(res)
}


#' Draw Venn plot
#'
#' @param venn_list A list of vectors for overlapping.
#' @param color A list of colors for different groups
#' @param ... Other parameters for ggvenn::ggvenn
#'
#' @returns
#' @export
#'
#' @examples
draw_Venn <- function(venn_list, color = NULL, ...) {
  if (is.null(color)) {
    color <- RColorBrewer::brewer.pal(length(venn_list), "Spectral")
  }
  return(
    ggvenn::ggvenn(
      data = venn_list, # 数据列表
      columns = NULL, # 对选中的列名绘图，最多选择4个，NULL为默认全选
      show_elements = F, # 当为TRUE时，显示具体的交集情况，而不是交集个数
      label_sep = "\n", # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
      show_percentage = T, # 显示每一组的百分比
      digits = 1, # 百分比的小数点位数
      fill_color = color, # 填充颜色
      fill_alpha = 0.8, # 填充透明度
      stroke_color = "white", # 边缘颜色
      stroke_alpha = 0.5, # 边缘透明度
      stroke_size = 0.5, # 边缘粗细
      stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
      set_name_color = "black", # 组名颜色
      set_name_size = 5, # 组名大小
      text_color = "black", # 交集个数颜色
      text_size = 4, # 交集个数文字大小
      ...
    )
  )
}
