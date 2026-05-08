#' Perform KEGG metabolite enrichment
#'
#' @param KEGGid A character vector of KEGG compound IDs (e.g., "C00160")
#' @param species A character string, KEGG species code (e.g., "hsa", "eco")
#' @param p.adjust.method Adjustment method for multiple testing, default is "BH"
#' @param population.id An integer denoting the number of all the compounds that is mapped by KEGG, used for targeted metabolomics
#'
#' @return A data.frame with enrichment results
#' @export
kegg_metabolite_enrichment <- function(KEGGid, species, p.adjust.method = "BH", population.id = NULL) {

  cache_dir <- file.path(tempdir(), "kegg_cache")
  dir.create(cache_dir, showWarnings = FALSE)

  # 缓存通路列表
  pathways_cache <- file.path(cache_dir, paste0(species, "_pathways.rds"))
  if (file.exists(pathways_cache)) {
    pathways <- readRDS(pathways_cache)
  } else {
    pathways <- KEGGREST::keggList("pathway", species)
    saveRDS(pathways, pathways_cache)
  }

  path_ids <- sub("path:", "", names(pathways))

  pathway2compound <- list()

  pb <- txtProgressBar(max = length(path_ids), style = 3)

  for (i in seq_along(path_ids)) {
    pid <- path_ids[i]

    pid_cache <- file.path(cache_dir, paste0(pid, ".rds"))
    entry <- if (file.exists(pid_cache)) {
      readRDS(pid_cache)
    } else {
      tryCatch({
        e <- KEGGREST::keggGet(pid)[[1]]
        saveRDS(e, pid_cache)
        e
      }, error = function(e) NULL)
    }

    if (!is.null(entry) && !is.null(entry$COMPOUND)) {
      compounds <- names(entry$COMPOUND)
      pathway2compound[[pid]] <- compounds
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  all_cmpds <- unique(unlist(pathway2compound))

  results_list <- list()  # 改为列表收集

  N <- length(KEGGid)
  M <- ifelse(is.null(population.id),length(all_cmpds), length(population.id))

  for (pid in names(pathway2compound)) {
    pw_cmpds <- pathway2compound[[pid]]
    overlap <- intersect(KEGGid, pw_cmpds)
    k <- length(overlap)
    n <- ifelse(is.null(population.id),length(pw_cmpds), length(intersect(pw_cmpds, population.id)))
    if (k * n > 0) {
      p <- phyper(k - 1, n, M - n, N, lower.tail = FALSE)
      enrich_ratio <- (k * M) / (n * N)
      matched <- paste(overlap, collapse = "; ")

      results_list[[length(results_list) + 1]] <- data.frame(
        ID = pid,
        Description = pathways[pid],
        GeneRatio = paste0(k, "//", N),
        BgRatio = paste0(n, "//", M),
        pvalue = p,
        EnrichmentRatio = enrich_ratio,
        geneID = matched,
        Count = k
      )
    }
  }

  results <- do.call(rbind, results_list)

  results$p.adjust <- p.adjust(results$pvalue, method = p.adjust.method)
  results <- results[order(results$p.adjust), ]
  return(results)
}


#' @keywords internal
convert_greek_letters <- function(text) {
  # 希腊字母映射表（Unicode → 英文名称）
  greek_map <- c(
    # 小写
    "α" = "alpha", "β" = "beta", "γ" = "gamma", "δ" = "delta", "ε" = "epsilon",
    "ζ" = "zeta", "η" = "eta", "θ" = "theta", "ι" = "iota", "κ" = "kappa",
    "λ" = "lambda", "μ" = "mu", "ν" = "nu", "ξ" = "xi", "ο" = "omicron",
    "π" = "pi", "ρ" = "rho", "σ" = "sigma", "τ" = "tau", "υ" = "upsilon",
    "φ" = "phi", "χ" = "chi", "ψ" = "psi", "ω" = "omega",
    # 大写
    "Α" = "Alpha", "Β" = "Beta", "Γ" = "Gamma", "Δ" = "Delta", "Ε" = "Epsilon",
    "Ζ" = "Zeta", "Η" = "Eta", "Θ" = "Theta", "Ι" = "Iota", "Κ" = "Kappa",
    "Λ" = "Lambda", "Μ" = "Mu", "Ν" = "Nu", "Ξ" = "Xi", "Ο" = "Omicron",
    "Π" = "Pi", "Ρ" = "Rho", "Σ" = "Sigma", "Τ" = "Tau", "Υ" = "Upsilon",
    "Φ" = "Phi", "Χ" = "Chi", "Ψ" = "Psi", "Ω" = "Omega"
  )

  # 批量替换（高效循环）
  for (greek in names(greek_map)) {
    text <- gsub(greek, greek_map[[greek]], text, fixed = TRUE)
  }

  return(text)
}
