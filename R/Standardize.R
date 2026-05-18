
#' Convert Gene ID Conveniently
#'
#' @param geneList A vector of gene ids or a matrix contains a column of gene ids named "Geneid".
#' @param species A string denotes the species. "Hs" for human and "Mm" for mouse.
#' @param from A string denotes the id type in the geneList.
#' @param to A string denotes the target id type.
#' @param single_mapping A logical value denotes whether one input gene id is mapped to only one target gene id, even though it has several alias. Default to be true.
#' @param to_rownames For matrix = TRUE only. A logical value denotes whether the converted ids are rownames of the resultant matrix. If rownames are duplicated, the resultant row takes the mean value.
#'
#' @returns A vector or matrix.
#' @export
#'
#' @examples
#' SYMBOL <- Quick_ID_conversion(ENTREZID, species = "Mm", input_datatype = "ENTREZID", target_datatype = "SYMBOL", matrix = F, single_mapping = T, to_rownames = F)
Quick_ID_conversion <- function(data, species, from, to, single_mapping = T, to_rownames = F, ...){

  args <- list(...)
  if ("input_datatype" %in% names(args) || "target_datatype" %in% names(args)) {
    lifecycle::deprecate_warn(...)
    if ("input_datatype" %in% names(args)) from <- args$input_datatype
    if ("target_datatype" %in% names(args)) to <- args$target_datatype
  }
    
    if (missing(from) || is.null(from) || missing(to) || is.null(to)) {
      stop("Argument 'from' & 'to' is required.")
    }
    
    df <- is.data.frame(data)
    if(df){
      if(!"Geneid" %in% colnames(data)) stop("The column named 'Geneid' is not found")
      ng <- Quick_ID_conversion(data$Geneid, species = species, from = from, to = to, single_mapping = single_mapping)
      data$Geneid <- ng
      if(to_rownames){
        data <- na.omit(data)
        if(anyDuplicated(ng) > 0){
          print("Duplication detected. Removing duplicated genes....")
          data <- remove_dup(data, index = "Geneid", method = "max")
        }
        row.names(data) <- data$Geneid
        data$Geneid <- NULL
      }
      return(data)
    }
    
    if(grepl("mm|mice|mouse|mus",stringr::str_to_lower(species))){
       if(single_mapping){
         gene_db_mmu <- geneSync::gene_db_mmu
         species = "mmu"
       }else{
         db = org.Mm.eg.db::org.Mm.eg.db
       } 
    }else if(grepl("hs|human|homo",stringr::str_to_lower(species))){
       if(single_mapping){
         gene_db_homo <- geneSync::gene_db_homo
         species = "homo"
       }else{
         db = org.Hs.eg.db::org.Hs.eg.db
       }
    }else{
      stop("This species is not supported")
    }

    from <- standardize_info(from, single_mapping)
    to <- standardize_info(to, single_mapping)
    if(from == to) stop("The from and to parameters are the same")
    
    
      if(single_mapping){
         data <- as.character(data)
          ta <- geneSync::gene_convert(data, species = species, from = from, to = to)
          if(to == "symbol") to <- "final_symbol"
          return(ta[match(data, ta$input), to])
      
      }else{
        conversion <- AnnotationDbi::select(db,keys = data,keytype = from,columns = to)
        agg <- stats::aggregate(
        conversion[, 2] ~ conversion[, 1],
        data = conversion,
        FUN = function(x) paste(unique(x), collapse = ";"))
        colnames(agg) <- c(from, to)
        res <- agg[match(data, agg[, 1]), 2]
        res[is.na(res)] <- data[is.na(res)]
        return(res)
      }
    

} 

#' @keywords internal
standardize_info <- function(arg, single_mapping){
  arg <- stringr::str_to_lower(arg)  
  if(grepl("ensem",arg)){
        return(ifelse(single_mapping, "ensembl_id", "ENSEMBL"))
    }else if(grepl("symb",arg)){
        return(ifelse(single_mapping, "symbol", "SYMBOL"))
    }else if(grepl("id",arg)){
        return(ifelse(single_mapping, "gene_id", "ENTREZID"))
    }else{
      stop("This gene id type is not supported")
    }
}


#' Remove duplicated gene in the expression matrix
#'
#' @param expr The expression matrix. 
#' @param index The column index of the gene. 
#' @param method The method to remove duplication. "max" takes the maximum for each value in duplicated row and combine them together. "mean" takes the mean. "kmax" keeps the row with the maximum sum (it will not change the raw data in the row).  
#'
#' @returns
#'
#' @export
#' @examples
remove_dup <- function(expr, index, method = "max") {
    if(!method %in% c("max", "mean", "kmax")){
      stop("Error: this method is not supported.")
    }
  
    if(is.numeric(index)){
      index <- colnames(expr)[index]
    } 
    
    g <- expr[, index, drop = TRUE]
    dpr <- g %in% unique(g[duplicated(g)])
    dpex <- expr[dpr,]
    if(nrow(dpex) == 0) return(expr)  
    expr <- expr[!dpr,]
    
    if(method == "kmax"){
      dpex <- dpex %>%
            dplyr::group_by(across(all_of(index))) %>%
            dplyr::slice_max(rowSums(across(where(is.numeric))), n = 1, with_ties = FALSE) %>%
            dplyr::ungroup()
    }else if(method == "max"){
       dpex <- dpex %>%
            dplyr::group_by(across(all_of(index))) %>%
            dplyr::summarise(across(everything(), ~if(is.numeric(.)) max(.) else first(.)), .groups = "drop")
    }else if(method == "mean"){
        dpex <- dpex %>%
            dplyr::group_by(across(all_of(index))) %>%
            dplyr::summarise(across(everything(), ~if(is.numeric(.)) mean(.) else first(.)), .groups = "drop")
    }
    
    dpex <- as.data.frame(dpex)
    rownames(dpex) <- NULL
    return(rbind(expr, dpex))
}