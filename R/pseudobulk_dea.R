source("R/helper.R")
library(stringr)

#' Prepare pseudo bulk data
#'
#' Generates bulked data and corresponding metadata for subsequent pseudo-bulk
#' differential expression analysis (pseudo-bulk DEA). Pseudo-bulk DEA works by
#' summarizing gene expression within a cluster and a condition. Instead of
#' operating on expression of single cells, summarized expression values of a
#' bulk of cells are used. In contrast to bulk RNA-Seq the bulk RNA counts are
#' generated in-silico. Therefore it is called pseudo bulk as it is known how
#' the expression values are generated.
#'
#' @param seurat_object a Seurat object.
#' @param rep.suffix the function expects a standardized way of naming the
#' biological replicates used in each condition. It is expected that each
#' sample name ends with a string passed to the function via rep.suffix (usually
#' _R or _REP) followed by a unique number within each condition. Usually, the
#' replicated are numbered incrementally from one in each condition. Default
#' "_R".
#' @param slot The slot of Seurat object containing the count data. Default
#' "counts".
#' @param assay The assay from which to select the slot. Default "RNA".
#' @param cluster_slot The variable of metadata in Seurat containing the
#' clusters.
#'
#' @return A list containing two lists "data" and "metadata". The former
#' containing the the bulked counts data frames of each cluster. The second
#' containing the corresponding bulked metadata.
#'
#' @export
#'
prepare_pseudobulk <- function(seurat_object,
                               rep.suffix="_R",
                               slot="counts",
                               assay="RNA",
                               cluster_slot="seurat_clusters"){
  
  pattern=paste(rep.suffix,"[123456789]*",sep="")
  
  ## 1 ##
  print("Extract Metadata...")
  colData <-
    seurat_object@meta.data %>%
    tibble::add_column(Cell=rownames(.), .before=1) %>%
     dplyr::arrange(.$Cell) %>%
    tibble::add_column(condition=gsub(pattern = pattern,
                                      replacement = "",
                                      .$orig.ident))
  
  ## 2 ##
  print("Extract counts...")
  counts <- Seurat::GetAssayData(seurat_object, slot = slot, assay=assay) %>% as.matrix()

  counts <- t(counts ) 
  print(counts)
  ## 3 ##
  print("Merge metadata and counts and prepare table for bulking...")
  counts <-
    counts %>% as.data.frame() %>%
    tibble::add_column(Cell=rownames(.), .before = 1) %>%
     dplyr::filter(.$Cell %in% colData$Cell) %>%
     dplyr::arrange(.$Cell) %>%
    tibble::add_column(clusters= 
                        (colData %>% dplyr::pull(as.symbol(cluster_slot))), 
                       .before = 1) %>%
    tibble::add_column(sample=colData$orig.ident, .before = 1) %>%
     dplyr::arrange(.$clusters, .$sample)
  
  ## 4 ##
  print("Transform single-cell signals to pseudo bulk signals...")
  
  counts.sum <-
    counts  %>% 
    as.data.frame() %>%
    dplyr::select(-3) %>%
    group_by(clusters, sample) 
  
  counts.sum <- 
    counts.sum %>%
    dplyr::summarise(
      dplyr::across(.cols = 
                      as.name(colnames(.)[3]):as.name(colnames(.)[ncol(.)]), 
                    .fns = sum))
  counts.sum <- counts.sum %>% as.data.frame()
  
  ## 5 ##
  print("Create colData dataframes for all clusters...")
  colData <-
    counts.sum %>% dplyr::select(1,2) %>%
    tibble::add_column(condition=gsub(pattern = pattern,
                                      replacement = "",
                                      .$sample))
  colDatas <-
    colData %>%
    as.data.frame() %>%
    dplyr::group_by(.$clusters) %>%
    dplyr::group_split()
  
  ## 6 ##
  print("Create counts dataframes for all clusters")
  counts.sum.df <-
    counts.sum %>% dplyr::group_by(.$clusters) %>% dplyr::group_split()
  all.counts <- vector("list", length(counts.sum.df))
  cluster_names <- unique(counts.sum$clusters)
  
  for(i in 1:length(counts.sum.df)){
    print(paste("Cluster ",i,sep=""))
    all.counts[[i]] <-
      counts.sum.df[[i]] %>%
      tibble::column_to_rownames("sample") %>%
      dplyr::select(-1) %>% t() %>%
      as.data.frame()
    
    rownames <- all.counts[[i]] %>% rownames()
    all.counts[[i]] <- all.counts[[i]] %>% sapply(as.numeric) %>% as.data.frame()
    rownames(all.counts[[i]]) <- rownames
  }
  names(all.counts) <- cluster_names
  
  ## 7 ##
  print("Pseudo bulk preparation finished. Results returned.")
  return(list(metadata=colDatas, data=all.counts))
}

#' DESeq2 wrapper
#'
#' A wrapper of the DESeq workflow with optional correction of the p-value
#' distribution and related adjusted p-values.
#'
#' @param data the count data.
#' @param colData the metadata describing the experiment.
#' @param as.DT whether to return a DT version of the adjusted result table.
#' @param add.fdrtool whether to add fdrtool-correction to the result table.
#' @param used.statistic statistic of fdrtool used to adjust p-values. supported
#' values are "normal" and "pvalue". The DEA-table shoudl contain a column
#' 'stat' or 'pvalue' respectively.
#' @param design a string variation of a design formula. Default: ~condition.
#' @param plot.fdrtool whether to plot the plots generated by fdrtool. Default:
#' F. Only works if add.fdrtool ist set to be true.
#' @importFrom dplyr %>%
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 DESeq
#' @importFrom stats na.omit
#' @importFrom dplyr arrange
#' @importFrom tibble add_column
#' @importFrom fdrtool fdrtool
#'
#'
#' @return  a list containing a DESeq2-object, after running DESeq analysis,
#' the results (optionally extended with an adjusted version of the results
#' containing corrected p-values using the r-package fdrtool), and (optionally)
#' a DT version of the latter results table.
#'
#' @export
#'
run_deseq <- function(data,
                      colData,
                      as.DT=F, 
                      add.fdrtool=T,
                      used.statistic="stat",
                      design="condition",
                      reference="WT",
                      treatment="CS",
                      factor="condition",
                      plot.fdrtool=F
                      ){
  
  d <- stats::as.formula(paste("~",design,sep=""))
  ### Generate DESeq object, based on data, metadata and design.
  DDS <- DESeq2::DESeqDataSetFromMatrix(countData=data,
                                        colData = colData,
                                        design=d)
  
  ### Get all used factors in the formula.
  factors <- 
    str_split(design, "\\+")[[1]] %>% 
    as.data.frame() %>% dplyr::rename(EV=colnames(.)[1])  %>% 
    filter(!grepl(pattern = "\\*", .$EV)) %>% unlist() %>% unname()
  
  ### Set the reference values for each factor.
  for(i in 1:length(factors)){
    DDS[[factors[i]]] <- relevel(DDS[[factors[i]]], reference[i])
  }
  
  ### Run DESeq analysis
  DDS <- DESeq(DDS)
  
  print(colData(DDS))
  
  ### Get analysis results.
  results <-  DESeq2::results(DDS, contrast = c("condition",treatment,reference)) %>%
              as.data.frame() %>%
              filter(!is.na(pvalue))%>%
              tibble::add_column(symbol=rownames(.), .before=1)
  
  ### Optionally perform correction of p-value.
  if(add.fdrtool){
    if(used.statistic=="stat"){
      corrected <-
        fdrtool::fdrtool(results$stat, statistic = "normal", plot=plot.fdrtool)
    } else if(used.statistic=="pvalue"){
      corrected <-
        fdrtool::fdrtool(results$pvalue, statistic="pvalue", plot=plot.fdrtool)
    }
    results <-
      results %>%
      tibble::add_column(padj.fdrtool=corrected$qval, 
                         pvalue.fdrtool=corrected$pval) %>%
      stats::na.omit() %>%
      dplyr::arrange(.$padj)
  }
  
  ### Optionally generate a DT table version of the results.
  if(as.DT){
    return(list(analysis=DDS,
                results=results,
                results.DT=(results %>% dt_table())))
  } else{
    return(list(analysis=DDS, results=results))
  }
}

#' Runs pseudo bulk RNA-Seq analysis
#'
#' Pseudo-bulk DEA works by summarizing gene expression within each cluster and
#' condition. Instead of operating on expression of single cells, summarized
#' expression values of a bulk of cells are used. In contrast to bulk RNA-Seq
#' the bulk RNA counts are generated in-silico. Therefore it is called pseudo
#' bulk as it is known how the expression values are generated. After generating
#' the pseudo bulk data, a standard differential expression analysis of RNA-Seq
#' data is performed. For now the function supports the standard DESeq2 workflow.
#'
#' @param seurat_object a Seurat object.
#' @param rep.suffix the function expects a standardized way of naming the
#' biological replicates used in each condition. It is expected that each
#' sample name ends with a string passed to the function via rep.suffix (usually
#' _R or _REP) followed by a unique number within each condition. Usually, the
#' replicated are numbered incrementally from one in each condition. Default
#' "_R".
#' @param slot The slot of Seurat object containing the count data. Default
#' "counts".
#' @param assay The assay from which to select the slot. Default "RNA".
#' @param cluster_slot The variable of metadata in Seurat containing the
#' clusters.
#' @param as.DT whether to return a DT version of the adjusted result table.
#' @param add.fdrtool whether to add fdrtool-correction to the result table.
#' @param used.statistic statistic of fdrtool used to adjust p-values. supported
#' values are "normal" and "pvalue". The DEA-table should contain a column
#' 'stat' or 'pvalue' respectively.
#' @param design a string variation of a design formula. Default: ~condition.
#'
#' @return a list containing pseudo bulk analysis results for all provided
#' clusters.
#' @export
#'
run_pseudobulk <- function(seurat_object,
                           rep.suffix="_R",
                           slot="counts",
                           assay="RNA",
                           cluster_slot="seurat_clusters",
                           as.DT=F,
                           add.fdrtool=T,
                           used.statistic="stat",
                           treatment="KO",
                           reference="WT",
                           design="condition"){
  prepared_pseudobulk <-  prepare_pseudobulk(seurat_object=seurat_object,
                                             rep.suffix = rep.suffix,
                                             slot = slot,
                                             assay = assay,
                                             cluster_slot=cluster_slot)
  
  print("finished")
  saveRDS(prepared_pseudobulk, file = "tmp/tmp.prepared.pseudobulk.rds")
  print("saved")
  print("getting clusters")
  clusters <-
    seurat_object@meta.data[, cluster_slot] %>%
    unique() %>% sort()
  print(clusters)
  saveRDS(clusters, file = "tmp/tmp.prepared.pseudobulk.rds")
  

  
  all.results <- vector("list", length(clusters))
  for(i in 1:length(clusters)){
    all.results[[i]] <-
      run_deseq(data=prepared_pseudobulk[["data"]][[i]] %>% na.omit(),
                colData=prepared_pseudobulk[["metadata"]][[i]],
                as.DT=as.DT, treatment = treatment, reference = reference,
                add.fdrtool=add.fdrtool,
                used.statistic=used.statistic,
                design=design
      )
  }
  names(all.results) <- clusters
  return(all.results)
}

#' This functions combines counts matrix and metadata of RNA-Seq experiment in 
#' one table containing all genes and samples row-wise, the expression values 
#' are also stored in one column along with corresponding metadata. This 
#' function provides a long format of the rnaseq counts, and is mainly used for 
#' as input for visualization with ggplot.
#'
#' @param counts the raw counts from RNA-Seq experiment.
#' @param metadata the corresponding metadata.
#' @param sample_name_col the name of the column containing the samples names.
#' @param colnames_replicate_infix the replicate pattern in the sample name.
#' @param db.key the key of the database used to fetch gene names.
#' @param species the name of the species as used in the database (e.g. 
#' mmusculus, hsapiens used when accessing ensembl with biomaRt etc.)
#' @param ensembl.version the version ensembl db to use (see useEnsembl from
#' biomaRt version `r eval(package.version("biomarRt))`
#' @param ensembl.mirror see function useEnsembl from biomaRt package.
#' @param ensembl.host see function useEnsembl from biomaRt package.
#' @param ensembl.gene.id  see function getBM from biomaRt package.
#' @param ensembl.gene.name  see function getBM from biomaRt package.
#'
#' @return a table containing all genes and their values for all samples as well
#' as the provided metadata in long format.
#' @export
#'
combine_data_and_metadata <- function(counts, 
                                      metadata, 
                                      sample_name_col = "sample", 
                                      colnames_replicate_infix = "_R",
                                      species = "mmusculus",
                                      db.key = "ensembl", 
                                      ensembl.version = NULL,
                                      ensembl.mirror = NULL,
                                      ensembl.host = NULL,
                                      ensembl.gene.id  = "enseml_gene_id",
                                      ensembl.gene.name ="external_gene_name"){
  
  counts.long <- 
    counts %>% 
    add_column(symbol=rownames(.)) %>% 
    tidyr::pivot_longer(cols=contains(colnames_replicate_infix), 
                        names_to="sample", 
                        values_to="counts")
  
  if(db.key=="ensembl"){
    ensmart <- 
      useEnsembl(biomart = "genes", 
                 dataset = paste(species,"_gene_ensembl",sep=""), 
                 version = ensembl.version,
                 mirror  = mirror,
                 host    = host)
    table <- 
      getBM(attributes=c(ensembl.gene.id, ensembl.gene.name), 
            mart=ensmart, 
            filters = c(ensembl.gene.id), 
            values=counts.long$ensembl_gene_id)
    colnames(table) <- c(ensembl.gene.id, "symbol")
    
    counts.long <- 
      counts.long %>% 
      merge(table, by=ensembl.gene.id) %>% 
      relocate(symbol, .before=ensembl.gene.id)
  }
  
  counts.long <- counts.long %>% merge(metadata, by=sample_name_col)
  return(counts.long)
}