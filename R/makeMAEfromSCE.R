makeMAEfromSCE <- function(sce_list, 
                           name_assays,
                           sampleID_name = "cell_id"){
  
  for(i in 1:length(sce_list)){
    colnames(colData(sce_list[[i]]))[which(colnames(colData(sce_list[[i]])) == sampleID_name)] <- "cell_id"
  }
  
  #sce_list must have a sampleID/cellID 
  se_list <- lapply(sce_list, function(x){
    SummarizedExperiment(
      assays = list(exprs = assay(x)),
      colData = colData(x),
      rowData = rowData(x)
    )
  })
  
  names(se_list) <- name_assays
  
  map_list <- lapply(se_list, function(y){
    data.frame(primary = as.character(y$cell_id),
               colname = as.character(y$cell_id), 
               stringsAsFactors = FALSE)
  })
  
  names(map_list) <- names(se_list)
  
  df_map <- listToMap(map_list)
  
  sam_list <- lapply(se_list, function(u){
    colData(u) %>% data.frame() %>% as.tibble()
  })
  
  sam_full <- sam_list %>% purrr::reduce(full_join, 
                                         by = "cell_id") %>% as.data.frame
  

  rownames(sam_full) <- sam_full$cell_id 
  
  mae <- MultiAssayExperiment(experiments = se_list,
                              colData = sam_full, 
                              sampleMap = df_map)
  return(mae)
  
}