







#' Title
#'
#' @param eset
#' @param pdata
#' @param id_pdata
#' @param BOR 
#' @param palette_heatmap 
#' @param path 
#'
#' @return
#' @export
#'
#' @examples
lira_heatmap<-function(eset, pdata = NULL, id_pdata = "ID", BOR = "BOR", palette_heatmap = 6, path = NULL){
  
  
  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("1-LIRA-heatmap")
  }
  
  if(!file.exists(file_store)) dir.create(file_store)
  abspath<-paste(getwd(),"/",file_store,"/",sep ="" )
  
  
  ###########################################
  if(!is.matrix(eset)) eset<-as.matrix(eset)
  ###########################################
  
  gene_type<- autoDetectGeneIdType(rownames(eset)[1])
  if(gene_type=="ensembl") {
    data(ref_count_op_ensembl)
    eset_ref<- ref_count_op_ensembl
  }
  if(gene_type=="symbol"){
    data(ref_count_op_symbol)
    eset_ref<- ref_count_op_symbol
  } 
  

  cat(crayon::bgYellow(">>>-- Drawing heatmap one by one...\n"))

  for (i in 1:length(colnames(eset))){
    pat_id <-  colnames(eset)[i]
    cat(crayon::underline(paste0(">>>---- Processing sample: ", colnames(eset)[i]), "\n"))
    message()
   if(colnames(eset)[i]%in%colnames(eset_ref)){
      colnames(eset_ref)[which(colnames(eset_ref)==colnames(eset)[i])] <- paste0(colnames(eset)[i], "_dup")
   }
   if(pat_id%in%ref_id){
     ref_id<- ref_id[!ref_id==pat_id]
     ref_pdata_heatmap <- ref_pdata_heatmap[ref_pdata_heatmap$ID%in%ref_id, ]
   }  
    
    
    eset_l<- as.matrix(eset[, i])
    rownames(eset_l) <- rownames(eset)
    colnames(eset_l)<- colnames(eset)[i]
    eset_m<- merge(eset_l, eset_ref, by = "row.names", all = FALSE)
    
    # print(eset_m[1:5,1:5])
    eset_m <- remove_duplicate_genes(eset = eset_m, column_of_symbol = "Row.names")
    eset_m <- count2tpm( countMat = eset_m, idType = gene_type)
    
    # id_cr<-  ref_lira_score[ref_lira_score$BOR=="CR"&ref_lira_score$PFS_months>20&ref_lira_score$score2=="low"&ref_lira_score$OS_status==0, ]$ID; id_cr
    # id_pr<-  ref_lira_score[ref_lira_score$BOR=="PR"&ref_lira_score$PFS_months>28&ref_lira_score$score2=="low"&ref_lira_score$OS_status==0, ]$ID; id_pr
    # id_sd<-  ref_lira_score[ref_lira_score$BOR=="SD"&ref_lira_score$PFS_months<2.5&ref_lira_score$score2=="high"&ref_lira_score$OS_status==1, ]$ID; id_sd
    # id_pd<-  ref_lira_score[ref_lira_score$BOR=="PD"&ref_lira_score$OS_months<3&ref_lira_score$score2=="high"&ref_lira_score$OS_status==1, ]$ID; id_pd
    # id_pd <- id_pd[5:10]; id_pd
    # ref_id <- c(id_cr, id_pr, id_sd, id_pd)
    # ref_pdata_heatmap <- ref_lira_score[ref_lira_score$ID%in%ref_id,]
    ##############################
    
    eset_m<- eset_m[, colnames(eset_m)%in%c(ref_id, pat_id) ]
    
    eset_m<- log2eset(eset_m)
    
    # print(eset_m[1:5, 1:5])
    ##############################
    
    
    if(is.null(pdata)){
      
      pdata2 <- data.frame(matrix(nrow = 1, ncol = ncol(ref_pdata_heatmap)))
      colnames(pdata2) <- colnames(ref_pdata_heatmap)
      pdata2[] <- NaN
      pdata2$ID <- pat_id
      pdata2$BOR <- "NE"
      
    }else{
      
      pdata1 <- data.frame(matrix(nrow = 1, ncol = ncol(ref_pdata_heatmap)))
      colnames(pdata1) <- colnames(ref_pdata_heatmap)
      pdata1[] <- NaN
      pdata1$ID <- pat_id
      colnames(pdata)[which(colnames(pdata)==id_pdata)] <- "ID"
      colnames(pdata)[which(colnames(pdata)==BOR)] <- "BOR"
      pdata1$BOR <- unique(pdata[pdata$ID==pat_id, BOR])
      pdata2 <- pdata1
    }
    
    pdata_heat <- rbind(pdata2, ref_pdata_heatmap)
    pdata_heat <- pdata_heat[pdata_heat$ID%in%colnames(eset_m),]
    
    # print(pdata_heat)
    ##############################
    
    
    # help("combine_pd_eset")
    input <- combine_pd_eset(eset = eset_m, 
                             pdata = pdata_heat,
                             feas = gsub(rf_feas, pattern = "_", replacement = "-"),
                             choose_who_when_duplicate = "eset")
    
    feas = gsub(rf_feas, pattern = "_", replacement = "-")
    input$BOR <- ifelse(input$BOR=="PR", 'dPR', input$BOR)
    input$BOR <- ifelse(input$BOR=="SD", 'eSD', input$BOR)
    
    sig_heatmap(input = input, 
                features = feas, 
                group = "BOR",
                scale = TRUE,
                palette = palette_heatmap, 
                width = 5.4, show_heatmap_col_name = TRUE, angle_col = 60,
                path = path, index = paste0(i, "-", pat_id) )
    
    save(input, file = paste0(abspath, i, "-", pat_id, "-heatmap-data.RData"))
    ############################
    message("")

  }
  
  return(input)
  
}
