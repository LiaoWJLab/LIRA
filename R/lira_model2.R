

#' Title
#'
#' @param eset
#' @param pdata
#' @param id_pdata
#' @param scale
#' @param plot
#' @param check_eset
#' @param ref 
#' @param loop 
#' @param return_all 
#' @param org 'hsa' or 'mus'
#'
#' @return
#' @export
#'
#' @examples
lira_model2<-function(eset, org = "hsa", pdata = NULL, id_pdata = "ID", scale = TRUE, plot = FALSE, ref = TRUE, loop = TRUE, check_eset = TRUE, return_all = FALSE){
  
  if(!is.matrix(eset)) eset<-as.matrix(eset)
  ###########################################
  
  if(ref){
    
    gene_type<- autoDetectGeneIdType(rownames(eset)[1])
    if(gene_type=="ensembl") {
      data(ref_count_op_ensembl)
      eset_ref<- ref_count_op_ensembl
    }
    if(gene_type=="symbol"){
      data(ref_count_op_symbol)
      eset_ref<- ref_count_op_symbol
    } 
    
   if(loop){
     
     
     cat(crayon::bgYellow(">>>-- Predicting new data one by one...\n"))
     res<- data.frame(NULL)
     for (i in 1:length(colnames(eset))){
       
       cat(crayon::underline(paste0(">>>---- Processing sample: ", colnames(eset)[i]), "\n"))
       message()
       if(colnames(eset)[i]%in%colnames(eset_ref)){
         colnames(eset_ref)[which(colnames(eset_ref)==colnames(eset)[i])] <- paste0(colnames(eset)[i], "_dup")
       }
       eset_l<- as.matrix(eset[, i])
       rownames(eset_l) <- rownames(eset)
       colnames(eset_l)<- colnames(eset)[i]
       eset_m<-merge(eset_l, eset_ref, by = "row.names", all = FALSE)
       eset_m <- remove_duplicate_genes(eset = eset_m, column_of_symbol = "Row.names")
       
       eset_m<-count2tpm( countMat = eset_m, idType = gene_type, org = org, source = "local")
       res_l<- lira_model(eset_m, pdata = pdata, id_pdata = id_pdata, scale = scale, plot = FALSE, check_eset = TRUE)
       res_l<- res_l$score
       res<-rbind(res, res_l[1, ])
       message("")
     }
     if(return_all){
       res<- rbind(res, res_l)
       res<- res[!duplicated(res$ID), ]
       res$cohort <- ifelse(res$ID%in%colnames(eset), "validation", "reference")
     }
     
   }else{
     
     cat(crayon::bgYellow(">>>-- Predicting new data with combined gene expression data...\n"))
     # change the duplicated names
     for (jj in 1:length(colnames(eset))) {
       if(colnames(eset)[jj]%in%colnames(eset_ref)){
         colnames(eset_ref)[which(colnames(eset_ref)==colnames(eset)[jj])] <- paste0(colnames(eset)[jj], "_dup")
       }
     }
     eset_m<-merge(eset, eset_ref, by = "row.names", all = FALSE)
     eset_m <- remove_duplicate_genes(eset = eset_m, column_of_symbol = "Row.names")
     
     eset_m<-count2tpm( countMat = eset_m, idType = gene_type, org = org, source = "local")
     res<- lira_model(eset_m, pdata = pdata, id_pdata = id_pdata, scale = scale, plot = FALSE, check_eset = TRUE)
     res<- res$score
     if(!return_all){
       res<- res[res$ID%in%colnames(eset), ]
     }else{
       res$cohort <- ifelse(res$ID%in%colnames(eset), "validation", "reference")
     }
     res<- res[!duplicated(res$ID), ]
   }
    
  }else{
    
    res<- lira_model(eset, pdata = pdata, id_pdata = id_pdata, scale = scale, plot = FALSE, check_eset = TRUE)
    res<- es$score
    
  }
  
  if(scale) res$riskscore <- c(as.numeric(scale(res$riskscore, center = T)) + 2.5) * 10
  return(res)
}
