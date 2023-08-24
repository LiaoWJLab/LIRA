

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
#' @param remove_batch
#' @param method
#'
#' @return
#' @export
#'
#' @examples
lira_model2<-function(eset, org = "hsa", pdata = NULL, id_pdata = "ID", scale = TRUE, plot = FALSE, ref = TRUE, loop = TRUE,
                      check_eset = TRUE, return_all = FALSE, remove_batch = FALSE, method = "tpm"){

  if(!is.matrix(eset)) eset<-as.matrix(eset)
  ###########################################

  if(ref){

    gene_type<- autoDetectGeneIdType(rownames(eset)[1])
    if(gene_type=="ensembl") {
      data(ref_count_op_ensembl)
      eset_ref<- ref_count_op_ensembl
    }
    if(gene_type=="symbol"){
      data(ref_count_op_symbol, package = "LIRA")
      eset_ref<- ref_count_op_symbol
    }

   if(loop){


     cat(crayon::bgYellow(">>>-- Predicting new data one by one...\n"))
     res<- data.frame(NULL)
     for (i in 1:length(colnames(eset))){

       pat_id <- colnames(eset)[i]
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
       ####################################

       if(remove_batch){
         if(method == "count"){
           eset_m2 <- eset_m
           eset_m2[,paste0(pat_id, "_b")] <-eset_m2[, pat_id]
           pd<- data.frame("ID" = colnames(eset_m2), "batch" = ifelse(colnames(eset_m2)%in%c(pat_id, paste0(pat_id, "_b")), "validation", "reference"))
           eset_m2 <- sva::ComBat_seq(eset_m2, batch = pd$batch)
           eset_m <- eset_m2[,-which(colnames(eset_m2)==paste0(pat_id, "_b"))]
         }
       }

       ####################################
       eset_m<-count2tpm( countMat = eset_m, idType = gene_type, org = org, source = "local")
       ####################################
       if(remove_batch){

         if(method == "tpm"){
           eset_m2 <- eset_m
           eset_m2[,paste0(pat_id, "_b")] <-eset_m2[, pat_id]
           eset_m2[,paste0(pat_id, "_b")] <- eset_m2[,paste0(pat_id, "_b")] * 1.2
           ESET1<- as.matrix(eset_m2[,!colnames(eset_m2)%in%c(pat_id, paste0(pat_id, "_b"))])
           ESET2 <- as.matrix(eset_m2[, colnames(eset_m2)%in%c(pat_id, paste0(pat_id, "_b"))])

           eset_m2 <- remove_batcheffect(eset1 = ESET2,
                                         eset2 = ESET1,
                                         path =  NULL,
                                         save_plot = TRUE)

           eset_m <- eset_m2[,-which(colnames(eset_m2)==paste0(pat_id, "_b"))]
           # print(colnames(eset_m))
         }
       }
       #######################################


       res_l<- lira_model(eset_m, pdata = pdata, id_pdata = id_pdata, scale = scale, check_eset = TRUE)
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
     res<- lira_model(eset_m, pdata = pdata, id_pdata = id_pdata, scale = scale, check_eset = TRUE)
     res<- res$score
     if(!return_all){
       res<- res[res$ID%in%colnames(eset), ]
     }else{
       res$cohort <- ifelse(res$ID%in%colnames(eset), "validation", "reference")
     }
     res<- res[!duplicated(res$ID), ]
   }

  }else{

    res<- lira_model(eset, pdata = pdata, id_pdata = id_pdata, scale = scale, check_eset = TRUE)
    res<- res$score

  }

  if(scale) res$riskscore <- c(as.numeric(scale(res$riskscore, center = T)) + 2.5) * 10
  return(res)
}
