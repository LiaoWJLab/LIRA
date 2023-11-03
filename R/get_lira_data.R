





#' Title
#'
#' @param id LIRA0001
#' @param data_type "tpm", "scale_tpm", "count"
#' @param path  "E:/03-NSCLC/19-NSCLC-LIRA/4-analysis/0-data/"
#'
#' @return
#' @export
#'
#' @examples
get_lira_data <- function(id = "LIRA0001", data_type = c("tpm", "scale_tpm", "count"), path = "E:/03-NSCLC/19-NSCLC-LIRA/4-analysis/0-data/"){

  if(data_type=="tpm"){
    (load(paste0(path, "Result-of-",id,"/4-2-eset-corrected-tpm-of-",id,".RData" )))
    eset <- eset_tpm_corrected
  }else if(data_type=="scale_tpm"){
    (load(paste0(path, "Result-of-",id,"/4-3-eset-corrected-tpm-scale-of-",id,".RData" )))
    eset <- eset2_tpm_scale2
  }else if(data_type=="count"){
    (load(paste0(path, "Result-of-",id,"/4-1-eset-corrected-count-of-",id,".RData" )))
    eset <- eset_count_corrected
  }
  ############################################
  # data("rf_feas", package = "LIRA")
  # feas <- gsub(rf_feas, pattern = "\\_", replacement = "-")
  # detect_miss(eset = eset, signature = feas)

  message(">>>=== Dimension of eset: ")
  print(dim(eset))
  ############################################

  return(eset)
}




#' Title Detect missing genes of signature gene list
#'
#' @param eset expression set
#' @param signature a gene list
#'
#' @return
#' @export
#'
#' @examples
detect_miss <- function(eset, signature){

  if(class(signature)=="list"){

    for (i in 1:length(signature)) {

      sigs <- signature[i]
      print(paste0(">>>== ", names(signature)[i]))

      # print("")
      freq1<-length(intersect(sigs, rownames(eset)))/length(sigs)
      if(freq1<0.5){
        msg1<- paste0(paste0(sprintf(">>>== Only %1.2f%%", 100*freq1)," of signature genes appear on gene matrix, interpret results with caution"))
        print(paste0(msg1))
        message("Missing genes: ")
        message(paste0(sigs[!sigs%in%rownames(eset)], collapse  = ", "))

      }else if(freq1>=0.5){
        message(paste0(paste0(sprintf(">>>==  %1.2f%%", 100*freq1)," of signature genes appear on gene matrix")))
        message("Missing genes: ")
        message(paste0(sigs[!sigs%in%rownames(eset)], collapse  = ", "))
      }
    }

  }else{
    sigs <- signature
    freq1<-length(intersect(sigs, rownames(eset)))/length(sigs)
    if(freq1<0.5){
      msg1<- paste0(paste0(sprintf(">>>== Only %1.2f%%", 100*freq1)," of signature genes appear on gene matrix, interpret results with caution"))
      print(paste0(msg1))
      message("Missing genes: ")
      message(paste0(sigs[!sigs%in%rownames(eset)], collapse  = ", "))
    }else if(freq1>=0.5){
      message(paste0(paste0(sprintf(">>>==  %1.2f%%", 100*freq1)," of signature genes appear on gene matrix")))
      message("Missing genes: ")
      message(paste0(sigs[!sigs%in%rownames(eset)], collapse  = ", "))
    }
  }


}





#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
get_pdata_op <- function(path = "E:/03-NSCLC/17-NSCLC-OP/3-data-analysis/1-model/7-final-model/best-model-type3-34639-pfs-0.001-ICN-good/4-riskscore-pdata.RData"){

  (load(path))
  colnames(pdata_score)[2] <- "lira_ref"
  pdata_score$BOR <- ifelse(pdata_score$BOR=="", "NE", pdata_score$BOR)
  pdata_score$ref_heatmap <- ifelse(pdata_score$ID%in%ref_pdata_heatmap$ID, "ref_heatmap", pdata_score$ID)
  print(head(pdata_score))
  print(table(pdata_score$treat_score, pdata_score$BOR))
  return(pdata_score)
}
