




#' Title
#'
#' @param sur_model
#' @param eset_new
#' @param feas
#' @param pdata_new
#' @param id_pdata_new
#' @param prefix
#' @param na_replace
#'
#' @return
#' @export
#'
#' @examples
predict_rf_model_cox<- function(sur_model,
                                eset_new,
                                pdata_new    = NULL,
                                feas         = NULL,
                                prefix       = c("-", "\\:"),
                                id_pdata_new = "ID",
                                na_replace   = 1){

  #######################################
  # library(tidyverse)
  library(randomForestSRC)
  ############################################
  #get features of model
  ############################################
  if(is.null(feas)){

    coef<-sur_model$importance
    feas<-c(colnames(sur_model$yvar), names(coef))
    coef<- coef %>% {
      data.frame(gene.name = names(.),
                 importance = .,
                 stringsAsFactors = FALSE)
    } %>%
      arrange(gene.name)
  }else{
    if(model_type=="rf"){
      coef<-sur_model$importance
      coef<- coef %>% {
        data.frame(gene.name = names(.),
                   importance = .,
                   stringsAsFactors = FALSE)
      } %>%
        arrange(gene.name)
    }
  }
  #################################################
  #################################################
  #make input data
  # eset_ref<-scale(t(eset_ref))
  # if(!is.null(prefix)){
  #   colnames(eset_ref)<- gsub(colnames(eset_ref), pattern = prefix[1], replacement = prefix[2])
  # }
  # eset_ref<-eset_ref[, colnames(eset_ref)%in%feas]

  eset_ref<- as.data.frame(matrix(0, nrow = 10, ncol = length(feas)))
  colnames(eset_ref)<-feas

  # normalization
  ##################################################
  eset_new<-scale(t(eset_new))

  if(!is.null(prefix)){
    for (dd in 1:length(prefix)) {
      colnames(eset_new)<- gsub(colnames(eset_new), pattern = prefix[dd], replacement = "_")
    }
  }
  ##################################################
  # message("For rf model, 'feas' must containe 'time', 'status' and features: c('time','status', genes_train) ")
  eset_new<-eset_new[ ,colnames(eset_new)%in% feas]
  eset_new<- IOBR:: assimilate_data(eset_ref, eset_new)
  # print(eset_new[1:5, 1:5])
  #################################################
  riskscore      <- predict(sur_model,
                            newdata    = eset_new,
                            na.action  = "na.impute",
                            importance = TRUE,
                            type = "risk")
  riskscore<-riskscore$predicted
  riskscore<-data.frame("ID" = rownames(eset_new), "riskscore" = riskscore)

  if(!is.null(pdata_new)){
    colnames(pdata_new)[which(colnames(pdata_new)==id_pdata_new)]<-"ID"
    riskscore<-merge(riskscore, pdata_new, by = "ID", all.x = TRUE, all.y = FALSE)
  }

  res<-list("score" = riskscore, "importance" = coef)
  return(res)
}
