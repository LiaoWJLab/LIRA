



#' Title
#'
#' @param eset processed gene expression data is reconmend
#' @param pdata
#' @param id_pdata
#' @param scale
#' @param plot
#' @param check_eset
#'
#' @return
#' @export
#'
#' @examples
lira_model<-function(eset, pdata = NULL, id_pdata = "ID", scale = FALSE, plot = FALSE, check_eset = FALSE){

  if(!is.matrix(eset)) eset<-as.matrix(eset)
  ###########################################

  ############################################
  if(scale){
    cat(crayon::green(">>>-- Scaling data...\n"))
    eset<-t(scale(t(eset)))
  }
  #############################################
  if(check_eset){
    cat(crayon::green(">>>-- Removing outlier genes...\n"))
    genes<- rownames(eset)
    genes<- IOBR::feature_manipulation(data = eset, feature = genes, is_matrix = TRUE, print_result = T)
    eset<-eset[rownames(eset)%in%genes, ]
  }
  #############################################
  data("rf_feas")
  #############################################
  # freq1<-length(intersect(rf_feas,rownames(eset)))/length(rf_feas)
  # if(freq1<0.5){
  #   msg1<- paste0(paste0(sprintf(">>>-- Only %1.2f%%", 100*freq1)," of model genes appear on gene matrix,\n interpret results with caution \n"))
  #   cat(crayon::bgRed(msg1))
  # }else if(freq1>=0.5){
  #   msg2<- paste0(paste0(sprintf(">>>-- %1.2f%%", 100*freq1)," of model genes appear on gene matrix\n"))
  #   cat(crayon::green(msg2))
  # }
  ###########################################
  data("rf_model")
  ############################################
  # https://github.com/r-lib/crayon
  cat(crayon::green(">>>-- Predicting new data with LIRA model...\n"))
  res<-      predict_rf_model_cox(sur_model    = rf_model,
                                  eset_new     = eset,
                                  pdata_new    = pdata,
                                  feas         = NULL,
                                  feas_genes   = rf_feas,
                                  prefix       = c("\\-","_"),
                                  id_pdata_new = id_pdata)
  cat(crayon::green(">>>-- DONE! \n"))

  return(res)

}
