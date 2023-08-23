



#' This function uses the LIRA model to predict new data and return the results.
#'
#' @param eset The input expression matrix.
#' @param pdata additional patient data to include in the analysis (default is NULL).
#' @param id_pdata The identifier column name in pdata (default is "ID").
#' @param scale A logical value indicating whether to scale the data (default is FALSE).
#' @param check_eset A logical value indicating whether to remove outlier genes from the input expression matrix (default is FALSE).
#' @param from_rbatch A logical value indicating whether the sample was processed by the rbatch pipeline (default is TRUE).
#' @param model option = 1, or 2
#'
#' @return
#' @export
#'
#' @examples
lira_model<-function(eset, pdata = NULL, id_pdata = "ID", scale = FALSE, check_eset = FALSE, from_rbatch = TRUE, model = 1){

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

  if(model==1){
    data("rf_model1", package = "LIRA")
    sur_model <- rf_model1
    data("rf_feas1", package = "LIRA")
    rf_feas <- rf_feas1
  }else if(model ==2 ){
    data("rf_model2", package = "LIRA")
    sur_model <- rf_model2
    data("rf_feas2", package = "LIRA")
    rf_feas <- rf_feas2
  }

  message(">>>=== Feature of model: ")
  print(rf_feas)
  #############################################
  #############################################
  freq1<-length(intersect(rf_feas,rownames(eset)))/length(rf_feas)
  if(freq1<0.5){
    msg1<- paste0(paste0(sprintf(">>>-- Only %1.2f%%", 100*freq1)," of model genes appear on gene matrix,\n interpret results with caution \n"))
    cat(crayon::bgRed(msg1))
  }else if(freq1>=0.5){
    msg2<- paste0(paste0(sprintf(">>>-- %1.2f%%", 100*freq1)," of model genes appear on gene matrix\n"))
    cat(crayon::green(msg2))
  }
  ###########################################

  # https://github.com/r-lib/crayon
  cat(crayon::green(">>>-- Predicting new data with LIRA model...\n"))
  res<-      predict_rf_model_cox(sur_model    = sur_model,
                                  eset_new     = eset,
                                  pdata_new    = pdata,
                                  feas         = NULL,
                                  feas_genes   = rf_feas,
                                  prefix       = c("\\-","_"),
                                  id_pdata_new = id_pdata)
  ##########################################
  if(from_rbatch){
    cat(crayon::red(">>>=== This sample was testing by RNAseq and was processessed by rbatch pipline... \n"))
    cat(crayon::red(">>>=== The LIRA score will be normalised to a range of 1-10 ... \n"))
    res$score$LIRA <- (res$score$LIRA - 100)/15

  }

  cat(crayon::green(">>>-- DONE! \n"))

  return(res)

}
