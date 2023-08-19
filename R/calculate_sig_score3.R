






#' Title
#'
#' @param eset
#' @param signature_zscore
#' @param signature_pca
#' @param signature_ssgsea
#'
#' @return
#' @export
#'
#' @examples
#'
calculate_sig_score3 <- function(eset, signature_zscore = feas_sig_op_zscore, signature_pca = feas_sig_op_pca, signature_ssgsea = feas_sig_op_ssgsea){


  eset1 <- calculate_sig_score(eset = eset, signature = signature_zscore, method = "zscore", mini_gene_count = 3)

  eset2 <- calculate_sig_score(eset = eset, signature = signature_pca, method = "pca", mini_gene_count = 3)

  eset3 <- calculate_sig_score(eset = eset, signature = signature_ssgsea, method = "ssgsea", mini_gene_count = 3)

  input <- inner_join(eset1, eset2, by = "ID")
  input <- inner_join(input, eset3, by = "ID")

  return(input)

}
