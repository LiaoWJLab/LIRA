


#' Title
#'
#' @param id 
#'
#' @return
#' @export
#'
#' @examples
autoDetectGeneIdType <- function(id){
  type <- NA
  if(grepl("^[Ee][Nn][Ss][A-Za-z]{0,3}[Gg][0-9]+", id)) type <- "ensembl"
  else if(grepl("^[0-9]+$", id)) type <- "entrez"
  else if(grepl("^[Yy][A-Za-z]{2}[0-9]{3}[A-Za-z]", id)) type <- "sgd"
  else if(grepl("^[Aa][Tt][0-9][A-Za-z][0-9]{5}", id)) type <- "tair"
  else type<- "symbol"
  return(type)
}




#' Title
#'
#' @param pkg 
#' @param type 
#'
#' @return
#' @export
#'
#' @examples
isAvailable <- function(pkg, type = ""){
  
  if(!(pkg %in% .packages(all.available=TRUE))) {
    message(paste0("Corresponding ", type,  " package not found: ",
                   pkg, "\nMake sure that you have it installed."))
    choice <- readline("Install it now? (y/n): ")
    if(choice == "y")
    {
      if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    }
    else stop(paste("Package", pkg, "is not available"))
  }
  require(pkg, character.only = TRUE)
}