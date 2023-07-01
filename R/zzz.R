



.onLoad <- function(libname, pkgname) {

  invisible(suppressPackageStartupMessages(library("ggplot2")))
  invisible(suppressPackageStartupMessages(library("IOBR")))

  invisible(suppressPackageStartupMessages(
    sapply(c("crayon", "randomForest", "randomForestSRC", "survminer", "ggplot2"),
           requireNamespace, quietly = TRUE)
    ))
}



##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("==========================================================================\n",
                " ", pkgname, " v", pkgVersion, "  ",

                "  For help: https://github.com/LiaoWJLab/LIRA/issues", "\n\n")

  citation <-paste0(" If you use ", pkgname, " in published research, please cite:\n",
                    " DQ Zeng, YR Fang, G Chen, …, WJ Liao*.\n",
                    " LIRA: Predicting NSCLC immunotherapy resistance based on Random Forest Model, Under Review \n ",
                    # " Journal of Thoracic Oncology. 12:687975,(2023). \n",
                    # " XXXX, 2020", "\n",
                    # " DOI: 10.3389/JTO.2023.687975\n",
                    # " PMID:  ","\n",
                    "==========================================================================")

  packageStartupMessage(paste0(msg, citation))
}



