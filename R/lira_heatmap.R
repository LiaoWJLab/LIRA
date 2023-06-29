







#' Title
#'
#' @param eset
#' @param pdata
#' @param id_pdata
#' @param BOR
#' @param palette_heatmap
#' @param path
#' @param anonymous
#' @param remove_batch
#'
#' @return
#' @export
#'
#' @examples
lira_heatmap<-function(eset, pdata = NULL, id_pdata = "ID", BOR = "BOR", palette_heatmap = 6, path = NULL, anonymous = TRUE, remove_batch = TRUE, method = "tpm"){


  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("1-LIRA-heatmap")
  }

  if(!file.exists(file_store)) dir.create(file_store)
  abspath<-paste(getwd(),"/",file_store,"/",sep ="" )

  print(paste0(">>>--- Result will be stored in ", file_store))
  ###########################################
  if(!is.matrix(eset)) eset<-as.matrix(eset)
  ###########################################

  print(paste0(">>>>--", rownames(eset)[1]))
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

    eset_l  <- as.matrix(eset[, i])
    rownames(eset_l) <- rownames(eset)
    colnames(eset_l)<- colnames(eset)[i]


    eset_m<- merge(eset_l, eset_ref, by = "row.names", all = FALSE)

    # print(eset_m[1:5,1:5])
    eset_m <- remove_duplicate_genes(eset = eset_m, column_of_symbol = "Row.names")
    feas_gene <- feature_manipulation(data = eset_m, feature = rownames(eset_m), is_matrix = TRUE)
    eset_m<- eset_m[rownames(eset_m)%in%feas_gene, ]
    ##############################

    if(remove_batch){
      if(method == "count"){
        eset_m2 <- eset_m
        eset_m2[,paste0(pat_id, "_b")] <-eset_m2[, pat_id]
        pd<- data.frame("ID" = colnames(eset_m2), "batch" = ifelse(colnames(eset_m2)%in%c(pat_id, paste0(pat_id, "_b")), "validation", "reference"))
        eset_m2 <- sva::ComBat_seq(eset_m2, batch = pd$batch)
        eset_m <- eset_m2[,-which(colnames(eset_m2)==paste0(pat_id, "_b"))]
      }
    }

    #############################

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
    ##############################
    eset_m<- log2eset(eset_m)

    ###############################

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

      }
    }


    #将基因表达值矩阵作个转置，使行为样本，列为基因
    gene <- t(eset_m)
    #我们使用 FactoMineR 包中的方法，实现 PCA 分析和聚类添加
    # library(FactoMineR)
    #样本中基因表达值的 PCA 分析
    gene.pca <-FactoMineR:: PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)
    p<-plot(gene.pca)  #PCA 简图
    print(p)
    ggsave(p, filename = paste0( i, "-", pat_id, "-PCA.pdf"), width = 6, height = 6, path = file_store)
    ##############################
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


    if(anonymous){

      for (dd in 1:dim(input)[1]) {
        if(input$ID[dd] ==pat_id){
          input$ID[dd] <- pat_id
        }else{
          input$ID[dd] <- paste0("Reference_", dd)
        }

      }

    }

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
