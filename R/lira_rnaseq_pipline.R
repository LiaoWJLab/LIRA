




#' Title
#'
#' @param pat_id
#' @param eset
#' @param data_type
#' @param palette
#' @param path1
#' @param path2
#' @param model
#' @param data_path
#'
#' @return
#' @export
#'
#' @examples
lira_rnaseq_pipline <- function(pat_id, eset = NULL, data_type = "tpm", model = 1, palette = "nrc", path1 = NULL, path2 = NULL, data_path = "E:/03-NSCLC/19-NSCLC-LIRA/4-analysis/0-data/"){


  if(is.null(path1)){
    if(is.null(path2)){
      file_name <- creat_folder(pat_id)
    }else{
      file_name <- creat_folder(path2)
    }
  }else{
    if(is.null(path2)){
      file_name <- creat_folder(path1, pat_id)
    }else{
      file_name <- creat_folder(path1, path2)
    }
  }
  ##################################

  pdata_op <- get_pdata_op()
  # table(pdata_op$BOR)
  ###################################
  if(is.null(eset)){
    eset <- get_lira_data(id = pat_id, data_type = data_type, path = data_path)  #也可选择 tpm_scale
    print(eset[1:5, 1:5])
  }

  ####################################
  res1 <- lira_model(eset        = eset,
                     check_eset  = TRUE,
                     scale       = TRUE,
                     pdata       = NULL,
                     id_pdata    = "ID",
                     model       = model,
                     from_rbatch = TRUE)
  score <- res1$score
  head(score)
  ####################################
  # (load("E:/03-NSCLC/17-NSCLC-OP/3-data-analysis/1-model/7-final-model/best-model-type3-34639-pfs-0.001-ICN-good/4-riskscore-pdata.RData"))
  # head(pdata_op)
  table(pdata_op$BOR)
  pdata_op$BOR <- ifelse(pdata_op$BOR=="", "NE", pdata_op$BOR)
  pdata_op$BOR <- ifelse(pdata_op$BOR=="CR", "CRPR", pdata_op$BOR)
  pdata_op$BOR <- ifelse(pdata_op$BOR=="PR", "CRPR", pdata_op$BOR)
  # colnames(pdata_op)[2] <- "LIRA_OP"
  ###################################
  score <- merge(score, pdata_op, by = "ID", all.x = T, all.y = TRUE)
  score <- score[!is.na(score$LIRA), ]

  #' 需要设置新数据的关键变量，否则后续画图会报错
  ####################################
  score[score$ID==pat_id, "BOR"] <- "NE"
  score[score$ID==pat_id, "ARM"] <- "IO"
  ####################################
  range(score$LIRA)
  ####################################
  input <- score[score$ARM=="IO", ]

  #####################################
  p2<-ggplot(input, aes(x= LIRA, fill= BOR)) +
    geom_histogram(bins = 30, colour = "grey", alpha = 0.5)+
    scale_fill_manual(values= palettes(palette = "jama"))+
    # geom_density(alpha=.2, fill="grey", weight = 1)+
    xlab("LIRA score (combine, TPM-scale)")+
    theme_light()+design_mytheme(legend.position = "bottom", axis_angle = 0, plot_title_size = 2)
  p2
  #####################################
  p1<- sig_box(data = input, signature = "LIRA", variable = "BOR")
  p<- p1|p2
  ggsave(p, filename = paste0("0-score-and-BOR-model-",model,".pdf"), width = 12, height = 6.3, path = file_name$folder_name)
  #####################################


  #' 寻找最佳的cutoff值
  bc1 <- best_cutoff2(pdata = input, variable = "LIRA", time = "PFS_months", status = "PFS_status" )
  bc2 <- best_cutoff2(pdata = input, variable = "LIRA", time = "OS_months", status = "OS_status" )

  sig_surv_plot(input_pdata =  input,
                signature   = "LIRA",
                time        = "OS_months",
                status      = "OS_status",
                save_path   = file_name$folder_name,
                mini_sig = "LIRA",
                index       = paste0("1-",pat_id,"-model-", model),
                project     = paste0("OS-model-",model))
  #####################################

  sig_surv_plot(input_pdata =  input,
                signature   = "LIRA",
                time        = "PFS_months",
                status      = "PFS_status",
                save_path   = file_name$folder_name,
                mini_sig = "LIRA",
                index       = paste0("2-",pat_id,"-model-", model),
                project     = paste0("PFS-model-",model))
  #####################################


  range(score$LIRA)
  head(score)
  score[score$ID==pat_id, ]
  #####################################

  lira_score_location(score        = input,
                      pat_id       = pat_id,
                      col_score    = "LIRA",
                      ref_score    = NULL,
                      best_cutoff  = round(bc1$best_cutoff,2),
                      palette      = palette,
                      panel        = "PFS",
                      cols         = NULL,
                      index        = paste0("3-PFS-model", model),
                      palette_line = "jama",
                      path         = file_name$folder_name)
  #################################

  lira_score_location(score        =  input,
                      pat_id       = pat_id,
                      col_score    = "LIRA",
                      ref_score    = NULL,
                      best_cutoff  = round(bc2$best_cutoff,2),
                      palette      = palette,
                      panel        = "OS",
                      cols         = NULL,
                      palette_line = "jama",
                      index        = paste0("3-OS-model", model),
                      path         = file_name$folder_name)

  lira_score_location(score        =  input,
                      pat_id       = pat_id,
                      col_score    = "LIRA",
                      ref_score    = NULL,
                      best_cutoff  = round(mean(input$LIRA),2),
                      palette      = palette,
                      panel        = "OS",
                      cols         = NULL,
                      palette_line = "jama",
                      index        = paste0("3-mean-model", model),
                      path         = file_name$folder_name)
  #################################

  lira_input <- input
  #################################
  index <- c(pat_id, ref_pdata_heatmap$ID);index

  #################################

  if(model==1){
    data("rf_feas_condiction1", package = "LIRA")
    rf_feas_condiction <- rf_feas_condiction1
    colnames(rf_feas_condiction)[1] <- "vars"
    data("rf_feas1", package = "LIRA")
    feas <- rf_feas1
  }else if(model==2){
    data("rf_feas_condiction2", package = "LIRA")
    rf_feas_condiction <- rf_feas_condiction2
    data("rf_feas2", package = "LIRA")
    feas <- rf_feas2
  }

  ################################
  feas <- gsub(feas, pattern = "_", replacement = "-")
  input <-IOBR:: combine_pd_eset(eset[,colnames(eset)%in%index], pdata = score[score$ARM=="IO", ], feas = feas)

  for (dd in 1:dim(input)[1]) {
    if(input$ID[dd] ==pat_id){
      input$ID[dd] <- pat_id
    }else{
      input$ID[dd] <- paste0("Reference_", dd)
    }
  }
  ######################################################

  sig_heatmap(input                 = input,
              features              = feas,
              condiction            = rf_feas_condiction,
              id_condiction         = "vars",
              col_condiction        = "condiction",
              group                 = "BOR",
              scale                 = TRUE,
              palette               = 6,
              width                 = 5,
              show_heatmap_col_name = TRUE,
              angle_col             = 60,
              path                  = file_name$folder_name,
              index                 = paste0("4-",pat_id, "-response-model", model),
              show_plot             = TRUE)
  ######################################################

  sig_heatmap(input                 = input,
              features              = feas,
              condiction            = rf_feas_condiction,
              id_condiction         = "vars",
              col_condiction        = "condiction2",
              group                 = "BOR",
              scale                 = TRUE,
              palette               = 6,
              width                 = 5,
              show_heatmap_col_name = TRUE,
              angle_col             = 60,
              path                  = file_name$folder_name,
              index                 =  paste0("4-",pat_id, "-os-model", model),
              show_plot             = TRUE)
  ######################################################



  sig_score <- calculate_sig_score3(eset = eset)
  head(sig_score)
  sig_score[,2:ncol(sig_score)] <- scale(sig_score[,2:ncol(sig_score)])
  input <- merge(score, sig_score, by = "ID")
  input <- input[input$ARM=="IO", ]
  print(colnames(input))
  #####################################################
  ######################################################

  if("Immune_Checkpoint"%in%colnames(input)){
    target <- 'Immune_Checkpoint'
    p<- sig_box(data = input, signature = target, variable = "BOR")

    ggsave(p, filename = paste0("5-",target,"-and-BOR.pdf"), width = 7, height = 6.3, path = file_name$folder_name)
    ######################################################
    bc1 <- best_cutoff2(pdata = input, variable = "Immune_Checkpoint", time = "PFS_months", status = "PFS_status" )
    lira_score_location(score        = input,
                        pat_id       = pat_id,
                        col_score    = "Immune_Checkpoint",
                        title        = "ICBscore",
                        ref_score    = NULL,

                        best_cutoff  = round(bc1$best_cutoff,2),
                        palette      = palette,
                        panel        = "PFS",
                        cols         = NULL,
                        index        = paste0("5-",pat_id, "-PFS-ICB-"),
                        palette_line = "jama",
                        path         = file_name$folder_name)
  }
  ########################################################

  if("IFNG_signature_Ayers_et_al"%in%colnames(input)){
    target <- 'IFNG_signature_Ayers_et_al'
    p<- sig_box(data = input, signature = target, variable = "BOR")

    ggsave(p, filename = paste0("5-",target,"-and-BOR.pdf"), width = 7, height = 6.3, path = file_name$folder_name)
    ######################################################
    bc1 <- best_cutoff2(pdata = input, variable = "IFNG_signature_Ayers_et_al", time = "PFS_months", status = "PFS_status" )
    lira_score_location(score        = input,
                        pat_id       = pat_id,
                        col_score    = "IFNG_signature_Ayers_et_al",
                        title        = "GEP",
                        ref_score    = NULL,

                        best_cutoff  = round(bc1$best_cutoff,2),
                        palette      = palette,
                        panel        = "PFS",
                        cols         = NULL,
                        index        = paste0("5-",pat_id, "-PFS-GEP-"),
                        palette_line = "jama",
                        path         = file_name$folder_name)
  }
  ########################################################


  input <- input[input$ID%in%index, ]
  ######################################################
  sig_heatmap(input                 = input,
              features              = names(feas_sig_op),
              condiction            = rf_sig_condiction,
              col_condiction        = "condiction",
              id_condiction         = "vars",
              group                 = "BOR",
              scale                 = FALSE,
              palette               = 6,
              width                 = 7,
              show_heatmap_col_name = TRUE,
              angle_col             = 60,
              path                  = file_name$folder_name,
              index                 = paste0("6-", pat_id, "-signature-response"),
              show_plot             = TRUE)

  sig_heatmap(input                 = input,
              features              = names(feas_sig_op),
              condiction            = rf_sig_condiction,
              col_condiction        = "condiction2",
              id_condiction         = "vars",
              group                 = "BOR",
              scale                 = FALSE,
              palette               = 6,
              width                 = 7,
              show_heatmap_col_name = TRUE,
              angle_col             = 60,
              path                  = file_name$folder_name,
              index                 = paste0("6-", pat_id, "-signature-os"),
              show_plot             = TRUE)
 return(lira_input)
}
