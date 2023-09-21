









#' Location of lira_score of each patient
#'
#' @param score  score
#' @param palette palette of response
#' @param showplot default is FALSE
#' @param path path to save result
#' @param palette_line palette of line
#' @param cols manual colors for objective response
#' @param panel option = `OS`, `PFS`
#' @param ref_score reference score of plot
#' @param id_score patient identifier
#' @param col_score
#' @param best_cutoff
#' @param bins_width
#' @param pat_id
#' @param index
#' @param title
#'
#' @return
#' @export
#'
#' @examples
lira_score_location<-function(score, pat_id, id_score = "ID", col_score = "riskscore", best_cutoff = NULL, ref_score = NULL,
                              palette = "nrc", cols = NULL, palette_line = "jama", showplot = TRUE, path = NULL, panel = "OS", bins_width = 33, index = NULL, title = "LIRA score"){


  score<-as.data.frame(score)
  if(!is.null(path)){
    if(!file.exists(path)) dir.create(path)
  }else{
    path<- "lira_score-location"
    if(!file.exists(path)) dir.create(path)
  }

  if(is.null(cols)){
    cols<-IOBR::palettes(category = "box",
                         palette = palette,
                         show_message = FALSE,
                         show_col = FALSE, alpha = 0.75)
  }else{
    cols<-cols
  }
  cols2<-IOBR::palettes(category = "box",
                        palette = palette_line,
                        show_message = FALSE,
                        show_col = FALSE,
                        alpha = 1)
  ###############################

  colnames(score)[which(colnames(score)==id_score)] <- "ID"
  colnames(score)[which(colnames(score)==col_score)] <- "riskscore"

  if(is.null(ref_score)){
    ref_score <- score
    print(ref_score[ref_score$ID==pat_id, ])
    ref_score[ref_score$ID==pat_id, "BOR"] <- "NE"
  }else{
    data("ref_lira_score", package = "LIRA")
    ref_score <- ref_lira_score
  }
  ###############################
  ###############################
  var<- title

  pat<-pat_id
  # print(paste0(">>> Processing patient: ", pat))

  if(is.null(best_cutoff)){
    if(panel == "PFS"){
      cutoff_all<-18.948
    }else if(panel == "OS"){
      cutoff_all<-19.464
    }
  }else{
    cutoff_all <- best_cutoff
  }

  pat_score<-score[score$ID==pat,]$riskscore
  pat_score<-round(pat_score, 3)
  message(paste0(">>> ", title, " of ", pat, " is ", pat_score))
  target<-sym("riskscore")

  pat_split<-unlist(stringr::str_split(pat, pattern = "_"))

  subt<-paste0("Sample name: ", pat_split[1])
  ###################################################################
  #参考链接：https://blog.csdn.net/weixin_45387324/article/details/109408376

  p<- ggplot(ref_score, aes(x= !!target, fill= BOR)) +
    geom_histogram(bins = bins_width, alpha = 0.76, colour = "grey")+    #aes(y=..density..),
    scale_fill_manual(values= cols)+
    geom_density(alpha=.2, fill="grey", weight = 1)+

    labs(title=  paste0(title, " = ", pat_score),
         subtitle= paste0(subt),
         caption = paste0(" Data of RNAseq: ",panel, ";  ","BC: best cutoff;   ", date()))+

    # xlab(paste0(target))+
    theme_light()+
    design_mytheme(legend.position = "bottom", axis_angle = 0, plot_title_size = 1.7, axis_title_size = 1.4)+
    xlab(paste0(title, " of 439 NSCLC patients treated with ICB")) +ylab("Count of patients")

  if(pat_score >= cutoff_all){
    cols_pat <- '#B24745FF'
  }else{
    cols_pat <- '#00A1D5FF'
  }
  p<-p+geom_vline(aes(xintercept = cutoff_all),
                  linetype="dashed",
                  color = cols2[1],
                  size = 0.45)+
    annotate(geom = "text", fontface = "plain", color= cols2[1],
             x = cutoff_all - 1.5, y= 15,hjust = 0,
             label = paste0('BC of ', panel, ' = ', cutoff_all), size=4.5, angle = 60)+


    geom_vline(aes(xintercept = pat_score),
               linetype="dashed",color = "black", size = 0.55)+
    annotate(geom = "text", fontface = "plain", color= cols_pat,
             x = pat_score - 1.5, y = 20, hjust = 0,
             label = paste0( var, ' = ', pat_score), size= 5.5, angle = 60)

  if(showplot) print(p)
  if(is.null(index)){
    ggsave(p,filename =paste0(pat,"-",var,".pdf"),
           width = 7.64,height = 5.76, path = path, dpi = 300)
  }else{
    ggsave(p,filename =paste0(index, "-", pat,"-",  var, ".pdf"),
           width = 7.64,height = 5.76, path = path, dpi = 300)
  }



}
