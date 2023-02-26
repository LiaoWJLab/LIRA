









#' Location of lira_score of each patient
#'
#' @param score  score
#' @param palette palette of response
#' @param showplot default is FALSE
#' @param path path to save result
#' @param palette_line palette of line
#' @param cols 
#' @param panel 
#' @param ref_score 
#'
#' @return
#' @export
#'
#' @examples
lira_score_location<-function(score, ref_score = ref_lira_score, palette = "nrc", cols = NULL, palette_line = "jama", showplot = TRUE, path = NULL, panel = "OS"){
  
  
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
  # data(ref_lira_score)
  ###############################
  pats<-as.character(score$ID)
  
  var<- "LIRA score"
  for (i in 1:length(pats)) {
    
    pat<-pats[i]
    # print(paste0(">>> Processing patient: ", pat))
    
    
    if(panel == "PFS"){

      cutoff_all<-18.948
      pat_score<-score[score$ID==pat,]$riskscore
        
    }else if(panel == "OS"){
      
      cutoff_all<-19.464
      pat_score<-score[score$ID==pat,]$riskscore
    }
    
    
    pat_score<-round(pat_score, 3)
    message(paste0(">>> ", "LIRA score", " of ", pat, " is ", pat_score))
    target<-sym("riskscore")

    pat_split<-unlist(stringr::str_split(pat, pattern = "_"))
    
    subt<-paste0("Sample name: ", pat_split[1])
    ###################################################################
    #参考链接：https://blog.csdn.net/weixin_45387324/article/details/109408376
    p<-ggplot(ref_score, aes(x= !!target, fill= BOR)) +
      geom_histogram(aes(y=..density..), bins = 30, colour = "grey", alpha = 0.5)+
      scale_fill_manual(values= cols)+
      geom_density(alpha=.2, fill="grey", weight = 1)+
      
      labs(title=  paste0("LIRA score", " = ", pat_score),
           subtitle= paste0(subt),
           caption = paste0(" Data of RNAseq: ",panel, ";  ","BC: best cutoff;   ", date()))+
      
      # xlab(paste0(target))+
      theme_light()+
      design_mytheme(legend.position = "bottom", axis_angle = 0, plot_title_size = 2)+
      xlab("LIRA score")
    
    p<-p+geom_vline(aes(xintercept = cutoff_all),
                    linetype="dashed", 
                    color = cols2[1],
                    size = 0.70)+
      annotate(geom = "text", fontface = "plain", color= cols2[1],
               x = 15, y= 0.15,hjust = 0,
               label = paste0('BC of ', panel, ' = ', cutoff_all), size=4.5)+
      
      
      geom_vline(aes(xintercept = pat_score),
                 linetype="dashed",color = "black", size = 0.70)+
      annotate(geom = "text", fontface = "plain", color= "black",
               x = pat_score - 8, y = 0.23,hjust = 0,
               label = paste0( var, ' of smaple = ', pat_score), size=4.5)
    
    if(showplot) print(p)
    ggsave(p,filename =paste0(i,"-",pat,"-",var,".pdf"),
           width = 7.64,height = 5.76, path = path, dpi = 300)
    
  }
  
}
