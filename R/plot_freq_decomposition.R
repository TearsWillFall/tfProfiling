#' Plot frequency decomposition for TFBS coverage data
#'
#' This function takes a DATA.FRAME/ TXT file with the frequency data of the TFBS coverage and
#' generates a plot for it.
#'
#' @param data LIST with data or Path to TXT file
#' @param output_dir Directory to output results.
#' @export
#' @import ggplot2

plot_freq_decomposition=function(data="",output_dir=""){

  if(!is.list(data)){
    cov_data=read.table(data,header=TRUE)
    name=ULPwgs::get_sample_name(data)
  }else{
    cov_data=data["COV_DATA"]
    name=data$STATS$TF
  }

  sep="/"

  if(output_dir==""){
    sep=""
  }


  df <- cov_data %>%dplyr::select(POSITION_RELATIVE_TO_TFBS,MEAN_DEPTH, HIGH, LOW) %>% dplyr::mutate(MEAN_DEPTH=MEAN_DEPTH/mean(MEAN_DEPTH)) %>%
  tidyr::gather(key = "FACTOR", value = "TYPE",-POSITION_RELATIVE_TO_TFBS) %>% dplyr::mutate(FACTOR=relevel(factor(FACTOR),"MEAN_DEPTH","HIGH","LOW")) %>% dplyr::mutate(SIZE=ifelse(FACTOR=="MEAN_DEPTH",0.1,ifelse(FACTOR=="HIGH",0.11,0.12)))

  p1=ggplot(cov_data,aes(x=POSITION_RELATIVE_TO_TFBS, y = MEAN_DEPTH/mean(MEAN_DEPTH))) +
  geom_line(col="red") +theme_classic() +labs(y="NORM_MEAN_DEPTH") + ggtitle("ORIGINAL")

  p2=ggplot(cov_data,aes(x=POSITION_RELATIVE_TO_TFBS, y = HIGH)) +
  geom_line(col="black") +theme_classic() + ggtitle("HIGH_FREQUENCY") +labs(y="NORM_MEAN_DEPTH")

  p3=ggplot(cov_data,aes(x=POSITION_RELATIVE_TO_TFBS, y = LOW)) +
  geom_line(col="blue") +theme_classic() + ggtitle("LOW_FREQUENCY") +labs(y="NORM_MEAN_DEPTH")


  p4=ggplot(df,aes(x=POSITION_RELATIVE_TO_TFBS, y = TYPE)) +
  geom_line(aes(color = FACTOR)) +
  scale_color_manual(name = "TYPE", labels = c("HIGH","LOW","ORIGINAL"),values = c("red","black", "blue"))+
  theme_classic()+
  theme(legend.position="bottom") +labs(y="NORM_MEAN_DEPTH")+
  ggtitle("ORIGINAL+HIGH_&_LOW_FREQUENCY") +caption(name)



  p=gridExtra::grid.arrange(
          grobs = list(p1,p2,p3,p4),
          widths = c(1, 1, 1,1),
          layout_matrix = rbind(c(1, 1, 1,1),
                                c(2, 2, 3,3),
                                c(4, 4, 4,4)),
        common.legend = TRUE, legend="right")


out_file=paste0(output_dir,sep,name,".",max(cov_data$TFBS_ANALYZED),"TFBS.S",abs(min(cov_data$POSITION_RELATIVE_TO_TFBS)),"-E",max(cov_data$POSITION_RELATIVE_TO_TFBS),".FREQUENCY.pdf")

pdf(out_file)
print(p)
dev.off()


}
