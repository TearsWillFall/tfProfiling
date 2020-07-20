#' Plot Mean Depth Coverage for TFBS
#'
#' This function takes a DATA.FRAME with Mean Depth Coverage data as input, and plots it
#'
#' @param data Data.frame with mean depth coverage data.
#' @param trend_line Draw trend line. Default TRUE
#' @param tf_name Transcription Factor name
#' @param sample_name Sample name
#' @param output_dir Directory to output results.
#' @export




plot_motif_coverage=function(data="",trend_line=TRUE,tf_name="",sample_name="",output_dir=""){
  sep="/"

  if(output_dir==""){
    sep=""
  }
  out_file=paste0(output_dir,sep,title,".pdf")
  pdf(out)
  p=ggplot2::ggplot(data,  ggplot2::aes(x=POSITION_RELATIVE_TO_TFBS,y=MEAN_DEPTH))+
  ggplot2::geom_ribbon(ggplot2::aes(ymin=CI95_LOWER_BOUND, ymax=CI95_UPPER_BOUND), fill="red", alpha=0.1) +
  ggplot2::ggtitle(paste(tf_name,"for sample",sample_name)) +
  ggplot2::theme_classic()
  if(trend){
      p=p+ggplot2::geom_line(col="red")+ ggplot2::geom_smooth(method = "loess", formula = y ~ x, size = 1)
  }
  print(p)
  dev.off()
}
