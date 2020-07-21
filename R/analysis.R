#' Analyze Mean Coverage Depth around TFBS
#'
#' This function takes a BED file with TFBS, a BAM file with an aligned sequence and a TXT file with the normalized
#' local coverage values for CNA, and outputs the corrected mean depth coverage values around TFBS, as well as
#' generates a pdf with a plot for them. For better understanding check: https://www.nature.com/articles/s41467-019-12714-4
#'
#'
#' @param bin_path Path to binary. Default tools/bedtools2/bin/bedtools
#' @param bin_path2 Path to secondary binary. Default tools/samtools/samtools
#' @param bed Path to BED file.
#' @param bam Path to BAM file.
#' @param tfbs_start Number of bases to analyze forward from TFBS central point. Default 1000
#' @param tfbs_end Number of bases to analyze  backward from TFBS central point. Default 1000
#' @param mean_cov Mean genome wide coverage. If not provided it will be estimated.
#' @param norm Path to TXT file with normalized local coverage
#' @param cov_limit Max base depth. Default 1000
#' @param max_regions Max number of TFBS to analyze. Default 100000
#' @param mapq Min quality of mapping reads. Default 0
#' @param threads Number of threads. Default 1
#' @param verbose Enables progress messages. Default False.
#' @param output_dir Directory to output results.
#' @export


analyze_tfbs_around_position=function(bin_path="tools/bedtools2/bin/bedtools",bin_path2="tools/samtools/samtools",bed="",bam="",tfbs_start=1000,tfbs_end=1000,mean_cov="",norm="",threads=1,cov_limit=1000,max_regions=100000,mapq=0,verbose=FALSE,output_dir=""){
  tictoc::tic("Analysis time")


  options(scipen=999)
  ref_data=read.table(bed)
  sample_name=ULPwgs::get_sample_name(bam)
  tf_name=ULPwgs::get_sample_name(bed)


  print(paste("Analyzing",tf_name,"Binding Sites for",sample_name))

  sep="/"

  if(output_dir==""){
    sep=""
  }

  output_dir=paste0(output_dir,sep,sample_name)

  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }

  # Calculate Mean Genome-wide Coverage

  if (mean_cov==""){
    if (!file.exists(paste0(output_dir,"/",sample_name,"_GENOME_COVERAGE/",sample_name,"_genome_coverage.txt"))){
      print("Calculating Genome-Wide Coverage")
      tictoc::tic("Calculation time")
      system.time(calculate_genowide_coverage(bin_path=bin_path,bam=bam,verbose=verbose,output_dir=output_dir))
      mean_cov=get_mean_coverage(file=paste0(output_dir,"/",sample_name,"_GENOME_COVERAGE/",sample_name,"_genome_coverage.txt"),output_dir=paste0(output_dir,"/",sample_name,"_GENOME_COVERAGE/"),sample_name=sample_name,save=TRUE)
      tictoc::toc()
    }else{
      print(paste("Genome-wide Coverage for sample",sample_name,"already exists"))
      print(paste("Loading Genome-wide Coverage"))
      mean_cov=get_mean_coverage(file=paste0(output_dir,"/",sample_name,"_GENOME_COVERAGE/",sample_name,"_genome_coverage.txt"),output_dir=paste0(output_dir,"/",sample_name,"_GENOME_COVERAGE/"),sample_name=sample_name,save=TRUE)
    }

  }
  print(paste0("Mean genome-wide coverage: ",mean_cov))




  output_dir=paste0(output_dir,"/",tf_name)


  if(!dir.exists(output_dir)){
      dir.create(output_dir)
  }



  print(paste(nrow(ref_data)," TFBS found in total"))

  if ((max_regions) !=0 & max_regions<nrow(ref_data)){
    print(paste("Total TFBS > Maximum Number of regions to analyze",paste0("(",max_regions,")")))
    print(paste("Random sampling",max_regions,"regions from all TFBS"))
    char2seed(tf_name,set=TRUE)
    ref_data=ref_data[sample(c(1:nrow(ref_data)),max_regions),]

  }



  # Read normalized local coverage

  norm_log2=read.table(norm,header=TRUE)
  colnames(norm_log2)=c("chr","start","end","log2")
  if(any(is.na(norm_log2))){
    warning("NAs found substituted with 0s")
    norm_log2[is.na(norm_log2)]=0
  }


  # Calculate Mean Depth Coverage Around TFBS

  print("Calculating Mean Depth Coverage Around TFBS")
  coverage_list=calculate_coverage_tfbs(bin_path=bin_path2,ref_data=ref_data,bam=bam,norm_log2=norm_log2,tfbs_start=tfbs_start,tfbs_end=tfbs_end,cov_limit=cov_limit,output_dir=output_dir,mapq=mapq,tf_name=tf_name,sample_name=sample_name,threads=threads,mean_cov=mean_cov)
  log_data=get_mean_and_conf_intervals(cov_data=coverage_list)

  out_file=paste0(output_dir,"/",tf_name,".tss")

  write.table(log_data,quote=FALSE,row.names=FALSE,out_file)

  print("Generating plots")
  tictoc::tic("Generation time")
  plot_motif_coverage(log_data,tf_name=tf_name,sample_name=sample_name,output_dir=output_dir)
  tictoc::toc()

  print(paste("Analysis finished for", tf_name))
  tictoc::toc()
  print("######################################################")
  gc()
}
file="~/SRR11742864/Androgen/Androgen.tss"



library(ggplot2)
accessibility_score=function(file="",output_dir=""){
    cov_data=read.table(file,header=TRUE)
    sample_name=ULPwgs::get_sample_name(file)
    cov_data$LOW<-get_low_signal(cov_data$MEAN_DEPTH)
    cov_data$HIGH<-get_high_signal(cov_data$MEAN_DEPTH,cov_data$LOW)

    stats=()
    n<-cov_data$TFBS_ANALYZED[1]
    range<-max(cov_data$HIGH) - min(cov_data$LOW)
    peaks<-find_peaks(cov_data$HIGH,m=20)
    peak_positions = cov_data$POSITION_RELATIVE_TO_TFBS[peaks]
    peak_distance = c(diff(peak_positions))
    mean_peak_distance = mean(peak_distance)
    median_peak_distance = median(peak_distance)

    df <- cov_data %>%
    dplyr::select(POSITION_RELATIVE_TO_TFBS,MEAN_DEPTH, HIGH, LOW) %>% dplyr::mutate(MEAN_DEPTH=MEAN_DEPTH/mean(MEAN_DEPTH)) %>%
    tidyr::gather(key = "variable", value = "TYPE",-POSITION_RELATIVE_TO_TFBS)

    p1=ggplot(cov_data,aes(x=POSITION_RELATIVE_TO_TFBS, y = MEAN_DEPTH/mean(MEAN_DEPTH))) +
    geom_line(col="black") +theme_classic() +labs(y="NORM_MEAN_DEPTH") + ggtitle("ORIGINAL")

    p2=ggplot(cov_data,aes(x=POSITION_RELATIVE_TO_TFBS, y = HIGH)) +
    geom_line(col="darkred") +theme_classic() + ggtitle("HIGH_FREQUENCY") +labs(y="NORM_MEAN_DEPTH")

    p3=ggplot(cov_data,aes(x=POSITION_RELATIVE_TO_TFBS, y = LOW)) +
    geom_line(col="steelblue") +theme_classic() + ggtitle("LOW_FREQUENCY") +labs(y="NORM_MEAN_DEPTH")


    p4=ggplot(df,aes(x=POSITION_RELATIVE_TO_TFBS, y = TYPE)) +
    geom_line(aes(color = variable)) +
    scale_color_manual(name = "TYPE", labels = c("HIGH","LOW","ORIGINAL"),values = c("darkred", "steelblue","black"))+
    theme_classic()+
    theme(legend.position="bottom") +labs(y="NORM_MEAN_DEPTH")+
    ggtitle("ORIGINAL+HIGH_&_LOW_FREQUENCY")


          p=gridExtra::grid.arrange(
                  grobs = list(p1,p2,p3,p4),
                  widths = c(1, 1, 1,1),
                  layout_matrix = rbind(c(1, 1, 1,1),
                                        c(2, 2, 3,3),
                                        c(4, 4, 4,4)),
                common.legend = TRUE, legend="right")

                    p=p+annotation_custom(arrow)


        fit <- loess(cov_data$HIGH ~ log(cov_data$POSITION_RELATIVE_TO_TFBS),data=cov_data)
    data$prior<-predict(fit,data)

    png(args[3])
    ggplot()+geom_point(aes(x=data$V2,y=data$prior),color="blue")+geom_point(aes(x=data$V2,y=data$V3),color=rgb(0,0,0,0.3))+ylim(0,1)
    dev.off()

    data$adjusted_range<-data$V3 / data$prior
    data$adjusted_range_normalized <- data$adjusted_range / data$adjusted_range[which(data$V1 == "CTCF.50perc_hg19.tss")]
    data$adjusted_range_rank <- rank(data$adjusted_range)
    #ggplot(filtered_data,aes(x=V2,y=adjusted_range))+geom_point()




  }
