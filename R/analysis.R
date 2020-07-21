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
    print(paste("Random sampling",max_regions,"regions from all TFBS with seed",char2seed(tf_name,set=FALSE)))
    char2seed(tf_name)
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

  out_file=paste0(output_dir,"/",sample_name,"_",tf_name,".",max(log_data$TFBS_ANALYZED),"TFBS.S",tfbs_start,"-E",tfbs_end,".tss")
  write.table(log_data,quote=FALSE,row.names=FALSE,out_file)

  print("Generating plots")
  tictoc::tic("Generation time")
  plot_motif_coverage(log_data,tf_name=tf_name,sample_name=sample_name,output_dir=output_dir)
  tictoc::toc()

  print(paste("Analysis finished for", tf_name))
  tictoc::toc()
  print("######################################################")
  rm(list = ls())
  gc()
}


#' Estimates the unranked Accesibility score for TFBS
#'
#' This function takes a DATA.FRAME/ TXT file with the TFBS coverage data and
#' estimates the Accesibility score as the range between the global maximum and minimum of
#' the high-frequency signal. In order, to estimate the high-frequency signal, and remove
#' the local biases the Savitzky-Golay filter is used.

#' @param data DATA.FRAME with data or Path to TXT file
#' @param output_dir Directory to output results.
#' @export



accessibility_score=function(data="",output_dir="",name=""){

    tictoc::tic("Calculation time")

    if(!is.data.frame(data)){
      cov_data=read.table(data,header=TRUE)
      name=ULPwgs::get_sample_name(data)
    }else{
      cov_data=data
    }


    sep="/"

    if(output_dir==""){
      sep=""
    }


    filter_length=max(cov_data$POSITION_RELATIVE_TO_TFBS)

    cov_data$LOW<-get_low_signal(cov_data$MEAN_DEPTH,ifelse(filter_length%%2==0,filter_length+1,filter_length))
    cov_data$HIGH<-get_high_signal(cov_data$MEAN_DEPTH,cov_data$LOW,ifelse(floor(filter_length/20)%%2==0,floor(filter_length/20)+1,floor(filter_length/20)))

    n=floor(mean(cov_data$TFBS_ANALYZED))
    range=max(cov_data$HIGH) - min(cov_data$LOW)
    peaks=find_peaks(cov_data$HIGH,m=20)
    peak_positions = cov_data$POSITION_RELATIVE_TO_TFBS[peaks]
    peak_distance = c(diff(peak_positions))
    mean_peak_distance = mean(peak_distance)
    median_peak_distance = median(peak_distance)


    stats=list(TF=name,MEAN_NUMBER_TFBS_ANALYZED=n,RANGE=range,MEAN_PEAK_DISTANCE=mean_peak_distance,MEDIAN_PEAK_DISTANCE=median_peak_distance,PEAKS=peaks,PEAK_POSITIONS=peak_positions,PEAK_DISTANCE =peak_distance)

    info=list(COV_DATA=cov_data,STATS=stats)
    tictoc::toc()

    out_file=paste0(output_dir,sep,name,".",max(cov_data$TFBS_ANALYZED),"TFBS.S",abs(min(cov_data$POSITION_RELATIVE_TO_TFBS)),"-E",max(cov_data$POSITION_RELATIVE_TO_TFBS),".FREQUENCY.txt")

    ## Generate LOG with data

    cat(paste(Sys.time(),"\n\n"),file=out_file,append=FALSE)
    cat(paste("## COVERAGE \n"),file=log_file,append=TRUE)
    write.table(cov_data,file=out_file,append=TRUE,sep="\t",quote=FALSE,row.names=FALSE)
    cat("\n",file=out_file,append=TRUE)
    cat(paste("## STATS \n"),file=log_file,append=TRUE)
    cat(stats,file=out_file,append=TRUE)

    return(info)
  }
