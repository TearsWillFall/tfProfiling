#' Analyze TFs for a sample
#'
#' This function takes the path to the directory with BED files for each TF, a BAM file with an aligned sequence and a TXT file with the normalized
#' local coverage values for CNA, and estimates the corrected mean depth coverage values around TFBS, as well as the accessibilty scores for them.
#' For better understanding check: https://www.nature.com/articles/s41467-019-12714-4
#'
#'
#' @param bin_path Path to binary. Default tools/bedtools2/bin/bedtools
#' @param bin_path2 Path to secondary binary. Default tools/samtools/samtools
#' @param bed_dir Path to directory with BED files for TFBS
#' @param tfs Number of TFs to analyze or a list with TFs to analyze. Default none will analyze all TFs in BED directory.
#' @param bam Path to BAM file
#' @param tfbs_start Number of bases to analyze forward from TFBS central point. Default 1000
#' @param tfbs_end Number of bases to analyze  backward from TFBS central point. Default 1000
#' @param mean_cov Mean genome wide coverage. If not provided it will be estimated.
#' @param norm Path to TXT file with normalized local coverage
#' @param cov_limit Max base depth. Default 1000
#' @param max_regions Max number of TFBS to analyze. Default 100000
#' @param mapq Min quality of mapping reads. Default 0
#' @param threads Number of threads. Default 1
#' @param verbose Enables progress messages. Default FALSE
#' @param output_dir Directory to output results. If not provided then outputs in current directory
#' @param plot Create plots. Default TRUE.
#' @return A DATA.FRAME with coverage data
#' @export


analyze_tfs=function(bin_path="tools/bedtools2/bin/bedtools",bin_path2="tools/samtools/samtools",bed_dir="",tfs="",bam="",tfbs_start=1000,tfbs_end=1000,mean_cov="",norm="",threads=1,cov_limit=1000,max_regions=100000,mapq=0,verbose=FALSE,output_dir="",plot=TRUE){

  sep="/"
  if(output_dir==""){
    sep=""
  }


tictoc::tic("Time")

sample_name=ULPwgs::get_sample_name(bam)

print(paste("Analyzing sample ",sample_name))


bed_files=list.files(bed_dir,full.names=TRUE)
bed_files=bed_files[grepl(".bed",bed_files)]
if (is.numeric(tfs)){
  if (tfs<length(bed_files)){
  print(paste(length(bed_files),"TFs in total"))
  print(paste("Total TFs > Number of TFs to analyze",paste0("(",tfs,")")))
  options(warn=-1)
  print(paste("Random sampling",tfs,"TFs from all TFs with seed",char2seed(bed_dir,set=FALSE)))
  char2seed(bed_dir)
  bed_files=sample(bed_files,tfs)
  options(warn=0)
}else{
  bed_files
}
}else if(is.vector(tfs)&& !tfs==""){
  bed_files=bed_files[grepl(paste(tfs,collapse="|"),bed_files)]
  if (length(bed_files)<length(tfs)){
    warning("Number of TFs to analyze > Number of TFs listed. Perhaps the TFs names provided are too unspecific?")
  }
}

print("=========================================================")

FUN=function(bed,bin_path,bin_path2,bam,tfbs_start,tfbs_end,mean_cov,norm,threads,cov_limit,max_regions,mapq,verbose,output_dir,plot){
  accessibility_score(analyze_tfbs_around_position(bin_path=bin_path,bin_path2=bin_path2,bam=bam,bed=bed,norm=norm,threads=threads,tfbs_start=tfbs_start,tfbs_end=tfbs_end,output_dir=output_dir,plot=plot,mapq=mapq,cov_limit=cov_limit,mean_cov=mean_cov,max_regions=max_regions,verbose=verbose),output_dir=output_dir,verbose=verbose)
}

all_stats=mapply(bed_files,FUN=FUN,SIMPLIFY=FALSE,bin_path=bin_path,bin_path2=bin_path2,bam=bam,norm=norm,threads=threads,tfbs_start=tfbs_start,tfbs_end=tfbs_end,output_dir=output_dir,plot=plot,mapq=mapq,cov_limit=cov_limit,mean_cov=mean_cov,max_regions=max_regions,verbose=verbose)

all_stats=suppressMessages(all_stats %>% dplyr::bind_rows())

output_dir=paste0(output_dir,sep,sample_name)

out_file=paste0(output_dir,"/",sample_name,".ALL.ANALYZED.TFS.STATS.txt")


if (file.exists(out_file)){
  write.table(all_stats,quote=FALSE,row.names=FALSE,out_file,append=TRUE,col.names=FALSE)

}else{
  write.table(all_stats,quote=FALSE,row.names=FALSE,out_file)
}

tictoc::toc()


}


#' Analyze Mean Coverage Depth around TFBS
#'
#' This function takes a BED file with TFBS, a BAM file with an aligned sequence and a TXT file with the normalized
#' local coverage values for CNA, and outputs the corrected mean depth coverage values around TFBS, as well as
#' generates a PDF with a plot for them. For better understanding check: https://www.nature.com/articles/s41467-019-12714-4
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
#' @param verbose Enables progress messages. Default FALSE
#' @param output_dir Directory to output results. If not provided then outputs in current directory
#' @param plot Create a plot with coverage data. Default TRUE.
#' @return A DATA.FRAME with coverage data
#' @export


### ////TO DO IMPLEMENT verbose to the rest of the functions

analyze_tfbs_around_position=function(bin_path="tools/bedtools2/bin/bedtools",bin_path2="tools/samtools/samtools",bed="",bam="",tfbs_start=1000,tfbs_end=1000,mean_cov="",norm="",threads=1,cov_limit=1000,max_regions=100000,mapq=0,verbose=FALSE,output_dir="",plot=TRUE){
  tictoc::tic("Analysis time")


  options(scipen=999)

  chr_check=system(paste(bin_path2,"view",bam," | head -n 1 | awk -F \"\t\" '{print $3}'"),intern=TRUE)

  ref_data=read.table(bed,comment.char="")
  if (!grepl("chr",chr_check)){
    ref_data[,1]=gsub("chr","",ref_data[,1])
  }


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
    warning("NAs found substituted with 1s")
    norm_log2[is.na(norm_log2)]=1
  }


  # Calculate Mean Depth Coverage Around TFBS

  print("Calculating Mean Depth Coverage Around TFBS")
  coverage_list=calculate_coverage_tfbs(bin_path=bin_path2,ref_data=ref_data,bam=bam,norm_log2=norm_log2,tfbs_start=tfbs_start,tfbs_end=tfbs_end,cov_limit=cov_limit,output_dir=output_dir,mapq=mapq,tf_name=tf_name,sample_name=sample_name,threads=threads,mean_cov=mean_cov)

  log_data=get_mean_and_conf_intervals(cov_data=coverage_list)
  log_data=cbind(TF=paste0(sample_name,"_",tf_name),log_data)

  out_file=paste0(output_dir,"/",sample_name,"_",tf_name,".",max(log_data$TFBS_ANALYZED),"TFBS.S",tfbs_start,"-E",tfbs_end,".tss")

  write.table(log_data,quote=FALSE,row.names=FALSE,out_file)
  if(plot){
    print("Generating plots")
    tictoc::tic("Generation time")
    plot_motif_coverage(log_data,tf_name=tf_name,sample_name=sample_name,output_dir=output_dir)
    tictoc::toc()
  }


  print(paste("Analysis finished for", tf_name))
  tictoc::toc()
  print("######################################################")

  return(log_data)
}


#' Estimates the unranked Accesibility score for TFBS
#'
#' This function takes a DATA.FRAME/ TXT file with the TFBS coverage data and
#' estimates the Accesibility score as the range between the global maximum and minimum of
#' the high-frequency signal. In order, to estimate the high-frequency signal, and remove
#' the local biases the Savitzky-Golay filter is used.

#' @param data DATA.FRAME with data or Path to TXT file
#' @param output_dir Directory to output results. If not provided then outputs in current directory
#' @param plot Create a plot with frequency data. Default TRUE.
#' @param verbose Enables progress messages. Default FALSE
#' @export



accessibility_score=function(data="",output_dir="",plot=TRUE,verbose=FALSE){


    tictoc::tic("Calculation time")

    if(!is.data.frame(data)){
      cov_data=read.table(data,header=TRUE)
      name=as.character(cov_data$TF[1])
    }else{
      cov_data=data
      name=as.character(cov_data$TF[1])
    }

    tf_name=strsplit(name,"_")[[1]][length(strsplit(name,"_")[[1]])]
    sample_name=paste0(strsplit(name,"_")[[1]][-length(strsplit(name,"_")[[1]])],collapse="_")


    print(paste("Estimating",tf_name,"Accesibility Score for",sample_name))

    sep="/"

    if(output_dir==""){
      sep=""
    }


    filter_length=max(cov_data$POSITION_RELATIVE_TO_TFBS)

    cov_data$LOW<-get_low_signal(cov_data$MEAN_DEPTH,ifelse(filter_length%%2==0,filter_length+1,filter_length))
    cov_data$HIGH<-get_high_signal(cov_data$MEAN_DEPTH,cov_data$LOW,ifelse(floor(filter_length/20)%%2==0,floor(filter_length/20)+1,floor(filter_length/20)))

    n=floor(mean(cov_data$TFBS_ANALYZED))
    range=max(cov_data$HIGH) - min(cov_data$HIGH)
    peaks=data.frame(PEAKS=find_peaks(cov_data$HIGH,m=20))
    peak_positions = data.frame(PEAK_POSITIONS=cov_data$POSITION_RELATIVE_TO_TFBS[peaks$PEAKS])
    peak_distance = data.frame(PEAK_DISTANCE=c(diff(peak_positions$PEAK_POSITIONS)))
    mean_peak_distance = mean(peak_distance$PEAK_DISTANCE)
    median_peak_distance = median(peak_distance$PEAK_DISTANCE)


    stats=data.frame(TF=cov_data$TF[1],MEAN_NUMBER_TFBS_ANALYZED=n,RANGE=range,MEAN_PEAK_DISTANCE=mean_peak_distance,MEDIAN_PEAK_DISTANCE=median_peak_distance)
    info=list(COV_DATA=cov_data,STATS=stats)


    output_dir=paste0(output_dir,sep,sample_name,"/",tf_name)
    out_file=paste0(output_dir,"/",name,".",max(cov_data$TFBS_ANALYZED),"TFBS.S",abs(min(cov_data$POSITION_RELATIVE_TO_TFBS)),"-E",max(cov_data$POSITION_RELATIVE_TO_TFBS),".FREQUENCY.txt")

    ## Generate LOG with data
    options(warn = -1)
    cat(paste(Sys.time(),"\n\n"),file=out_file,append=FALSE)
    cat(paste("# COVERAGE \n"),file=out_file,append=TRUE)
    write.table(cov_data,file=out_file,append=TRUE,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
    cat("\n",file=out_file,append=TRUE)
    cat(paste("## STATS \n"),file=out_file,append=TRUE)
    write.table(stats,file=out_file,append=TRUE,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
    cat("\n",file=out_file,append=TRUE)
    cat(paste("### PEAKS \n"),file=out_file,append=TRUE)
    write.table(peaks,file=out_file,append=TRUE,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
    cat("\n",file=out_file,append=TRUE)
    cat(paste("### PEAK_POSITIONS \n"),file=out_file,append=TRUE)
    write.table(peak_positions,file=out_file,append=TRUE,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
    cat("\n",file=out_file,append=TRUE)
    cat(paste("#### PEAK_DISTANCE \n"),file=out_file,append=TRUE)
    write.table(peak_distance,file=out_file,append=TRUE,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
    options(warn = 0)
    if(plot){
      print("Generating plots")
      tictoc::tic("Generation time")
      plot_freq_decomposition(info,output_dir)
      tictoc::tic()
    }
  tictoc::toc()
  print("-------------------------------------------------")
  return(stats)
  }


#' Ranks the Accessibility Score for all the TFs
#'
#' This function takes a path to a TXT file with the accessibility scores of all the TFs analyzed and ranks them.

#' @param data Path to TXT file
#' @param output_dir Directory to output results. If not provided then outputs in current directory
#' @param verbose Enables progress messages. Default FALSE
#' @export


rank_accessibility=function(data="",output_dir="",verbose=FALSE){

  stats_data=read.table(data,header=TRUE,comment.char="")
  sample_name=ULPwgs::get_sample_name(data)


  print(paste("Ranking Accesibility Score for",sample_name))

  sep="/"

  if(output_dir==""){
    sep=""
  }

  results=data.frame(TF=stats_data$TF,RANGE=stats_data$RANGE,RANK=rank(stats_data$RANGE)/length(rank(stats_data$RANGE)))

  output_dir=paste0(output_dir,sep,sample_name)

  out_file=paste0(output_dir,"/",sample_name,".RANKED.ACCESSIBILITY.SCORE.txt")



  write.table(results,quote=FALSE,row.names=FALSE,out_file)








  }
