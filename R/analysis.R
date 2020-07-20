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
#' @param tss_start Number of bases to analyze forward from TFBS central point. Default 1000
#' @param tss_end Number of bases to analyze  backward from TFBS central point. Default 1000
#' @param mean_cov Mean genome wide coverage. If not provided it will be estimated.
#' @param cov_limit Max base depth. Default 1000
#' @param max_regions Max number of TFBS to analyze. Default 100000
#' @param mapq Min quality of mapping reads. Default 0
#' @param threads Number of threads. Default 1
#' @param verbose Enables progress messages. Default False.
#' @export


analyze_tss_around_position=function(bin_path="tools/bedtools2/bin/bedtools",bin_path2="tools/samtools/samtools",bed="",bam="",tss_start=1000,tss_end=1000,mean_cov="",norm="",threads=1,cov_limit=1000,max_regions=100000,mapq=0,verbose=FALSE){
  tic("Analysis time")


  options(scipen=999)
  ref_data=read.table(bed)
  sample_name=ULPwgs::get_sample_name(bam)
  tf_name=ULPwgs::get_sample_name(bed)

  sep="/"

  if(output_dir==""){
    sep=""
  }

  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }

  output_dir=paste0(output_dir,sep,sample_name)

  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }

  output_dir=paste0(output_dir,"/",tf_name)


  if(!dir.exists(output_dir)){
      dir.create(output_dir)
  }


  print(paste("Analyzing",tf_name,"Binding Sites for",sample_name))
  print(paste(nrow(ref_data)," TSS found in total."))

  if ((max_regions) !=0 & max_regions<nrow(ref_data)){
    print(paste("Total TSS > Maximum Number of regions to analyze",paste0("(",max_regions,")")))
    print(paste("Random sampling",max_regions,"regions from all TSS."))
    char2seed(tf_name,FALSE)
    ref_data=sample(ref_data)

  }

  # Calculate Mean Genome-wide Coverage

  if (!mean_cov==""){
    print("Calculating Genome-Wide Coverage")
    tictoc::tic("")
    system.time(calculate_genowide_coverage(bin_path=bin_path,bam=bam,verbose=verbose))
    mean_cov=get_mean_coverage(file=paste0(sample_name,"_GENOME_COVERAGE/",sample_name,"_genome_coverage.txt"),output_dir=paste0(sample_name,"_GENOME_COVERAGE/"),sample_name=sample_name,save=TRUE)
    tictoc::toc("Calculation time")
  }
  print(paste0("Mean genome-wide coverage:",mean_cov))

  # Read normalized local coverage

  if(!norm==""){
    norm_log2=read.table(norm,header=TRUE)
    colnames(norm_log2)=c("chr","start","end","log2")
    if(any(is.na(norm_log2))){
      warning("NAs found substituted with 0s.")
      norm_log2[is.na(norm_log2)]=0
    }
  }

  # Calculate Mean Depth Coverage Around TFBS

  print("Calculating Mean Depth Coverage Around TFBS")
  tic("Calculation time")
  coverage_list=calculate_coverage_tss(bin_path=bin_path2,ref_data=ref_data,bam=bam,norm_log2=norm_log2,tss_start=tss_start,tss_end=tss_end,cov_limit=cov_limit,output_dir=output_dir,mapq=mapq,tf_name=tf_name,sample_name=sample_name)
  log_data=get_mean_and_conf_intervals(cov_data=coverage_list)
  toc()
  print("Generating plots")
  tic("Generation time")
  plot_motif_coverage(log_data,tf_name=tf_name,sample_name=sample_name)
  toc()

  print("Analysis finished for", tf_name)
  toc()
  print("######################################################")
}
