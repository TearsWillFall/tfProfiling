#' Calculate Genome-wide coverage for a BAM file
#'
#' This function takes a BAM file and generates a TXT file with its genome-wide coverage
#'
#' @param bin_path Path to binary. Default tools/bedtools2/bin/bedtools
#' @param bam Path to BAM file.
#' @param verbose Enables progress messages. Default False.
#' @param output_dir Directory to output results.
#' @export



calculate_genowide_coverage=function(bin_path="tools/bedtools2/bin/bedtools",bam="",verbose=FALSE,output_dir=""){

  sep="/"
  sample_name=ULPwgs::get_sample_name(bam)

  if(output_dir==""){
    sep=""
  }

  out_file=paste0(output_dir,sep,sample_name,"_GENOME_COVERAGE")

  if (!dir.exists(out_file)){
      dir.create(out_file)
    }

  out_file=paste0(out_file,"/",sample_name,"_genome_coverage.txt")

  if(verbose){
    print(paste(bin_path,"genomecov -ibam",bam,">",out_file))
  }
  system(paste(bin_path,"genomecov -ibam",bam,">",out_file))

}

#' Get normalized local coverage from segmentation data (CNA) for a single base position
#'
#' This function takes normalized data for local coverage as input,
#' as well as a valid genomic position for a single base and returns a log2 corrected value for said position
#'
#' @param pos Genomic position
#' @param chr Chromosome
#' @param norm_log2 Data.frame with normalized local coverage
#' @export


get_norm_local_coverage=function(pos="",chr="",norm_log2=""){
   norm=norm_log2[norm_log2$chr==chr & pos>norm_log2$start & pos<norm_log2$end,]$log2
   if(!length(norm)==0){
     return(norm^2)
   }else{
     return (1)
   }

}

#' Calculate mean and confidence intervals for positions relative to TFBSs
#'
#' This function takes a LIST of DATA.FRAMES as input(each DATA.FRAME for each TFBS), and returns a DATA.FRAME
#' with the mean and confidence intervals for the positions relative to them
#'
#' @param data LIST of DATA.FRAMES with the data
#' @param CI Confidence Interval
#' @export

get_mean_and_conf_intervals=function(cov_data="",CI=0.95){
  FUN=function(x,cov_data){
    dat=cov_data[[x]]
    return(dat[,"norm_cor_cov",drop=FALSE])
  }
  norm_cor_cov_list=lapply(seq(1:length(cov_data)),FUN=FUN,cov_data=cov_data)
  norm_cor_cov_df=suppressMessages(norm_cor_cov_list %>% dplyr::bind_cols())
  norm_cor_cov_means=rowMeans(norm_cor_cov_df,na.rm = TRUE)
  numb_analyz_tss=norm_cor_cov_df %>% is.na() %>% `!` %>% rowSums()
  if (all(numb_analyz_tss>3)){
      norm_cor_cov_mes=apply(norm_cor_cov_df,1,FUN=function(x){qt(CI,sum(!is.na(x))-1)*sd(x[!is.na(x)])/sqrt(sum(!is.na(x)))})
        results=data.frame(POSITION_RELATIVE_TO_TSS=cov_data[[1]]$pos_relative_to_tss,MEAN_DEPTH=norm_cor_cov_means,CI95_LOWER_BOUND=norm_cor_cov_means-norm_cor_cov_mes,CI95_UPPER_BOUND=norm_cor_cov_means+norm_cor_cov_mes,TSS_ANALYZED=numb_analyz_tss)
  }else{
  results=data.frame(POSITION_RELATIVE_TO_TSS=cov_data[[1]]$pos_relative_to_tss,MEAN_DEPTH=norm_cor_cov_means,TSS_ANALYZED=numb_analyz_tss)
  }
  return (results)
}

#' Calculate mean genome coverage
#'
#' This function takes a TXT file with genome-wide ocverage as input and returns the mean genome coverage
#'
#' @param file Path to file with genome-wide coverage
#' @param output_dir Directory to output results
#' @param region Region for which to get mean coverage. Default genome
#' @param sample_name Sample name
#' @param save Save as TXT. Default TRUE
#' @export


get_mean_coverage=function(file="",output_dir="",region="genome",sample_name="",save=TRUE){


  data=read.table(file)
  data=data %>% dplyr::mutate(tot_bases=(V2*V3)) %>% dplyr::group_by(V1) %>% dplyr::summarise(cnt=dplyr::n(),avg_cov=sum(tot_bases)/max(V4), .groups = 'drop') %>% as.data.frame()

  if (save){
    sep="/"


    if(output_dir==""){
      sep=""
    }
    out_file=paste0(output_dir,sep,sample_name,"_mean_genome_coverage.txt")
    write.table(data.frame(data$V1,data$avg_cov),file=out_file,sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

  }

  return(data[data$V1==region,]$avg_cov)
}

#' Calculate Mean Coverage Depth around TFBS
#'
#' This function takes a DATA.FRAME with TFBS positions, a BAM file with the sequence to analyze, a mean genome wide coverage,
#' local coverage values for CNA, and outputs the corrected mean depth coverage values around TFBS and
#' a TSS file with them. For better understanding check: https://www.nature.com/articles/s41467-019-12714-4
#'
#'
#' @param bin_path Path to binary. Default tools/samtools/samtools
#' @param ref_data Data from BED file as Data.frame
#' @param bam Path to BAM file.
#' @param sample_name Sample name
#' @param tf_name Transcription Factor name
#' @param mean_cov Mean genome wide coverage.
#' @param norm_log2 Data.frame with normalized local coverage
#' @param tss_start Number of bases to analyze forward from TFBS central point. Default 1000
#' @param tss_end Number of bases to analyze  backward from TFBS central point. Default 1000
#' @param cov_limit Max base depth. Default 1000
#' @param max_regions Max number of TFBS to analyze. Default 100000
#' @param mapq Min quality of mapping reads. Default 0
#' @param threads Number of threads. Default 1
#' @param output_dir Directory to output results
#' @export
#' @import pbapply

calculate_coverage_tss=function(bin_path="tools/samtools/samtools",ref_data="",bam="",sample_name="",tf_name="",mean_cov="",norm_log2="",tss_start=1000,tss_end=1000,cov_limit=1000,mapq=0,threads=1,output_dir=""){
  if(ncol(ref_data)<6){
    ref_data=data.frame(chr=ref_data[,1],start=ref_data[,2],end=ref_data[,3],strand="+",pos=as.integer((ref_data[,3]+ref_data[,2])/2))

  }else{
    ref_data=data.frame(chr=ref_data[,1],start=ref_data[,2],end=ref_data[,3],strand=ref_data[,which(ref_data[1,]=="+" |ref_data[1,]=="-")],pos=as.integer((ref_data[,3]+ref_data[,2])/2))
  }

  tss_to_analyze=ref_data

  ## Filter for overlapping TSS or duplicated TSS to save time

  tss_to_analyze=tss_to_analyze %>% dplyr::filter(!grepl("_",chr),pos-tss_start>1) %>% dplyr::distinct(chr, pos, .keep_all = TRUE)

  FUN=function(x,bin_path,norm_log2,start,end,mean_cov,bam,cov_limit,mapq){
  tss_data=t(x)

  ## IF running sambamba
  ## cov_data=read.csv(text=system(paste(bin_path,"depth base -t 3 -z -L ",paste0(tss_data$chr,":",as.numeric(tss_data$pos)-start,"-",as.numeric(tss_data$pos)+end),bam),intern=TRUE),header=TRUE,sep="\t")
  ## cov_data=cbind(cov_data[,1:3],strand=tss_data$strand)
  ## names(cov_data)=c("chr","pos","cov","strand")

  cov_data=read.csv(text=system(paste(bin_path,"depth -aa -Q",mapq, "-r",paste0(tss_data[1],":",as.numeric(tss_data[5])-start,"-",as.numeric(tss_data[5])+end),bam),intern=TRUE),header=FALSE,sep="\t")
  colnames(cov_data)=c("chr","pos","cov")
  norm_cov=get_norm_local_coverage(pos=tss_data[5],chr=tss_data[1],norm_log2=norm_log2)
  if (norm_cov==0){
    norm_cov=0.001
  }

  cov_data=cbind(cov_data,strand=tss_data[5])

  return(cov_data %>% dplyr::mutate(cor_cov=cov/as.numeric(mean_cov))  %>% dplyr::mutate(norm_cor_cov=ifelse(cor_cov<cov_limit,cor_cov/norm_cov,NA),pos_relative_to_tss=dplyr::if_else(strand=="+",pos-as.numeric(tss_data[5]),-(pos-as.numeric(tss_data[5])))) %>% dplyr::arrange(pos_relative_to_tss))
  }
  cl=parallel::makeCluster(threads)
  coverage_list=pbapply(X=tss_to_analyze,1,FUN=FUN,bin_path=bin_path,norm_log2=norm_log2,start=tss_start,end=tss_end,mean_cov=mean_cov,bam=bam,cov_limit=cov_limit,mapq=mapq,cl=cl)
  ## coverage_list=pbapply::pblapply(seq(1,nrow(tss_to_analyze),1),FUN=FUN,tss_data=tss_to_analyze,bin_path=bin_path,norm_log2=norm_log2,start=tss_start,end=tss_end,mean_cov=mean_cov,bam=bam,cov_limit=cov_limit,mapq=mapq,cl=cl)
  on.exit(parallel::stopCluster(cl))
  print(paste("TSS analyzed:",nrow(tss_to_analyze)))
  print(paste("TSS skipped:",nrow(ref_data)-nrow(tss_to_analyze)))


  return(coverage_list)
  }
