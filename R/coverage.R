#' Calculate Genome-wide coverage for a BAM file
#'
#' This function takes a BAM file and generates a TXT file with its genome-wide coverage
#'
#' @param bin_path Path to binary. Default tools/bedtools2/bin/bedtools
#' @param bin_path2 Path to binary. Default tools/samtools/samtools
#' @param bam Path to BAM file.
#' @param verbose Enables progress messages. Default FALSE.
#' @param bed Path to BED with regions to filter.
#' @param genome Genome file FA format. Not required. Only to estimate coverage for specific regions.
#' @param output_dir Directory to output results. If not provided then outputs in current directory
#' @export

calculate_genowide_coverage=function(bin_path="tools/bedtools2/bin/bedtools",bin_path2="tools/samtools/samtools",bam="",verbose=FALSE,output_dir="",genome="",bed=""){

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

  if (genome!=""){
    genome=paste("-g",genome)
  }

  if (bed==""){
    if(verbose){
      print(paste(bin_path,"genomecov -ibam",bam,genome,">",out_file))
    }
    system(paste(bin_path,"genomecov -ibam",bam,genome,">",out_file))
  }else{
    if(verbose){
      print(paste(bin_path2," view -b -L",bed, bam,"|",bin_path,"genomecov -ibam stdin",genome,">",out_file))
    }
    system(paste(bin_path2," view -b -L",bed, bam,"|",bin_path,"genomecov -ibam stdin",genome,">",out_file))
  }
}

#' Get normalized local coverage from segmentation data (CNA) for a single base position
#'
#' This function takes normalized data for local coverage as input,
#' as well as a valid genomic position for a single base and returns a log2 corrected value for said position
#'
#' @param pos Genomic position
#' @param chr Chromosome
#' @param norm_log2 DATA.FRAME with normalized local coverage
#' @return A double with the local coverage
#' @export


get_norm_local_coverage=function(pos="",chr="",norm_log2=""){
   norm_log2$log2=ifelse(is.na(norm_log2$log2),0,norm_log2$log2)
   norm=as.numeric(norm_log2[norm_log2$chr==chr & pos>norm_log2$start & pos<norm_log2$end,]$log2)
   if(!length(norm)==0){
     norm=2**norm
     if (norm==0){
       return(0.001)
    }
     else{
       return(norm)
     }
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
#' @return A DATA.FRAME with the mean coverage values and CI
#' @export

get_mean_and_conf_intervals=function(cov_data="",CI=0.95){
  FUN=function(x,cov_data){
    dat=cov_data[[x]]
    return(dat[,"norm_cor_cov",drop=FALSE])
  }
  norm_cor_cov_list=lapply(seq(1:length(cov_data)),FUN=FUN,cov_data=cov_data)
  norm_cor_cov_df=suppressMessages(norm_cor_cov_list %>% dplyr::bind_cols())
  norm_cor_cov_means=rowMeans(norm_cor_cov_df,na.rm = TRUE)
  numb_analyz_tfbs=norm_cor_cov_df %>% is.na() %>% `!` %>% rowSums()
  if (all(numb_analyz_tfbs>3)){
      norm_cor_cov_mes=apply(norm_cor_cov_df,1,FUN=function(x){qt(CI,sum(!is.na(x))-1)*sd(x[!is.na(x)])/sqrt(sum(!is.na(x)))})
        results=data.frame(POSITION_RELATIVE_TO_TFBS=cov_data[[1]]$pos_relative_to_tfbs,MEAN_DEPTH=norm_cor_cov_means,CI95_LOWER_BOUND=norm_cor_cov_means-norm_cor_cov_mes,CI95_UPPER_BOUND=norm_cor_cov_means+norm_cor_cov_mes,TFBS_ANALYZED=numb_analyz_tfbs)
  }else{
  results=data.frame(POSITION_RELATIVE_TO_TFBS=cov_data[[1]]$pos_relative_to_tfbs,MEAN_DEPTH=norm_cor_cov_means,TFBS_ANALYZED=numb_analyz_tfbs)
  }
  return (results)
}



#' Get coverage for target regions
#'
#' This function takes a  BAM and a BED file with target regions and estimate per base coverage for them
#' using samtools depth
#'
#' @param bin_path Path to samtools. Default tools/samtools/samtools
#' @param bam Path to BAM file
#' @param bed Path to BED file
#' @param mapq Minimum mapping quality
#' @param threads Number of threads to use. Default 3
#' @param output_dir Directory to output results
#' @return A double with the mean genome coverage
#' @export

get_coverage_bed=function(bin_path="tools/samtools/samtools",bam="",bed="",mapq=0,threads=3,output_dir=""){
  region_data=read.table(bed,comment.char="$")
  region_data=region_data[,c(1,2,3,4)]
  colnames(region_data)=c("chr","start","end")
  cov_results=parallel::mclapply(1:nrow(region_data),FUN=function(x){
    cov_data=read.csv(text=system(paste(bin_path,"depth -a -Q",mapq, "-r",paste0(region_data[x,1],":",region_data[x,2],"-",region_data[x,3]),bam),intern=TRUE),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses=c("character"));
    cov_data$Original=paste0(region_data[x,1],":",region_data[x,2],"-",region_data[x,3]);
    cov_data$Id=region_data[x,4];
    cov_data$Position=as.integer((as.numeric(region_data[x,2])+as.numeric(region_data[x,3]))/2);
    cov_data$Relative_Position=cov_data$Position-as.numeric(cov_data$V3);
    return(cov_data)
  },mc.cores=threads)
  cov_results=dplyr::bind_rows(cov_results)
  return(cov_results)
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
#' @return A double with the mean genome coverage
#' @export


get_mean_coverage=function(file="",output_dir="",region="genome",sample_name="",save=TRUE){
  data=read.table(file)
  data=data %>% dplyr::mutate(tot_bases=(as.numeric(V2)*as.numeric(V3))) %>% dplyr::group_by(V1) %>% dplyr::summarise(cnt=dplyr::n(),avg_cov=sum(tot_bases)/max(V4), .groups = 'drop') %>% as.data.frame()

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
#' local coverage values for CNA, and outputs the corrected depth coverage values around TFBS and generates
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
#' @param tfbs_start Number of bases to analyze forward from TFBS central point. Default 1000
#' @param tfbs_end Number of bases to analyze  backward from TFBS central point. Default 1000
#' @param cov_limit Max base depth. Default 1000
#' @param max_regions Max number of TFBS to analyze. Default 100000
#' @param mapq Min quality of mapping reads. Default 0
#' @param threads Number of threads. Default 1
#' @param output_dir Directory to output results. If not provided then outputs in current directory
#' @return A list of DATA.FRAME with the coverage per base of each TFBS
#' @export
#' @import pbapply

calculate_coverage_tfbs=function(bin_path="tools/samtools/samtools",ref_data="",bam="",sample_name="",tf_name="",mean_cov="",norm_log2="",tfbs_start=1000,tfbs_end=1000,cov_limit=1000,mapq=0,threads=1,output_dir=""){

  ref_data=data.frame(chr=ref_data[,1],start=ref_data[,2],end=ref_data[,3],strand=ref_data[,which(ref_data[1,]=="+" | ref_data[1,]=="-")],pos=as.integer((ref_data[,3]+ref_data[,2])/2))
  tfbs_to_analyze=ref_data

  ## Filter for overlapping TFBS or duplicated TFBS to save time
  tfbs_to_analyze=tfbs_to_analyze %>% dplyr::filter(!grepl("_",chr),(pos-start)<1) %>% dplyr::distinct(chr, pos, .keep_all = TRUE)

  FUN=function(x,bin_path,norm_log2,start,end,mean_cov,bam,cov_limit,mapq){
    tfbs_data=t(x)

    ## IF running sambamba
    ## cov_data=read.csv(text=system(paste(bin_path,"depth base -t 3 -z -L ",paste0(tfbs_data$chr,":",as.numeric(tfbs_data$pos)-start,"-",as.numeric(tfbs_data$pos)+end),bam),intern=TRUE),header=TRUE,sep="\t")
    ## cov_data=cbind(cov_data[,1:3],strand=tfbs_data$strand)
    ## names(cov_data)=c("chr","pos","cov","strand")

    cov_data=read.csv(text=system(paste(bin_path,"depth -a -Q",mapq, "-r",paste0(tfbs_data[1],":",as.numeric(tfbs_data[5])-start,"-",as.numeric(tfbs_data[5])+end),bam),intern=TRUE),header=FALSE,sep="\t")
    colnames(cov_data)=c("chr","pos","cov")
    cov_data$chr=as.character(cov_data$chr)
    norm_cov=get_norm_local_coverage(pos=as.numeric(tfbs_data[5]),chr=tfbs_data[1],norm_log2=norm_log2)

    if(!(nrow(cov_data)==(start+end+1))){
      fix=(as.numeric(tfbs_data[5])-start):(as.numeric(tfbs_data[5])+end)
      fix=data.frame(chr=tfbs_data[1],pos=fix[!fix==cov_data$pos])
      cov_data=dplyr::bind_rows(cov_data,fix) %>% dplyr::arrange(pos)
    }

    cov_data=cbind(cov_data,strand=tfbs_data[4])
    cov_data=cov_data %>% dplyr::mutate(cor_cov=as.numeric(cov)/as.numeric(mean_cov))  %>% dplyr::mutate(norm_cor_cov=cor_cov/as.numeric(norm_cov),pos_relative_to_tfbs=dplyr::if_else(strand=="+",pos-as.numeric(tfbs_data[5]),-(pos-as.numeric(tfbs_data[5])))) %>% dplyr::arrange(pos_relative_to_tfbs)

    return(cov_data)
  }
  cl=parallel::makeCluster(threads)
  coverage_list=pbapply(X=tfbs_to_analyze,1,FUN=FUN,bin_path=bin_path,norm_log2=norm_log2,start=tfbs_start,end=tfbs_end,mean_cov=mean_cov,bam=bam,cov_limit=cov_limit,mapq=mapq,cl=cl)
  ## coverage_list=pbapply::pblapply(seq(1,nrow(tfbs_to_analyze),1),FUN=FUN,tfbs_data=tfbs_to_analyze,bin_path=bin_path,norm_log2=norm_log2,start=tfbs_start,end=tfbs_end,mean_cov=mean_cov,bam=bam,cov_limit=cov_limit,mapq=mapq,cl=cl)

  on.exit(parallel::stopCluster(cl))
  print(paste("TFBS analyzed:",nrow(tfbs_to_analyze)))
  print(paste("TFBS skipped:",nrow(ref_data)-nrow(tfbs_to_analyze)))


  return(coverage_list)
  }

#' Calculate Coverage around Genomic Position
#'
#' This function takes chromosome location, as well as strand information, and generates a coverage report around
#' for it. This coverage can be normalized if local segmentation values for the region are given and further corrected by the
#' the local coverage in the BAM
#'
#'
#' @param bin_path Path to binary. Default tools/samtools/samtools
#' @param chr Chromosome to analyze.
#' @param bam Path to BAM file.
#' @param position Genomic position within the chromosome
#' @param strand Strand direction. If not available asumes + strand.
#' @param norm_log2 Local segment log2
#' @param start Downstream distance from genomic position to estimate coverage
#' @param end Upstream distance from genomic position to estimate coverage
#' @param start_bin Downstream distance from genomic position to limit bin size of central region. Only if binned mode is selected.
#' @param end_bin Upstream distance from genomic position to limit bin size of central region. Only if binned mode is selected.
#' @param mean_cov Mean genome wide coverage.
#' @param method Method to use to estimate coverage. default/binned
#' @param mapq Min quality of mapping reads. Default 0
#' @return A DATA.FRAME with per base coverage por each genomic position
#' @import tidyverse
#' @export

calculate_coverage_around_gp=function(bin_path="tools/samtools/samtools",chr="",position="",strand="",bam="",norm_log2=1,start=1000,end=1000,mean_cov=1,mapq=0,method="default",start_bin=75,end_bin=75){
    if (method=="default"){
      cov_data=read.csv(text=system(paste(bin_path,"depth -a -Q",mapq, "-r",paste0(chr,":",as.numeric(position)-start,"-",as.numeric(position)+end),bam),intern=TRUE),header=FALSE,sep="\t",stringsAsFactors=FALSE)
      colnames(cov_data)=c("chr","pos","cov")
      cov_data$strand=strand
      cov_data=cov_data %>% dplyr::mutate(cor_cov=as.numeric(cov)/sum(as.numeric(cov)))  %>% dplyr::mutate(norm_cor_cov=cor_cov/as.numeric(norm_log2),pos_relative_to_tfbs=dplyr::if_else(strand=="+"|strand=="",pos-as.numeric(position),-(pos-as.numeric(position)))) %>% dplyr::arrange(pos_relative_to_tfbs)
      return(cov_data)

    }else if (method=="binned"){
      cov_data_central=read.csv(text=system(paste(bin_path,"depth -a -Q",mapq, "-r",paste0(chr,":",as.numeric(position)-start_bin,"-",as.numeric(position)+end_bin),bam),intern=TRUE),header=FALSE,sep="\t",stringsAsFactors=FALSE)
      colnames(cov_data_central)=c("chr","pos","cov")
      cov_data_central$bin="CENTRAL"
      cov_data_central$bin_pos=paste0(chr,":",as.numeric(position)-start_bin,"-",as.numeric(position)+end_bin)
      cov_data=cov_data_central %>% dplyr::group_by(bin) %>% dplyr::summarise(cov=cov)
      cov_data_left_flank=read.csv(text=system(paste(bin_path,"depth -a -Q",mapq, "-r",paste0(chr,":",as.numeric(position)-start,"-",as.numeric(position)-start_bin),bam),intern=TRUE),header=FALSE,sep="\t",stringsAsFactors=FALSE)
      colnames(cov_data_left_flank)=c("chr","pos","cov")
      cov_data_left_flank$bin="LEFT_FLANK"
      cov_data_left_flank$bin_pos=paste0(chr,":",as.numeric(position)-start,"-",as.numeric(position)-start_bin)
      cov_data_left_flank=cov_data_left_flank %>% dplyr::group_by(bin,bin_pos) %>% dplyr::summarise(cov=cov)
      cov_data_right_flank=read.csv(text=system(paste(bin_path,"depth -a -Q",mapq, "-r",paste0(chr,":",as.numeric(position)+end_bin,"-",as.numeric(position)+end),bam),intern=TRUE),header=FALSE,sep="\t",stringsAsFactors=FALSE)
      colnames(cov_data_right_flank)=c("chr","pos","cov")
      cov_data_right_flank$bin="RIGHT_FLANK"
      cov_data_right_flank$bin_pos=paste0(chr,":",as.numeric(position)+end_bin,"-",as.numeric(position)+end)
      cov_data_right_flank=cov_data_right_flank %>% dplyr::group_by(bin,bin_pos) %>% dplyr::summarise(cov=cov)
      cov_data=rbind(cov_data_central,cov_data_left_flank,cov_data_right_flank)
      cov_data=cov_data %>% dplyr::mutate(cor_cov=as.numeric(cov)/mean_cov)  %>% dplyr::mutate(norm_cor_cov=cor_cov/as.numeric(norm_log2),bin=dplyr::if_else(strand=="+"|strand=="",bin,dplyr::if_else(bin=="RIGHT_FLANK","LEFT_FLANK","RIGHT_FLANK"))) %>% dplyr::arrange(bin)
      return(cov_data)
    }
}
