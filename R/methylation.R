#' Calculate Mean Methylation Ratio around TFBS
#'
#' This function takes a DATA.FRAME with TFBS positions, a BAM file with the sequence to analyze
#' and generates the values of mean Methylation Ratios per base around TFBS which then outputs into
#' a TSS file. For better understanding check: https://www.nature.com/articles/s41467-019-12714-4
#'
#'
#' @param bin_path Path to binary. Default tools/samtools/samtools
#' @param ref_data Data from BED file as Data.frame
#' @param bam Path to BAM file.
#' @param ref_genome Path to reference genome FA file.
#' @param sample_name Sample name
#' @param tf_name Transcription Factor name
#' @param tfbs_start Number of bases to analyze forward from TFBS central point. Default 1000
#' @param tfbs_end Number of bases to analyze  backward from TFBS central point. Default 1000
#' @param mapq Min quality of mapping reads. Default 10
#' @param phred Min phred base quality. Default 5
#' @param keep_strand Use strand information from BED files if available. Default TRUE.
#' @param bin_width Width of the the bins in which to group methylation data. Default 50.
#' @param verbose Enables progress messages. Default FALSE
#' @param output_dir Directory to output results. If not provided then outputs in current directory
#' @return A list of DATA.FRAMEs with the methylation ratio per base and per bin_width with default bin_width 50
#' @export
#' @import pbapply

calculate_MR_tfbs=function(bin_path="tools/PileOMeth/output/MethylDackel",ref_data="",bam="",sample_name="",tf_name="",ref_genome="",tfbs_start=1000,tfbs_end=1000,mapq=10,phred=5,output_dir="",keep_strand=TRUE,bin_width=50,verbose=FALSE){

	sep="/"
  if(output_dir==""){
    sep=""
  }

	out_file=paste0(output_dir,sep,sample_name,"_",tf_name)

  if (!dir.exists(out_file)){
      dir.create(out_file)
    }

  ref_data=data.frame(chr=ref_data[,1],start=ref_data[,2],end=ref_data[,3],strand=ref_data[,which(ref_data[1,]=="+" | ref_data[1,]=="-")],pos=as.integer((ref_data[,3]+ref_data[,2])/2))
	numb_tfbs=nrow(ref_data)
	tfbs_to_analyze=ref_data

  ## Filter for overlapping TFBS or duplicated TFBS to save time

  tfbs_to_analyze=tfbs_to_analyze %>% dplyr::filter(!grepl("_",chr),(pos-start)<1) %>% dplyr::distinct(chr, pos, .keep_all = TRUE)


	tfbs_to_analyze=tfbs_to_analyze %>% dplyr::mutate(start=pos-tfbs_start,end=pos+tfbs_end)
	ref_data=ref_data %>% dplyr::rowwise() %>% dplyr::do(data.frame(chr = .$chr,pos_relative_to_tfbs =seq(-tfbs_start,tfbs_end,by=1), pos = seq(.$start-tfbs_start, .$start+tfbs_end, by = 1)))

	## Create a temporary bed file with regions +-1000 pb from TFBS to feed to MethylDackel

  tmp_bed(chr=tfbs_to_analyze$chr,start=tfbs_to_analyze$start,end=tfbs_to_analyze$end,strand=tfbs_to_analyze$strand,name=paste0(output_dir,sep,sample_name,"_",tf_name))

	strand=""
	if (keep_strand){
		strand="--keepStrand"
	}
	if (verbose){
		print(paste(bin_path," extract ",ref_genome, bam,"-l ",paste0(output_dir,sep,sample_name,"_",tf_name,".bed.tmp")," -o ",out_file,strand," -q ",mapq," -p ",phred))
	}
	system(paste(bin_path," extract ",ref_genome, bam,"-l ",paste0(output_dir,sep,sample_name,"_",tf_name,".bed.tmp")," -o ",out_file,strand," -q ",mapq," -p ",phred))
	system(paste0("rm ",paste0(output_dir,sep,sample_name,"_",tf_name,".bed.tmp")))
	tfbs=read.table(paste0(out_file,"_CpG.bedGraph"),skip=1)
	tfbs$pos=as.integer((tfbs$V2+tfbs$V3)/2)
	colnames(tfbs)=c("chr","start","end","MR","nC","nT","pos")

	# Generate per base mean methylation data across all TFBS

	merg_tfbs1=dplyr::left_join(ref_data,tfbs,by=c("chr","pos"))%>% dplyr::group_by(pos_relative_to_tfbs) %>% dplyr::mutate( x_bins = ifelse(is.na(cut(pos_relative_to_tfbs, breaks = seq(-tfbs_end,tfbs_start,1),include.lowest=FALSE,labels=FALSE)),0,cut(pos_relative_to_tfbs, breaks = seq(-tfbs_end,tfbs_start,1),include.lowest=FALSE,labels=FALSE)))%>% dplyr::group_by(x_bins) %>%
	dplyr::mutate(x_bins=as.integer((max(pos_relative_to_tfbs)+min(pos_relative_to_tfbs))/2)) %>%
	dplyr::summarise(MEAN_MR=mean(MR,na.rm=TRUE),CI=qt(0.95,(sum(!is.na(MR))-1))*sd(MR,na.rm=TRUE)/sqrt(sum(!is.na(MR))),DATA_POINTS_ANALYZED=sum(!is.na(MR)))
	merg_tfbs1= dplyr::rename(merg_tfbs1,POSITION_RELATIVE_TO_TFBS=x_bins) %>% dplyr::mutate(CI95_UPPER_BOUND=MEAN_MR+CI,CI95_LOWER_BOUND=MEAN_MR-CI,TFBS_ANALYZED=nrow(tfbs_to_analyze),BIN_WIDTH=1)

	# Generate per bin_width mean methylation data across all TFBS. Default bin_width 50. This is done because
	# methylation ratio is scarcely distributed across all TFBS, so even though we analyze 1000 TFBS not all of them return methylation info.

	merg_tfbs2=dplyr::left_join(ref_data,tfbs,by=c("chr","pos"))%>% dplyr::group_by(pos_relative_to_tfbs) %>% dplyr::mutate( x_bins = ifelse(is.na(cut(pos_relative_to_tfbs, breaks = seq(-tfbs_end,tfbs_start,bin_width),include.lowest=FALSE,labels=FALSE)),0,cut(pos_relative_to_tfbs, breaks = seq(-tfbs_end,tfbs_start,bin_width),include.lowest=FALSE,labels=FALSE)))%>% dplyr::group_by(x_bins) %>%
	dplyr::mutate(x_bins=as.integer((max(pos_relative_to_tfbs)+min(pos_relative_to_tfbs))/2)) %>%
	dplyr::summarise(MEAN_MR=mean(MR,na.rm=TRUE),CI=qt(0.95,(sum(!is.na(MR))-1))*sd(MR,na.rm=TRUE)/sqrt(sum(!is.na(MR))),DATA_POINTS_ANALYZED=sum(!is.na(MR)))
	merg_tfbs2= dplyr::rename(merg_tfbs2,POSITION_RELATIVE_TO_TFBS=x_bins) %>% dplyr::mutate(CI95_UPPER_BOUND=MEAN_MR+CI,CI95_LOWER_BOUND=MEAN_MR-CI,TFBS_ANALYZED=nrow(tfbs_to_analyze),BIN_WIDTH=bin_width)

  print(paste("TFBS analyzed:",nrow(tfbs_to_analyze)))
  print(paste("TFBS skipped:",nrow(numb_tfbs)-nrow(tfbs_to_analyze)))


	return(list(merg_tfbs1,merg_tfbs2))
  }
