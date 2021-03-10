#' Generate seed from non-symbolic characters
#'
#' This function takes a string as input and transforms non-symbolic characters into a randomizer seed
#'
#' @param x A string
#' @return A number
#' @export



char2seed <- function(x,set=TRUE){

	tmp <- c(0:9,0:25,0:25)
	names(tmp) <- c(0:9,letters,LETTERS)

	x <- gsub("[^0-9a-zA-Z]","",as.character(x))

	xsplit <- tmp[ strsplit(x,'')[[1]] ]

	seed <- sum(rev( 7^(seq(along=xsplit)-1) ) * xsplit)
        seed <- as.integer( seed %% (2^31-1) )

	if(set){
		set.seed(seed)
		return(invisible(seed))
	} else {
		return(seed)
	}
}

#' Generates a bed file with temp suffix
#'
#'
#'
#' @param chr A vector with chromosome names
#' @param start A vector with starting position
#' @param end A vector with end position
#' @param strand A vector with strand information
#' @param name A string with file name
#' @return
#' @export

tmp_bed=function(chr="",start="",end="",strand="",name="File"){
	dat=data.frame(chr=chr,start=start,end=end,fill1=".",fill2=".",strand=strand)
	write.table(paste0(name,".bed.tmp"),quote=FALSE,col.names=FALSE,row.names=FALSE)
}
