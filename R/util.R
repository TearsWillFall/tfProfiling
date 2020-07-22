#' Generate seed from non-symbolic characters
#'
#' This function takes a string as input and transforms non-symbolic characters into a number
#'
#' @param x A string
#' @return A number
#' @export



char2seed <- function(x){

	tmp <- c(0:9,0:25,0:25)
	names(tmp) <- c(0:9,letters,LETTERS)

	x <- gsub("[^0-9a-zA-Z]","",as.character(x))

	xsplit <- tmp[ strsplit(x,'')[[1]] ]

	seed <- sum(rev( 7^(seq(along=xsplit)-1) ) * xsplit)
        seed <- as.integer( seed %% (2^31-1) )

	return(seed)

}
