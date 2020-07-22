#' Filters low signal frequency in TFBS coverage data
#'
#' This function takes a list with the coverage data and normalizes and filters it leaving only the low frequency signal
#'
#' @param coverage_signal A list with coverage data
#' @param m Filter length
#' @return A list with low frequency signal values
#' @export


get_low_signal <- function(coverage_signal,m=1001) {
   normalized <- coverage_signal / mean(coverage_signal)
   low<-signal::sgolayfilt(normalized,3,m)
   return(low)
}


#' Filters high signal frequency in TFBS coverage data
#'
#' This function takes a list with the coverage data and normalizes and filters it leaving only the high frequency signal
#'
#' @param coverage_signal A list with coverage data
#' @param low_signal A list with the low frequency signal values
#' @param m Filter length
#' @return A list with high frequency signal values
#' @export


get_high_signal <- function(coverage_signal,low_signal,m=51) {
   normalized <- coverage_signal / mean(coverage_signal)
   high<-signal::sgolayfilt(normalized,3,m)
   high_adjusted<-(high/low_signal)
   return(high_adjusted)
}

#' Calls maximum and minimum peaks
#'
#' This function takes a list with a series of values and calls the minimum and maximum
#'
#' @param x A list with values
#' @param m Filter length
#' @return A list with peak values
#' @export




find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
       z <- i - m + 1
       z <- ifelse(z > 0, z, 1)
       w <- i + m + 1
       w <- ifelse(w < length(x), w, length(x))
       if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
     pks <- unlist(pks)
     pks
}
