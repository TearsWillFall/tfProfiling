
get_low_signal <- function(coverage_signal) {
   normalized <- coverage_signal / mean(coverage_signal)
   low<-sgolayfilt(normalized,3,1001)
   return(low)
}

getHighSignal <- function(coverage_signal,low_signal) {
   normalized <- coverage_signal / mean(coverage_signal)
   high<-sgolayfilt(normalized,3,51)
   high_adjusted<-(high/low_signal)
   return(high_adjusted)
}

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
