Mode <- function(x) {
  #if ( length(x) <= 2 ) return(x[1])
  #if ( anyNA(x) ) x = x[!is.na(x)]
  #ux <- unique(x)
  #ux[which.max(tabulate(match(x, ux)))]
  
  
  blah2=density(x)
  blah2$x[which.max(blah2$y)]
}