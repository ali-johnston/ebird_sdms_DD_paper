sci_notation <- function(x, digits = 1) {
  if (length(x) > 1) {
    return(append(sci_notation(x[1]), sci_notation(x[-1])))
  }
  if (is.na(x)) {
    return(0)
  } else if (x == 0) {
    return(NA)
  } else if (x %in% c(1, 10)) {
    return(x)
  }
  exponent <- floor(log10(x))
  base <- round(x / 10^exponent, digits)
  if (base == 1) {
    e <- as.expression(substitute(10^exponent, list(exponent = exponent)))
  } else {
    e <- as.expression(substitute(base%*%10^exponent, 
                                  list(base = base, exponent = exponent)))
  }
  return(e) 
}