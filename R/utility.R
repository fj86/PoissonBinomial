check.args.GPB <- function(x, probs, val_p, val_q, wts, method, log.p = FALSE){
  # check if 'x' contains only integers
  if(!is.null(x) && any(x - round(x) != 0)){
    warning("'x' should contain integers only! Using rounded off values.")
    x <- floor(x)
  }
  
  # check if 'probs' contains only probabilities
  if(is.null(probs) || any(is.na(probs) | probs < 0 | probs > 1))
    stop("'probs' must contain real numbers between 0 and 1!")
  
  # number of probabilities
  n <- length(probs)
  
  # check if 'val_p' and 'val_q' have the same length as 'probs'
  if(length(val_p) != n || length(val_q) != n) stop("'probs', 'val_p' and 'val_q' must have the same length!")
  
  if(!is.null(wts) && length(wts) != n)
    stop("'probs' and 'wts' (if not NULL) must have the same length!")
  
  # check if 'val_p' contains only integers
  if(!is.null(val_p) && any(val_p - round(val_p) != 0)){
    warning("'val_p' should contain integers only! Using rounded off values.")
    val_p <- floor(val_p)
  }
  
  # check if 'val_q' contains only integers
  if(!is.null(val_q) && any(val_q - round(val_q) != 0)){
    warning("'val_q' should contain integers only! Using rounded off values.")
    val_q <- floor(val_q)
  }
  
  # check if 'wts' contains only integers (zeros are allowed)
  if(!is.null(wts) && any(is.na(wts) | wts < 0 | abs(wts - round(wts)) > 1e-07))
    stop("'wts' must contain non-negative integers!")
  
  # make sure that the value of 'method' matches one of the implemented procedures
  method <- match.arg(method, c("DivideFFT", "Convolve", "Characteristic", "Normal", "RefinedNormal"))
  
  # if all checks were successful, return matched 'method'
  return(method)
}


transformGPB <- function(x, probs, val_p, val_q, wts){
  # number of probabilities
  n <- length(probs)
  
  ## expand 'probs', 'val_p' and 'val_q' according to the counts in 'wts'
  # if 'wts' is NULL, set it to be a vector of ones
  if(is.null(wts))
    wts <- rep(1, n)
  
  # expand 'probs', 'val_p', 'val_q'
  probs <- rep(probs, wts)
  val_p <- rep(val_p, wts)
  val_q <- rep(val_q, wts)
  
  # reorder 'val_p' and 'val_q' so that values in 'val_p' are always greater
  val_gr <- pmax(val_p, val_q)
  val_lo <- pmin(val_p, val_q)
  probs[val_gr > val_p] <- 1 - probs[val_gr > val_p]
  
  # re-compute length of 'probs' (= sum of 'wts')
  n <- sum(wts)
  
  ## determine relevant range of observations
  # determine minimum and maximum possible observations
  sum_min <- sum(val_lo)
  sum_max <- sum(val_gr)
  
  # which probabilities are 0 or 1, which val_p and val_q are equal
  idx.0 <- which(probs == 0)
  idx.1 <- which(probs == 1)
  idx.v <- which(val_gr == val_lo & probs > 0 & probs < 1)
  idx.r <- setdiff(1:n, union(union(idx.0, idx.1), idx.v))
  
  # guaranteed
  val_gr_sure <- val_gr[idx.1]
  val_lo_sure <- val_lo[idx.0]
  vals_equal <- val_gr[idx.v]# equal to val_lo[idx.v]
  sum_sure <- sum(val_gr_sure, val_lo_sure, vals_equal)
  
  # limit 'probs', 'val_p' and 'val_q' to relevant range
  np <- length(idx.r)
  if(np){
    probs <- probs[idx.r]
    val_gr <- val_gr[idx.r]
    val_lo <- val_lo[idx.r]
  }else{
    probs <- 1
    val_gr <- 0
    val_lo <- 0
  }
  
  # compute differences and their GCD
  diffs <- val_gr - val_lo
  
  # bounds of relevant observations
  sum_min_in <- sum(val_lo) + sum_sure
  sum_max_in <- sum(val_gr) + sum_sure
  
  return(list(probs = probs, val_p = val_gr, val_q = val_lo, compl.range = sum_min:sum_max, inner.range = sum_min_in:sum_max_in, inner.size = sum_max_in - sum_min_in + 1, n = np, diffs = diffs))
}
