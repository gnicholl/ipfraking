
#' Summary of Raked Weights
#'
#' Implements generic `summary` method for output of `raking` function.
#'
#' @export
summary.ipfraking = function(object,...) {
  print("this is a summary")

  # inflation table?

  # quantile table?
}


#' Histogram of raked weights
#'
#' Implements generic `hist` method for output of `raking` function.
#'
#' @export
hist.ipfraking = function(object,...) {

  # extract initial and raked weights
  wgt0 = object$initwgt
  wgt1 = object$weights

  # scale weights relative to their mean
  wgt0 = wgt0 / mean(wgt0)
  wgt1 = wgt1 / mean(wgt1)

  # plot
  xlim=c(min(wgt0,wgt1),max(wgt0,wgt1))
  tmp = hist(wgt1)
  final = hist( wgt0, col=rgb(0,0,1,0.1), xlim=xlim, breaks = tmp$breaks,
                main="Raked wgts (red) vs. init wgts (blue)",
                xlab="weight (relative to mean)")
  final = hist( wgt1, col=rgb(1,0,0,0.1), xlim=xlim, breaks = tmp$breaks, add=TRUE)

}










