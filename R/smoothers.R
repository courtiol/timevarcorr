#' Smoothing by kernel regression
#'
#' The function perform the smoothing of a time-series by non-parametric kernel regression.
#'
#' The function is essentially a wrapper that calls different underlying
#' functions depending on the kernel that is selected:
#' - [`lpridge::lpepa`] for "epanechnikov".
#' - [`stats::ksmooth`] for "normal" and "box".
#' The argument `param_smoother` can be used to pass additional arguments to
#' these functions.
#'
#' @param x a numeric vector of the series to be smoothed.
#' @param t a (numeric or Date) vector of time points. If missing, observations
#'   are considered to correspond to sequential time steps (i.e. 1, 2 ...).
#' @param h a scalar indicating the bandwidth used by the smoothing function.
#' @param t.for.pred a (numeric or Date) vector of time points at which to
#'   evaluate the smoothed fit. If missing, `t` is used.
#' @param kernel a character string indicating which kernel to use: "epanechnikov"
#'    (the default), "box", or "normal" (abbreviations also work).
#' @param param_smoother a list of additional parameters to provide to the
#'   internal smoothing function (see **Details**).
#'
#' @return a dataframe of time points (`t.for.pred`) and corresponding fitted
#'   values.
#'
#' @seealso [`tcor`]
#'
#' @export
#'
#' @references
#' A short post I found useful: \url{http://users.stat.umn.edu/~helwig/notes/smooth-notes.html}
#'
#' @examples
#'
#' ## Smooth 10 first values of a vector
#'
#' kern_smooth(stockprice$DAX[1:20], h = 5)
#'
#'
#' ## Prediction at time step 2 and 3
#'
#' kern_smooth(stockprice$DAX, h = 1, t.for.pred = c(2, 3))
#'
#'
#' ## Smoothing using a vector of dates for time
#'
#' kern_smooth(x = stockprice$DAX[1:10], t = stockprice$DateID[1:10], h = 5)
#'
#'
#' ## Smoothing conserves original order
#'
#' kern_smooth(x = stockprice$DAX[10:1], t = stockprice$DateID[10:1], h = 5)
#'
#'
#' ## Effect of the bandwidth
#'
#' plot(stockprice$DAX[1:100] ~ stockprice$DateID[1:100],
#'      las = 1, ylab = "DAX index", xlab = "Date")
#' points(kern_smooth(stockprice$DAX[1:100], stockprice$DateID[1:100], h = 1),
#'        type = "l", col = "grey")
#' points(kern_smooth(stockprice$DAX[1:100], stockprice$DateID[1:100], h = 3),
#'        type = "l", col = "blue")
#' points(kern_smooth(stockprice$DAX[1:100], stockprice$DateID[1:100], h = 10),
#'        type = "l", col = "red")
#' legend("topright", fill = c("grey", "blue", "red"),
#'        legend = c("1", "3", "10"), bty = "n", title = "Bandwidth (h)")
#'
#'
#' ## Effect of the kernel
#'
#' plot(stockprice$DAX[1:100] ~ stockprice$DateID[1:100],
#'      las = 1, ylab = "DAX index", xlab = "Date")
#' points(kern_smooth(stockprice$DAX[1:100], stockprice$DateID[1:100], h = 10),
#'        type = "l", col = "orange")
#' points(kern_smooth(stockprice$DAX[1:100], stockprice$DateID[1:100], h = 10, kernel = "box"),
#'        type = "l", col = "blue")
#' points(kern_smooth(stockprice$DAX[1:100], stockprice$DateID[1:100], h = 10, kernel = "norm"),
#'        type = "l", col = "red")
#' legend("topright", fill = c("orange", "blue", "red"),
#'        legend = c("epanechnikov", "box", "normal"), bty = "n", title = "Kernel method")
#'
kern_smooth <- function(x, t = seq_along(x), h, t.for.pred = t,
                        kernel = c("epanechnikov", "box", "normal"), param_smoother = list()) {

  ## method selection
  kernel <- match.arg(kernel)

  ## handling missing values
  xx <- x[!is.na(x)]
  tt <- t[!is.na(x)]

  ## smoothing
  if (kernel %in% c("normal", "box")) {

    param_list <- c(param_smoother, list(x = tt, y = xx, kernel = kernel, bandwidth = h, x.points = t.for.pred))
    res <- do.call(stats::ksmooth, param_list)
    original_order <- order(t.for.pred)
    return(data.frame(t = res$x[original_order], x = res$y[original_order]))

  } else if (kernel == "epanechnikov") {

    param_list <- c(param_smoother, list(x = tt, y = xx, bandwidth = h, x.out = t.for.pred))
    res <- do.call(lpridge::lpepa, param_list)
    original_order <- order(t.for.pred)
    new.t <- res$x.out[original_order]
    class(new.t) <- class(t) # restore class, otherwise lost
    return(data.frame(t = new.t, x = res$est[original_order]))
  }

  stop("Argument `kernel` unknown.")
}
