#' @describeIn tcor Internal function selecting the optimal bandwidth parameter `h`.
#'
#' See **Bandwidth selection** for details.
#'
#' @order 3
#'
#'
#' @export
#'
select_h <- function(x, y, t = seq_along(x), cor.method = c("pearson", "spearman"),
                     kernel = c("epanechnikov", "box", "normal"), param_smoother = list(),
                     verbose = FALSE) {

  CV_error <- NA # default value for export if not computed

  ## check nb.cores input
  if (is.null(options("mc.cores")[[1]])) {
    nb.cores <- 1
  } else {
    nb.cores <- options("mc.cores")[[1]]
  }

  if (Sys.info()[['sysname']] != "Windows" && parallel::detectCores() < nb.cores) {
    stop(paste("\nYour computer does not allow such a large value for the argument `nb.cores`. The maximum value you may consider is", parallel::detectCores(), ".\n"))
  }

  if (Sys.info()[['sysname']] != "Windows" && parallel::detectCores() > 1 && nb.cores == 1) {
    message("\nYou may set `nb.cores` to a number higher than 1 for faster computation.\n")
  }

  if (Sys.info()[['sysname']] == "Windows" && nb.cores > 1) {
    message("\nThe argument `nb.cores` is ignore on Windows-based infrastructure.\n")
  }

  time <- system.time({

    h_selection <- "LOO-CV"

    if (length(x) > 500) message("Bandwidth selection using LOO-CV... (may take a while)")

    ## define function for performing leave-one-out cross-validation ## TODO: extract code into standalone function
    CV <- function(h) {
      CVi <- mclapply(seq_along(t), function(oob) {
        obj <- calc_rho(x = x[-oob], y = y[-oob], t = t[-oob], h = h, t.for.pred = t[oob],
                        cor.method = cor.method, kernel = kernel, param_smoother = param_smoother)
        (obj$rho_smoothed - ((x[oob] - obj$x_smoothed) * (y[oob] - obj$y_smoothed)) / (obj$sd_x_smoothed * obj$sd_y_smoothed))^2
      })
      CVi_num <- as.numeric(CVi)
      failing <- is.na(CVi_num)
      res <- mean(CVi_num[!failing])
      if (verbose) {
        print(paste("h =", round(h, digits = 1L), "    # failing =", sum(failing),"    CV =", res))
      }
      res
    }

    ## estimate h by cross-validation
    h_max <- 8*sqrt(length(x))
    opt <- stats::optimize(CV, interval = c(1, h_max), tol = 1) ## TODO: check if other optimizer would be best

    ## if h_max is best, use elbow criterion instead
    h <- opt$minimum
    CV_error <- opt$objective
    print(paste("h selected using LOO-CV =", round(h, digits = 1L)))
    CV_bound <- CV(h_max)

    if (CV_bound <= opt$objective) {
      CV_error <- NA # remove stored value since not used in this case
      h_selection <- "elbow criterion"

      if (length(x) > 500) message("Bandwidth selection using elbow criterion... (may take a while)")

      RMSE_Lh0 <- function(h0) {
        d <- data.frame(h = 1:h0)
        d$CV_h <- sapply(h, CV)
        fit_Lh0 <- stats::lm(CV_h ~ h, data = d)
        sqrt(mean(stats::residuals(fit_Lh0)^2))
      }

      RMSE_Rh0 <- function(h0) {
        d <- data.frame(h = (h0 + 1):h_max)
        d$CV_h <- sapply(h, CV)
        fit_Rh0 <- stats::lm(CV_h ~ h, data = d)
        sqrt(mean(stats::residuals(fit_Rh0)^2))
      }

      RMSE_h0 <- function(h0) {
        (h0 - 1)/(h_max - 1)*RMSE_Lh0(h0) + (h_max - h0)/(h_max - 1)*RMSE_Rh0(h0)
      }

      opt <- stats::optimize(RMSE_h0, interval = c(3, h_max - 1), tol = 1)
      h <- opt$minimum
      print(paste("h selected using elbow criterion =", round(h, digits = 1L)))
    }
  })
  list(h = h, h_selection = h_selection, CV_error = CV_error, time = time)
}
