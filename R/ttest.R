#' Compute equality test between correlation coefficient estimates at two time points
#'
#' This function tests whether smoothed correlation values at two time points are equal (H0) or not.
#' The test is described page 341 in Choi & Shin, 2021.
#'
#' Two different test statistics can be used, one is asymptotically Student-t distributed under H0 and one is chi-square distributed.
#' In practice, it seems to give very similar results.
#'
#' @param tcor_obj the output of a call to [`tcor`] with `CI = TRUE`.
#' @param t1 the first time point used by the test (by default, the first time point in the time series).
#' @param t2 the second time point used by the test (by default, the last time point in the time series).
#' @param test a character string indicating which test to use ("student", the default; or "chi2").
#'
#' @return a data.frame with the result of the test, including the effect size (`delta_r = r[t2] - r[t1]`).
#'
#' @export
#'
#' @seealso [`tcor`]
#'
#' @examples
#' res <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = 50, CI = TRUE))
#' equality_test(res)
#'
#' ## Chi2 instead of Student's t-test
#'
#' equality_test(res, test = "chi2")
#'
#'
#' ## time point can be dates or indices (mixing possible)
#'
#' equality_test(res, t1 = "2000-04-04", t2 = 1000)
#'
equality_test <- function(tcor_obj, t1 = 1, t2 = nrow(tcor_obj), test = c("student", "chi2")) {

  ## check and prepare inputs
  test <- match.arg(test)

  if (!attr(tcor_obj, "CI")) {
    stop("For running an equality test, you must produce `tcor_obj` using the argument `CI = TRUE`.")
  }

  if (attr(tcor_obj, "h_selection") == "fixed by user") {
    message("The bandwidth stored in `tcor_obj` and used by `equality_test()` has not been automatically selected but was set by the user. Rerun `tcor()` without specifying `h` if you want otherwise.")
  }

  ## turn non default date into index
  if (t1 != 1 || t2 != nrow(tcor_obj)) {
    if (!is.numeric(t1)) {
      t1 <- which(tcor_obj$t == t1)
      if (length(t1) == 0) stop("`t1` not found in `tcor_obj`.")
    }
    if (!is.numeric(t2)) {
      t2 <- which(tcor_obj$t == t2)
      if (length(t2) == 0) stop("`t2` not found in `tcor_obj`.")
    }
  }

  ## compute t-test statistics
  delta_r <- diff(tcor_obj$r[c(t1, t2)])

  if (is.na(delta_r)) {
    stop("`t1` and/or `t2` correspond to time points where `r` is not defined. Please check `tcor_obj` and adjust `t1` and/or `t2` before testing.")
  }

  SE_delta_r <- sqrt(sum(tcor_obj$SE[c(t1, t2)]^2))
  Tstat <- delta_r/SE_delta_r

  ## run and output tests
  if (test == "student") {

    df <- 2*nrow(tcor_obj[!is.na(tcor_obj$r), ]) - 2
    pv_student <- 2*stats::pt(abs(Tstat), df = df, lower.tail = FALSE)
    return(data.frame(delta_r = delta_r, SE_delta_r = SE_delta_r, T_stat = Tstat, df_student = df, pv_student = pv_student))

  } else if (test == "chi2") {

    pv_chi2 <- stats::pchisq(Tstat^2, df = 1, lower.tail = FALSE)
    return(data.frame(delta_r = delta_r, SE_delta_r = SE_delta_r, chi2_stat = Tstat^2, df_chi2 = 1, pv_chi2 = pv_chi2))

  }
  stop("`test` unknown.")
}
