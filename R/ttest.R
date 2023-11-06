#' Compute equality test between correlation coefficient estimates at two time points
#'
#' This function tests whether smoothed correlation values at two time points are equal (H0) or not.
#' The test is described page 341 in Choi & Shin (2021).
#'
#' Two different test statistics can be used, one is asymptotically Student-t distributed under H0 and one is chi-square distributed.
#' In practice, it seems to give very similar results.
#'
#' @param tcor_obj the output of a call to [`tcor()`] with `CI = TRUE`.
#' @param t1 the first time point used by the test (by default, the first time point in the time series).
#' @param t2 the second time point used by the test (by default, the last time point in the time series).
#' @param test a character string indicating which test to use ("student", the default; or "chi2").
#'
#' @return a data.frame with the result of the test, including the effect size (`delta_r = r[t2] - r[t1]`).
#'
#' @export
#'
#' @seealso [`test_ref()`], [`tcor()`]
#'
#' @examples
#' ## Simple example
#'
#' res <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = 50, CI = TRUE))
#' test_equality(res)
#'
#' ## Chi2 instead of Student's t-test
#'
#' test_equality(res, test = "chi2")
#'
#'
#' ## Time point can be dates or indices (mixing possible) but output as in input data
#'
#' test_equality(res, t1 = "2000-04-04", t2 = 1000)
#' res[1000, "t"] ## t2 matches with date in `res`
#' stockprice[1000, "DateID"] ## t2 does not match with date `stockprice` due to missing values
#'
#'
#' ## It could be useful to use `keep.missing = TRUE` for index to match original data despite NAs
#'
#' res2 <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID,
#'                               h = 50, CI = TRUE, keep.missing = TRUE))
#' test_equality(res2, t1 = "2000-04-04", t2 = 1000)
#' res[1000, "t"] ## t2 matches with date in `res`
#' stockprice[1000, "DateID"] ## t2 does match with date `stockprice` despite missing values
#'
test_equality <- function(tcor_obj, t1 = 1, t2 = nrow(tcor_obj), test = c("student", "chi2")) {

  ## check and prepare inputs
  test <- match.arg(test)

  if (!attr(tcor_obj, "CI")) {
    stop("For running an equality test, you must produce `tcor_obj` using the argument `CI = TRUE`.")
  }

  if (attr(tcor_obj, "h_selection") == "fixed by user") {
    message("The bandwidth stored in `tcor_obj` and used by `test_equality()` has not been automatically selected but was set by the user. Rerun `tcor()` without specifying `h` if you want otherwise.")
  }

  if (length(t1) != 1 || length(t2) != 1) {
    stop("`t1` and `t2` should be of length 1.")
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

  ## turn index into dates
  t1_date <- tcor_obj$t[t1]
  t2_date <- tcor_obj$t[t2]
  if ("Date" %in% class(tcor_obj$t)) {
    t1_date <- as.Date(t1_date)
    t2_date <- as.Date(t2_date)
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
    return(data.frame(t1 = t1_date, r1 = tcor_obj$r[t1], t2 = t2_date, r2 = tcor_obj$r[t2], delta_r = delta_r, SE_delta_r = SE_delta_r, T_stat = Tstat, df = df, p = pv_student))

  } else if (test == "chi2") {

    pv_chi2 <- stats::pchisq(Tstat^2, df = 1, lower.tail = FALSE)
    return(data.frame(t1 = t1_date, r1 = tcor_obj$r[t1], t2 = t2_date, r2 = tcor_obj$r[t2], delta_r = delta_r, SE_delta_r = SE_delta_r, chi2_stat = Tstat^2, df = 1, p = pv_chi2))

  }
  stop("`test` unknown.")
}


#' Test difference between correlation coefficient estimates and a value of reference
#'
#' This function tests whether smoothed correlation values are equal (H0) or not to a reference value (default = `0`).
#' The test is not described in Choi & Shin, 2021, but it is based on the idea behind [`test_equality()`].
#'
#' Two different test statistics can be used, one is asymptotically Student-t distributed under H0 and one is chi-square distributed.
#' In practice, it seems to give very similar results.
#'
#' @inheritParams test_equality
#' @param t a vector of time point(s) used by the test (by default, all time points are considered).
#' @param r_ref a scalar indicating the reference value for the correlation coefficient to be used in the test (default = `0`).
#' @param p.adjust.methods a character string indicating the method used to adjust p-values for multiple testing (see [`p.adjust()`]; default = "none").
#'
#' @return a data.frame with the result of the test, including the effect size (`delta_r = r[t] - r_ref`).
#'
#' @export
#'
#' @seealso [`test_equality()`], [`tcor()`]
#'
#' @examples
#' ## Comparison of all correlation values to reference of 0.5
#'
#' res <- with(stockprice, tcor(x = SP500, y = FTSE100, t = DateID, h = 300, CI = TRUE))
#' ref <- 0.5
#' test_against_ref <- test_ref(res, r_ref = ref)
#' head(test_against_ref)
#'
#'
#' ## Plot to illustrate the correspondance with confidence intervals
#'
#' plot(res$r ~ res$t, type = "l", ylim = c(0, 1), col = NULL)
#' abline(v = test_against_ref$t[test_against_ref$p > 0.05], col = "lightgreen")
#' abline(v = test_against_ref$t[test_against_ref$p < 0.05], col = "red")
#' points(res$r ~ res$t, type = "l")
#' points(res$upr ~ res$t, type = "l", lty = 2)
#' points(res$lwr ~ res$t, type = "l", lty = 2)
#' abline(h = ref, col = "blue")
#'
#'
#' ## Test correlation of 0 a specific time points (using index or dates)
#'
#' test_ref(res, t = c(100, 150))
#' test_ref(res, t = c("2000-08-18", "2000-10-27"))
#'
test_ref <- function(tcor_obj, t = tcor_obj$t, r_ref = 0, test = c("student", "chi2"), p.adjust.methods = c("none", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr")) {

  ## check and prepare inputs
  if (length(r_ref) != 1) {
    stop("`r_ref` must be of length 1.")
  }

  if (r_ref < -1 || r_ref > 1) {
    stop("`r_ref must belong to [-1, 1].`")
  }

  test <- match.arg(test)

  p.adjust.methods <- match.arg(p.adjust.methods)


  if (!attr(tcor_obj, "CI")) {
    stop("For running an equality test, you must produce `tcor_obj` using the argument `CI = TRUE`.")
  }

  if (attr(tcor_obj, "h_selection") == "fixed by user") {
    message("The bandwidth stored in `tcor_obj` and used by `test_ref()` has not been automatically selected but was set by the user. Rerun `tcor()` without specifying `h` if you want otherwise.")
  }


  ## turn non default date into index
  t_ori <- t_index_full <- t

  if (!is.numeric(t)) {
    if (any("Date" %in% class(tcor_obj$t))) {
      t <- as.Date(t) ## match() is very picky, so the class need to be correct
    }
    t <- t_index_full <- match(t, tcor_obj$t)
    t <- t[!is.na(t)]
  }

  ## compute t-test statistics
  delta_r <- tcor_obj$r[t] - r_ref

  if (all(is.na(delta_r))) {
    stop("`t` correspond to time points where `r` is not defined. Please check `tcor_obj` and adjust `t` before testing.")
  }

  SE_delta_r <- tcor_obj$SE[t]
  Tstat <- delta_r/SE_delta_r

  ## run and output tests
  if (test == "student") {

    df <- 2*nrow(tcor_obj[!is.na(tcor_obj$r), ]) - 2
    pv_student <- 2*stats::pt(abs(Tstat), df = df, lower.tail = FALSE)
    pv_student_corrected <- stats::p.adjust(pv_student, method = p.adjust.methods)

    d <- data.frame(r = tcor_obj$r[t], r_ref = r_ref, t_index = t, delta_r = delta_r, SE_delta_r = SE_delta_r, T_stat = Tstat, df = df, p = pv_student_corrected, p_adjustment = p.adjust.methods)
    d <- merge(d, data.frame(t_index = t_index_full, t = t_ori), all.y = TRUE)
    return(cbind(t = d$t, d[, !colnames(d) %in% c("t", "t_index")]))

  } else if (test == "chi2") {

    pv_chi2 <- stats::pchisq(Tstat^2, df = 1, lower.tail = FALSE)
    pv_chi2_corrected <- stats::p.adjust(pv_chi2, method = p.adjust.methods)

    d <- data.frame(r = tcor_obj$r[t], r_ref = r_ref, t_index = t, delta_r = delta_r, SE_delta_r = SE_delta_r, chi2_stat = Tstat^2, df = 1, p = pv_chi2_corrected, p_adjustment = p.adjust.methods)
    d <- merge(d, data.frame(t_index = t_index_full, t = t_ori), all.y = TRUE)
    return(cbind(t = d$t, d[, !colnames(d) %in% c("t", "t_index")]))
  }
  stop("`test` unknown.")
}
