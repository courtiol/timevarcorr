test_that("dates are not lost", {
  job <- with(stockprice[1:6, ], calc_rho(x = SP500, y = FTSE100, t = DateID, h = 20))
  expect_true(inherits(job$t, "Date"))
})

test_that("order of time points do not alter the result", {
  ordered <- with(stockprice[1:6, ], calc_rho(x = SP500, y = FTSE100, t = DateID, h = 20))
  disorderd <- with(stockprice[6:1, ], calc_rho(x = SP500, y = FTSE100, t = DateID, h = 20))
  reordered <- disorderd[6:1, ]
  rownames(reordered) <- NULL
  expect_equal(ordered, reordered)
})

test_that("correlations under `h = Inf` match time-invariant correlations", {
  rho_inf <- unique(with(stockprice[1:6, ], calc_rho(x = SP500, y = FTSE100, t = DateID, h = Inf, kernel = "box"))$rho)
  rho_inv <- cor(stockprice[1:6, "SP500"], stockprice[1:6, "FTSE100"])[[1]]
  expect_equal(rho_inf, rho_inv)
})
