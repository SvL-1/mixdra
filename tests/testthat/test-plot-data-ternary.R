mock_res <- function(A1 = 0.7, A2 = -0.3, A3 = -1.8, A4 = 1,
                     ec50 = c(2, 5, 0.5), max = 100) {
  list(reference = "CA", response = "continuous",
       base = c(max = max, slope1 = 3, slope2 = 3, slope3 = 3,
                ec50_1 = ec50[1], ec50_2 = ec50[2], ec50_3 = ec50[3]),
       pairwise = c(A1 = A1, A2 = A2, A3 = A3), A4_overall = A4,
       individual = data.frame(ratio = c("a", "b"), C1 = c(1/3, 0.6),
                               C2 = c(1/3, 0.2), C3 = c(1/3, 0.2),
                               A4 = c(0.5, -0.3), n = c(5L, 5L),
                               stringsAsFactors = FALSE))
}

test_that("isoplane_plot_data stacks SA, ASA and markers with a series label", {
  res <- mock_res()
  d <- isoplane_plot_data(res, df = NULL, n = 6)
  expect_true(all(c("C1", "C2", "C3", "series") %in% names(d)))
  expect_setequal(unique(d$series), c("CA+S/A", "CA+S/A+S/A", "EC50"))
  expect_equal(sum(d$series == "EC50"), nrow(res$individual))
  expect_equal(sum(d$series == "CA+S/A"),     (6 + 1) * (6 + 2) / 2)
  expect_equal(sum(d$series == "CA+S/A+S/A"), (6 + 1) * (6 + 2) / 2)
})

test_that("sigma_tu_plot_data stacks SA and ASA per-chemical curves", {
  res <- mock_res()
  d <- sigma_tu_plot_data(res, n = 11)
  expect_true(all(c("chem", "z", "sigma_tu", "series") %in% names(d)))
  expect_setequal(unique(d$series), c("CA+S/A", "CA+S/A+S/A"))
  expect_setequal(unique(d$chem), c("C1", "C2", "C3"))
  expect_equal(sum(d$series == "CA+S/A"),     3 * 11)
  expect_equal(sum(d$series == "CA+S/A+S/A"), 3 * 11)
})

test_that("isoplane_plot_data handles empty individual (no ternary ratios)", {
  res <- mock_res()
  res$individual <- res$individual[0, ]   # 0-row, same columns
  d <- isoplane_plot_data(res, df = NULL, n = 6)
  # no EC50 rows, but both isoplane series present and no error
  expect_setequal(unique(d$series), c("CA+S/A", "CA+S/A+S/A"))
  expect_equal(sum(d$series == "EC50"), 0)
})
