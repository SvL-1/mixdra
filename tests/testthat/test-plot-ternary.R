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

test_that("plot_isoplane returns a plotly object with three traces", {
  skip_if_not_installed("plotly")
  p <- plot_isoplane(mock_res(), df = NULL, n = 8)
  expect_s3_class(p, "plotly")
  expect_equal(length(plotly::plotly_build(p)$x$data), 3)  # SA, ASA, markers
})

test_that("plot_sigma_tu returns a plotly object with the reference line", {
  skip_if_not_installed("plotly")
  p <- plot_sigma_tu(mock_res(), n = 11)
  expect_s3_class(p, "plotly")
  # 3 chems x 2 models = 6 line traces + 1 reference line = 7
  expect_equal(length(plotly::plotly_build(p)$x$data), 7)
})

test_that("plot_isoplane works with empty individual (skips the EC50 trace)", {
  skip_if_not_installed("plotly")
  res <- mock_res()
  res$individual <- res$individual[0, ]
  p <- plot_isoplane(res, df = NULL, n = 8)
  expect_s3_class(p, "plotly")
  expect_equal(length(plotly::plotly_build(p)$x$data), 2)  # SA + ASA, no markers
})
