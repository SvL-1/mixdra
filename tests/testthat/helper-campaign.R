# Shared campaign test fixture. Used by test-app-io.R and test-app-modules.R.
# testthat auto-sources helper-*.R before running tests.

campaign_fixture <- function() {
  data.frame(
    C1  = c(0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 1),
    C2  = c(0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 1, 2, 1, 2, 1),
    C3  = c(0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 1, 2, 1),
    Res = c(100, 90, 80, 70, 92, 84, 76, 95, 88, 80, 60, 40, 62, 44, 30)
  )
}
