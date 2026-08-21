# Shared campaign test fixture. Used by test-app-io.R and test-app-modules.R.
# testthat auto-sources helper-*.R before running tests.

# A deterministic simulated campaign written to a temp CSV; returns the path.
# The campaign_server tests that must exercise a REAL upload (a changed
# `input$file`, and a quantal campaign, neither of which the bundled continuous
# example can provide) upload one of these.
campaign_csv <- function(par, response = "continuous", ...) {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(simulate_mixture(par, "CA", "SA", response, cv = 0, ...),
                   path, row.names = FALSE)
  path
}

campaign_fixture <- function() {
  data.frame(
    C1  = c(0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 1),
    C2  = c(0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 1, 2, 1, 2, 1),
    C3  = c(0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 1, 2, 1),
    Res = c(100, 90, 80, 70, 92, 84, 76, 95, 88, 80, 60, 40, 62, 44, 30)
  )
}
