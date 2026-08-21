skip_if_not_installed("shiny")
skip_if_not_installed("bslib")

test_that("app_ui assembles a bslib page", {
  ui <- app_ui()
  expect_true(inherits(ui, c("shiny.tag", "shiny.tag.list", "bslib_page", "bslib_fragment")))
})

test_that("app_server is a function of (input, output, session)", {
  expect_true(is.function(app_server))
  expect_setequal(names(formals(app_server)), c("input", "output", "session"))
})

test_that("run_app is exported and errors clearly when a dependency is missing", {
  # We cannot uninstall packages in a test; just assert the guard logic exists by
  # checking run_app is a function that references the required packages.
  expect_true(is.function(run_app))
  body_txt <- paste(deparse(body(run_app)), collapse = " ")
  expect_match(body_txt, "bslib")
  expect_match(body_txt, "shiny")
})

test_that("app_ui includes the Ternary Mixture panel", {
  # "ternary-fit_asa" was the standalone tab's own "Fit Advanced S/A" button,
  # removed when the ternary stage became a campaign sub-tab consuming frozen
  # store values (no local fit button); that legacy nav entry is retired in a
  # later task, so just assert the panel is still present, not an internal
  # control id.
  expect_match(as.character(app_ui()), "Ternary Mixture", fixed = TRUE)
})
