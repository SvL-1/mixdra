skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("plotly")

test_that("potentiation_ui has the reference-series selector", {
  ui <- potentiation_ui("p")
  expect_match(paste(as.character(ui), collapse = " "), "Reference series",
               fixed = TRUE)
})

test_that("potentiation_server analyses the bundled example until a file is uploaded", {
  shiny::testServer(potentiation_server, args = list(meta = shiny::reactiveValues()), {
    expect_equal(result()$comparisons$series, c("Level 1", "Level 2", "Level 3"))
    expect_match(controls_note(), "the 20 control observations of Stressor 1",
                 all = FALSE, fixed = TRUE)
    expect_match(controls_note(), "used for both curves of the comparison with Level 1",
                 all = FALSE, fixed = TRUE)
  })
})

test_that("the plot title uses the typed stressor names, with defaults when blank", {
  shiny::testServer(potentiation_server, args = list(meta = shiny::reactiveValues()), {
    title <- function()
      plotly::plotly_build(plot_obj())$x$layout$title$text
    expect_equal(title(), "Stressor 1 with and without Stressor 2")
    session$setInputs(name1 = "Chem", name2 = "Heat")
    expect_equal(title(), "Chem with and without Heat")
  })
})

test_that("an upload that is not a CSV gives a message, not a raw error", {
  f <- tempfile(fileext = ".csv")
  writeBin(as.raw(c(0, 1, 2, 255, 254, 0, 0, 7)), f)
  shiny::testServer(potentiation_server, args = list(meta = shiny::reactiveValues()), {
    session$setInputs(file = list(datapath = f, name = "x.csv"))
    html <- as.character(output$errors$html)
    expect_match(html, "text-danger", fixed = TRUE)
  })
})

test_that("an upload whose analysis errors shows the message in the errors area", {
  # No small input slips past validation and still breaks the engine, so the
  # engine failure is simulated.
  testthat::local_mocked_bindings(
    analyse_potentiation = function(...) stop("optimiser blew up"))
  f <- tempfile(fileext = ".csv")
  file.copy(system.file("extdata", "potentiation_example_continuous.csv",
                        package = "mixdra"), f)
  shiny::testServer(potentiation_server, args = list(meta = shiny::reactiveValues()), {
    session$setInputs(file = list(datapath = f, name = "x.csv"))
    html <- as.character(output$errors$html)
    expect_match(html, "text-danger", fixed = TRUE)
    expect_match(html, "optimiser blew up", fixed = TRUE)
    expect_error(result(), class = "shiny.silent.error")
  })
})
