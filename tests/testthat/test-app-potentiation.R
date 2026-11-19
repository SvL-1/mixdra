skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("plotly")

ex_df <- function()
  read.csv(system.file("extdata", "potentiation_example_continuous.csv",
                       package = "mixdra"))

# The example analysis, fitted once for the tests that only read it.
ex_result <- local({
  res <- NULL
  function() {
    if (is.null(res)) res <<- analyse_potentiation(ex_df())
    res
  }
})

test_that("potentiation_ui has the reference-series selector", {
  ui <- potentiation_ui("p")
  expect_match(paste(as.character(ui), collapse = " "), "Reference series",
               fixed = TRUE)
})

test_that("the sidebar explains quantal counts and why series order matters", {
  txt <- paste(as.character(potentiation_ui("p")), collapse = " ")
  expect_match(txt, "still showing the response (e.g. alive)", fixed = TRUE)
  expect_match(txt, "so the proportion declines with dose", fixed = TRUE)
  expect_match(txt, "increasing stressor-2 level", fixed = TRUE)
})

test_that("controls note names series without controls when nothing could be borrowed", {
  # Level 2 as reference: it has no controls, so Level 1 and 3 cannot borrow.
  note <- potentiation_controls_note(analyse_potentiation(ex_df(),
                                                          reference = "Level 2"))
  expect_match(note, "Level 1, Level 2 and Level 3 have no control observations",
               all = FALSE, fixed = TRUE)
  expect_match(note, "so their max is estimated from the dosed rows only",
               all = FALSE, fixed = TRUE)
  expect_false(any(grepl("each series uses its own", note, fixed = TRUE)))
  expect_false(any(grepl("counted twice", note, fixed = TRUE)))
})

test_that("controls note: one series without controls is named in the singular", {
  d <- ex_df()
  d <- d[d$Series %in% c("Stressor 1", "Level 2"), ]
  note <- potentiation_controls_note(analyse_potentiation(d, reference = "Level 2"))
  expect_equal(note, paste0("Controls: Level 2 has no control observations, so ",
                            "its max is estimated from the dosed rows only."))
})

test_that("controls note says 'own controls' only when every series has them", {
  d <- ex_df()
  ctl <- d[d$Series == "Stressor 1" & d$C1 == 0, ]
  for (s in c("Level 1", "Level 2", "Level 3"))
    d <- rbind(d, transform(ctl, Series = s))
  note <- potentiation_controls_note(analyse_potentiation(d))
  expect_equal(note, "Controls: each series uses its own control observations.")
})

test_that("the tests table shows SS to 2 decimals", {
  tb <- potentiation_tests_table(ex_result())
  expect_equal(names(tb), c("Series", "n", "SS run 1", "SS run 2", "X²", "p",
                            "Verdict"))
  cmp <- ex_result()$comparisons
  expect_equal(tb[["SS run 1"]], round(cmp$obj1, 2))
  expect_equal(tb[["SS run 2"]], round(cmp$obj2, 2))
  expect_equal(tb[["SS run 1"]][[1]], 301550.97)        # not signif(, 4) = 301600
  expect_equal(tb[["SS run 2"]][[1]], 341687.92)
  expect_equal(tb$Verdict, ex_result()$comparisons$direction)
})

test_that("typing a stressor-1 name updates the conclusion without refitting", {
  real <- analyse_potentiation
  calls <- 0
  local_mocked_bindings(analyse_potentiation = function(...) {
    calls <<- calls + 1
    real(...)
  })
  shiny::testServer(potentiation_server, args = list(meta = shiny::reactiveValues()), {
    expect_match(as.character(output$conclusion$html),
                 "Stressor 1 is significantly less toxic", fixed = TRUE)
    n <- calls
    session$setInputs(name1 = "Chem")
    expect_match(as.character(output$conclusion$html),
                 "Chem is significantly less toxic", fixed = TRUE)
    expect_equal(calls, n)
  })
})

test_that("an upload with missing values is rejected with the row numbers", {
  d <- ex_df(); d$Res[12] <- NA
  f <- tempfile(fileext = ".csv")
  utils::write.csv(d, f, row.names = FALSE)
  shiny::testServer(potentiation_server, args = list(meta = shiny::reactiveValues()), {
    session$setInputs(file = list(datapath = f, name = "x.csv"))
    html <- as.character(output$errors$html)
    expect_match(html, "Row(s) 12 have missing values", fixed = TRUE)
    expect_error(result(), class = "shiny.silent.error")
  })
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
    # read.csv warns about the junk bytes; only the message matters here.
    suppressWarnings(session$setInputs(file = list(datapath = f, name = "x.csv")))
    html <- suppressWarnings(as.character(output$errors$html))
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
