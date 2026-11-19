ex <- function() read.csv(test_path("fixtures", "potentiation",
                                    "potentiation_example_continuous.csv"))

# Every analysis result produced below is kept so the obj1 <= obj2 invariant
# can be checked across all of them.
pot_results <- new.env()
keep <- function(name, r) { assign(name, r, envir = pot_results); r }

# Replicated continuous series from one log-logistic curve; one seed per
# replicate. `controls = FALSE` drops the C1 == 0 rows (the series then borrows
# the reference's controls).
sim_series <- function(label, ec50, seeds, controls = TRUE) {
  d <- do.call(rbind, lapply(seeds, function(s)
    simulate_single(c(max = 100, slope = 2, ec50 = ec50),
                    conc = c(0, 2^(-3:3)), cv = 0.05, seed = s)))
  if (!controls) d <- d[d$C1 > 0, , drop = FALSE]
  data.frame(Series = label, d)
}

recovery_data <- function() {
  rbind(sim_series("reference", 1,   1:4),
        sim_series("same",      1,   11:14, controls = FALSE),
        sim_series("potent",    0.5, 21:24, controls = FALSE),
        sim_series("weak",      2,   31:34, controls = FALSE))
}

test_that("workbook level 1: separate vs shared EC50 reproduces Sam's sheet", {
  r <- keep("example", analyse_potentiation(ex()))
  c1 <- r$comparisons[r$comparisons$series == "Level 1", ]
  expect_equal(c1$n, 90)
  expect_equal(c1$obj1, 301550.972, tolerance = 1e-6)
  expect_equal(c1$obj2, 341687.917, tolerance = 1e-6)
  expect_equal(c1$chi, 11.2463, tolerance = 1e-4)
  expect_equal(c1$p, 0.000798, tolerance = 1e-2)
  expect_equal(c1$ec50_series, 0.13417, tolerance = 1e-3)
  expect_true(c1$significant); expect_equal(c1$direction, "less toxic")
  expect_true(c1$controls_borrowed)
})

test_that("workbook level 2: no influence, X2 floored at ~0", {
  c2 <- analyse_potentiation(ex())$comparisons
  c2 <- c2[c2$series == "Level 2", ]
  expect_equal(c2$obj1, 393465.224, tolerance = 1e-6)
  expect_gte(c2$chi, 0); expect_lt(c2$chi, 1e-3)
  expect_false(c2$significant); expect_equal(c2$direction, "no significant influence")
})

test_that("workbook level 3: corrected shared fit, extrapolated EC50", {
  c3 <- analyse_potentiation(ex())$comparisons
  c3 <- c3[c3$series == "Level 3", ]
  expect_equal(c3$obj1, 379874.993, tolerance = 1e-6)
  expect_equal(c3$obj2, 411358.37, tolerance = 1e-6)   # NOT the workbook's 399798.70 (cell error)
  expect_equal(c3$chi, 7.166, tolerance = 1e-3)
  expect_equal(c3$ec50_series, 0.3607, tolerance = 1e-3)
  expect_true(c3$significant); expect_equal(c3$direction, "less toxic")
  expect_true(c3$ec50_extrapolated)                    # max C1 in Level 3 is 0.1
})

test_that("series and comparisons keep file order", {
  r <- analyse_potentiation(ex())
  expect_equal(r$series, c("Stressor 1", "Level 1", "Level 2", "Level 3"))
  expect_equal(r$comparisons$series, c("Level 1", "Level 2", "Level 3"))
  expect_equal(r$reference, "Stressor 1")
  expect_equal(names(r$curves), r$series)
  expect_equal(names(r$shared), r$comparisons$series)
  expect_identical(r$summary, NA_character_)
})

test_that("recovery: no shift is not significant, a shift is detected with the right direction", {
  r <- keep("recovery", analyse_potentiation(recovery_data()))
  cmp <- r$comparisons
  rownames(cmp) <- cmp$series
  expect_equal(r$reference, "reference")
  expect_false(cmp["same", "significant"])
  expect_equal(cmp["same", "direction"], "no significant influence")
  expect_true(cmp["potent", "significant"])
  expect_equal(cmp["potent", "direction"], "more toxic")
  expect_lt(cmp["potent", "ratio"], 1)
  expect_equal(cmp["weak", "direction"], "less toxic")
  expect_true(all(cmp$controls_borrowed))
  expect_equal(r$curves$potent$par[["ec50"]], 0.5, tolerance = 0.1)
  expect_equal(r$curves$weak$par[["ec50"]], 2, tolerance = 0.1)
})

test_that("a series with its own controls does not borrow", {
  d <- recovery_data()
  own_ctl <- data.frame(Series = "potent", C1 = 0,
                        Res = c(98, 101, 103, 99, 100, 102))
  d <- rbind(d, own_ctl)
  r <- keep("own_controls", analyse_potentiation(d))
  cmp <- r$comparisons[r$comparisons$series == "potent", ]
  n_potent <- sum(d$Series == "potent")
  n_ref <- sum(d$Series == "reference")
  expect_false(cmp$controls_borrowed)
  expect_false(r$curves$potent$controls_borrowed)
  expect_equal(cmp$n, n_potent + n_ref)
  expect_equal(nrow(r$curves$potent$data), n_potent)
  # the other series still borrow the reference's 4 controls
  expect_equal(nrow(r$curves$weak$data), sum(d$Series == "weak") + 4)
})

test_that("quantal path uses deviance", {
  sim_q <- function(label, ec50, seeds, controls = TRUE) {
    d <- do.call(rbind, lapply(seeds, function(s)
      simulate_single(c(max = 0.9, slope = 2, ec50 = ec50),
                      conc = c(0, 2^(-3:3)), response = "binary",
                      group_size = 20, seed = s)))
    if (!controls) d <- d[d$C1 > 0, , drop = FALSE]
    data.frame(Series = label, d)
  }
  d <- rbind(sim_q("alone", 1, 1:3), sim_q("with", 0.4, 4:6, controls = FALSE))
  r <- keep("quantal", analyse_potentiation(d, response = "quantal"))
  cmp <- r$comparisons
  expect_identical(r$response, "quantal")
  expect_identical(cmp$chi, pmax(0, cmp$obj2 - cmp$obj1))
  expect_equal(cmp$direction, "more toxic")
  for (cv in r$curves) expect_lte(cv$par[["max"]], 1)
  expect_lte(r$shared$with$par[["max_reference"]], 1)
  expect_lte(r$shared$with$par[["max_series"]], 1)
})

test_that("reference can be any series", {
  r <- keep("ref_level2", analyse_potentiation(ex(), reference = "Level 2"))
  expect_equal(r$reference, "Level 2")
  expect_equal(r$comparisons$series, c("Stressor 1", "Level 1", "Level 3"))
  expect_equal(r$series, c("Stressor 1", "Level 1", "Level 2", "Level 3"))
})

test_that("no controls anywhere still fits", {
  d <- ex()
  d <- d[d$C1 > 0, ]
  r <- keep("no_controls", expect_no_error(analyse_potentiation(d)))
  expect_false(any(r$comparisons$controls_borrowed))
  expect_false(any(vapply(r$curves, `[[`, logical(1), "controls_borrowed")))
  expect_equal(r$comparisons$n, c(50, 50, 50))
})

test_that("flat series does not error and is flagged", {
  d <- ex()
  flat <- data.frame(Series = "Flat", C1 = c(rep(0, 5), d$C1[d$Series == "Level 1"]),
                     Res = 100)
  d <- rbind(d, flat)
  r <- keep("flat", expect_no_error(analyse_potentiation(d)))
  cf <- r$comparisons[r$comparisons$series == "Flat", ]
  expect_type(cf$significant, "logical")
  expect_false(is.na(cf$significant))
  expect_true(cf$ec50_extrapolated)
  expect_true(r$curves$Flat$ec50_extrapolated)
})

test_that("obj1 <= obj2 always", {
  results <- as.list(pot_results)
  expect_gte(length(results), 7)
  for (nm in names(results)) {
    cmp <- results[[nm]]$comparisons
    expect_true(all(cmp$obj1 <= cmp$obj2), label = nm)
    expect_true(all(cmp$chi >= 0), label = nm)
  }
})

test_that("the exact gradient matches a central finite difference", {
  num_grad <- function(f, p, h = 1e-6) vapply(seq_along(p), function(i) {
    e <- replace(numeric(length(p)), i, h * abs(p[[i]]))
    (f(p + e) - f(p - e)) / (2 * e[[i]])
  }, numeric(1))
  dc <- data.frame(C1 = c(0, 0, 0.05, 0.1, 0.2, 0.4), Res = c(98, 102, 90, 55, 20, 3))
  dq <- data.frame(C1 = c(0, 0.5, 1, 2, 4), Exposed = 20,
                   Affected = c(18, 15, 9, 4, 0))
  for (case in list(list(d = dc, resp = "continuous", p = c(100, 3, 0.12)),
                    list(d = dq, resp = "quantal",    p = c(0.9, 1.7, 1.1)))) {
    f <- function(p) pot_objective(case$d, p[1], p[2], p[3], case$resp)
    expect_equal(unname(pot_gradient(case$d, case$p[1], case$p[2], case$p[3],
                                     case$resp)),
                 num_grad(f, case$p), tolerance = 1e-5, label = case$resp)
  }
})

test_that("errors on bad input", {
  expect_error(analyse_potentiation(ex(), reference = "Nope"), "Nope")
  one <- ex(); one <- one[one$Series == "Stressor 1", ]
  expect_error(analyse_potentiation(one), "two series")
  tiny <- ex(); tiny <- tiny[tiny$C1 > 0, ]
  tiny <- rbind(tiny, data.frame(Series = "Tiny", C1 = c(0.1, 0.2, 0.4),
                                 Res = c(500, 300, 10)))
  expect_error(analyse_potentiation(tiny), "Tiny")
})
