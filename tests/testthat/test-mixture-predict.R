# Equivalence test: the unified `mix_response()` must reproduce every one of
# the 16 verbatim-ported predictors (CA/IA x reference/SA/DR/DL x binary/ternary)
# across a dense concentration grid, in both slope-sign regimes (all-positive =
# decreasing response, all-negative = increasing response).
#
# The old functions are the ground truth. CA models are solved by bisection, so
# they are compared at a loose tolerance; IA models are closed-form and compared
# tightly.
#
# Binary DR in the original carries a single interaction param `b` applied to z1
# only: exp((a + b*z1) * z1 * z2). The unified form takes a per-chemical b-vector,
# so the binary mapping is b -> c(b, 0).

CA_TOL <- 1e-4
IA_TOL <- 1e-8

# --- parameter sets -----------------------------------------------------------
# Increasing response = all slopes negative; decreasing = all slopes positive.
pars_bi  <- list(max = 800, ec50s = c(0.08, 50),       a = 0.5)
pars_tri <- list(max = 800, ec50s = c(0.08, 50, 10),   a = 0.5)
slopes_bi_dec  <- c(6, 0.4)
slopes_tri_dec <- c(6, 0.4, 2)

# binary deviation params
b_dr_bi <- 0.3            # old single-b form
b_dl_bi <- 0.2
# ternary deviation params
b_dr_tri <- c(0.3, 0.1, 0.2)
b_dl_tri <- 0.2

# --- concentration grids ------------------------------------------------------
grid_bi <- expand.grid(c1 = c(0, 0.05, 0.2), c2 = c(0, 10, 60))
grid_tri <- expand.grid(c1 = c(0, 0.05, 0.2),
                        c2 = c(0, 10, 60),
                        c3 = c(0, 5, 30))

expect_close <- function(actual, expected, tol, info) {
  # Skip rows the old model leaves undefined (NULL/NaN); the unified function is
  # only required to agree where the original returns a finite number.
  if (is.null(expected) || length(expected) == 0 ||
      !is.finite(expected)) return(invisible())
  expect_equal(actual, expected, tolerance = tol, info = info)
}

# ------------------------------------------------------------------------------
# BINARY
# ------------------------------------------------------------------------------
test_that("mixture_predict reproduces all binary CA models", {
  for (sgn in c(1, -1)) {
    sl <- sgn * slopes_bi_dec
    for (i in seq_len(nrow(grid_bi))) {
      c1 <- grid_bi$c1[i]; c2 <- grid_bi$c2[i]
      concs <- c(c1, c2)

      expect_close(
        mix_response(concs, pars_bi$max, sl, pars_bi$ec50s, "CA", "reference"),
        ca_bi(c1, c2, pars_bi$max, sl[1], sl[2], pars_bi$ec50s[1], pars_bi$ec50s[2]),
        CA_TOL, info = sprintf("CA ref bi sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, pars_bi$max, sl, pars_bi$ec50s, "CA", "SA", a = pars_bi$a),
        ca_sa_bi(c1, c2, pars_bi$max, sl[1], sl[2], pars_bi$ec50s[1], pars_bi$ec50s[2], pars_bi$a),
        CA_TOL, info = sprintf("CA SA bi sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, pars_bi$max, sl, pars_bi$ec50s, "CA", "DR",
                        a = pars_bi$a, b = c(b_dr_bi, 0)),
        ca_dr_bi(c1, c2, pars_bi$max, sl[1], sl[2], pars_bi$ec50s[1], pars_bi$ec50s[2], pars_bi$a, b_dr_bi),
        CA_TOL, info = sprintf("CA DR bi sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, pars_bi$max, sl, pars_bi$ec50s, "CA", "DL",
                        a = pars_bi$a, b = b_dl_bi),
        ca_dl_bi(c1, c2, pars_bi$max, sl[1], sl[2], pars_bi$ec50s[1], pars_bi$ec50s[2], pars_bi$a, b_dl_bi),
        CA_TOL, info = sprintf("CA DL bi sgn=%d row=%d", sgn, i))
    }
  }
})

test_that("mixture_predict reproduces all binary IA models", {
  for (sgn in c(1, -1)) {
    sl <- sgn * slopes_bi_dec
    for (i in seq_len(nrow(grid_bi))) {
      c1 <- grid_bi$c1[i]; c2 <- grid_bi$c2[i]
      concs <- c(c1, c2)

      expect_close(
        mix_response(concs, pars_bi$max, sl, pars_bi$ec50s, "IA", "reference"),
        ia_bi(c1, c2, pars_bi$max, sl[1], sl[2], pars_bi$ec50s[1], pars_bi$ec50s[2]),
        IA_TOL, info = sprintf("IA ref bi sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, pars_bi$max, sl, pars_bi$ec50s, "IA", "SA", a = pars_bi$a),
        ia_sa_bi(c1, c2, pars_bi$max, sl[1], sl[2], pars_bi$ec50s[1], pars_bi$ec50s[2], pars_bi$a),
        IA_TOL, info = sprintf("IA SA bi sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, pars_bi$max, sl, pars_bi$ec50s, "IA", "DR",
                        a = pars_bi$a, b = c(b_dr_bi, 0)),
        ia_dr_bi(c1, c2, pars_bi$max, sl[1], sl[2], pars_bi$ec50s[1], pars_bi$ec50s[2], pars_bi$a, b_dr_bi),
        IA_TOL, info = sprintf("IA DR bi sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, pars_bi$max, sl, pars_bi$ec50s, "IA", "DL",
                        a = pars_bi$a, b = b_dl_bi),
        ia_dl_bi(c1, c2, pars_bi$max, sl[1], sl[2], pars_bi$ec50s[1], pars_bi$ec50s[2], pars_bi$a, b_dl_bi),
        IA_TOL, info = sprintf("IA DL bi sgn=%d row=%d", sgn, i))
    }
  }
})

# ------------------------------------------------------------------------------
# TERNARY
# ------------------------------------------------------------------------------
test_that("mixture_predict reproduces all ternary CA models", {
  for (sgn in c(1, -1)) {
    sl <- sgn * slopes_tri_dec
    e <- pars_tri$ec50s; m <- pars_tri$max
    for (i in seq_len(nrow(grid_tri))) {
      c1 <- grid_tri$c1[i]; c2 <- grid_tri$c2[i]; c3 <- grid_tri$c3[i]
      concs <- c(c1, c2, c3)

      expect_close(
        mix_response(concs, m, sl, e, "CA", "reference"),
        ca_tri(c1, c2, c3, m, sl[1], sl[2], sl[3], e[1], e[2], e[3]),
        CA_TOL, info = sprintf("CA ref tri sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, m, sl, e, "CA", "SA", a = pars_tri$a),
        ca_sa_tri(c1, c2, c3, m, sl[1], sl[2], sl[3], e[1], e[2], e[3], pars_tri$a),
        CA_TOL, info = sprintf("CA SA tri sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, m, sl, e, "CA", "DR", a = pars_tri$a, b = b_dr_tri),
        ca_dr_tri(c1, c2, c3, m, sl[1], sl[2], sl[3], e[1], e[2], e[3],
                  pars_tri$a, b_dr_tri[1], b_dr_tri[2], b_dr_tri[3]),
        CA_TOL, info = sprintf("CA DR tri sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, m, sl, e, "CA", "DL", a = pars_tri$a, b = b_dl_tri),
        ca_dl_tri(c1, c2, c3, m, sl[1], sl[2], sl[3], e[1], e[2], e[3], pars_tri$a, b_dl_tri),
        CA_TOL, info = sprintf("CA DL tri sgn=%d row=%d", sgn, i))
    }
  }
})

test_that("mixture_predict reproduces all ternary IA models", {
  for (sgn in c(1, -1)) {
    sl <- sgn * slopes_tri_dec
    e <- pars_tri$ec50s; m <- pars_tri$max
    for (i in seq_len(nrow(grid_tri))) {
      c1 <- grid_tri$c1[i]; c2 <- grid_tri$c2[i]; c3 <- grid_tri$c3[i]
      concs <- c(c1, c2, c3)

      expect_close(
        mix_response(concs, m, sl, e, "IA", "reference"),
        ia_tri(c1, c2, c3, m, sl[1], sl[2], sl[3], e[1], e[2], e[3]),
        IA_TOL, info = sprintf("IA ref tri sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, m, sl, e, "IA", "SA", a = pars_tri$a),
        ia_sa_tri(c1, c2, c3, m, sl[1], sl[2], sl[3], e[1], e[2], e[3], pars_tri$a),
        IA_TOL, info = sprintf("IA SA tri sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, m, sl, e, "IA", "DR", a = pars_tri$a, b = b_dr_tri),
        ia_dr_tri(c1, c2, c3, m, sl[1], sl[2], sl[3], e[1], e[2], e[3],
                  pars_tri$a, b_dr_tri[1], b_dr_tri[2], b_dr_tri[3]),
        IA_TOL, info = sprintf("IA DR tri sgn=%d row=%d", sgn, i))

      expect_close(
        mix_response(concs, m, sl, e, "IA", "DL", a = pars_tri$a, b = b_dl_tri),
        ia_dl_tri(c1, c2, c3, m, sl[1], sl[2], sl[3], e[1], e[2], e[3], pars_tri$a, b_dl_tri),
        IA_TOL, info = sprintf("IA DL tri sgn=%d row=%d", sgn, i))
    }
  }
})
