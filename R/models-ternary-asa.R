# Hand-written (NOT a verbatim port): ternary "Advanced S/A" CA predictor with
# per-pair (A1,A2,A3) and three-way (A4) interaction coefficients. Generalises
# ca_sa_tri (single shared `a`); the single-`a` model is the special case
# A1=A2=A3=A4. Deviation function (candidate; verified against the workbook):
#   F = exp(A1*z1*z2 + A2*z1*z3 + A3*z2*z3 + A4*z1*z2*z3)
# with zi = TUi/sum(TU over present chems), TUi = Ci/EC50i. Binary branches
# activate only their pair's term; the ternary branch activates all four.

#' ca_asa_tri: ternary Advanced S/A CA predictor
#'
#' @param c1,c2,c3 Concentrations.
#' @param max,slope1,slope2,slope3,ec50_1,ec50_2,ec50_3 Base curve parameters.
#' @param A1,A2,A3 Pairwise interaction coefficients (A1: C1-C2, A2: C1-C3,
#'   A3: C2-C3). @param A4 Three-way interaction coefficient.
#' @return Predicted response (scalar).
#' @keywords internal
ca_asa_tri <- function(c1, c2, c3, max, slope1, slope2, slope3,
                       ec50_1, ec50_2, ec50_3, A1, A2, A3, A4) {
  .fn <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3,
                  Ec50_1, Ec50_2, Ec50_3, A1, A2, A3, A4) {
    L <- 0; U <- Max; Y <- 0
    if (C1 == 0 & C2 == 0 & C3 == 0) {
      if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) return(Max)
      else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) return(Max)
      else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) return(0)
      else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) return(0)
    }
    if (C1 > 0 & C2 == 0 & C3 == 0) return(Max / (1 + (C1 / Ec50_1)^Slope1))
    if (C1 == 0 & C2 > 0 & C3 == 0) return(Max / (1 + (C2 / Ec50_2)^Slope2))
    if (C1 == 0 & C2 == 0 & C3 > 0) return(Max / (1 + (C3 / Ec50_3)^Slope3))

    bisect <- function(F4_of, Smask) {
      L <- 0; U <- Max; Y <- 0
      while (abs(U - L) > 1e-06) {
        Y <- (L + U) / 2
        G <- F4_of(Y)
        if (Smask == "pos") { if (G < 0) L <- Y else U <- Y }
        else                { if (G > 0) L <- Y else U <- Y }
      }
      Y
    }

    if (C1 > 0 & C2 > 0 & C3 == 0) {
      TU1 <- C1 / Ec50_1; TU2 <- C2 / Ec50_2
      z1 <- TU1 / (TU1 + TU2); z2 <- TU2 / (TU1 + TU2)
      Smask <- if (Slope1 > 0 & Slope2 > 0) "pos" else "neg"
      return(bisect(function(Y) {
        ec1 <- Ec50_1 * ((Max - Y) / Y)^(1 / Slope1)
        ec2 <- Ec50_2 * ((Max - Y) / Y)^(1 / Slope2)
        (C1 / ec1) + (C2 / ec2) - exp(A1 * z1 * z2)
      }, Smask))
    }
    if (C1 > 0 & C2 == 0 & C3 > 0) {
      TU1 <- C1 / Ec50_1; TU3 <- C3 / Ec50_3
      z1 <- TU1 / (TU1 + TU3); z3 <- TU3 / (TU1 + TU3)
      Smask <- if (Slope1 > 0 & Slope3 > 0) "pos" else "neg"
      return(bisect(function(Y) {
        ec1 <- Ec50_1 * ((Max - Y) / Y)^(1 / Slope1)
        ec3 <- Ec50_3 * ((Max - Y) / Y)^(1 / Slope3)
        (C1 / ec1) + (C3 / ec3) - exp(A2 * z1 * z3)
      }, Smask))
    }
    if (C1 == 0 & C2 > 0 & C3 > 0) {
      TU2 <- C2 / Ec50_2; TU3 <- C3 / Ec50_3
      z2 <- TU2 / (TU2 + TU3); z3 <- TU3 / (TU2 + TU3)
      Smask <- if (Slope2 > 0 & Slope3 > 0) "pos" else "neg"
      return(bisect(function(Y) {
        ec2 <- Ec50_2 * ((Max - Y) / Y)^(1 / Slope2)
        ec3 <- Ec50_3 * ((Max - Y) / Y)^(1 / Slope3)
        (C2 / ec2) + (C3 / ec3) - exp(A3 * z2 * z3)
      }, Smask))
    }
    if (C1 > 0 & C2 > 0 & C3 > 0) {
      TU1 <- C1 / Ec50_1; TU2 <- C2 / Ec50_2; TU3 <- C3 / Ec50_3
      s <- TU1 + TU2 + TU3
      z1 <- TU1 / s; z2 <- TU2 / s; z3 <- TU3 / s
      Smask <- if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) "pos" else "neg"
      return(bisect(function(Y) {
        ec1 <- Ec50_1 * ((Max - Y) / Y)^(1 / Slope1)
        ec2 <- Ec50_2 * ((Max - Y) / Y)^(1 / Slope2)
        ec3 <- Ec50_3 * ((Max - Y) / Y)^(1 / Slope3)
        F4 <- exp(A1 * z1 * z2 + A2 * z1 * z3 + A3 * z2 * z3 + A4 * z1 * z2 * z3)
        (C1 / ec1) + (C2 / ec2) + (C3 / ec3) - F4
      }, Smask))
    }
  }
  .fn(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3,
      A1, A2, A3, A4)
}

# Vectorised over the concentration vectors (matches ca_*_tri_vec convention).
ca_asa_tri_vec <- Vectorize(ca_asa_tri, vectorize.args = c("c1", "c2", "c3"))
