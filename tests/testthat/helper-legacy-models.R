# Legacy verbatim mixture predictors (CA/IA x reference/SA/DR/DL x binary/ternary).
# RELOCATED from R/models-binary.R + R/models-ternary.R. These are NO LONGER
# part of the package: production resolves mix_response() via model_spec().
# They are retained here solely as the equivalence oracle for
# test-mixture-predict.R (and the subjects of test-models-binary/ternary.R).
# testthat auto-sources helper-*.R before running tests.

# Binary (2-chemical) mixture model functions.
#
# Ported from MixTox_shiny_v2/functions/model_functions.R. Each source
# function already takes explicit named arguments; the only change here is
# renaming them to the package's lower-case convention. The bodies are
# reproduced verbatim (via deparse of the source), so the maths is identical
# (cross-checked numerically by tools/port_models.R: maxdiff == 0).

#' ca_bi: binary (2-chemical) mixture predictor (ported verbatim from CA_bi)
#'
#' @param c1,c2 Concentrations.
#' @param max,slope1,slope2,ec501,ec502 Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ca_bi <- function(c1, c2, max, slope1, slope2, ec501, ec502) {
  .fn <- function(C1, C2, max, slope1, slope2, ec501, ec502) {
    if (length(C1) > 1 || length(C2) > 1 || length(max) > 1 || 
        length(slope1) > 1 || length(slope2) > 1 || length(ec501) > 
        1 || length(ec502) > 1) {
        stop("CA_bi received vector inputs, expected scalars")
    }
    if (C1 == 0 && C2 == 0) {
        if (slope1 > 0 & slope2 > 0) {
            return(max)
        }
        if (slope1 < 0 && slope2 < 0) {
            return(0)
        }
    }
    if (C1 > 0 && C2 == 0) {
        return(max/(1 + (C1/ec501)^slope1))
    }
    if (C1 == 0 && C2 > 0) {
        return(max/(1 + (C2/ec502)^slope2))
    }
    L <- 0
    U <- max
    for (i in 1:50) {
        Y <- (L + U)/2
        if (Y <= 0 || Y >= max) {
            return(NaN)
        }
        ec1 <- ec501 * ((max - Y)/Y)^(1/slope1)
        ec2 <- ec502 * ((max - Y)/Y)^(1/slope2)
        G <- (C1/ec1) + (C2/ec2) - 1
        if (slope1 > 0 && slope2 > 0) {
            if (G < 0) {
                L <- Y
            }
            else {
                U <- Y
            }
        }
        if (slope1 < 0 && slope2 < 0) {
            if (G > 0) {
                L <- Y
            }
            else {
                U <- Y
            }
        }
    }
    return(Y)
}
  .fn(c1, c2, max, slope1, slope2, ec501, ec502)
}

#' ia_bi: binary (2-chemical) mixture predictor (ported verbatim from IA_bi)
#'
#' @param c1,c2 Concentrations.
#' @param max,slope1,slope2,ec501,ec502 Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ia_bi <- function(c1, c2, max, slope1, slope2, ec501, ec502) {
  .fn <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2) {
    if (C1 == 0 & C2 == 0) {
        if (Slope1 > 0 & Slope2 > 0) {
            return(Max)
        }
        if (Slope1 < 0 & Slope2 < 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
    F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
    if (Slope1 > 0 & Slope2 > 0) {
        return(Max * F1 * F2)
    }
    if (Slope1 < 0 & Slope2 < 0) {
        return(Max * (F1 + F2 - F1 * F2))
    }
}
  .fn(c1, c2, max, slope1, slope2, ec501, ec502)
}

#' ca_sa_bi: binary (2-chemical) mixture predictor (ported verbatim from CA_SA_bi)
#'
#' @param c1,c2 Concentrations.
#' @param max,slope1,slope2,ec501,ec502,a Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ca_sa_bi <- function(c1, c2, max, slope1, slope2, ec501, ec502, a) {
  .fn <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a) {
    TU1 <- C1/Ec50_1
    TU2 <- C2/Ec50_2
    z1 <- TU1/(TU1 + TU2)
    z2 <- TU2/(TU1 + TU2)
    if (C1 == 0 & C2 == 0) {
        if (Slope1 > 0 & Slope2 > 0) {
            return(Max)
        }
        if (Slope1 < 0 & Slope2 < 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    L <- 0
    U <- Max
    for (i in 1:50) {
        Y <- (L + U)/2
        ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
        ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
        F1 <- exp(a * z1 * z2)
        G <- (C1/ec1) + (C2/ec2) - F1
        if (Slope1 > 0 & Slope2 > 0) {
            if (G < 0) {
                L <- Y
            }
            else {
                U <- Y
            }
        }
        if (Slope1 < 0 & Slope2 < 0) {
            if (G > 0) {
                L <- Y
            }
            else {
                U <- Y
            }
        }
    }
    return(Y)
}
  .fn(c1, c2, max, slope1, slope2, ec501, ec502, a)
}

#' ca_dr_bi: binary (2-chemical) mixture predictor (ported verbatim from CA_DR_bi)
#'
#' @param c1,c2 Concentrations.
#' @param max,slope1,slope2,ec501,ec502,a,b Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ca_dr_bi <- function(c1, c2, max, slope1, slope2, ec501, ec502, a, b) {
  .fn <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a, b) {
    TU1 <- C1/Ec50_1
    TU2 <- C2/Ec50_2
    z1 <- TU1/(TU1 + TU2)
    z2 <- TU2/(TU1 + TU2)
    if (C1 == 0 & C2 == 0) {
        if (Slope1 > 0 & Slope2 > 0) {
            return(Max)
        }
        if (Slope1 < 0 & Slope2 < 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    L <- 0
    U <- Max
    for (i in 1:50) {
        Y <- (L + U)/2
        ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
        ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
        F1 <- exp((a + b * z1) * z1 * z2)
        G <- (C1/ec1) + (C2/ec2) - F1
        if (Slope1 > 0 & Slope2 > 0) {
            if (G < 0) {
                L <- Y
            }
            else {
                U <- Y
            }
        }
        if (Slope1 < 0 & Slope2 < 0) {
            if (G > 0) {
                L <- Y
            }
            else {
                U <- Y
            }
        }
    }
    return(Y)
}
  .fn(c1, c2, max, slope1, slope2, ec501, ec502, a, b)
}

#' ca_dl_bi: binary (2-chemical) mixture predictor (ported verbatim from CA_DL_bi)
#'
#' @param c1,c2 Concentrations.
#' @param max,slope1,slope2,ec501,ec502,a,b Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ca_dl_bi <- function(c1, c2, max, slope1, slope2, ec501, ec502, a, b) {
  .fn <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a, b) {
    if (C1 == 0 & C2 == 0) {
        if (Slope1 > 0 & Slope2 > 0) {
            return(Max)
        }
        if (Slope1 < 0 & Slope2 < 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    TU1 <- C1/Ec50_1
    TU2 <- C2/Ec50_2
    z1 <- TU1/(TU1 + TU2)
    z2 <- TU2/(TU1 + TU2)
    L <- 0
    U <- Max
    for (i in 1:50) {
        Y <- (L + U)/2
        ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
        ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
        F1 <- exp((a * (1 - (b * (TU1 + TU2)))) * z1 * z2)
        G <- (C1/ec1) + (C2/ec2) - F1
        if (Slope1 > 0 & Slope2 > 0) {
            if (G < 0) {
                L <- Y
            }
            else {
                U <- Y
            }
        }
        if (Slope1 < 0 & Slope2 < 0) {
            if (G > 0) {
                L <- Y
            }
            else {
                U <- Y
            }
        }
    }
    return(Y)
}
  .fn(c1, c2, max, slope1, slope2, ec501, ec502, a, b)
}

#' ia_sa_bi: binary (2-chemical) mixture predictor (ported verbatim from IA_SA_bi)
#'
#' @param c1,c2 Concentrations.
#' @param max,slope1,slope2,ec501,ec502,a Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ia_sa_bi <- function(c1, c2, max, slope1, slope2, ec501, ec502, a) {
  .fn <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a) {
    if (C1 == 0 & C2 == 0) {
        if (Slope1 > 0 & Slope2 > 0) {
            return(Max)
        }
        if (Slope1 < 0 & Slope2 < 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    TU1 <- C1/Ec50_1
    TU2 <- C2/Ec50_2
    z1 <- TU1/(TU1 + TU2)
    z2 <- TU2/(TU1 + TU2)
    F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
    F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
    F3 <- a * z1 * z2
    if (Slope1 > 0 & Slope2 > 0) {
        Trans <- qnorm(F1 * F2)
        P <- pnorm(Trans + F3)
    }
    if (Slope1 < 0 & Slope2 < 0) {
        Trans <- qnorm(F1 + F2 - F1 * F2)
        P <- pnorm(Trans - F3)
    }
    return(Max * P)
}
  .fn(c1, c2, max, slope1, slope2, ec501, ec502, a)
}

#' ia_dr_bi: binary (2-chemical) mixture predictor (ported verbatim from IA_DR_bi)
#'
#' @param c1,c2 Concentrations.
#' @param max,slope1,slope2,ec501,ec502,a,b Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ia_dr_bi <- function(c1, c2, max, slope1, slope2, ec501, ec502, a, b) {
  .fn <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a, b) {
    if (C1 == 0 & C2 == 0) {
        if (Slope1 > 0 & Slope2 > 0) {
            return(Max)
        }
        if (Slope1 < 0 & Slope2 < 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    TU1 <- C1/Ec50_1
    TU2 <- C2/Ec50_2
    z1 <- TU1/(TU1 + TU2)
    z2 <- TU2/(TU1 + TU2)
    F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
    F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
    F3 <- (a + b * z1) * z1 * z2
    if (Slope1 > 0 & Slope2 > 0) {
        Trans <- qnorm(F1 * F2)
        P <- pnorm(Trans + F3)
    }
    if (Slope1 < 0 & Slope2 < 0) {
        Trans <- qnorm(F1 + F2 - F1 * F2)
        P <- pnorm(Trans - F3)
    }
    return(Max * P)
}
  .fn(c1, c2, max, slope1, slope2, ec501, ec502, a, b)
}

#' ia_dl_bi: binary (2-chemical) mixture predictor (ported verbatim from IA_DL_bi)
#'
#' @param c1,c2 Concentrations.
#' @param max,slope1,slope2,ec501,ec502,a,b Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ia_dl_bi <- function(c1, c2, max, slope1, slope2, ec501, ec502, a, b) {
  .fn <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a, b) {
    if (C1 == 0 & C2 == 0) {
        if (Slope1 > 0 & Slope2 > 0) {
            return(Max)
        }
        if (Slope1 < 0 & Slope2 < 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    TU1 <- C1/Ec50_1
    TU2 <- C2/Ec50_2
    z1 <- TU1/(TU1 + TU2)
    z2 <- TU2/(TU1 + TU2)
    F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
    F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
    if (Slope1 > 0 & Slope2 > 0) {
        P50 <- 1 - (F1 * F2)
        F3 <- (a * (1 - (b * P50))) * z1 * z2
        Trans <- qnorm(F1 * F2)
        P <- pnorm(Trans + F3)
    }
    if (Slope1 < 0 & Slope2 < 0) {
        P50 <- (F1 + F2 - F1 * F2)
        F3 <- (a * (1 - (b * P50))) * z1 * z2
        Trans <- qnorm(F1 + F2 - F1 * F2)
        P <- pnorm(Trans - F3)
    }
    return(Max * P)
}
  .fn(c1, c2, max, slope1, slope2, ec501, ec502, a, b)
}

# Vectorised wrappers over the c1/c2 concentration vectors.
ca_bi_vec <- Vectorize(ca_bi, vectorize.args = c("c1", "c2"))
ia_bi_vec <- Vectorize(ia_bi, vectorize.args = c("c1", "c2"))
ca_sa_bi_vec <- Vectorize(ca_sa_bi, vectorize.args = c("c1", "c2"))
ca_dr_bi_vec <- Vectorize(ca_dr_bi, vectorize.args = c("c1", "c2"))
ca_dl_bi_vec <- Vectorize(ca_dl_bi, vectorize.args = c("c1", "c2"))
ia_sa_bi_vec <- Vectorize(ia_sa_bi, vectorize.args = c("c1", "c2"))
ia_dr_bi_vec <- Vectorize(ia_dr_bi, vectorize.args = c("c1", "c2"))
ia_dl_bi_vec <- Vectorize(ia_dl_bi, vectorize.args = c("c1", "c2"))


# Ternary (3-chemical) mixture model functions.
#
# Ported from MixTox_shiny_v2/functions/model_functions.R. Each source
# function already takes explicit named arguments; the only change here is
# renaming them to the package's lower-case convention. The bodies are
# reproduced verbatim (via deparse of the source), so the maths is identical
# (cross-checked numerically by tools/port_models.R: maxdiff == 0).

#' ca_tri: ternary (3-chemical) mixture predictor (ported verbatim from CA)
#'
#' @param c1,c2,c3 Concentrations.
#' @param max,slope1,slope2,slope3,ec50_1,ec50_2,ec50_3 Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ca_tri <- function(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3) {
  .fn <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3) {
    L <- 0
    U <- Max
    Y <- 0
    ec1 <- 0
    ec2 <- 0
    ec3 <- 0
    G <- 0
    if (C1 == 0 & C2 == 0 & C3 == 0) {
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            return(Max)
        }
        else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
            return(Max)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            return(0)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0 & C3 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0 & C3 == 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    if (C1 == 0 & C2 == 0 & C3 > 0) {
        return(Max/(1 + (C3/Ec50_3)^Slope3))
    }
    if (C1 > 0 & C2 > 0 & C3 == 0) {
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            G <- (C1/ec1) + (C2/ec2) - 1
            if (Slope1 > 0 & Slope2 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope2 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 > 0 & C2 == 0 & C3 > 0) {
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            G <- (C1/ec1) + (C3/ec3) - 1
            if (Slope1 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 == 0 & C2 > 0 & C3 > 0) {
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            G <- (C2/ec2) + (C3/ec3) - 1
            if (Slope2 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope2 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 > 0 & C2 > 0 & C3 > 0) {
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            G <- (C1/ec1) + (C2/ec2) + (C3/ec3) - 1
            if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
}
  .fn(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3)
}

#' ia_tri: ternary (3-chemical) mixture predictor (ported verbatim from IA)
#'
#' @param c1,c2,c3 Concentrations.
#' @param max,slope1,slope2,slope3,ec50_1,ec50_2,ec50_3 Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ia_tri <- function(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3) {
  .fn <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3) {
    F1 <- 0
    F2 <- 0
    F3 <- 0
    if (C1 == 0 & C2 == 0 & C3 == 0) {
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            return(Max)
        }
        else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
            return(Max)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            return(0)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0 & C3 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0 & C3 == 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    if (C1 == 0 & C2 == 0 & C3 > 0) {
        return(Max/(1 + (C3/Ec50_3)^Slope3))
    }
    if (C1 > 0 & C2 > 0 & C3 == 0) {
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        if (Slope1 > 0 & Slope2 > 0) {
            return(Max * F1 * F2)
        }
        else if (Slope1 < 0 & Slope2 < 0) {
            return(Max * (F1 + F2 - F1 * F2))
        }
    }
    if (C1 > 0 & C2 == 0 & C3 > 0) {
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        if (Slope1 > 0 & Slope3 > 0) {
            return(Max * F1 * F3)
        }
        else if (Slope1 < 0 & Slope3 < 0) {
            return(Max * (F1 + F3 - F1 * F3))
        }
    }
    if (C1 == 0 & C2 > 0 & C3 > 0) {
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        if (Slope2 > 0 & Slope3 > 0) {
            return(Max * F2 * F3)
        }
        else if (Slope2 < 0 & Slope3 < 0) {
            return(Max * (F2 + F3 - F2 * F3))
        }
    }
    if (C1 > 0 & C2 > 0 & C3 > 0) {
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            return(Max * F1 * F2 * F3)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            return(Max * (F1 + F2 + F3 - (F1 * F2) - (F1 * F3) - 
                (F2 * F3) + (F1 * F2 * F3)))
        }
    }
}
  .fn(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3)
}

#' ca_sa_tri: ternary (3-chemical) mixture predictor (ported verbatim from CA_SA)
#'
#' @param c1,c2,c3 Concentrations.
#' @param max,slope1,slope2,slope3,ec50_1,ec50_2,ec50_3,a Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ca_sa_tri <- function(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a) {
  .fn <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a) {
    F1 <- 0
    F2 <- 0
    F3 <- 0
    F4 <- 0
    z1 <- 0
    z2 <- 0
    z3 <- 0
    L <- 0
    U <- Max
    ec1 <- ec2 <- ec3 <- G <- 0
    Y <- 0
    if (C1 == 0 & C2 == 0 & C3 == 0) {
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            return(Max)
        }
        else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
            return(Max)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            return(0)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0 & C3 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0 & C3 == 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    if (C1 == 0 & C2 == 0 & C3 > 0) {
        return(Max/(1 + (C3/Ec50_3)^Slope3))
    }
    if (C1 > 0 & C2 > 0 & C3 == 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        z1 <- TU1/(TU1 + TU2)
        z2 <- TU2/(TU1 + TU2)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            F4 <- exp(a * z1 * z2)
            G <- (C1/ec1) + (C2/ec2) - F4
            if (Slope1 > 0 & Slope2 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope2 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 > 0 & C2 == 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU3)
        z3 <- TU3/(TU1 + TU3)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            F4 <- exp(a * z1 * z3)
            G <- (C1/ec1) + (C3/ec3) - F4
            if (Slope1 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 == 0 & C2 > 0 & C3 > 0) {
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z2 <- TU2/(TU2 + TU3)
        z3 <- TU3/(TU2 + TU3)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            F4 <- exp(a * z2 * z3)
            G <- (C2/ec2) + (C3/ec3) - F4
            if (Slope2 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope2 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 > 0 & C2 > 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU2 + TU3)
        z2 <- TU2/(TU1 + TU2 + TU3)
        z3 <- TU3/(TU1 + TU2 + TU3)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            F4 <- exp(a * z1 * z2 * z3)
            G <- (C1/ec1) + (C2/ec2) + (C3/ec3) - F4
            if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
}
  .fn(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a)
}

#' ca_dr_tri: ternary (3-chemical) mixture predictor (ported verbatim from CA_DR)
#'
#' @param c1,c2,c3 Concentrations.
#' @param max,slope1,slope2,slope3,ec50_1,ec50_2,ec50_3,a,b1,b2,b3 Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ca_dr_tri <- function(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a, b1, b2, b3) {
  .fn <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a, b1, b2, b3) {
    L <- 0
    U <- Max
    Y <- 0
    ec1 <- 0
    ec2 <- 0
    ec3 <- 0
    F1 <- 0
    G <- 0
    TU1 <- 0
    TU2 <- 0
    TU3 <- 0
    z1 <- 0
    z2 <- 0
    z3 <- 0
    if (C1 == 0 & C2 == 0 & C3 == 0) {
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            return(Max)
        }
        else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
            return(Max)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            return(0)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0 & C3 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0 & C3 == 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    if (C1 == 0 & C2 == 0 & C3 > 0) {
        return(Max/(1 + (C3/Ec50_3)^Slope3))
    }
    if (C1 > 0 & C2 > 0 & C3 == 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU2 + TU3)
        z2 <- TU2/(TU1 + TU2 + TU3)
        z3 <- TU3/(TU1 + TU2 + TU3)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            F1 <- exp((a + (b1 * z1) + (b2 * z2)) * z1 * z2)
            G <- (C1/ec1) + (C2/ec2) - F1
            if (Slope1 > 0 & Slope2 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope2 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 > 0 & C2 == 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU2 + TU3)
        z2 <- TU2/(TU1 + TU2 + TU3)
        z3 <- TU3/(TU1 + TU2 + TU3)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            F1 <- exp((a + (b1 * z1) + (b3 * z3)) * z1 * z3)
            G <- (C1/ec1) + (C3/ec3) - F1
            if (Slope1 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 == 0 & C2 > 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU2 + TU3)
        z2 <- TU2/(TU1 + TU2 + TU3)
        z3 <- TU3/(TU1 + TU2 + TU3)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            F1 <- exp((a + (b2 * z2) + (b3 * z3)) * z2 * z3)
            G <- (C2/ec2) + (C3/ec3) - F1
            if (Slope2 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope2 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 > 0 & C2 > 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU2 + TU3)
        z2 <- TU2/(TU1 + TU2 + TU3)
        z3 <- TU3/(TU1 + TU2 + TU3)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            F1 <- exp((a + (b1 * z1) + (b2 * z2) + (b3 * z3)) * 
                z1 * z2 * z3)
            G <- (C1/ec1) + (C2/ec2) + (C3/ec3) - F1
            if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
}
  .fn(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a, b1, b2, b3)
}

#' ca_dl_tri: ternary (3-chemical) mixture predictor (ported verbatim from CA_DL)
#'
#' @param c1,c2,c3 Concentrations.
#' @param max,slope1,slope2,slope3,ec50_1,ec50_2,ec50_3,a,b Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ca_dl_tri <- function(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a, b) {
  .fn <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a, b) {
    L <- 0
    U <- Max
    Y <- 0
    ec1 <- 0
    ec2 <- 0
    ec3 <- 0
    F1 <- 0
    G <- 0
    TU1 <- 0
    TU2 <- 0
    TU3 <- 0
    z1 <- 0
    z2 <- 0
    z3 <- 0
    if (C1 == 0 & C2 == 0 & C3 == 0) {
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            return(Max)
        }
        else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
            return(Max)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            return(0)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0 & C3 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0 & C3 == 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    if (C1 == 0 & C2 == 0 & C3 > 0) {
        return(Max/(1 + (C3/Ec50_3)^Slope3))
    }
    if (C1 > 0 & C2 > 0 & C3 == 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        z1 <- TU1/(TU1 + TU2)
        z2 <- TU2/(TU1 + TU2)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            F1 <- exp((a * (1 - (b * (TU1 + TU2)))) * z1 * z2)
            G <- (C1/ec1) + (C2/ec2) - F1
            if (Slope1 > 0 & Slope2 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope2 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 > 0 & C2 == 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU3)
        z3 <- TU3/(TU1 + TU3)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            F1 <- exp((a * (1 - (b * (TU1 + TU3)))) * z1 * z3)
            G <- (C1/ec1) + (C3/ec3) - F1
            if (Slope1 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 == 0 & C2 > 0 & C3 > 0) {
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z2 <- TU2/(TU2 + TU3)
        z3 <- TU3/(TU2 + TU3)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            F1 <- exp((a * (1 - (b * (TU2 + TU3)))) * z2 * z3)
            G <- (C2/ec2) + (C3/ec3) - F1
            if (Slope2 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope2 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
    if (C1 > 0 & C2 > 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU2 + TU3)
        z2 <- TU2/(TU1 + TU2 + TU3)
        z3 <- TU3/(TU1 + TU2 + TU3)
        L <- 0
        U <- Max
        while (abs(U - L) > 1e-06) {
            Y <- (L + U)/2
            ec1 <- Ec50_1 * ((Max - Y)/Y)^(1/Slope1)
            ec2 <- Ec50_2 * ((Max - Y)/Y)^(1/Slope2)
            ec3 <- Ec50_3 * ((Max - Y)/Y)^(1/Slope3)
            F1 <- exp((a * (1 - (b * (TU1 + TU2 + TU3)))) * z1 * 
                z2 * z3)
            G <- (C1/ec1) + (C2/ec2) + (C3/ec3) - F1
            if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
                if (G < 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
            if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
                if (G > 0) {
                  L <- Y
                }
                else {
                  U <- Y
                }
            }
        }
        return(Y)
    }
}
  .fn(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a, b)
}

#' ia_sa_tri: ternary (3-chemical) mixture predictor (ported verbatim from IA_SA)
#'
#' @param c1,c2,c3 Concentrations.
#' @param max,slope1,slope2,slope3,ec50_1,ec50_2,ec50_3,a Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ia_sa_tri <- function(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a) {
  .fn <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a) {
    F1 <- 0
    F2 <- 0
    F3 <- 0
    F4 <- 0
    z1 <- 0
    z2 <- 0
    z3 <- 0
    Trans <- 0
    P <- 0
    if (C1 == 0 & C2 == 0 & C3 == 0) {
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            return(Max)
        }
        else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
            return(Max)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            return(0)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0 & C3 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0 & C3 == 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    if (C1 == 0 & C2 == 0 & C3 > 0) {
        return(Max/(1 + (C3/Ec50_3)^Slope3))
    }
    if (C1 > 0 & C2 > 0 & C3 == 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        z1 <- TU1/(TU1 + TU2)
        z2 <- TU2/(TU1 + TU2)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F4 <- a * z1 * z2
        if (Slope1 > 0 & Slope2 > 0) {
            Trans <- qnorm(F1 * F2)
            P <- pnorm(Trans + F4)
        }
        if (Slope1 < 0 & Slope2 < 0) {
            Trans <- qnorm(F1 + F2 - F1 * F2)
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
    if (C1 > 0 & C2 == 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU3)
        z3 <- TU3/(TU1 + TU3)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        F4 <- a * z1 * z3
        if (Slope1 > 0 & Slope3 > 0) {
            Trans <- qnorm(F1 * F3)
            P <- pnorm(Trans + F4)
        }
        if (Slope1 < 0 & Slope3 < 0) {
            Trans <- qnorm(F1 + F3 - F1 * F3)
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
    if (C1 == 0 & C2 > 0 & C3 > 0) {
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z2 <- TU2/(TU2 + TU3)
        z3 <- TU3/(TU2 + TU3)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        F4 <- a * z2 * z3
        if (Slope2 > 0 & Slope3 > 0) {
            Trans <- qnorm(F2 * F3)
            P <- pnorm(Trans + F4)
        }
        if (Slope2 < 0 & Slope3 < 0) {
            Trans <- qnorm(F2 + F3 - F2 * F3)
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
    if (C1 > 0 & C2 > 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU2 + TU3)
        z2 <- TU2/(TU1 + TU2 + TU3)
        z3 <- TU3/(TU1 + TU2 + TU3)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        F4 <- a * z1 * z2 * z3
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            Trans <- qnorm(F1 * F2 * F3)
            P <- pnorm(Trans + F4)
        }
        if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            Trans <- qnorm(F1 + F2 + F3 - (F1 * F2) - (F1 * F3) - 
                (F2 * F3) + (F1 * F2 * F3))
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
}
  .fn(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a)
}

#' ia_dr_tri: ternary (3-chemical) mixture predictor (ported verbatim from IA_DR)
#'
#' @param c1,c2,c3 Concentrations.
#' @param max,slope1,slope2,slope3,ec50_1,ec50_2,ec50_3,a,b1,b2,b3 Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ia_dr_tri <- function(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a, b1, b2, b3) {
  .fn <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a, b1, b2, b3) {
    F1 <- 0
    F2 <- 0
    F3 <- 0
    F4 <- 0
    z1 <- 0
    z2 <- 0
    z3 <- 0
    Trans <- 0
    P <- 0
    if (C1 == 0 & C2 == 0 & C3 == 0) {
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            return(Max)
        }
        else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
            return(Max)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            return(0)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0 & C3 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0 & C3 == 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    if (C1 == 0 & C2 == 0 & C3 > 0) {
        return(Max/(1 + (C3/Ec50_3)^Slope3))
    }
    if (C1 > 0 & C2 > 0 & C3 == 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        z1 <- TU1/(TU1 + TU2)
        z2 <- TU2/(TU1 + TU2)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F4 <- (a + (b1 * z1) + (b2 * z2)) * z1 * z2
        if (Slope1 > 0 & Slope2 > 0) {
            Trans <- qnorm(F1 * F2)
            P <- pnorm(Trans + F4)
        }
        if (Slope1 < 0 & Slope2 < 0) {
            Trans <- qnorm(F1 + F2 - F1 * F2)
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
    if (C1 > 0 & C2 == 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU3)
        z3 <- TU3/(TU1 + TU3)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        F4 <- (a + (b1 * z1) + (b3 * z3)) * z1 * z3
        if (Slope1 > 0 & Slope3 > 0) {
            Trans <- qnorm(F1 * F3)
            P <- pnorm(Trans + F4)
        }
        if (Slope1 < 0 & Slope3 < 0) {
            Trans <- qnorm(F1 + F3 - F1 * F3)
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
    if (C1 == 0 & C2 > 0 & C3 > 0) {
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z2 <- TU2/(TU2 + TU3)
        z3 <- TU3/(TU2 + TU3)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        F4 <- (a + (b2 * z2) + (b3 * z3)) * z2 * z3
        if (Slope2 > 0 & Slope3 > 0) {
            Trans <- qnorm(F2 * F3)
            P <- pnorm(Trans + F4)
        }
        if (Slope2 < 0 & Slope3 < 0) {
            Trans <- qnorm(F2 + F3 - F2 * F3)
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
    if (C1 > 0 & C2 > 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU2 + TU3)
        z2 <- TU2/(TU1 + TU2 + TU3)
        z3 <- TU3/(TU1 + TU2 + TU3)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        F4 <- (a + (b1 * z1) + (b2 * z2) + (b3 * z3)) * z1 * 
            z2 * z3
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            Trans <- qnorm(F1 * F2 * F3)
            P <- pnorm(Trans + F4)
        }
        if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            Trans <- qnorm(F1 + F2 + F3 - (F1 * F2) - (F1 * F3) - 
                (F2 * F3) + (F1 * F2 * F3))
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
}
  .fn(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a, b1, b2, b3)
}

#' ia_dl_tri: ternary (3-chemical) mixture predictor (ported from IA_DL)
#'
#' NOT verbatim: the original port used the deviation term `(a + b*P50)` for the
#' C1+C3, C2+C3 and full-triple subsets, which contradicts Jonker et al. (2005)
#' Eq. 13 (`G = a*(1 - b_DL*P)*prod(z)`) and the C1+C2 subset / binary `ia_dl_bi`.
#' Those six occurrences were corrected to `(a*(1 - b*P50))`. Pinned by
#' tests/testthat/test-mixture-predict.R against the unified `mix_response()`.
#'
#' @param c1,c2,c3 Concentrations.
#' @param max,slope1,slope2,slope3,ec50_1,ec50_2,ec50_3,a,b Model parameters.
#' @return Predicted response (scalar).
#' @keywords internal
ia_dl_tri <- function(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a, b) {
  .fn <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a, b) {
    F1 <- 0
    F2 <- 0
    F3 <- 0
    F4 <- 0
    z1 <- 0
    z2 <- 0
    z3 <- 0
    Trans <- 0
    P <- 0
    P50 <- 0
    if (C1 == 0 & C2 == 0 & C3 == 0) {
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            return(Max)
        }
        else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
            return(Max)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            return(0)
        }
        else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
            return(0)
        }
    }
    if (C1 > 0 & C2 == 0 & C3 == 0) {
        return(Max/(1 + (C1/Ec50_1)^Slope1))
    }
    if (C1 == 0 & C2 > 0 & C3 == 0) {
        return(Max/(1 + (C2/Ec50_2)^Slope2))
    }
    if (C1 == 0 & C2 == 0 & C3 > 0) {
        return(Max/(1 + (C3/Ec50_3)^Slope3))
    }
    if (C1 > 0 & C2 > 0 & C3 == 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        z1 <- TU1/(TU1 + TU2)
        z2 <- TU2/(TU1 + TU2)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        P50 <- 1 - (F1 * F2)
        F4 <- (a * (1 - (b * P50))) * z1 * z2
        if (Slope1 > 0 & Slope2 > 0) {
            Trans <- qnorm(F1 * F2)
            P <- pnorm(Trans + F4)
        }
        if (Slope1 < 0 & Slope2 < 0) {
            P50 <- (F1 + F2 - F1 * F2)
            F4 <- (a * (1 - (b * P50))) * z1 * z2
            Trans <- qnorm(F1 + F2 - F1 * F2)
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
    if (C1 > 0 & C2 == 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU3)
        z3 <- TU3/(TU1 + TU3)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        P50 <- 1 - (F1 * F3)
        F4 <- (a * (1 - (b * P50))) * z1 * z3
        if (Slope1 > 0 & Slope3 > 0) {
            Trans <- qnorm(F1 * F3)
            P <- pnorm(Trans + F4)
        }
        if (Slope1 < 0 & Slope3 < 0) {
            P50 <- (F1 + F3 - F1 * F3)
            F4 <- (a * (1 - (b * P50))) * z1 * z3
            Trans <- qnorm(F1 + F3 - F1 * F3)
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
    if (C1 == 0 & C2 > 0 & C3 > 0) {
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z2 <- TU2/(TU2 + TU3)
        z3 <- TU3/(TU2 + TU3)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        P50 <- 1 - (F2 * F3)
        F4 <- (a * (1 - (b * P50))) * z2 * z3
        if (Slope2 > 0 & Slope3 > 0) {
            Trans <- qnorm(F2 * F3)
            P <- pnorm(Trans + F4)
        }
        if (Slope2 < 0 & Slope3 < 0) {
            P50 <- (F2 + F3 - F2 * F3)
            F4 <- (a * (1 - (b * P50))) * z2 * z3
            Trans <- qnorm(F2 + F3 - F2 * F3)
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
    if (C1 > 0 & C2 > 0 & C3 > 0) {
        TU1 <- C1/Ec50_1
        TU2 <- C2/Ec50_2
        TU3 <- C3/Ec50_3
        z1 <- TU1/(TU1 + TU2 + TU3)
        z2 <- TU2/(TU1 + TU2 + TU3)
        z3 <- TU3/(TU1 + TU2 + TU3)
        F1 <- 1/(1 + (C1/Ec50_1)^Slope1)
        F2 <- 1/(1 + (C2/Ec50_2)^Slope2)
        F3 <- 1/(1 + (C3/Ec50_3)^Slope3)
        P50 <- 1 - (F1 * F2 * F3)
        F4 <- (a * (1 - (b * P50))) * z1 * z2 * z3
        if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
            Trans <- qnorm(F1 * F2 * F3)
            P <- pnorm(Trans + F4)
        }
        if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
            P50 <- (F1 + F2 + F3 - (F1 * F2) - (F1 * F3) - (F2 * 
                F3) + (F1 * F2 * F3))
            F4 <- (a * (1 - (b * P50))) * z1 * z2 * z3
            Trans <- qnorm(F1 + F2 + F3 - (F1 * F2) - (F1 * F3) - 
                (F2 * F3) + (F1 * F2 * F3))
            P <- pnorm(Trans - F4)
        }
        return(Max * P)
    }
}
  .fn(c1, c2, c3, max, slope1, slope2, slope3, ec50_1, ec50_2, ec50_3, a, b)
}

# Vectorised wrappers over the c1/c2/c3 concentration vectors.
ca_tri_vec <- Vectorize(ca_tri, vectorize.args = c("c1", "c2", "c3"))
ia_tri_vec <- Vectorize(ia_tri, vectorize.args = c("c1", "c2", "c3"))
ca_sa_tri_vec <- Vectorize(ca_sa_tri, vectorize.args = c("c1", "c2", "c3"))
ca_dr_tri_vec <- Vectorize(ca_dr_tri, vectorize.args = c("c1", "c2", "c3"))
ca_dl_tri_vec <- Vectorize(ca_dl_tri, vectorize.args = c("c1", "c2", "c3"))
ia_sa_tri_vec <- Vectorize(ia_sa_tri, vectorize.args = c("c1", "c2", "c3"))
ia_dr_tri_vec <- Vectorize(ia_dr_tri, vectorize.args = c("c1", "c2", "c3"))
ia_dl_tri_vec <- Vectorize(ia_dl_tri, vectorize.args = c("c1", "c2", "c3"))

