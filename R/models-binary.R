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

