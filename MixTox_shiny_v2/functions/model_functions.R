# Load needed packages
library(plotly)
library(readxl)
library(dplyr)
library(ggplot2)
library(Cairo)
library(grid)

##### Functions to add identifier columns #######

############# Run this user defined functions ###########
###### CA (concentration addition) calculation #####
CA <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3) {
  
  # Initialize variables
  L <- 0
  U <- Max
  Y <- 0
  ec1 <- 0
  ec2 <- 0
  ec3 <- 0
  G <- 0
  
  # Case when all concentrations are 0
  if (C1 == 0 & C2 == 0 & C3 == 0) {
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      return(Max)
    } else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
      return(Max)
    } 
    else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      return(0)
    } else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0 & C3 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0 & C3 == 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when only C3 > 0
  if (C1 == 0 & C2 == 0 & C3 > 0) {
    return(Max / (1 + (C3 / Ec50_3) ^ Slope3))
  }
  
  # Case when C1 > 0 and C2 > 0, but C3 == 0
  if (C1 > 0 & C2 > 0 & C3 == 0) {
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      G <- (C1 / ec1) + (C2 / ec2) - 1
      if (Slope1 > 0 & Slope2 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope2 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C1 > 0 and C3 > 0, but C2 == 0
  if (C1 > 0 & C2 == 0 & C3 > 0) {
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      G <- (C1 / ec1) + (C3 / ec3) - 1
      if (Slope1 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C2 > 0 and C3 > 0, but C1 == 0
  if (C1 == 0 & C2 > 0 & C3 > 0) {
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      G <- (C2 / ec2) + (C3 / ec3) - 1
      if (Slope2 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope2 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C1 > 0, C2 > 0, and C3 > 0
  if (C1 > 0 & C2 > 0 & C3 > 0) {
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      G <- (C1 / ec1) + (C2 / ec2) + (C3 / ec3) - 1
      if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
}




###### CA_bi #####
CA_bi <- function(C1, C2, max, slope1, slope2, ec501, ec502) {
  if (length(C1) > 1 || length(C2) > 1 || length(max) > 1 || 
      length(slope1) > 1 || length(slope2) > 1 || 
      length(ec501) > 1 || length(ec502) > 1) {
    stop("CA_bi received vector inputs, expected scalars")
  }
  # Case when both concentrations are 0
  if (C1 == 0 && C2 == 0) {
    if (slope1 > 0 & slope2 > 0) {
      return(max)
    }
    if (slope1 < 0 && slope2 < 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 && C2 == 0) {
    return(max / (1 + (C1 / ec501) ^ slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 && C2 > 0) {
    return(max / (1 + (C2 / ec502) ^ slope2))
  }
  
  # Case when both C1 > 0 and C2 > 0
  L <- 0
  U <- max
  for (i in 1:50) {  # Perform 50 iterations as in the VBA code
    Y <- (L + U) / 2
    if (Y <= 0 || Y >= max) {
      return(NaN)
    }
    ec1 <- ec501 * ((max - Y) / Y) ^ (1 / slope1)
    ec2 <- ec502 * ((max - Y) / Y) ^ (1 / slope2)
    G <- (C1 / ec1) + (C2 / ec2) - 1
    
    if (slope1 > 0 && slope2 > 0) {
      if (G < 0) {
        L <- Y
      } else {
        U <- Y
      }
    }
    
    if (slope1 < 0 && slope2 < 0) {
      if (G > 0) {
        L <- Y
      } else {
        U <- Y
      }
    }
  }
  
  return(Y)
}

###### CA_SA ######
CA_SA <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a) {
  
  # Initialize variables
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
  
  # Case when all concentrations are 0
  if (C1 == 0 & C2 == 0 & C3 == 0) {
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      return(Max)
    } else if(Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
      return(Max)
    }
    else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      return(0)
    } else if(Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0 & C3 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0 & C3 == 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when only C3 > 0
  if (C1 == 0 & C2 == 0 & C3 > 0) {
    return(Max / (1 + (C3 / Ec50_3) ^ Slope3))
  }
  
  # Case when C1 > 0 and C2 > 0, but C3 == 0
  if (C1 > 0 & C2 > 0 & C3 == 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    z1 <- TU1 / (TU1 + TU2)
    z2 <- TU2 / (TU1 + TU2)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    
    L <- 0
    U <- Max
    
    # Use binary search to converge to the correct value of Y
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      F4 <- exp(a * z1 * z2)
      G <- (C1 / ec1) + (C2 / ec2) - F4
      if (Slope1 > 0 & Slope2 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope2 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C1 > 0 and C3 > 0, but C2 == 0
  if (C1 > 0 & C2 == 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU3)
    z3 <- TU3 / (TU1 + TU3)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    
    L <- 0
    U <- Max
    
    # Use binary search to converge to the correct value of Y
    while (abs(U - L) > 1e-6) {
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F4 <- exp(a * z1 * z3)
      G <- (C1 / ec1) + (C3 / ec3) - F4
      if (Slope1 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C2 > 0 and C3 > 0, but C1 == 0
  if (C1 == 0 & C2 > 0 & C3 > 0) {
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z2 <- TU2 / (TU2 + TU3)
    z3 <- TU3 / (TU2 + TU3)
    
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    
    L <- 0
    U <- Max
    
    # Use binary search to converge to the correct value of Y
    while (abs(U - L) > 1e-6) {
      Y <- (L + U) / 2
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F4 <- exp(a * z2 * z3)
      G <- (C2 / ec2) + (C3 / ec3) - F4
      if (Slope2 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope2 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C1 > 0, C2 > 0, and C3 > 0
  if (C1 > 0 & C2 > 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU2 + TU3)
    z2 <- TU2 / (TU1 + TU2 + TU3)
    z3 <- TU3 / (TU1 + TU2 + TU3)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    
    L <- 0
    U <- Max
    
    # Use binary search to converge to the correct value of Y
    while (abs(U - L) > 1e-6) {
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F4 <- exp(a * z1 * z2 * z3)
      G <- (C1 / ec1) + (C2 / ec2) + (C3 / ec3) - F4
      if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
}


###### CA_SA_bi #####
CA_SA_bi <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a) {
  # Initialize variables
  TU1 <- C1 / Ec50_1
  TU2 <- C2 / Ec50_2
  z1 <- TU1 / (TU1 + TU2)
  z2 <- TU2 / (TU1 + TU2)
  
  # Case when both concentrations are 0
  if (C1 == 0 & C2 == 0) {
    if (Slope1 > 0 & Slope2 > 0) {
      return(Max)
    }
    if (Slope1 < 0 & Slope2 < 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when both C1 > 0 and C2 > 0
  L <- 0
  U <- Max
  for (i in 1:50) {  # Perform 50 iterations as in the VBA code
    Y <- (L + U) / 2
    ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
    ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
    F1 <- exp(a * z1 * z2)
    G <- (C1 / ec1) + (C2 / ec2) - F1
    
    if (Slope1 > 0 & Slope2 > 0) {
      if (G < 0) {
        L <- Y
      } else {
        U <- Y
      }
    }
    
    if (Slope1 < 0 & Slope2 < 0) {
      if (G > 0) {
        L <- Y
      } else {
        U <- Y
      }
    }
  }
  
  return(Y)
}


###### CA_Advanced S/A & Advanced S/A + S/A #####
CA_Tri_SA_Adv <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a1, a2, a3, a4) {
  
  # Initialize variables
  L <- 0
  U <- 0
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
  
  # Case when all concentrations are 0
  if (C1 == 0 & C2 == 0 & C3 == 0) {
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      return(Max)
    } else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
      return(Max)
    } 
    else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      return(0)
    } else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0 & C3 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0 & C3 == 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when only C3 > 0
  if (C1 == 0 & C2 == 0 & C3 > 0) {
    return(Max / (1 + (C3 / Ec50_3) ^ Slope3))
  }
  
  # Case when C1 > 0 and C2 > 0, but C3 == 0
  if (C1 > 0 & C2 > 0 & C3 == 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    z1 <- TU1 / (TU1 + TU2)
    z2 <- TU2 / (TU1 + TU2)
    
    L <- 0
    U <- Max
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      F1 <- exp(a1 * z1 * z2)
      G <- (C1 / ec1) + (C2 / ec2) - F1
      if (Slope1 > 0 & Slope2 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope2 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C1 > 0 and C3 > 0, but C2 == 0
  if (C1 > 0 & C2 == 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU3)
    z3 <- TU3 / (TU1 + TU3)
    
    L <- 0
    U <- Max
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F1 <- exp(a2 * z1 * z3)
      G <- (C1 / ec1) + (C3 / ec3) - F1
      if (Slope1 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C2 > 0 and C3 > 0, but C1 == 0
  if (C1 == 0 & C2 > 0 & C3 > 0) {
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z2 <- TU2 / (TU2 + TU3)
    z3 <- TU3 / (TU2 + TU3)
    
    L <- 0
    U <- Max
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F1 <- exp(a3 * z2 * z3)
      G <- (C2 / ec2) + (C3 / ec3) - F1
      if (Slope2 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope2 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C1 > 0, C2 > 0, and C3 > 0
  if (C1 > 0 & C2 > 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU2 + TU3)
    z2 <- TU2 / (TU1 + TU2 + TU3)
    z3 <- TU3 / (TU1 + TU2 + TU3)
    
    L <- 0
    U <- Max
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F1 <- exp(((((z1 * z2) / 0.25) * a1) + (((z1 * z3) / 0.25) * a2) + (((z2 * z3) / 0.25) * a3) + a4) * z1 * z2 * z3)
      G <- (C1 / ec1) + (C2 / ec2) + (C3 / ec3) - F1
      if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
}

###### CA_DR #####
CA_DR <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a, b1, b2, b3) {
  
  # Initialize variables
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
  
  # Case when all concentrations are 0
  if (C1 == 0 & C2 == 0 & C3 == 0) {
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      return(Max)
    } else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
      return(Max)
    }
    else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      return(0)
    } else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0 & C3 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0 & C3 == 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when only C3 > 0
  if (C1 == 0 & C2 == 0 & C3 > 0) {
    return(Max / (1 + (C3 / Ec50_3) ^ Slope3))
  }
  
  # Case when C1 > 0 and C2 > 0, but C3 == 0
  if (C1 > 0 & C2 > 0 & C3 == 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU2 + TU3)
    z2 <- TU2 / (TU1 + TU2 + TU3)
    z3 <- TU3 / (TU1 + TU2 + TU3)
    
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      F1 <- exp((a + (b1 * z1) + (b2 * z2)) * z1 * z2)
      G <- (C1 / ec1) + (C2 / ec2) - F1
      if (Slope1 > 0 & Slope2 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope2 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C1 > 0 and C3 > 0, but C2 == 0
  if (C1 > 0 & C2 == 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU2 + TU3)
    z2 <- TU2 / (TU1 + TU2 + TU3)
    z3 <- TU3 / (TU1 + TU2 + TU3)
    
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F1 <- exp((a + (b1 * z1) + (b3 * z3)) * z1 * z3)
      G <- (C1 / ec1) + (C3 / ec3) - F1
      if (Slope1 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C2 > 0 and C3 > 0, but C1 == 0
  if (C1 == 0 & C2 > 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU2 + TU3)
    z2 <- TU2 / (TU1 + TU2 + TU3)
    z3 <- TU3 / (TU1 + TU2 + TU3)
    
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F1 <- exp((a + (b2 * z2) + (b3 * z3)) * z2 * z3)
      G <- (C2 / ec2) + (C3 / ec3) - F1
      if (Slope2 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope2 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C1 > 0, C2 > 0, and C3 > 0
  if (C1 > 0 & C2 > 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU2 + TU3)
    z2 <- TU2 / (TU1 + TU2 + TU3)
    z3 <- TU3 / (TU1 + TU2 + TU3)
    
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F1 <- exp((a + (b1 * z1) + (b2 * z2) + (b3 * z3)) * z1 * z2 * z3)
      G <- (C1 / ec1) + (C2 / ec2) + (C3 / ec3) - F1
      if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
}

###### CA_DR_binary #####
CA_DR_bi <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a, b) {
  # Initialize variables
  TU1 <- C1 / Ec50_1
  TU2 <- C2 / Ec50_2
  z1 <- TU1 / (TU1 + TU2)
  z2 <- TU2 / (TU1 + TU2)
  
  # Case when both concentrations are 0
  if (C1 == 0 & C2 == 0) {
    if (Slope1 > 0 & Slope2 > 0) {
      return(Max)
    }
    if (Slope1 < 0 & Slope2 < 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when both C1 > 0 and C2 > 0
  L <- 0
  U <- Max
  for (i in 1:50) {  # Perform 50 iterations as in the VBA code
    Y <- (L + U) / 2
    ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
    ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
    F1 <- exp((a + b * z1) * z1 * z2)
    G <- (C1 / ec1) + (C2 / ec2) - F1
    
    if (Slope1 > 0 & Slope2 > 0) {
      if (G < 0) {
        L <- Y
      } else {
        U <- Y
      }
    }
    
    if (Slope1 < 0 & Slope2 < 0) {
      if (G > 0) {
        L <- Y
      } else {
        U <- Y
      }
    }
  }
  
  return(Y)
}
###### CA_DL #####
CA_DL <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a, b) {
  
  # Initialize variables
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
  
  # Case when all concentrations are 0
  if (C1 == 0 & C2 == 0 & C3 == 0) {
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      return(Max)
    } else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
      return(Max)
    }
    else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      return(0)
    } else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0 & C3 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0 & C3 == 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when only C3 > 0
  if (C1 == 0 & C2 == 0 & C3 > 0) {
    return(Max / (1 + (C3 / Ec50_3) ^ Slope3))
  }
  
  # Case when C1 > 0 and C2 > 0, but C3 == 0
  if (C1 > 0 & C2 > 0 & C3 == 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    z1 <- TU1 / (TU1 + TU2)
    z2 <- TU2 / (TU1 + TU2)
    
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      F1 <- exp((a * (1 - (b * (TU1 + TU2)))) * z1 * z2)
      G <- (C1 / ec1) + (C2 / ec2) - F1
      if (Slope1 > 0 & Slope2 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope2 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C1 > 0 and C3 > 0, but C2 == 0
  if (C1 > 0 & C2 == 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU3)
    z3 <- TU3 / (TU1 + TU3)
    
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F1 <- exp((a * (1 - (b * (TU1 + TU3)))) * z1 * z3)
      G <- (C1 / ec1) + (C3 / ec3) - F1
      if (Slope1 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C2 > 0 and C3 > 0, but C1 == 0
  if (C1 == 0 & C2 > 0 & C3 > 0) {
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z2 <- TU2 / (TU2 + TU3)
    z3 <- TU3 / (TU2 + TU3)
    
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F1 <- exp((a * (1 - (b * (TU2 + TU3)))) * z2 * z3)
      G <- (C2 / ec2) + (C3 / ec3) - F1
      if (Slope2 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope2 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
  
  # Case when C1 > 0, C2 > 0, and C3 > 0
  if (C1 > 0 & C2 > 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU2 + TU3)
    z2 <- TU2 / (TU1 + TU2 + TU3)
    z3 <- TU3 / (TU1 + TU2 + TU3)
    
    L <- 0
    U <- Max
    # Dynamically adapt loop based on convergence criteria
    while (abs(U - L) > 1e-6) {  # Stop when the difference between L and U is sufficiently small
      Y <- (L + U) / 2
      ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
      ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
      ec3 <- Ec50_3 * ((Max - Y) / Y) ^ (1 / Slope3)
      F1 <- exp((a * (1 - (b * (TU1 + TU2 + TU3)))) * z1 * z2 * z3)
      G <- (C1 / ec1) + (C2 / ec2) + (C3 / ec3) - F1
      if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
        if (G < 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
      if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
        if (G > 0) {
          L <- Y
        } else {
          U <- Y
        }
      }
    }
    return(Y)
  }
}



###### CA_DL_binary #####
CA_DL_bi <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a, b) {
  # Case when both concentrations are 0
  if (C1 == 0 & C2 == 0) {
    if (Slope1 > 0 & Slope2 > 0) {
      return(Max)
    }
    if (Slope1 < 0 & Slope2 < 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when both C1 > 0 and C2 > 0
  TU1 <- C1 / Ec50_1
  TU2 <- C2 / Ec50_2
  z1 <- TU1 / (TU1 + TU2)
  z2 <- TU2 / (TU1 + TU2)
  
  L <- 0
  U <- Max
  for (i in 1:50) {  # Perform 50 iterations as in the VBA code
    Y <- (L + U) / 2
    ec1 <- Ec50_1 * ((Max - Y) / Y) ^ (1 / Slope1)
    ec2 <- Ec50_2 * ((Max - Y) / Y) ^ (1 / Slope2)
    F1 <- exp((a * (1 - (b * (TU1 + TU2)))) * z1 * z2)
    G <- (C1 / ec1) + (C2 / ec2) - F1
    
    if (Slope1 > 0 & Slope2 > 0) {
      if (G < 0) {
        L <- Y
      } else {
        U <- Y
      }
    }
    
    if (Slope1 < 0 & Slope2 < 0) {
      if (G > 0) {
        L <- Y
      } else {
        U <- Y
      }
    }
  }
  
  return(Y)
}


###### IA #####
IA <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3) {
  
  # Initialize F1, F2, F3
  F1 <- 0
  F2 <- 0
  F3 <- 0
  
  # Case when all concentrations are 0
  if (C1 == 0 & C2 == 0 & C3 == 0) {
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      return(Max)
    } else if(Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
      return(Max)
    }
    else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      return(0)
    } else if (Slope1 < 0 & Slope2 <0 & Slope3 == 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0 & C3 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0 & C3 == 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when only C3 > 0
  if (C1 == 0 & C2 == 0 & C3 > 0) {
    return(Max / (1 + (C3 / Ec50_3) ^ Slope3))
  }
  
  # Case when C1 > 0 and C2 > 0, but C3 == 0
  if (C1 > 0 & C2 > 0 & C3 == 0) {
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    if (Slope1 > 0 & Slope2 > 0) {
      return(Max * F1 * F2)
    } else if (Slope1 < 0 & Slope2 < 0) {
      return(Max * (F1 + F2 - F1 * F2))
    }
  }
  
  # Case when C1 > 0 and C3 > 0, but C2 == 0
  if (C1 > 0 & C2 == 0 & C3 > 0) {
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    if (Slope1 > 0 & Slope3 > 0) {
      return(Max * F1 * F3)
    } else if (Slope1 < 0 & Slope3 < 0) {
      return(Max * (F1 + F3 - F1 * F3))
    }
  }
  
  # Case when C2 > 0 and C3 > 0, but C1 == 0
  if (C1 == 0 & C2 > 0 & C3 > 0) {
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    if (Slope2 > 0 & Slope3 > 0) {
      return(Max * F2 * F3)
    } else if (Slope2 < 0 & Slope3 < 0) {
      return(Max * (F2 + F3 - F2 * F3))
    }
  }
  
  # Case when C1 > 0, C2 > 0, and C3 > 0
  if (C1 > 0 & C2 > 0 & C3 > 0) {
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      return(Max * F1 * F2 * F3)
    } else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      return(Max * (F1 + F2 + F3 - (F1 * F2) - (F1 * F3) - (F2 * F3) + (F1 * F2 * F3)))
    }
  }
}

###### IA_bi ######
IA_bi <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2) {
  # Case when both concentrations are 0
  if (C1 == 0 & C2 == 0) {
    if (Slope1 > 0 & Slope2 > 0) {
      return(Max)
    }
    if (Slope1 < 0 & Slope2 < 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when both C1 > 0 and C2 > 0
  F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
  F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
  
  if (Slope1 > 0 & Slope2 > 0) {
    return(Max * F1 * F2)
  }
  
  if (Slope1 < 0 & Slope2 < 0) {
    return(Max * (F1 + F2 - F1 * F2))
  }
}

###### IA_SA #####
IA_SA <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a) {
  
  # Initialize F1, F2, F3, F4, z1, z2, z3
  F1 <- 0
  F2 <- 0
  F3 <- 0
  F4 <- 0
  z1 <- 0
  z2 <- 0
  z3 <- 0
  Trans <- 0
  P <- 0
  
  # Case when all concentrations are 0
  if (C1 == 0 & C2 == 0 & C3 == 0) {
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      return(Max)
    } else if (Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
      return(Max)
    }
    else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      return(0)
    } else if (Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0 & C3 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0 & C3 == 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when only C3 > 0
  if (C1 == 0 & C2 == 0 & C3 > 0) {
    return(Max / (1 + (C3 / Ec50_3) ^ Slope3))
  }
  
  # Case when C1 > 0 and C2 > 0, but C3 == 0
  if (C1 > 0 & C2 > 0 & C3 == 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    z1 <- TU1 / (TU1 + TU2)
    z2 <- TU2 / (TU1 + TU2)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F4 <- a * z1 * z2
    
    if (Slope1 > 0 & Slope2 > 0) {
      Trans <- qnorm(F1 * F2)  # Equivalent to NormSInv in Excel
      P <- pnorm(Trans + F4)   # Equivalent to NormSDist in Excel
    }
    if (Slope1 < 0 & Slope2 < 0) {
      Trans <- qnorm(F1 + F2 - F1 * F2)
      P <- pnorm(Trans - F4)
    }
    
    return(Max * P)
  }
  
  # Case when C1 > 0 and C3 > 0, but C2 == 0
  if (C1 > 0 & C2 == 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU3)
    z3 <- TU3 / (TU1 + TU3)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
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
  
  # Case when C2 > 0 and C3 > 0, but C1 == 0
  if (C1 == 0 & C2 > 0 & C3 > 0) {
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z2 <- TU2 / (TU2 + TU3)
    z3 <- TU3 / (TU2 + TU3)
    
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
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
  
  # Case when C1 > 0, C2 > 0, and C3 > 0
  if (C1 > 0 & C2 > 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU2 + TU3)
    z2 <- TU2 / (TU1 + TU2 + TU3)
    z3 <- TU3 / (TU1 + TU2 + TU3)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    F4 <- a * z1 * z2 * z3
    
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      Trans <- qnorm(F1 * F2 * F3)
      P <- pnorm(Trans + F4)
    }
    if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      Trans <- qnorm(F1 + F2 + F3 - (F1 * F2) - (F1 * F3) - (F2 * F3) + (F1 * F2 * F3))
      P <- pnorm(Trans - F4)
    }
    
    return(Max * P)
  }
}
###### IA_SA_bi ######
IA_SA_bi <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a) {
  # Case when both concentrations are 0
  if (C1 == 0 & C2 == 0) {
    if (Slope1 > 0 & Slope2 > 0) {
      return(Max)
    }
    if (Slope1 < 0 & Slope2 < 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when both C1 > 0 and C2 > 0
  TU1 <- C1 / Ec50_1
  TU2 <- C2 / Ec50_2
  z1 <- TU1 / (TU1 + TU2)
  z2 <- TU2 / (TU1 + TU2)
  
  F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
  F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
  F3 <- a * z1 * z2
  
  if (Slope1 > 0 & Slope2 > 0) {
    Trans <- qnorm(F1 * F2)  # Normal inverse (qnorm)
    P <- pnorm(Trans + F3)   # Normal cumulative distribution (pnorm)
  }
  
  if (Slope1 < 0 & Slope2 < 0) {
    Trans <- qnorm(F1 + F2 - F1 * F2)  # Normal inverse (qnorm)
    P <- pnorm(Trans - F3)   # Normal cumulative distribution (pnorm)
  }
  
  return(Max * P)
}

###### IA_DR #####
IA_DR <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a, b1, b2, b3) {
  
  # Initialize variables
  F1 <- 0
  F2 <- 0
  F3 <- 0
  F4 <- 0
  z1 <- 0
  z2 <- 0
  z3 <- 0
  Trans <- 0
  P <- 0
  
  # Case when all concentrations are 0
  if (C1 == 0 & C2 == 0 & C3 == 0) {
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      return(Max)
    } else if(Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
      return(Max)
    }
    else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      return(0)
    } else if(Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0 & C3 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0 & C3 == 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when only C3 > 0
  if (C1 == 0 & C2 == 0 & C3 > 0) {
    return(Max / (1 + (C3 / Ec50_3) ^ Slope3))
  }
  
  # Case when C1 > 0 and C2 > 0, but C3 == 0
  if (C1 > 0 & C2 > 0 & C3 == 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    z1 <- TU1 / (TU1 + TU2)
    z2 <- TU2 / (TU1 + TU2)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F4 <- (a + (b1 * z1) + (b2 * z2)) * z1 * z2
    
    if (Slope1 > 0 & Slope2 > 0) {
      Trans <- qnorm(F1 * F2)  # NormSInv in Excel
      P <- pnorm(Trans + F4)   # NormSDist in Excel
    }
    if (Slope1 < 0 & Slope2 < 0) {
      Trans <- qnorm(F1 + F2 - F1 * F2)
      P <- pnorm(Trans - F4)
    }
    
    return(Max * P)
  }
  
  # Case when C1 > 0 and C3 > 0, but C2 == 0
  if (C1 > 0 & C2 == 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU3)
    z3 <- TU3 / (TU1 + TU3)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
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
  
  # Case when C2 > 0 and C3 > 0, but C1 == 0
  if (C1 == 0 & C2 > 0 & C3 > 0) {
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z2 <- TU2 / (TU2 + TU3)
    z3 <- TU3 / (TU2 + TU3)
    
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
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
  
  # Case when C1 > 0, C2 > 0, and C3 > 0
  if (C1 > 0 & C2 > 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU2 + TU3)
    z2 <- TU2 / (TU1 + TU2 + TU3)
    z3 <- TU3 / (TU1 + TU2 + TU3)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    F4 <- (a + (b1 * z1) + (b2 * z2) + (b3 * z3)) * z1 * z2 * z3
    
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      Trans <- qnorm(F1 * F2 * F3)
      P <- pnorm(Trans + F4)
    }
    if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      Trans <- qnorm(F1 + F2 + F3 - (F1 * F2) - (F1 * F3) - (F2 * F3) + (F1 * F2 * F3))
      P <- pnorm(Trans - F4)
    }
    
    return(Max * P)
  }
}

###### IA_DR_bi ######
IA_DR_bi <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a, b) {
  # Case when both concentrations are 0
  if (C1 == 0 & C2 == 0) {
    if (Slope1 > 0 & Slope2 > 0) {
      return(Max)
    }
    if (Slope1 < 0 & Slope2 < 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when both C1 > 0 and C2 > 0
  TU1 <- C1 / Ec50_1
  TU2 <- C2 / Ec50_2
  z1 <- TU1 / (TU1 + TU2)
  z2 <- TU2 / (TU1 + TU2)
  
  F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
  F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
  F3 <- (a + b * z1) * z1 * z2
  
  if (Slope1 > 0 & Slope2 > 0) {
    Trans <- qnorm(F1 * F2)  # Normal inverse (qnorm)
    P <- pnorm(Trans + F3)   # Normal cumulative distribution (pnorm)
  }
  
  if (Slope1 < 0 & Slope2 < 0) {
    Trans <- qnorm(F1 + F2 - F1 * F2)  # Normal inverse (qnorm)
    P <- pnorm(Trans - F3)   # Normal cumulative distribution (pnorm)
  }
  
  return(Max * P)
}

###### IA_DL #####
IA_DL <- function(C1, C2, C3, Max, Slope1, Slope2, Slope3, Ec50_1, Ec50_2, Ec50_3, a, b) {
  
  # Initialize variables
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
  
  # Case when all concentrations are 0
  if (C1 == 0 & C2 == 0 & C3 == 0) {
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      return(Max)
    } else if(Slope1 > 0 & Slope2 > 0 & Slope3 == 0) {
      return(Max)
    }
    else if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      return(0)
    } else if(Slope1 < 0 & Slope2 < 0 & Slope3 == 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0 & C3 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0 & C3 == 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when only C3 > 0
  if (C1 == 0 & C2 == 0 & C3 > 0) {
    return(Max / (1 + (C3 / Ec50_3) ^ Slope3))
  }
  
  # Case when C1 > 0 and C2 > 0, but C3 == 0
  if (C1 > 0 & C2 > 0 & C3 == 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    z1 <- TU1 / (TU1 + TU2)
    z2 <- TU2 / (TU1 + TU2)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    
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
  
  # Case when C1 > 0 and C3 > 0, but C2 == 0
  if (C1 > 0 & C2 == 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU3)
    z3 <- TU3 / (TU1 + TU3)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    
    P50 <- 1 - (F1 * F3)
    F4 <- (a + (b * P50)) * z1 * z3
    
    if (Slope1 > 0 & Slope3 > 0) {
      Trans <- qnorm(F1 * F3)
      P <- pnorm(Trans + F4)
    }
    if (Slope1 < 0 & Slope3 < 0) {
      P50 <- (F1 + F3 - F1 * F3)
      F4 <- (a + (b * P50)) * z1 * z3
      Trans <- qnorm(F1 + F3 - F1 * F3)
      P <- pnorm(Trans - F4)
    }
    
    return(Max * P)
  }
  
  # Case when C2 > 0 and C3 > 0, but C1 == 0
  if (C1 == 0 & C2 > 0 & C3 > 0) {
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z2 <- TU2 / (TU2 + TU3)
    z3 <- TU3 / (TU2 + TU3)
    
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    
    P50 <- 1 - (F2 * F3)
    F4 <- (a + (b * P50)) * z2 * z3
    
    if (Slope2 > 0 & Slope3 > 0) {
      Trans <- qnorm(F2 * F3)
      P <- pnorm(Trans + F4)
    }
    if (Slope2 < 0 & Slope3 < 0) {
      P50 <- (F2 + F3 - F2 * F3)
      F4 <- (a + (b * P50)) * z2 * z3
      Trans <- qnorm(F2 + F3 - F2 * F3)
      P <- pnorm(Trans - F4)
    }
    
    return(Max * P)
  }
  
  # Case when C1 > 0, C2 > 0, and C3 > 0
  if (C1 > 0 & C2 > 0 & C3 > 0) {
    TU1 <- C1 / Ec50_1
    TU2 <- C2 / Ec50_2
    TU3 <- C3 / Ec50_3
    z1 <- TU1 / (TU1 + TU2 + TU3)
    z2 <- TU2 / (TU1 + TU2 + TU3)
    z3 <- TU3 / (TU1 + TU2 + TU3)
    
    F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
    F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
    F3 <- 1 / (1 + (C3 / Ec50_3) ^ Slope3)
    
    P50 <- 1 - (F1 * F2 * F3)
    F4 <- (a + (b * P50)) * z1 * z2 * z3
    
    if (Slope1 > 0 & Slope2 > 0 & Slope3 > 0) {
      Trans <- qnorm(F1 * F2 * F3)
      P <- pnorm(Trans + F4)
    }
    if (Slope1 < 0 & Slope2 < 0 & Slope3 < 0) {
      P50 <- (F1 + F2 + F3 - (F1 * F2) - (F1 * F3) - (F2 * F3) + (F1 * F2 * F3))
      F4 <- (a + (b * P50)) * z1 * z2 * z3
      Trans <- qnorm(F1 + F2 + F3 - (F1 * F2) - (F1 * F3) - (F2 * F3) + (F1 * F2 * F3))
      P <- pnorm(Trans - F4)
    }
    
    return(Max * P)
  }
}


###### IA_DL_bi ######
IA_DL_bi <- function(C1, C2, Max, Slope1, Slope2, Ec50_1, Ec50_2, a, b) {
  # Case when both concentrations are 0
  if (C1 == 0 & C2 == 0) {
    if (Slope1 > 0 & Slope2 > 0) {
      return(Max)
    }
    if (Slope1 < 0 & Slope2 < 0) {
      return(0)
    }
  }
  
  # Case when only C1 > 0
  if (C1 > 0 & C2 == 0) {
    return(Max / (1 + (C1 / Ec50_1) ^ Slope1))
  }
  
  # Case when only C2 > 0
  if (C1 == 0 & C2 > 0) {
    return(Max / (1 + (C2 / Ec50_2) ^ Slope2))
  }
  
  # Case when both C1 > 0 and C2 > 0
  TU1 <- C1 / Ec50_1
  TU2 <- C2 / Ec50_2
  z1 <- TU1 / (TU1 + TU2)
  z2 <- TU2 / (TU1 + TU2)
  
  F1 <- 1 / (1 + (C1 / Ec50_1) ^ Slope1)
  F2 <- 1 / (1 + (C2 / Ec50_2) ^ Slope2)
  
  if (Slope1 > 0 & Slope2 > 0) {
    P50 <- 1 - (F1 * F2)
    F3 <- (a * (1 - (b * P50))) * z1 * z2
    Trans <- qnorm(F1 * F2)  # Normal inverse (qnorm)
    P <- pnorm(Trans + F3)   # Normal cumulative distribution (pnorm)
  }
  
  if (Slope1 < 0 & Slope2 < 0) {
    P50 <- (F1 + F2 - F1 * F2)
    F3 <- (a * (1 - (b * P50))) * z1 * z2
    Trans <- qnorm(F1 + F2 - F1 * F2)  # Normal inverse (qnorm)
    P <- pnorm(Trans - F3)   # Normal cumulative distribution (pnorm)
  }
  
  return(Max * P)
}

