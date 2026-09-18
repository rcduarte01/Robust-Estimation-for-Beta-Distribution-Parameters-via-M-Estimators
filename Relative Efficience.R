library(MASS)

MSE <- function(theta, theta_estimado){
  
  mean((theta - theta_estimado)^2, na.rm = TRUE)
  
}
eficiencia_MCSE <- function(theta, est_MLE, est_rob){
  ok <- complete.cases(est_MLE, est_rob)
  
  est_MLE <- est_MLE[ok]
  est_rob <- est_rob[ok]
  
  R <- length(est_MLE)
  
  if(R < 2){
    return(c(
      MSE_MLE = NA,
      MSE_Robust = NA,
      RE = NA,
      MCSE = NA,
      R = R
    ))
  }
  
  
  q_MLE <- (est_MLE - theta)^2
  q_rob <- (est_rob - theta)^2
  
  
  mse_MLE <- mean(q_MLE)
  mse_rob <- mean(q_rob)
  
  
  RE <- mse_MLE / mse_rob
  
  
  S <- cov(cbind(q_MLE, q_rob))
  
  grad <- c(
    1 / mse_rob,
    -mse_MLE / mse_rob^2
  )
  

  MCSE <- sqrt(
    drop(t(grad) %*% S %*% grad) / R
  )
  
  
  return(c(
    MSE_MLE = mse_MLE,
    MSE_Robust = mse_rob,
    RE = RE,
    MCSE = MCSE,
    R = R
  ))
}



# ============================================================
# Simulation of Relative Efficience
# ============================================================

eficienciaRelativa <- function(n, N, a, b, rob = TRUE, alpha = 0){
  
  MIEst  <- matrix(NA_real_, nrow = N, ncol = 2)
  MVEst  <- matrix(NA_real_, nrow = N, ncol = 2)
  MDPEst <- matrix(NA_real_, nrow = N, ncol = 2)
  
  colnames(MIEst)  <- c("a", "b")
  colnames(MVEst)  <- c("a", "b")
  colnames(MDPEst) <- c("a", "b")
  
  
  # ----------------------------------------------------------
  # Monte Carlo
  # ----------------------------------------------------------
  
  for(i in seq_len(N)){
    
    skip_to_next <- FALSE
    
    tryCatch({
      mm <- rbeta(
        n,
        shape1 = a,
        shape2 = b
      )
      
      if(rob){
        
        xbar <- median(mm)
        s2   <- mad(mm)^2
        
      }else{
        
        xbar <- mean(mm)
        s2   <- var(mm)
        
      }
      
      
      alpha0 <- -xbar * (-xbar + xbar^2 + s2) / s2
      
      beta0 <- (alpha0 - alpha0*xbar) / xbar
      
      
      # ------------------------------------------------------
      # MI Estimator
      # ------------------------------------------------------
      
      MIEst[i, ] <- MI_estimador_beta(
        mm,
        rob
      )
      
      
      # ------------------------------------------------------
      # MDP Estimator
      # ------------------------------------------------------
      
      MDPEst[i, ] <- MDPDE_estimador_beta(
        mm,
        alpha,
        rob
      )
      
      
      # ------------------------------------------------------
      # MLE estimator
      # ------------------------------------------------------
      
      MVEst[i, ] <- fitdistr(
        mm,
        densfun = "beta",
        start = list(
          shape1 = alpha0,
          shape2 = beta0
        )
      )$estimate
      
      
    },
    
    error = function(e){
      skip_to_next <<- TRUE
    },
    
    message = function(m){
      skip_to_next <<- TRUE
    })
    
    
    if(skip_to_next){
      next
    }
    
  }
  

  
  MI_a <- eficiencia_MCSE(
    theta   = a,
    est_MLE = MVEst[, "a"],
    est_rob = MIEst[, "a"]
  )
  
  
  MI_b <- eficiencia_MCSE(
    theta   = b,
    est_MLE = MVEst[, "b"],
    est_rob = MIEst[, "b"]
  )
  
  
  
  MDP_a <- eficiencia_MCSE(
    theta   = a,
    est_MLE = MVEst[, "a"],
    est_rob = MDPEst[, "a"]
  )
  
  MDP_b <- eficiencia_MCSE(
    theta   = b,
    est_MLE = MVEst[, "b"],
    est_rob = MDPEst[, "b"]
  )
  
  resultados <- data.frame(
    
    n = rep(n, 4),
    
    a = rep(a, 4),
    
    b = rep(b, 4),
    
    Estimator = c(
      "MI",
      "MI",
      "MDP",
      "MDP"
    ),
    
    Parameter = c(
      "a",
      "b",
      "a",
      "b"
    ),
    
    RelativeEfficiency = round(c(
      MI_a["RE"],
      MI_b["RE"],
      MDP_a["RE"],
      MDP_b["RE"]
    ),4),
    
    MCSE = round(c(
      MI_a["MCSE"],
      MI_b["MCSE"],
      MDP_a["MCSE"],
      MDP_b["MCSE"]
    ),4),
    
    EffectiveReplications = c(
      MI_a["R"],
      MI_b["R"],
      MDP_a["R"],
      MDP_b["R"]
    )
    
  )
  
  

  rownames(resultados) <- NULL
  
  return(resultados)
  
}

#### SIMULATIONS
set.seed(2022)
eficienciaRelativa(
  n = 500,
  N = 1000,
  a = 1,
  b = 2,
  rob = FALSE,
  alpha = 0.4
)









