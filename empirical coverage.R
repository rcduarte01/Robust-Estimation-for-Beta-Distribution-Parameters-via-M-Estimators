cub_empirico<- function(a,b, x0=0, cont=0, alpha=0){
  #x0: contamination
  #cont: percent of contamination
  MVa <-  0
  MVb <-  0
  
  MIa <-  0
  MIb <-  0
  
  MDPa <- 0
  MDPb <- 0
  
  for(i in 1:1000){
    skip_to_next <- FALSE
    tryCatch({
    muestra <- rbeta(100, shape1 = a, shape2 = b)
    if(cont!=0){
      muestra[1:cont]<- rep(x0,cont)
    }
    
    
    ajusteMV <- confint(fitdist(muestra,"beta"))
    ajusteMI <- Ajuste_MI_estimador_beta(muestra)[,3:4]
    ajusteMDP<- Ajuste_MDPE_estimador_beta(muestra, alpha)[,3:4]
    
    }, error = function(e) {skip_to_next <<- TRUE },
    message = function(m) {skip_to_next <<- TRUE })
    
    if(skip_to_next) { next }
    
    MVa <- MVa + ifelse(ajusteMV[1,1]<=a & ajusteMV[1,2]>=a,1,0)
    MVb <- MVb + ifelse(ajusteMV[2,1]<=b & ajusteMV[2,2]>=b,1,0)
    
    MIa <- MIa + ifelse(ajusteMI[1,1]<=a & ajusteMI[1,2]>=a,1,0)
    MIb <- MIb + ifelse(ajusteMI[2,1]<=b & ajusteMI[2,2]>=b,1,0)
    
    MDPa <- MDPa + ifelse(ajusteMDP[1,1]<=a & ajusteMDP[1,2]>=a,1,0)
    MDPb <- MDPb + ifelse(ajusteMDP[2,1]<=b & ajusteMDP[2,2]>=b,1,0)
  }
  
  return(c(MVa,MVb,MIa,MIb,MDPa,MDPb)/1000)
  
}

set.seed(2022)
cub_empirico(a=2,
             b=6,
             x0=0.9,
             cont =10,
             alpha = 0.3)



