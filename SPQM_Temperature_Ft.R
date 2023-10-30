SPQM_Temperature_Ft <- function(observations,simulations) {

  library(zoo)
  library(fGarch)
  
  ##### prepare observation data
  tmaxOb = observations$tmax
  tminOb = observations$tmin
  
  
  ex_dateH = observations$date
  date.dfH <- as.data.frame(ex_dateH)
  toIndObs=tidyr::separate(data = date.dfH, col = ex_dateH, into = c("year","month","day"), sep = "-")
  toIndObs[,4] = (1:length(tmaxOb))
  
  
  ##### prepare simulation data  
  tmaxSim = simulations$tmax
  tminSim = simulations$tmin
  
  # tmaxSim = tmeanSim+trangeSim/2
  # tminSim = tmeanSim-trangeSim/2
  
  ex_date <- simulations$date
  date.df <- as.data.frame(ex_date)
  toInd=tidyr::separate(data = date.df, col = ex_date, into = c("year","month","day"), sep = "-")
  toInd[,4] = (1:length(tmaxSim))
  yr = unique(toInd[,1])
  
  
  
 
  
  #### quantile mapping starts - on monthly basis
  #################### tmax
  biasCor = matrix(NA, nrow = length(tmaxSim),ncol = 1)
  for (m in 1:12) {
    
    if (m < 10) {
      idx = which(toInd[,2]==paste('0',m,sep = ""))
      idxO = which(toIndObs[,2]==paste('0',m,sep = ""))
      
    } else {
      idx = which(toInd[,2]==m)
      idxO = which(toIndObs[,2]==m)
      
    }
    
    chk1 = tmaxSim[idx]
    t = toInd[idx,4]
    linRes=lm(as.array(chk1)~t)
    
    tt = t-t[1]
    temL1 = (coef(linRes)[[2]]*tt)
    
    chk2 = (chk1 - temL1)
    
    tem=rank(chk2)/(length(chk2)+1)
    
    
    oM = tmaxOb[idxO]
    tO = toIndObs[idxO,4]
    linResO=lm(as.array(oM)~tO)
    
    ttO = tO-tO[1]
    temL1O = (coef(linResO)[[2]]*ttO)
    
    oM1 = (oM - temL1O)
    snormPar = snormFit(oM1)
    
    tem1=qsnorm(tem,snormPar$par[1],snormPar$par[2],snormPar$par[3])
    
    tem2 = tem1 + temL1
    
    biasCor[idx]=tem2
    
    
  }
  biasCor1 = (na.approx(biasCor, rule = 2))
  
  x <- (1:length(tmaxSim))
  linResChk=lm(as.array(biasCor1)~x)
  y <- predict(linResChk,data.frame(x))
  
  xO <- (1:length(tmaxOb))
  linResChkO=lm(as.array(tmaxOb)~xO)
  yO <- predict(linResChkO,data.frame(xO))
  
  deltaTr = as.numeric(yO[length(yO)]) - y[1]
  
  biasCorTmax = biasCor1 + deltaTr
  
  
  #################### tmin
  biasCor = matrix(NA, nrow = length(tminSim),ncol = 1)
  for (m in 1:12) {
    
    if (m < 10) {
      idx = which(toInd[,2]==paste('0',m,sep = ""))
      idxO = which(toIndObs[,2]==paste('0',m,sep = ""))
      
    } else {
      idx = which(toInd[,2]==m)
      idxO = which(toIndObs[,2]==m)
      
    }
    
    chk1 = tminSim[idx]
    t = toInd[idx,4]
    linRes=lm(as.array(chk1)~t)
    
    tt = t-t[1]
    temL1 = (coef(linRes)[[2]]*tt)
    
    chk2 = (chk1 - temL1)
    
    tem=rank(chk2)/(length(chk2)+1)
    
    
    oM = tminOb[idxO]
    tO = toIndObs[idxO,4]
    linResO=lm(as.array(oM)~tO)
    
    ttO = tO-tO[1]
    temL1O = (coef(linResO)[[2]]*ttO)
    
    oM1 = (oM - temL1O)
    snormPar = snormFit(oM1)
    
    tem1=qsnorm(tem,snormPar$par[1],snormPar$par[2],snormPar$par[3])
    
    tem2 = tem1 + temL1
    
    biasCor[idx]=tem2
    
    
  }
  biasCor1 = (na.approx(biasCor, rule = 2))
  
  x <- (1:length(tminSim))
  linResChk=lm(as.array(biasCor1)~x)
  y <- predict(linResChk,data.frame(x))
  
  xO <- (1:length(tminOb))
  linResChkO=lm(as.array(tminOb)~xO)
  yO <- predict(linResChkO,data.frame(xO))
  
  deltaTr = as.numeric(yO[length(yO)]) - y[1]
  biasCorTmin = biasCor1 + deltaTr
  
  biasCorTmean = (biasCorTmin + biasCorTmax)/2
  biasCorTrange = abs(biasCorTmax - biasCorTmin)
  
  output = data.frame(date = ex_date, tmean = biasCorTmean, trange = biasCorTrange)
  return(output)
}

