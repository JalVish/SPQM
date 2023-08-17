precip_BC_threshBinarySel_function <- function(observations,simulations,thresh) {

  library(fitdistrplus)
  library(Kendall)
  ##### prepare observation data
  pr = observations$data
  pr[pr<=thresh] = 0
  
  pr1 = pr
  pr1[pr1!=0] = (pr[pr!=0] - thresh)
  
  
  ex_dateH = observations$date
  date.dfH <- as.data.frame(ex_dateH)
  toIndH=tidyr::separate(data = date.dfH, col = ex_dateH, into = c("year","month","day"), sep = "-")
  toIndH[,4]= seq(1,length(ex_dateH), by=1)
  
  
  ##### prepare simulation data  
  chk = simulations$data
  chk[chk<0] = sqrt(sqrt(.Machine$double.eps))
  
  ex_date <- simulations$date
  date.df <- as.data.frame(ex_date)
  toInd=tidyr::separate(data = date.df, col = ex_date, into = c("year","month","day"), sep = "-")
  toInd[,4] = ((length(ex_dateH)+1):(length(ex_dateH)+dim(date.df)[1]))
  yr = unique(toInd[,1])
  
  
  #### quantile mapping starts - on monthly basis
  m1 = matrix(NA,nrow = length(chk), ncol = 1)
  neglected = matrix(NA,nrow = 12, ncol = 1)
  for (m in 1:12) {
    
    if (m < 10) {
      idx=which(toInd[,2]==paste('0',m,sep = ""))
      idxH=which(toIndH[,2]==paste('0',m,sep = ""))
    } else {
      idx=which(toInd[,2]==m)
      idxH=which(toIndH[,2]==m)
    }
    
    
    ####### match binary timeseries and get threshold
    chk2=pr[idxH]
    Bi_H = chk2
    Bi_H[Bi_H>0.2] = 1
    Bi_H[Bi_H<=0.2] = 0
    toIndH1=toIndH[idxH,4]
    
    pdryObs = length(chk2[chk2==0])/length(chk2)
    #remove(chk2)
    
    if (pdryObs==1 | (length(chk2)-length(chk2[chk2==0]))<30) {
      m1[idx] = NaN
    } else {
      
    
    
    
    sl = lm(Bi_H~toIndH1)
    val_end = sl$coefficients[2]*toIndH1[length(toIndH1)]+sl$coefficients[1]
    lineH = sl$coefficients[2]*toIndH1+sl$coefficients[1]
    remove(sl)
    
    chk2=chk[idx]
    val = seq(min(chk2),max(chk2),by=0.1)
    pval = matrix(NA,ncol = 1,nrow = length(val))
    for (t in 1:length(val)) {
      pval[t] = (length(chk2)-sum(chk2>val[t]))/length(chk2)
    }
    chk3 = chk2
    tem = approx(pval,val,pdryObs)
    chk3[chk3<tem$y] = 0
    chk31 = chk3[chk3>0]
    
    chk4 = chk3
    chk4[chk4>tem$y] = 1
    
    dateF = ex_date[idx]
    toInd1 = toInd[idx,4]
    
    sl = lm(chk4~toInd1)
    val_start = sl$coefficients[2]*toInd1[1]+sl$coefficients[1]
    
    er = abs(val_end-val_start)
    
    treshVal = data.frame(matrix(ncol = 2, nrow = 17))
    colnames(treshVal) <- c("t","e")
    treshVal$t[1] = tem$y
    treshVal$e[1] = er
    
    ii <- 1
    while (er >= 0.01) {
      new_mean = (mean(chk4)-(val_start-val_end))
      tem = approx(pval,val,(1-new_mean))
      chk3 = chk2
      chk3[chk3<tem$y] = 0
      
      chk31 = chk3[chk3>0]
      
      chk4 = chk3
      chk4[chk4>tem$y] = 1
      
      dateF = ex_date[idx]
      toInd1 = toInd[idx,4]
      
      sl = lm(chk4~toInd1)
      val_start = sl$coefficients[2]*toInd1[1]+sl$coefficients[1]
      er = abs(round(val_end,3)-round(val_start,3))
      ii=ii+1
      treshVal$t[ii] = tem$y
      treshVal$e[ii] = er
      if (ii == 16) {
        break
      }
    } 
    tv = na.omit(treshVal)
    threshF = tv$t[tv$e == min(tv$e)]
    
    ####### remove trend (if significant) in observations and fit distribution
    
    prJ = pr1[idxH]
    prJNZ = prJ[prJ!=0]
    
    mkTest = MannKendall(prJNZ)   #MK Test for significance of the trend
    if (mkTest$sl > 0.05) {
      paramG = tryCatch(
        {
          paramG = fitdist(prJNZ, "gamma")
          
        },
        error=function(e) {
          paramG = fitdist(prJNZ, "gamma",method = "mse")
          return(paramG)
        }
      )
      
      paramG = fitdist(prJNZ, "gamma",method = "mse")
    } else {
      xH = toIndH[idxH,4]
      xH1 = xH[prJ!=0]
      slH = lm(prJNZ ~ xH1)
      lineHist = xH1*slH$coefficients[2]+slH$coefficients[1]
      prJNZ_nt = prJNZ/lineHist
      
      neglected = sum(prJNZ_nt<0)
      
      if ((neglected>0)) {
        prJNZ_nt1 = ifelse(lineHist < 0, prJNZ, prJNZ / lineHist)
        paramG = fitdist(prJNZ_nt1, "gamma") 
        
      } else {
        paramG = fitdist(prJNZ_nt, "gamma")
      }
    }
    
    ######## remove trend in the projections
    chkk2 = chk[idx]
    xF = toInd[idx,4]
    slF = lm(chkk2 ~ xF)
    mkTest = MannKendall(chkk2)   #MK Test for significance of the trend
    if (mkTest$sl > 0.05) {
      
      lineFut = matrix(1,nrow = length(chkk2))
      chk2 = chkk2/lineFut
      
      lineFut1 = matrix(0,nrow = length(chkk2))
      
      chk3 = chk2[chk2>threshF] ## nonzero values
      
      loc <- match(sort(chk3), chk2, nomatch = 0) ### location of non-zero values
    } else {
      
      
      lineFut = predict(slF)
      chk2 = chkk2/lineFut
      
      lineFut1 = lineFut-lineFut[1]
      
      chk3 = chk2[chk2>threshF] ## nonzero values
      
      loc <- match(sort(chk3), chk2, nomatch = 0) ### location of non-zero values
      
    }
    
    
    ###### SPQM
    tem = sort(rank(chk3)/(length(chk3)+1))
    tem1 = qgamma(tem,paramG$estimate[1],paramG$estimate[2])+0.2
    tem2 = matrix(0,nrow = length(chk2),ncol = 1)
    tem2[loc] = tem1
    
    
    slop_obs = lm(pr[idxH]~ex_dateH[idxH])
    slop_pre = lm(tem2~ex_date[idx])
    
    valObs = predict(slop_obs)
    valPre = predict(slop_pre)
    
    ratio = valObs[length(valObs)]/valPre[1]
    ratio1 = ratio+lineFut1
    
    tem21 = tem2
    tem21[tem21!=0] = (tem2[tem2!=0])*ratio1[tem2!=0]
    
  
    m1[idx] = tem21
    
    }
    
  }
  
  ### adjustment for negative values and unusually high values
  m1[m1<0]=0
  t = max(simulations$data)-median(simulations$data)
  maxV = which(m1>10*t) # based on MonteCarlo Simulations (MonteCarloSimulations_maxThreshold.R)
  
  if (length(maxV)>0){
    m2=m1
    m2[m2>10*t] = NA
    m1 = (na.approx(m2, rule = 2))
  }
  
  
  #output = list(dat = data.frame(date = ex_date, data = m1), neglected=neglected)
  output = data.frame(date = ex_date, data = m1)
  return(output)
}

