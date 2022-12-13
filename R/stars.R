
#' flr2stars()
#' @param object of class FLStockR  
#' @param quantities default is 95CIs as c(0.025,0.975)
#' @return STARS list with $timeseris and $refpts
#' @export
flr2stars <- function(object,quantiles = c(0.025,0.975)){
  
  
    x= stockMedians(object)
    endyr = dims(x)$maxyear
    refpts = round(rbind(x@refpts,FLPar(
                     Fcur = an(fbar(simplify(x))[,ac(endyr)]),
                     Bcur=an(unitSums(ssb(x)[,ac(endyr)])),
                     B0.33 = quantile(an(unitSums(ssb(x))),0.33,na.rm=T),
                     B0.66 = quantile(an(unitSums(ssb(x))),0.66,na.rm=T))),3)
    
    refpts = data.frame(RefPoint=row.names(refpts),Value=as.data.frame(refpts[,1])$data)
    rownames(refpts) = 1:nrow(refpts)
  
    if(dim(object)[6]==1){  

    timeseries =  data.frame(Year=dims(x)$minyear:dims(x)$maxyear,
                             Rec_lower=NA,
                             Rec=round(an(rec(x)),0),Rec_upper=NA,
                             SSB_lower=NA,
                             SSB=round(an(ssb(x)),1),
                             SSB_upper=NA,
                             
                             TSB_lower=NA,
                             TSB=round(an(computeStock(x)),1),
                             TSB_upper=NA,
                             
                             Catches=round(an(catch(x)),2),
                             Landings = round(an(landings(x)),2),
                             Discards = round(an(discards(x)),2),
                             
                             F_lower=NA,
                             F = round(an(fbar(x)),3),
                             F_upper=NA,
                             
                             Fishing2_lower=NA,
                             Fishing2=NA,
                             Fishing2_upper=NA)
    } 
    
    if(dim(object)[6]>1){
      
      timeseries =  data.frame(Year=dims(object)$minyear:dims(object)$maxyear,
                               
                               Rec_lower=round(an(quantile(rec(object),quantiles[1])),0),
                               Rec=round(an(quantile(rec(object),0.5)),0),
                               Rec_upper=round(an(quantile(rec(object),quantiles[2])),0),
                               SSB_lower=round(an(quantile(ssb(object),quantiles[1])),1),
                               SSB=round(an(quantile(ssb(object),0.5)),1),
                               SSB_upper=round(an(quantile(ssb(object),quantiles[2])),1),
                               
                               TSB_lower=an(round(quantile(computeStock(object),quantiles[1]),1)),
                               TSB=an(round(quantile(computeStock(object),0.5),1)),
                               TSB_upper=an(round(quantile(computeStock(object),quantiles[2]),1)),
                               
                               Catches=round(an(quantile(catch(object),0.5)),2),
                               Landings = an(round(quantile(landings(object),0.5),2)),
                               Discards = an(round(quantile(discards(object),0.5),2)),
                               
                               F_lower=an(round(quantile(fbar(object),quantiles[1]),3)),
                               F = round(an(quantile(fbar(object),0.5)),3),
                               F_upper=an(round(quantile(fbar(object),quantiles[2]),3)),
                               
                               Fishing2_lower=NA,
                               Fishing2=NA,
                               Fishing2_upper=NA)   
      } 
    
    
    return(list(timeseries=timeseries,refpts=refpts))
  } 
  


#' ss2stars()
#' @param mvln output of ssmvln() 
#' @param output choice c("iters","mle")[1]
#' @param quantiles default is 95CIs as c(0.025,0.975)
#' @return STARS list with $timeseris and $refpts
#' @export
ss2stars <- function(mvln,output=c("iters","mle")[1],quantiles = c(0.025,0.975)){
  kbinp = FALSE
  if(is.null(mvln$kb)){
    if(output=="mle")
      stop("Output option mle requires to load the mvln object, not only $kb")
    kbinp=TRUE
    kb = mvln
    
  }  else {
    kb = mvln$kb
    mle = mvln$mle
  }
  
  
  if(output=="iters"){
    quants=c("SSB","F","Catch","stock","harvest","Recr")
    
    mu = aggregate(cbind(stock,harvest,SSB,F,Catch,Recr)~year+run,kb,
                   quantile,c(0.5,quantiles))
    
    timeseries =  data.frame(year=mu$year,
                             Rec_lower=mu[,quants[6]][,2],
                             Rec=mu[,quants[6]][,1],
                             Rec_upper=mu[,quants[6]][,3],
                             SSB_lower=mu[,quants[1]][,2],
                             SSB=mu[,quants[1]][,1],
                             SSB_upper=mu[,quants[1]][,3],
                             
                             Bratio_lower=mu[,quants[4]][,2],
                             Bratio=mu[,quants[4]][,1],
                             Bratio_upper=mu[,quants[4]][,3],
                             
                             Catches=mu[,quants[3]][,1],
                             Landings=NA,
                             Discards=NA,
                             
                             F_lower=mu[,quants[2]][,2],
                             F=mu[,quants[2]][,1],
                             F_upper=mu[,quants[2]][,3],
                             
                             Fratio_lower=mu[,quants[5]][,2],
                             Fratio=mu[,quants[5]][,1],
                             Fratio_upper=mu[,quants[5]][,3]
    )
    timeseries[-1] =     round(timeseries[-1],3)  
    endyr = which(mu$year==max(mu$year))
    
    
    refpts = round(t(data.frame(
      Ftgt = median(kb$F/kb$harvest),
      Btgt = median(kb$SSB/kb$stock),
      Bthr= NA,
      Blim= NA,
      Fcur  = timeseries$F[endyr],
      Bcur  = timeseries$SSB[endyr],
      B0.33= quantile(timeseries$SSB,0.33),
      B0.66= quantile(timeseries$SSB,0.66))),3)
    
  }
  if(output=="mle"){
    
    
    timeseries =  data.frame(year=mle$year,
                             Rec_lower=NA,
                             Rec=mle$Recr,
                             Rec_upper=NA,
                             SSB_lower=,
                             SSB=mle$SSB,
                             SSB_upper=NA,
                             
                             Bratio_lower=NA,
                             Bratio=mle$stock,
                             Bratio_upper=NA,
                             
                             Catches=mle$Catch,
                             Landings=NA,
                             Discards=NA,
                             
                             F_lower=NA,
                             F=mle$F,
                             F_upper=NA,
                             
                             Fratio_lower=NA,
                             Fratio=mle$harvest,
                             Fratio_upper=NA
    )
    
    
    timeseries[-1] =     round(timeseries[-1],3)  
    endyr = which(mle$year==max(mle$year))                  
    
    refpts = round(t(data.frame(
      Ftgt = mvln$refpts[1,2],
      Btgt = mvln$refpts[2,2],
      Bthr= NA,
      Blim= NA,
      Fcur  = timeseries$F[endyr],
      Bcur  = timeseries$SSB[endyr],
      B0.33= quantile(timeseries$SSB,0.33,na.rm=T),
      B0.66= quantile(timeseries$SSB,0.66,na.rm=T))),3)
  }
  
  refpts = data.frame(RefPoint=row.names(refpts),Value=refpts[,1])
  rownames(refpts) = 1:nrow(refpts)
  
  return(list(timeseries=timeseries,refpts=refpts))
}


  
  #' jabba2stars()
  #' @param jabba fit from JABBA fit_jabba() or jabba$kbtrj 
  #' @param quantiles default is 95CIs as c(0.025,0.975)
  #' @param bfrac biomass fraction of Bmsy, default 0.3Bmsy (ICES) 
  #' @return STARS list with $timeseris and $refpts
  #' @export
  jabba2stars <- function(jabba,quantiles = c(0.025,0.975),
                          bfrac = 0.3){
    kbinp = FALSE
    if(is.null(jabba$assessment)){
      kbinp=TRUE
      kb = jabba$kbtrj
    }  else {
      kb = jabba$kbtrj
    }
    quants=c("B","H","Catch","stock","harvest")
    mu = aggregate(cbind(stock,harvest,B,H,Catch)~year+run,kb,
                    quantile,c(0.5,quantiles))
    
    timeseries =  data.frame(year=mu$year,
                             Rec_lower=NA,Rec=NA,Rec_upper=NA,
                             Biomass_lower=mu[,quants[1]][,2],
                             Biomass=mu[,quants[1]][,1],
                             Biomass_upper=mu[,quants[1]][,3],
                             
                             Bratio_lower=mu[,quants[4]][,2],
                             Bratio=mu[,quants[4]][,1],
                             Bratio_upper=mu[,quants[4]][,3],
                             
                             Catches=mu[,quants[3]][,1],
                             Landings=NA,
                             Discards=NA,
                             
                            F_lower=mu[,quants[2]][,2],
                            F=mu[,quants[2]][,1],
                            F_upper=mu[,quants[2]][,3],
                             
                            Fratio_lower=mu[,quants[5]][,2],
                            Fratio=mu[,quants[5]][,1],
                            Fratio_upper=mu[,quants[5]][,3]
                            )
                            
    timeseries[-1] =     round(timeseries[-1],3)                     
    
     endyr = which(mu$year==max(mu$year))
   
   
      refpts = round(t(data.frame(
        Ftgt = median(kb$H/kb$harvest),
        Btgt = median(kb$B/kb$stock),
        Bthr= NA,
        Blim= median(bfrac*kb$B/kb$stock),
        Fcur  = timeseries$F[endyr],
        Bcur  = timeseries$Biomass[endyr],
        B0.33= quantile(timeseries$Biomass,0.33),
        B0.66= quantile(timeseries$Biomass,0.66))),3)
      
      refpts = data.frame(RefPoint=row.names(refpts),Value=refpts[,1])
     rownames(refpts) = 1:nrow(refpts)
      
      return(list(timeseries=timeseries,refpts=refpts))
  }
  
  