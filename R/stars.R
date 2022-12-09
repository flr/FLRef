


  
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
                             Stock1_lower=mu[,quants[1]][,2],
                             Stock1=mu[,quants[1]][,1],
                             Stock1_upper=mu[,quants[1]][,3],
                             
                             Stock2_lower=mu[,quants[4]][,2],
                             Stock2=mu[,quants[4]][,1],
                             Stock2_upper=mu[,quants[4]][,3],
                             
                             Catches=mu[,quants[3]][,1],
                             Landings=NA,
                             Discards=NA,
                             
                            Fishing1_lower=mu[,quants[2]][,2],
                            Fishing1=mu[,quants[2]][,1],
                            Fishing1_upper=mu[,quants[2]][,3],
                             
                            Fishing2_lower=mu[,quants[5]][,2],
                            Fishing2=mu[,quants[5]][,1],
                            Fishing2_upper=mu[,quants[5]][,3]
                            )
                            
    timeseries[-1] =     round(timeseries[-1],3)                     
    
     endyr = which(mu$year==max(mu$year))
   
   
      refpts = round(t(data.frame(
        Ftgt = median(kb$H/kb$harvest),
        Btgt = median(kb$B/kb$stock),
        Bthr= NA,
        Blim= median(bfrac*kb$B/kb$stock),
        Fcur  = timeseries$Fishing1[endyr],
        Bcur  = timeseries$Stock1[endyr],
        B0.33= quantile(timeseries$Stock1,0.33),
        B0.66= quantile(timeseries$Stock1,0.66))),3)
      
      refpts = data.frame(RefPoint=row.names(refpts),Value=refpts[,1])
     rownames(refpts) = 1:nrow(refpts)
      
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
                             Stock1_lower=mu[,quants[1]][,2],
                             Stock1=mu[,quants[1]][,1],
                             Stock1_upper=mu[,quants[1]][,3],
                             
                             Stock2_lower=mu[,quants[4]][,2],
                             Stock2=mu[,quants[4]][,1],
                             Stock2_upper=mu[,quants[4]][,3],
                             
                             Catches=mu[,quants[3]][,1],
                             Landings=NA,
                             Discards=NA,
                             
                             Fishing1_lower=mu[,quants[2]][,2],
                             Fishing1=mu[,quants[2]][,1],
                             Fishing1_upper=mu[,quants[2]][,3],
                             
                             Fishing2_lower=mu[,quants[5]][,2],
                             Fishing2=mu[,quants[5]][,1],
                             Fishing2_upper=mu[,quants[5]][,3]
    )
    timeseries[-1] =     round(timeseries[-1],3)  
    endyr = which(mu$year==max(mu$year))
    
    
    refpts = round(t(data.frame(
      Ftgt = median(kb$F/kb$harvest),
      Btgt = median(kb$SSB/kb$stock),
      Bthr= NA,
      Blim= NA,
      Fcur  = timeseries$Fishing1[endyr],
      Bcur  = timeseries$Stock1[endyr],
      B0.33= quantile(timeseries$Stock1,0.33),
      B0.66= quantile(timeseries$Stock1,0.66))),3)
    
    }
    if(output=="mle"){
      
      
      timeseries =  data.frame(year=mle$year,
                               Rec_lower=NA,
                               Rec=mle$Recr,
                               Rec_upper=NA,
                               Stock1_lower=,
                               Stock1=mle$SSB,
                               Stock1_upper=NA,
                               
                               Stock2_lower=NA,
                               Stock2=mle$stock,
                               Stock2_upper=NA,
                               
                               Catches=mle$Catch,
                               Landings=NA,
                               Discards=NA,
                               
                               Fishing1_lower=NA,
                               Fishing1=mle$F,
                               Fishing1_upper=NA,
                               
                               Fishing2_lower=NA,
                               Fishing2=mle$harvest,
                               Fishing2_upper=NA
      )
      
    
      timeseries[-1] =     round(timeseries[-1],3)  
      endyr = which(mle$year==max(mle$year))                  
    
    refpts = round(t(data.frame(
      Ftgt = mvln$refpts[1,2],
      Btgt = mvln$refpts[2,2],
      Bthr= NA,
      Blim= NA,
      Fcur  = timeseries$Fishing1[endyr],
      Bcur  = timeseries$Stock1[endyr],
      B0.33= quantile(timeseries$Stock1,0.33),
      B0.66= quantile(timeseries$Stock1,0.66))),3)
    }
    
    refpts = data.frame(RefPoint=row.names(refpts),Value=refpts[,1])
    rownames(refpts) = 1:nrow(refpts)
    
    return(list(timeseries=timeseries,refpts=refpts))
  }
  