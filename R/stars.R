
#' ss2FLStockR()
#' @param mvln output from ssmvln() 
#' @param output expected outputs presented as "mle" or median of "iters"
#' @param thin thinnig rate of retained iters
#' @return FLStockR with refpts
#' @export
ss2FLStockR <- function(mvln,thin=10, output=NULL){
  kbinp = FALSE
  if(!is.null(mvln$mle)){
    if(is.null(output)) output ="mle"
  } else {
    output = "iters"
  }
  
  
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
    df = kb
    df = df[seq(1,nrow(df),thin),]
    year = unique(df$year)
    for(i in 1:length(year)){
      df[df$year%in%year[i],]$iter = 1:nrow(df[df$year%in%year[i],]) 
    }
    df = df[order(df$year),]
    
    
    
    N = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                              season="all",area="unique",iter=df$iter,data=df$Recr))
    C = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                              season="all",area="unique",iter=df$iter,data=df$Catch))
    Mat = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                                season="all",area="unique",iter=df$iter,data=df$SSB/df$Recr))
    H = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                              season="all",area="unique",iter=df$iter,data=df$F))
    year = unique(df$year)
  } else {
    N = as.FLQuant(data.frame(age=1,year=mle$year,unit="unique",
                              season="all",area="unique",iter=1,data=mle$Recr))
    C = as.FLQuant(data.frame(age=1,year=mle$year,unit="unique",
                              season="all",area="unique",iter=1,data=mle$Catch))
    Mat = as.FLQuant(data.frame(age=1,year=mle$year,unit="unique",
                                season="all",area="unique",iter=1,data=mle$SSB/mle$Recr))
    H = as.FLQuant(data.frame(age=1,year=mle$year,unit="unique",
                              season="all",area="unique",iter=1,data=mle$F)) 
    year = unique(mle$year)
    
  }
  
  stk = FLStockR(
    stock.n=N,
    catch.n = C,
    landings.n = C,
    discards.n = FLQuant(0, dimnames=list(age="1", year = year)),
    stock.wt=FLQuant(1, dimnames=list(age="1", year = year )),
    landings.wt=FLQuant(1, dimnames=list(age="1", year = (year))),
    discards.wt=FLQuant(1, dimnames=list(age="1", year = (year))),
    catch.wt=FLQuant(1, dimnames=list(age="1", year = (year))),
    mat=Mat,
    m=FLQuant(0.0001, dimnames=list(age="1", year = (year))),
    harvest = H,
    m.spwn = FLQuant(0, dimnames=list(age="1", year = (year))),
    harvest.spwn = FLQuant(0.0, dimnames=list(age="1", year = (year)))
  )
  units(stk) = standardUnits(stk)
  stk@catch = computeCatch(stk)
  stk@landings = computeLandings(stk)
  stk@discards = computeStock(stk)
  stk@stock = computeStock(stk)
  
  if(kbinp){
    stk@refpts = FLPar(
      Ftgt = median(kb$F/kb$harvest),
      Btgt = median(kb$SSB/kb$stock)
    )
  } else {
    stk@refpts = FLPar(
      Ftgt = mvln$refpts[1,2],
      Btgt = mvln$refpts[2,2],
      MSY = mvln$refpts[3,2],
      B0 = mvln$refpts[4,2],
      R0 = mvln$refpts[5,2]
    )
  }
  return(stk)
}


#' Function to summarise forecast results
#' @param object *FLStocks* with list of *FLStockR* objects 
#' @param eval.yrs evaluation years of forecast 
#' @param rel if TRUE ratios B/Btgt and F/Ftgt are shown
#' @return data.frame
#' @export
fwd2stars <- function(object,eval.yrs=NULL, rel=FALSE){
  if(!class(object)=="FLStocks"){
    object = FLStocks(forecast=object)
  }
  if(is.null(eval.yrs)){
    eval.yrs = an(range(object[[1]])["maxyear"])
  }
  if(!rel){
    df = do.call(rbind,Map(function(x,y){
      stk = window(x,start=min(eval.yrs),end=max(eval.yrs))
      flqs = FLQuants(
        Cy = round(landings(x)[,ac(eval.yrs)],3),
        By = round(ssb(x)[,ac(eval.yrs)],3),
        Fy = round(fbar(x)[,ac(eval.yrs)],3)
      )
      out = as.data.frame(flqs)
      data.frame(scenario=y, t(as.matrix(out$data)))
    },x=object,y=names(object))
    )
    names(df) = c("scenario",paste0("C",eval.yrs),
                  paste0("B",eval.yrs),
                  paste0("F",eval.yrs))
  }
  if(rel){
    
    df = do.call(rbind,Map(function(x,y){
      if(!class(x)=="FLStockR") 
        stop("input must be FLStockR object with @ref
            pts")
      
      stk = window(x,start=min(eval.yrs),end=max(eval.yrs))
      flqs = FLQuants(
        Cy = round(landings(x)[,ac(eval.yrs)],3),
        By = round(ssb(x)[,ac(eval.yrs)],3)/stk@refpts[[2]],
        Fy = round(fbar(x)[,ac(eval.yrs)]/stk@refpts[[1]],3)
      )
      out = as.data.frame(flqs)
      data.frame(scenario=y, t(as.matrix(out$data)))
    },x=object,y=names(object))
    )
    refn = rownames(object[[1]]@refpts)[1:2]
    names(df) = c("scenario",paste0("C",eval.yrs),
                  paste0("B",eval.yrs,"/",refn[2]),
                  paste0("F",eval.yrs,"/",refn[1]))
  }    
  rownames(df) = 1:nrow(df)
  return(df)  
} # End of function



#' jabba2FLStockR()
#' @param jabba fit from JABBA fit_jabba() or jabba$kbtrj 
#' @param blim biomass limit reference point as fraction of Bmsy
#' @param bpa biomass precautionary reference point as fraction of Bmsy
#' @param thin thinnig rate of retained iters 
#' @param rel if TRUE ratios BBmsy and FFmsy are stored
#' @return FLStockR with refpts
#' @export
jabba2FLStockR <- function(jabba,blim=0.3,bthr=0.5,thin=10,rel=FALSE){
  kbinp = FALSE
  if(is.null(jabba$assessment)){
    kbinp=TRUE
    kb = jabba
  }  else {
    kb = jabba$kbtrj
  }
  
  df = kb
  df = df[seq(1,nrow(df),thin),]
  year = unique(df$year)
  for(i in 1:length(year)){
    df[df$year%in%year[i],]$iter = 1:nrow(df[df$year%in%year[i],]) 
  }
  
  
  N = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                            season="all",area="unique",iter=df$iter,data=exp(df$Bdev)))
  C = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                            season="all",area="unique",iter=df$iter,data=df$Catch))
  
  if(!rel){
    Mat = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                                season="all",area="unique",iter=df$iter,data=df$B/exp(df$Bdev)))
    
    H = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                              season="all",area="unique",iter=df$iter,data=df$H))
    B= as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                             season="all",area="unique",iter=df$iter,data=df$B))
  } else {
    Mat = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                                season="all",area="unique",iter=df$iter,data=df$stock/exp(df$Bdev)))
    
    H = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                              season="all",area="unique",iter=df$iter,data=df$harvest))
    B= as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                             season="all",area="unique",iter=df$iter,data=df$stock))
    
  }
  
  
  stk = FLStockR(
    stock.n=N,
    catch.n = C,
    landings.n = C,
    discards.n = FLQuant(0, dimnames=list(age="1", year = (year))),
    stock.wt=FLQuant(1, dimnames=list(age="1", year = (year))),
    landings.wt=FLQuant(1, dimnames=list(age="1", year = year)),
    discards.wt=FLQuant(1, dimnames=list(age="1", year = year)),
    catch.wt=FLQuant(1, dimnames=list(age="1", year = year)),
    mat=Mat,
    m=FLQuant(0.0001, dimnames=list(age="1", year = year)),
    harvest = H,
    m.spwn = FLQuant(0, dimnames=list(age="1", year = year)),
    harvest.spwn = FLQuant(0.0, dimnames=list(age="1", year = year))
  )
  units(stk) = standardUnits(stk)
  stk@catch = computeCatch(stk)
  stk@landings = computeLandings(stk)
  stk@discards = computeStock(stk)
  stk@stock = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                                    season="all",area="unique",iter=df$iter,data=df$B))
  
  if(kbinp){
    stk@refpts = FLPar(
      Fmsy = median(kb$H)/median(kb$harvest),
      Bmsy = median(kb$B)/median(kb$stock),
      MSY = NA,
      Bthr= median(bthr*kb$B)/median(kb$stock),
      Blim= median(blim*kb$B)/median(kb$stock),
      B0 = median(kb$B)/median(kb$BB0),
    )
  } else {
    stk@refpts = FLPar(
      Fmsy = median(kb$H)/median(kb$harvest),
      Bmsy = median(kb$B)/median(kb$stock),
      MSY = jabba$refpts$msy[1],
      Blim= median(blim*kb$B)/median(kb$stock),
      Bpa= median(bthr*kb$B)/median(kb$stock),
      B0 = median(kb$B)/median(kb$BB0),
    )
  }
  if(rel){
    stk@refpts[1:2] = 1
    stk@refpts["Blim"] = blim
    stk@refpts["Bthr"] = bthr
    stk@refpts["B0"] = median(kb$B/kb$BB0)/median(kb$B/kb$stock)
  }
  stk@desc = "spm"
  return(stk)
}


#' spict2FLQuant()
#' @param x fit from SPICT
#' @param osa add one-step-ahead forecast
#' @param forecast TRUE/FALSE
#' @return FLQuant  
#' @author adopted from Laurie Kell (biodyn)
#' @export
spict2FLQuant <- function(x,metric=c("ssb","fbar","catch","stock","harvest")[1],osa=FALSE,forecast=F){
  
    
  vals=c("logB","logFnotS","logCpred","logBBmsy","logFFmsynotS")
  metrics=c("ssb","fbar","catch","stock","harvest")
  
  val = vals[which(metrics%in%metric)]
  
  if(val=="logB"){
    vec = x$par.random         
  } else {
    vec = x$value
  }  
  
  quant= an(exp(vec[which(names(vec)==val)]))
  if(val=="logCpred"){
    season = x$inp$dtc[1]
  } else {
    season = round(1/x$inp$dteuler)
  }
  
  timerange = c(x$inp$timerange[1],x$inp$timerange[2])
  if(forecast) timerange[2] = x$inp$maninterval[2]
  
  year.obs  =c(rep(seq(x$inp$timerange[1],x$inp$timerange[2]),each=season))
  year  =c(rep(seq(timerange[1],timerange[2]),each=season))
  seas=c(rep(seq(season),timerange[2]-timerange[1]+1))
  
  if(forecast == TRUE){
    year = year[1:length(quant)]
    seas = seas[1:length(quant)]
    seas[length(seas)] = season
  } else {
     quant = quant[1:length(year)]
   }         
    
  if(forecast == TRUE & val=="logCpred"){
    if(length(year.obs)<length(year))
            quant[(length(year.obs)+1):length(year)] = NA
   }


  dat=data.frame(age=1,year  =year,
                 season=seas,
                 data  =quant)
  if(val=="logCpred"){
    out= seasonSums(as.FLQuant(dat))} else {
      out= as.FLQuant(dat)[,,,season]
      dimnames(out)$season="all"
    }
  if(!forecast & !osa) out = window(out,end=floor(max(x$inp$timeC)))
  
  return(out)
} # End



#' spict2FLStockR()
#' @param res fit from SPICT  
#' @param blim biomass limit reference point as fraction of Bmsy
#' @param bthr biomass precautionary reference point as fraction of Bmsy
#' @param rel if TRUE ratios BBmsy and FFmsy are stored
#' @param osa add one-step-ahead forecast 
#' @param forecast extract forecast TRUE/FALSE
#' @return FLStockR with refpts
#' @export
spict2FLStockR <- function(res,blim=0.3,bthr=0.5,rel=FALSE,osa=FALSE,forecast=NULL){
  
  if(!is.null(forecast)){
    fw = forecast
  } else {
  if(is.null(res$man)) fw = FALSE
  if(!is.null(res$man)) fw = TRUE
  }
  if(is.null(res$value)){
    runs = ref
  }
  if(!is.null(res$value) & is.null(res$man) & is.null(res$retro)){
     runs = list(run=res)
  }
  if(!is.null(res$man)){
    runs = res$man
  }
  
  if(!is.null(res$retro)){
    runs = res$retro
  }
  
  # start 
  stks =  FLStocks(lapply(runs,function(res){
  
  
  b = spict2FLQuant(res,metric="ssb",forecast=fw,osa=osa)
  f = spict2FLQuant(res,metric="fbar",forecast=fw,osa=osa)
  
  if(!rel){
    B = spict2FLQuant(res,metric="ssb",forecast=fw,osa=osa)
    H = spict2FLQuant(res,metric="fbar",forecast=fw,osa=osa)
  } else {
    B = spict2FLQuant(res,metric="stock",forecast=fw,osa=osa)
    H = spict2FLQuant(res,metric="harvest",forecast=fw,osa=osa)
  }
  
  endyr = max(an(dimnames(B)$year))
  C = spict2FLQuant(res,metric="catch",forecast=fw,osa=osa)
  C =window(C,end=endyr)
  if(fw){
  intyr = dims(C[,!is.na(C)])$maxyear
  finyr = dims(C)$maxyear
  C[,ac(intyr:finyr)] = b[,ac(intyr:finyr)]*f[,ac(intyr:finyr)]
  }
  df = as.data.frame(B)
  year = unique(df$year)
  N = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                            season="all",area="unique",iter=1,data=1))
  Mat = B
  
  
  
  stk = FLStockR(
    stock.n=N,
    catch.n = C,
    landings.n = C,
    discards.n = FLQuant(0, dimnames=list(age="1", year = (year))),
    stock.wt=FLQuant(1, dimnames=list(age="1", year = (year))),
    landings.wt=FLQuant(1, dimnames=list(age="1", year = year)),
    discards.wt=FLQuant(1, dimnames=list(age="1", year = year)),
    catch.wt=FLQuant(1, dimnames=list(age="1", year = year)),
    mat=Mat,
    m=FLQuant(0.0001, dimnames=list(age="1", year = year)),
    harvest = H,
    m.spwn = FLQuant(0, dimnames=list(age="1", year = year)),
    harvest.spwn = FLQuant(0.0, dimnames=list(age="1", year = year))
  )
  units(stk) = standardUnits(stk)
  stk@catch = computeCatch(stk)
  stk@landings = computeLandings(stk)
  stk@discards = computeStock(stk)
  stk@stock = B
  
  stk@refpts = FLPar(
    Fmsy = res$report$Fmsy,
    Bmsy = res$report$Bmsy,
    MSY = res$report$MSY,
    Blim= res$report$Bmsy*blim,
    Bthr= res$report$Bmsy*bthr,
    B0 = res$value["K"],
  )
  
  if(rel){
    stk@refpts[1:2] =1 
    stk@refpts["Blim"] = blim
    stk@refpts["Bthr"] = bthr
    stk@refpts["B0"] = res$value["K"]/res$report$Bmsy
  }
   stk@desc = "spm"
   stk
   }))
  
  stks@desc ="spm"
  if(length(stks)==1) stks = stks[[1]] 
  
  return(stks)
}



#' flr2stars()
#' @param object of class FLStockR  
#' @param quantities default is 90CIs as c(0.05,0.95)
#' @return STARS list with $timeseris and $refpts
#' @export
flr2stars <- function(object,quantiles = c(0.05,0.95)){
  
  
  
    x= stockMedians(object)
    if((dims(object)$season+dims(object)$unit+dims(object)$area)>3){
      y = simplify(x)
    } else {
      y= x
    }
    
    endyr = dims(x)$maxyear
    refpts = round(rbind(x@refpts,FLPar(
                     Fcur = an(fbar(y)[,ac(endyr)]),
                     Bcur=an(unitSums(ssb(x)[,ac(endyr)])),
                     B0.33 = quantile(an(unitSums(ssb(x))),0.33,na.rm=T),
                     B0.66 = quantile(an(unitSums(ssb(x))),0.66,na.rm=T))),3)
    
    refpts = data.frame(RefPoint=row.names(refpts),Value=as.data.frame(refpts[,1])$data)
    rownames(refpts) = 1:nrow(refpts)
  
    if(dim(object)[6]==1){  

    timeseries =  data.frame(Year=dims(x)$minyear:dims(x)$maxyear,
                             Rec_lower=NA,
                             Rec=round(an(rec(x)),1),Rec_upper=NA,
                             SSB_lower=NA,
                             SSB=round(an(ssb(x)),1),
                             SSB_upper=NA,
                             
                             Bratio_lower=NA,
                             Bratio=round(an(ssb(x))/x@refpts[[2]],3),
                             Bratio_upper=NA,
                             
                             Catches=round(an(catch(x)),2),
                             Landings = round(an(landings(x)),2),
                             Discards = round(an(discards(x)),2),
                             
                             F_lower=NA,
                             F = round(an(fbar(x)),3),
                             F_upper=NA,
                             
                             Fratio_lower=NA,
                             Fratio=round(an(fbar(x))/x@refpts[[1]],3),
                             Fratio_upper=NA)
    } 
    
    if(dim(object)[6]>1){
      
      timeseries =  data.frame(Year=dims(object)$minyear:dims(object)$maxyear,
                               
                               Rec_lower=round(an(quantile(rec(object),quantiles[1])),0),
                               Rec=round(an(quantile(rec(object),0.5)),1),
                               Rec_upper=round(an(quantile(rec(object),quantiles[2])),0),
                               SSB_lower=round(an(quantile(ssb(object),quantiles[1])),1),
                               SSB=round(an(quantile(ssb(object),0.5)),1),
                               SSB_upper=round(an(quantile(ssb(object),quantiles[2])),1),
                               
                               Bratio_lower=round(an(quantile(ssb(object),quantiles[1]))/object@refpts[[2]],3),
                               Bratio=round(an(quantile(ssb(object),0.5))/object@refpts[[2]],3),
                               Bratio_upper=round(an(quantile(ssb(object),quantiles[2]))/object@refpts[[2]],3),
                               
                               Catches=round(an(quantile(catch(object),0.5)),2),
                               Landings = an(round(quantile(landings(object),0.5),2)),
                               Discards = an(round(quantile(discards(object),0.5),2)),
                               
                               F_lower=an(round(quantile(fbar(object),quantiles[1]),3)),
                               F = round(an(quantile(fbar(object),0.5)),3),
                               F_upper=an(round(quantile(fbar(object),quantiles[2]),3)),
                               
                               Fratio_lower=an(round(quantile(fbar(object),quantiles[1])/object@refpts[[1]],3)),
                               Fratio=round(an(quantile(fbar(object),0.5)/object@refpts[[1]]),3),
                               Fratio_upper=an(round(quantile(fbar(object),quantiles[2])/object@refpts[[1]],3)))   
      } 
    
    
    return(list(timeseries=timeseries,refpts=refpts))
  } 
  

#' updstars()
#' @param star output of star list 
#' @param newrefpts manually adjusted reference points
#' @return STARS list with $timeseris and $refpts
#' @export

updstars <- function(star,newrefpts){
  newts = star$timeseries
  bmsyr = star$refpts[2,2]/newrefpts[2,2]
  fmsyr = star$refpts[1,2]/newrefpts[1,2]
  newts[,8:10] = newts[,8:10]*bmsyr 
  newts[,8:10] = newts[,17:19]*fmsyr 
  out = star
  out$timeseries = newts
  out$refpts = newrefpts
  return(out)
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
                             Landings=mu[,quants[3]][,1],
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
      Ftgt = median(kb$F)/median(kb$harvest),
      Btgt = median(kb$SSB)/median(kb$stock),
      Bthr= NA,
      Blim= NA,
      Fcur  = timeseries$F[endyr],
      Bcur  = timeseries$SSB[endyr],
      B0.33= quantile(timeseries$SSB,0.33),
      B0.66= quantile(timeseries$SSB,0.66))),3)
    
  }
  if(output=="mle"){
    quants=c("SSB","F","Catch","stock","harvest","Recr")
    
    mu = aggregate(cbind(stock,harvest,SSB,F,Catch,Recr)~year+run,kb,
                   quantile,c(0.5,quantiles))
    
    timeseries =  data.frame(year=mle$year,
                             Rec_lower=mu[,quants[6]][,2],
                             Rec=mle$Recr,
                             Rec_upper=mu[,quants[6]][,3],
                             SSB_lower=mu[,quants[1]][,2],
                             SSB=mle$SSB,
                             SSB_upper=mu[,quants[1]][,3],
                             
                             Bratio_lower=mu[,quants[4]][,2],
                             Bratio=mle$stock,
                             Bratio_upper=mu[,quants[6]][,3],
                             
                             Catches=mle$Catch,
                             Landings=NA,
                             Discards=NA,
                             
                             F_lower=mu[,quants[2]][,2],
                             F=mle$F,
                             F_upper=mu[,quants[2]][,3],
                             
                             Fratio_lower=mu[,quants[5]][,2],
                             Fratio=mle$harvest,
                             Fratio_upper=mu[,quants[5]][,3]
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
  #' @param quantiles default is 90CIs as c(0.05,0.95)
  #' @param blim biomass limit point as fraction of Bmsy, default 0.3Bmsy (ICES) 
  #' @param bthr biomass precautionary point as fraction of Bmsy, default 0.5Bmsy (ICES) 
  #' @return STARS list with $timeseris and $refpts
  #' @export
  jabba2stars <- function(jabba,quantiles = c(0.05,0.95),
                          blim = 0.3,bthr=0.5){
    kbinp = FALSE
    if(is.null(jabba$assessment)){
      kbinp=TRUE
      kb = jabba
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
        Ftgt = median(kb$H)/median(kb$harvest),
        Btgt = median(kb$B)/median(kb$stock),
        Bthr= median(bthr*kb$B)/median(kb$stock),
        Blim= median(blim*kb$B)/median(kb$stock),
        Fcur  = timeseries$F[endyr],
        Bcur  = timeseries$Biomass[endyr],
        B0.33= quantile(timeseries$Biomass,0.33),
        B0.66= quantile(timeseries$Biomass,0.66))),3)
      
      refpts = data.frame(RefPoint=row.names(refpts),Value=refpts[,1])
     rownames(refpts) = 1:nrow(refpts)
      
      return(list(timeseries=timeseries,refpts=refpts))
  }
  
  #' spict2stars()
  #' @param spict fit from  fit.spict()  
  #' @param blim biomass limit point as fraction of Bmsy, default 0.3Bmsy (ICES) 
  #' @param bthr biomass precautionary point as fraction of Bmsy, default 0.5Bmsy (ICES) 
  #' @return STARS list with $timeseris and $refpts
  #' @export
  spict2stars <- function(spict,
                          blim = 0.3,bthr=0.5){
    
    stk = spict2FLStockR(spict)
    stkr = spict2FLStockR(spict,rel=TRUE)
    yrs = an(dimnames(stk)$year)
    # TODO add uncertainty
    timeseries =  data.frame(year=yrs,
                             Rec_lower=NA,Rec=NA,Rec_upper=NA,
                             Biomass_lower=NA,
                             Biomass=an(ssb(stk)),
                             Biomass_upper=NA,
                             
                             Bratio_lower=NA,
                             Bratio=an(ssb(stkr)),
                             Bratio_upper=NA,
                             
                             Catches=an(catch(stk)),
                             Landings=an(landings(stk)),
                             Discards=NA,
                             
                             F_lower=NA,
                             F=an(fbar(stk)),
                             F_upper=NA,
                             
                             Fratio_lower=NA,
                             Fratio=an(fbar(stkr)),
                             Fratio_upper=NA
    )
    
    timeseries[-1] =     round(timeseries[-1],3)                     
    endyr = which(max(yrs)==yrs)
   # TODO change system to relative
    refpts = round(t(data.frame(
      Ftgt = stk@refpts[[1]],
      Btgt = stk@refpts[[2]],
      Bthr= bthr*stk@refpts[[2]],
      Blim= blim*stk@refpts[[2]],
      Fcur  = timeseries$F[endyr],
      Bcur  = timeseries$Biomass[endyr],
      B0.33= quantile(timeseries$Biomass,0.33),
      B0.66= quantile(timeseries$Biomass,0.66))),3)
    
    refpts = data.frame(RefPoint=row.names(refpts),Value=refpts[,1])
    rownames(refpts) = 1:nrow(refpts)
    
    return(list(timeseries=timeseries,refpts=refpts))
  }
  
  
  #' stock2ratios()
  #' @param object of class *FLStockR*  
  #' @return FLStockR with ratios F/Ftgt and B/Btgt
  #' @export
  stock2ratios <- function(object){
      stk = object  
      B = as.FLQuant(ssb(object)/object@refpts[[2]])
      H = fbar(object)/object@refpts[[1]]
      C = computeLandings(object)
      dimnames(B)$age = "1"
      dimnames(C)$age = "1"
      dimnames(H)$age = "1"
      
      df = as.data.frame(B)
       year = unique(df$year)
       iters = an(unique(df$iter))
      
      
      
      stk = FLStockR(stock.n=FLQuant(1, dimnames=list(age="1", year = (year),iter=iters)),
      catch.n = C,
      landings.n = C,
      discards.n = FLQuant(0, dimnames=list(age="1", year = (year),iter=iters)),
      stock.wt=FLQuant(1, dimnames=list(age="1", year = (year),iter=iters)),
      landings.wt=FLQuant(1, dimnames=list(age="1", year = year,iter=iters)),
      discards.wt=FLQuant(1, dimnames=list(age="1", year = year,iter=iters)),
      catch.wt=FLQuant(1, dimnames=list(age="1", year = year,iter=iters)),
      mat=B,
      m=FLQuant(0.0001, dimnames=list(age="1", year = year)),
      harvest = H,
      m.spwn = FLQuant(0, dimnames=list(age="1", year = year)),
      harvest.spwn = FLQuant(0.0, dimnames=list(age="1", year = year))
    )
    units(stk) = standardUnits(stk)
    stk@catch = computeCatch(stk)
    stk@landings = computeLandings(stk)
    stk@discards = computeStock(stk)
    stk@stock = B
    stk@refpts = FLPar(Ftgt=1,Btgt=1,MSY=object@refpts[[3]])
    row.names(stk@refpts)[1:2] = row.names(object@refpts)[1:2] 
    stk@desc = object@desc
    
    return(stk)
  }


