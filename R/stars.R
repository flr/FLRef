
#' ss2FLStockR()
#' @param mvln output from ssmvln() 
#' @param output expected outputs presented as "mle" or median of "iters"
#' @param thin thinnig rate of retained iters
#' @param rel if TRUE ratios B/Btgt and F/Ftgt are shown
#' @return FLStockR with refpts
#' @export
ss2FLStockR <- function(mvln,thin=1, output=NULL,rel=FALSE){
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
    
    if(rel){
      df$b = df$stock  
      df$f = df$harvest  
    } else {
      df$b = df$SSB  
      df$f = df$F
    }
    
    N = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                              season="all",area="unique",iter=df$iter,data=df$Recr))
    C = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                              season="all",area="unique",iter=df$iter,data=df$Catch))
    Mat = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                                season="all",area="unique",iter=df$iter,data=df$b/df$Recr))
    
    H = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                              season="all",area="unique",iter=df$iter,data=df$f))
    year = unique(df$year)
  } else {
    
    if(rel){
      mle$b = mle$stock  
      mle$f = mle$harvest  
    } else {
      mle$b = mle$SSB  
      mle$f = mle$F
    }
    N = as.FLQuant(data.frame(age=1,year=mle$year,unit="unique",
                              season="all",area="unique",iter=1,data=mle$Recr))
    C = as.FLQuant(data.frame(age=1,year=mle$year,unit="unique",
                              season="all",area="unique",iter=1,data=mle$Catch))
    Mat = as.FLQuant(data.frame(age=1,year=mle$year,unit="unique",
                                season="all",area="unique",iter=1,data=mle$b/mle$Recr))
    H = as.FLQuant(data.frame(age=1,year=mle$year,unit="unique",
                              season="all",area="unique",iter=1,data=mle$f)) 
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
    kbfit = kb[kb$type=="fit",]
    stk@refpts = FLPar(
      Ftgt = median(kbfit$F)/(kbfit$harvest),
      Btgt = median(kbfit$SSB)/(kbfit$stock)
    )
    
    if(rel){
      stk@refpts[1:2][] = 1
    }
      
      
  } else {
    stk@refpts = FLPar(
      Ftgt = mvln$refpts[1,2],
      Btgt = mvln$refpts[2,2],
      MSY = mvln$refpts[3,2],
      B0 = mvln$refpts[4,2],
      R0 = mvln$refpts[5,2]
    )
    
    if(rel){
      stk@refpts["B0"] = stk@refpts["B0"]/stk@refpts[2]
      stk@refpts[1:2][] = 1
    }
  }
  return(stk)
}


#' Function to summarise forecast results
#' @param object *FLStocks* with list of *FLStockR* objects 
#' @param eval.yrs evaluation years of forecast 
#' @param rel if TRUE ratios B/Btgt and F/Ftgt are shown
#' @param dB computes change in percentage biomass to refyr
#' @return data.frame
#' @export
fwd2stars <- function(object,eval.yrs=NULL, rel=NULL,dB=NULL,refyr=NULL){

  
  if(!class(object)=="FLStocks"){
    object = FLStocks(forecast=object)
  }
  if(any(c("Btgt","Bmsy")%in%rownames(object[[1]]@refpts))){
    if(is.null(dB)) dB = FALSE
    if(is.null(rel)) rel = TRUE
  } else {
    if(is.null(dB)) dB =TRUE
    if(is.null(rel)) rel = FALSE
    
  }
  
  if(is.null(refyr)){
    refyr = an(range(object[[1]])["minyear"])
  }
    
  if(is.null(eval.yrs)){
    eval.yrs = an(range(object[[1]])["maxyear"])
  }
  if(!rel){
    df = do.call(rbind,Map(function(x,y){
      stk = window(x,start=min(eval.yrs),end=max(eval.yrs))
      flqs = FLQuants(
        Cy = round(catch(x)[,ac(eval.yrs)],1),
        Fy = round(fbar(x)[,ac(eval.yrs)],3),
        By = round(ssb(x)[,ac(eval.yrs)],1)
      )
      out = as.data.frame(flqs)
      data.frame(scenario=y, t(as.matrix(out$data)))
    },x=object,y=names(object))
    )
    names(df) = c("scenario",paste0("C",eval.yrs),
                  paste0("F",eval.yrs),
                  paste0("B",eval.yrs))
  }
  
  if(rel){
    
  
    
    df = do.call(rbind,Map(function(x,y){
      if(!class(x)=="FLStockR") 
        stop("input must be FLStockR object with @ref
            pts")
      if(!any(c("Btgt","Bmsy")%in%rownames(object[[1]]@refpts))){
        stop("No B target (Btgt) provided to compute ratio of B/Btgt")
      }
      stk = window(x,start=min(eval.yrs),end=max(eval.yrs))
      flqs = FLQuants(
        Cy = round(catch(x)[,ac(eval.yrs)],1),
        Fy = round(fbar(x)[,ac(eval.yrs)]/stk@refpts[[1]],3),
        By = round(ssb(x)[,ac(eval.yrs)]/stk@refpts[[2]],3))
    
      if(any(c("Bmsy","Btgt")%in%rownames(object[[1]]@refpts))){
        flqs = FLQuants(c(flqs[-3],FLQuants(Btgt=round(ssb(x)[,ac(eval.yrs)]/an(object[[1]]@refpts[2]),3))))
      }
      
      
      if("Bpa"%in%rownames(object[[1]]@refpts)){
        flqs = FLQuants(c(flqs,FLQuants(Bpa=round(ssb(x)[,ac(eval.yrs)]/an(object[[1]]@refpts["Bpa"]),3))))
      }
      if("Bthr"%in%rownames(object[[1]]@refpts)){
        flqs = FLQuants(c(flqs,FLQuants(Bthr=round(ssb(x)[,ac(eval.yrs)]/an(object[[1]]@refpts["Bthr"]),3))))
      }
      if("Blim"%in%rownames(object[[1]]@refpts)){
        flqs = FLQuants(c(flqs,FLQuants(Blim=round(ssb(x)[,ac(eval.yrs)]/an(object[[1]]@refpts["Blim"]),3))))
      }
      
      out = as.data.frame(flqs)
      data.frame(scenario=y, t(as.matrix(out$data)))
    },x=object,y=names(object))
    )
    
    refn = rownames(object[[1]]@refpts)[1:2]
    
    nam = c("scenario",paste0("C",eval.yrs),
                  paste0("F",eval.yrs,"/",refn[1]),
                  paste0("B",eval.yrs,"/",refn[2])
                  )
    
    #if(any(c("Bmsy","Btgt")%in%rownames(object[[1]]@refpts))){
    #  nam = c(nam[-4], paste0("B",eval.yrs,"/",rownames(object[[1]]@refpts)[2]))
    #}
    
    if("Bpa"%in%rownames(object[[1]]@refpts)){
      nam = c(nam, paste0("B",eval.yrs,"/","Bpa"))
    }
    if("Bthr"%in%rownames(object[[1]]@refpts)){
      nam = c(nam, paste0("B",eval.yrs,"/","Bthr"))
    }
    if("Blim"%in%rownames(object[[1]]@refpts)){
      nam = c(nam, paste0("B",eval.yrs,"/","Blim"))
    }
    names(df) = nam  
  }    
  
  if(dB){
    dBs =(do.call(rbind,lapply(fstks,function(x){
      round(matrix(an(100*(ssb(x)[,ac(eval.yrs)]/ssb(x)[,ac(rep(refyr,length(eval.yrs)))]-1)),ncol=length(eval.yrs)),2)
    })))
  
  nam = c(names(df),paste0("dB%",eval.yrs))
  df = cbind(df,dBs)
  names(df) = nam
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
      Bthr= median(bthr*kb$B)/median(kb$stock),
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
#' @param what mle or log.sd
#' @return FLQuant  
#' @author adopted from Laurie Kell (biodyn)
#' @export
spict2FLQuant <- function(x,metric=c("ssb","fbar","catch","stock","harvest")[1],osa=FALSE,forecast=F,what=c("mle")){
  
    
  vals=c("logB","logFnotS","logCpred","logBBmsy","logFFmsynotS")
  metrics=c("ssb","fbar","catch","stock","harvest")
  val = vals[which(metrics%in%metric)]
  
  if(what=="mle"){
   if(val=="logB"){
    vec = x$par.random
    vn = names(x$par.random)
   } else {
    vec = x$value
    vn = names(x$value)
   }  
  } else {
    if(val=="logB"){
      vec = x$diag.cov.random     
      vn = names(x$par.random)
    } else {
      vec = x$sd
      vn = names(x$value)
    } 
  }
  
  
  quant= an((vec[which(vn==val)]))
  if(what=="mle") quant = exp(quant)
   
  
     
  if(val=="logCpred"){
    season = x$inp$dtc[1]
  } else {
    season = round(1/x$inp$dteuler)
  }
 
 
  timerange = c(x$inp$timerange[1],x$inp$timerange[2])
  if(forecast) timerange[2] = x$inp$maninterval[2]
  
  year.obs  =c(rep(seq(x$inp$timerange[1],x$inp$timerange[2]),each=season))
  year  =floor(c(rep(seq(timerange[1],timerange[2]),each=season)))
  seas=c(rep(seq(season),length(year)/season))
  #seas=c(rep(seq(season),timerange[2]-timerange[1]+1))
  min.timeI = min(do.call(c,lapply(x$inp$timeI,function(t)min(t))))
  
  
  if(val=="logCpred"& min(x$inp$timeC)>min(min.timeI)){
    ctinp = floor(x$inp$timeC)
    #ctinp = c(ctinp,max(ctinp)+1) 
    cadj = rep(NA,length(ctinp))
    qs =quant[1:length(year)] 
    cadj[year%in%ctinp] = qs[!is.na(qs)][1:length(ctinp)]
    quant=cadj
  }
    
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
#' @param itsCI number of iterations to depict uncertainty in plots
#' @return FLStockR with refpts
#' @export
spict2FLStockR <- function(res,blim=0.3,bthr=0.5,rel=FALSE,osa=FALSE,forecast=NULL,itsCI=1){
  
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
  if(itsCI>1){
    sdb = spict2FLQuant(res,metric="ssb",forecast=fw,osa=osa,what="sd")
    b =rlnorm(itsCI,log(b),sdb)
    f = propagate(f,itsCI)
  }
  
  
  if(!rel){
    B = spict2FLQuant(res,metric="ssb",forecast=fw,osa=osa)
    H = spict2FLQuant(res,metric="fbar",forecast=fw,osa=osa)
     if(itsCI>1){
     sdB = spict2FLQuant(res,metric="ssb",forecast=fw,osa=osa,what="sd")
     sdH = spict2FLQuant(res,metric="fbar",forecast=fw,osa=osa,what="sd")
     B =rlnorm(itsCI,mean=log(B),sd=sdB)
     H = rlnorm(itsCI,log(H),sdB)
     }
    
  } else {
    B = spict2FLQuant(res,metric="stock",forecast=fw,osa=osa)
    H = spict2FLQuant(res,metric="harvest",forecast=fw,osa=osa)
     if(itsCI>1){
      sdB = spict2FLQuant(res,metric="stock",forecast=fw,osa=osa,what="sd")
      sdH = spict2FLQuant(res,metric="harvest",forecast=fw,osa=osa,what="sd")
      B =rlnorm(itsCI,log(B),sdB)
      H = rlnorm(itsCI,log(H),sdH)
     }
  }
  
  endyr = max(an(dimnames(B)$year))
  C = spict2FLQuant(res,metric="catch",forecast=fw,osa=osa)
  C =window(C,end=endyr)
  
  
  if(fw){
  intyr = dims(C[,!is.na(C)])$maxyear
  finyr = dims(C)$maxyear
  }
  
  
  if(itsCI>1){
    C = propagate(C,itsCI)
    sdC = spict2FLQuant(res,metric="catch",forecast=fw,osa=osa,what="sd")
    
    C[,(dimnames(sdC)$year)[!is.na(sdC)]] =rlnorm(itsCI,log(C)[,(dimnames(sdC)$year)[!is.na(sdC)]],sdC[,(dimnames(sdC)$year)[!is.na(sdC)]])
  }  
  
  
  if(fw){
    C[,ac(intyr:finyr)] = b[,ac(intyr:finyr)]*f[,ac(intyr:finyr)]
  }
  
  df = as.data.frame(B)
  year = unique(df$year)
  N = as.FLQuant(data.frame(age=1,year=year,unit="unique",
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
  stk = propagate(stk,itsCI)
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


#' sam2FLQuant()
#' @param object fit from FLSAM
#' @param forecast TRUE/FALSE
#' @return FLQuant  
#' @author adopted from Laurie Kell (biodyn)
#' @export
sam2FLQuant <- function(object,metric=c("ssb","fbar","catch","rec")[1],what=c("mle")){
  
  
  vals=c("logssb","logfbar","logCatch","logR","logLand")
  metrics=c("ssb","fbar","catch","rec","landings")
  val = vals[which(metrics%in%metric)]
  
  res =  subset(object@params,name%in%val)[,c("name","value","std.dev")]
  
  res = cbind(year=seq(object@range["minyear"],
                 object@range["maxyear"]),res)
  
  if(what=="mle"){
  
  dat=data.frame(age=1,year  =res$year,
                 season="all",
                 data  =exp(res$value))
  } else {
    
  dat=data.frame(age=1,year  =res$year,
                   season="all",
                   data  =res$std.dev)  
    
  }
  
  out= as.FLQuant(dat)
  
  return(out)
} # End


#' sam2FLStockR()
#' @param res fit from FLSAM  
#' @param itsCI number of iterations to depict uncertainty in plots
#' @return FLStockR with refpts
#' @export
sam2FLStockR <- function(res,itsCI=1000){
  
  # start 
    B = sam2FLQuant(res,metric="ssb")
    H = sam2FLQuant(res,metric="fbar")
    C = sam2FLQuant(res,metric="catch")
    R = sam2FLQuant(res,metric="rec")
    
    if(itsCI>1){
      sdb = sam2FLQuant(res,metric="ssb",what="sd")
      sdf = sam2FLQuant(res,metric="fbar",what="sd")
      sdc = sam2FLQuant(res,metric="catch",what="sd")
      sdr = sam2FLQuant(res,metric="rec",what="sd")
      B =rlnorm(itsCI,log(B),sdb)
      H = rlnorm(itsCI,log(H),sdf)
      C = rlnorm(itsCI,log(C),sdc)
      R = rlnorm(itsCI,log(R),sdr)
    }
    
    
    
    
    
    df = as.data.frame(B)
    year = unique(df$year)
    N = R
    Mat = B/R
    
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
    
    
    
    stk@desc = "aspm"
    stk

  
  
  return(stk)
}


#' flr2stars()
#' @param object of class FLStockR (MLE)
#' @param uncertainty of class FLStock with iters
#' @param quantities default is 90CIs as c(0.05,0.95)
#' @return STARS list with $timeseris and $refpts
#' @export
flr2stars <- function(object,uncertainty=NULL,quantiles = c(0.05,0.95)){
  
  
  
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
    
    refpts = data.frame(RefPoint=rownames(refpts),Value=as.data.frame(refpts[,1])$data)
    rownames(refpts) = 1:nrow(refpts)
  
   if(dim(object)[6]==1){  

    timeseries =  data.frame(Year=dims(x)$minyear:dims(x)$maxyear,
                             Rec_lower=NA,
                             Rec=round(an(rec(x)),1),
                             Rec_upper=NA,
                             
                             SSB_lower=NA,
                             SSB=round(an(ssb(x)),1),
                             SSB_upper=NA,
                             
                             Bratio_lower=NA,
                             Bratio=NA,
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
    
                             if(any(c("Bmsy","Btgt")%in%rownames(object@refpts))){
                             timeseries$Bratio=round(an(ssb(x))/x@refpts[[2]],3)
                             }
       } 
    
    if(dim(object)[6]>1){
      uncertainty = object
      timeseries =  data.frame(Year=dims(uncertainty)$minyear:dims(uncertainty)$maxyear,
                               
                               Rec_lower=round(an(quantile(rec(uncertainty),quantiles[1])),0),
                               Rec=round(an(quantile(rec(uncertainty),0.5)),1),
                               Rec_upper=round(an(quantile(rec(uncertainty),quantiles[2])),0),
                               SSB_lower=round(an(quantile(ssb(uncertainty),quantiles[1])),1),
                               SSB=round(an(quantile(ssb(uncertainty),0.5)),1),
                               SSB_upper=round(an(quantile(ssb(uncertainty),quantiles[2])),1),
                               
                               Bratio_lower=NA,
                               Bratio=NA,
                               Bratio_upper=NA,
                               
                               Catches=round(an(quantile(catch(uncertainty),0.5)),2),
                               Landings = an(round(quantile(landings(uncertainty),0.5),2)),
                               Discards = an(round(quantile(discards(uncertainty),0.5),2)),
                               
                               F_lower=an(round(quantile(fbar(uncertainty),quantiles[1]),3)),
                               F = round(an(quantile(fbar(uncertainty),0.5)),3),
                               F_upper=an(round(quantile(fbar(uncertainty),quantiles[2]),3)),
                               
                               Fratio_lower=an(round(quantile(fbar(uncertainty)/uncertainty@refpts[[1]],quantiles[1]),3)),
                               Fratio=round(an(quantile(fbar(uncertainty)/uncertainty@refpts[[1]],0.5)),3),
                               Fratio_upper=an(round(quantile(fbar(uncertainty)/object@refpts[[1]],quantiles[2]),3)))   
    
      
                              # Add Bratio
                              if(any(c("Bmsy","Btgt")%in%rownames(object@refpts))){
                                timeseries$Bratio_lower=round(an(quantile(ssb(uncertainty),quantiles[1]))/object@refpts[[2]],3)
                                timeseries$Bratio=round(an(quantile(ssb(uncertainty)/uncertainty@refpts[[2]],0.5)),3)
                                timeseries$Bratio_upper=round(an(quantile(ssb(uncertainty),quantiles[2]))/object@refpts[[2]],3)
                              }
      } 
    
      
    
    if(!is.null(uncertainty)){
      
      timeseries =  data.frame(Year=dims(uncertainty)$minyear:dims(uncertainty)$maxyear,
                               
                               Rec_lower=round(an(quantile(rec(uncertainty),quantiles[1])),0),
                               Rec=round(an(quantile(rec(object),0.5)),1),
                               Rec_upper=round(an(quantile(rec(uncertainty),quantiles[2])),0),
                               SSB_lower=round(an(quantile(ssb(uncertainty),quantiles[1])),1),
                               SSB=round(an(quantile(ssb(object),0.5)),1),
                               SSB_upper=round(an(quantile(ssb(uncertainty),quantiles[2])),1),
                              
                        
                               Bratio_lower=NA,
                               Bratio=NA,
                               Bratio_upper=NA,
                               
                               Catches=round(an(quantile(catch(object),0.5)),2),
                               Landings = an(round(quantile(landings(object),0.5),2)),
                               Discards = an(round(quantile(discards(object),0.5),2)),
                               
                               F_lower=an(round(quantile(fbar(uncertainty),quantiles[1]),3)),
                               F = round(an(quantile(fbar(object),0.5)),3),
                               F_upper=an(round(quantile(fbar(uncertainty),quantiles[2]),3)),
                               
                               Fratio_lower=an(round(quantile(fbar(uncertainty)/object@refpts[[1]],quantiles[1]),3)),
                               Fratio=round(an(quantile(fbar(object)/object@refpts[[1]],0.5)),3),
                               Fratio_upper=an(round(quantile(fbar(uncertainty)/object@refpts[[1]],quantiles[2]),3)))   
      
                              # Add Bratio
                              if(any(c("Bmsy","Btgt")%in%rownames(object@refpts))){
                                timeseries$Bratio_lower=round(an(quantile(ssb(uncertainty),quantiles[1]))/object@refpts[[2]],3)
                                timeseries$Bratio=round(an(quantile(ssb(object)/object@refpts[[2]],0.5)),3)
                                timeseries$Bratio_upper=round(an(quantile(ssb(uncertainty),quantiles[2]))/object@refpts[[2]],3)
                              }
    } 
    
    
    
    
    return(list(timeseries=timeseries,refpts=refpts))
  } 
  

#' updstars()
#' @param star output of star list with new refpoints 
#' @param newrefpts FLPar manually adjusted reference points
#' @return STARS list with $timeseris and $refpts
#' @export

updstars <- function(star,newrefpts){
  newts = star$timeseries
  bmsyr = star$refpts[2,2]/newrefpts[[2]]
  fmsyr = star$refpts[1,2]/newrefpts[[1]]
  
  star$refpts[1,2] = newrefpts[[1]]
  star$refpts[2,2] = newrefpts[[2]]
  if(any(c("Bpa","Bthr")%in%rownames(newrefpts))){
    star$refpts[3,2] =  c(newrefpts[c("Bpa","Bthr")[c("Bpa","Bthr")%in%rownames(newrefpts)]])
  }
  if(any(c("Blim")%in%rownames(newrefpts))){
    star$refpts[4,2] =  c(newrefpts[c("Blim")[c("Blim")%in%rownames(newrefpts)]])
  }
  
  newts[,8:10] = newts[,8:10]*bmsyr 
  newts[,17:19] = newts[,17:19]*fmsyr 
  out = star
  out$timeseries = newts
  out$refpts = star$refpts
  return(out)
}


#' ss2stars()
#' @param mvln output of ssmvln() 
#' @param output choice c("iters","mle")[1]
#' @param quantiles default is 95CIs as c(0.025,0.975)
#' @return STARS list with $timeseris and $refpts
#' @export
ss2stars <- function(mvln,output=c("iters","mle")[2],quantiles = c(0.05,0.95)){
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
                             Bratio_upper=mu[,quants[4]][,3],
                             
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
  #' @param bthr biomass precautionary point as fraction of Bmsy, default 0.5Bmsy (ICES) 
  #' @return STARS list with $timeseris and $refpts
  #' @export
  spict2stars <- function(spict,
                          blim = 0.3,bthr=0.5,quantiles = c(0.05,0.95)){
    
    stka = spict2FLStockR(spict,bthr=bthr,blim = blim)
    stk = spict2FLStockR(spict,rel=TRUE,bthr=bthr,blim = blim)
    stki = spict2FLStockR(spict,bthr=bthr,blim = blim,rel=T,itsCI=5000)
    stkai = spict2FLStockR(spict,bthr=bthr,blim = blim,rel=F,itsCI=5000)
    
     
    yrs = an(dimnames(stk)$year)
  
  
    # Added uncertainty
    timeseries =  data.frame(year=yrs,
                             Rec_lower=NA,Rec=NA,Rec_upper=NA,
                             Biomass_lower=round(an(quantile(ssb(stkai),quantiles[1])),1),
                             Biomass=round(an(ssb(stka)),1),
                             Biomass_upper=round(an(quantile(ssb(stkai),quantiles[2])),1),
  
                             Bratio_lower=round(an(quantile(ssb(stki),quantiles[1])),3),
                             Bratio=round(an(ssb(stk)),3),
                             Bratio_upper=round(an(quantile(ssb(stki),quantiles[2])),3),
                             
                             Catches=an(catch(stk)),
                             Landings=an(landings(stk)),
                             Discards=NA,
                             
                             F_lower=round(an(quantile(fbar(stkai),quantiles[1])),3),
                             F=round(an(fbar(stka)),3),
                             F_upper=round(an(quantile(fbar(stkai),quantiles[2])),3),
                             
                             Fratio_lower=round(an(quantile(fbar(stki),quantiles[1])),3),
                             Fratio=round(an(fbar(stk)),3),
                             Fratio_upper=round(an(quantile(fbar(stki),quantiles[2])),3)
    )
    
    endyr = which(max(yrs)==yrs)
    stka@refpts[[1]] = (timeseries$F/timeseries$Fratio)[endyr] # Adjust for time-varying
    stka@refpts[[2]] = (timeseries$Biomass/timeseries$Bratio)[endyr] # Adjust for time-varying
    
   # TODO change system to relative
    refpts = round(t(data.frame(
      Ftgt = stka@refpts[[1]],
      Btgt = stka@refpts[[2]],
      Bthr= bthr*stka@refpts[[2]],
      Blim= blim*stka@refpts[[2]],
      Fcur  = timeseries$F[endyr],
      Bcur  = timeseries$Biomass[endyr],
      B0.33= quantile(timeseries$Biomass,0.33),
      B0.66= quantile(timeseries$Biomass,0.66))),3)
    
    refpts = data.frame(RefPoint=row.names(refpts),Value=refpts[,1])
    rownames(refpts) = 1:nrow(refpts)
    
    return(list(timeseries=timeseries,refpts=refpts))
  }
  
  #' sam2stars()
  #' @param sam fit from FLSAM  
  #' @param refpts optional FLPar object with reference points 
  #' @return STARS list with $timeseris and $refpts
  #' @export
  sam2stars <- function(sam,
                          refpts = NULL,quantiles = c(0.05,0.95)){
    rps = refpts
    stka = sam2FLStockR(sam,itsCI=1)
    stkai = sam2FLStockR(sam,itsCI=5000)
    
    
    yrs = an(dimnames(stka)$year)
    
    Bratio_lower=Bratio=Bratio_upper = NA
    Fratio_lower=Fratio=Fratio_upper = NA
    Landings = sam2FLQuant(sam,metric="landings")
    
    if(!is.null(rps)){
      ftgt = an(rps[rownames(rps)[grep("F",rownames(rps))]])[1]
      rps.tgt = rps[!rownames(rps)%in%c("Blim","Bpa","Bthr")]
      btgt = an(rps.tgt[rownames(rps.tgt)[grep("B",rownames(rps.tgt))]])[1]
      Bratio_lower=round(an(quantile(ssb(stkai),quantiles[1])),3)/btgt
      Bratio=round(an(ssb(stka)),3)/btgt
      Bratio_upper=round(an(quantile(ssb(stkai),quantiles[2])),3)/btgt
      
      Fratio_lower=round(an(quantile(fbar(stkai),quantiles[1])),3)/ftgt
      Fratio=round(an(fbar(stka)),3)/ftgt
      Fratio_upper=round(an(quantile(fbar(stkai),quantiles[2])),3)/ftgt
   }
    
    
    
    
    # Added uncertainty
    timeseries =  data.frame(year=yrs,
                             Rec_lower=round(an(quantile(rec(stkai),quantiles[1])),1),
                             Rec=round(an(rec(stka)),1),
                             Rec_upper=round(an(quantile(rec(stkai),quantiles[2])),1),
                             
                             Biomass_lower=round(an(quantile(ssb(stkai),quantiles[1])),1),
                             Biomass=round(an(ssb(stka)),1),
                             Biomass_upper=round(an(quantile(ssb(stkai),quantiles[2])),1),
                             
                             Bratio_lower=Bratio_lower,
                             Bratio=Bratio,
                             Bratio_upper=Bratio_upper,
                             
                               Catches=round(an(catch(stka)),1),
                             Landings=round(an(Landings),1),
                             Discards=NA,
                             
                             F_lower=round(an(quantile(fbar(stkai),quantiles[1])),3),
                             F=round(an(fbar(stka)),3),
                             F_upper=round(an(quantile(fbar(stkai),quantiles[2])),3),
                             
                             Fratio_lower=Fratio_lower,
                             Fratio=Fratio,
                             Fratio_upper=Fratio_upper
    )
    
    if(sum(round(an(catch(stka)),1)-round(an(Landings),1))>1){
      timeseries$Discards = round(an(catch(stka)),1)-round(an(Landings),1)
      
    }
    
    if(is.null(rps)){
      
    }
    
    endyr = which(max(yrs)==yrs)
    
    
    
    # TODO change system to relative
    refpts = round(t(data.frame(
      Ftgt = NA,
      Btgt = NA,
      Bthr= NA,
      Blim= NA,
      Fcur  = timeseries$F[endyr],
      Bcur  = timeseries$Biomass[endyr],
      B0.33= quantile(timeseries$Biomass,0.33),
      B0.66= quantile(timeseries$Biomass,0.66))),3)
    
    refpts = data.frame(RefPoint=row.names(refpts),Value=refpts[,1])
    rownames(refpts) = 1:nrow(refpts)
    
    
     if(!is.null(rps)){
      rps.tgt = rps[!rownames(rps)%in%c("Blim","Bpa","Bthr")]
      refpts[1,2] =  an(rps[rownames(rps)[grep("F",rownames(rps))]])[1]
      refpts[2,2] = an(rps.tgt[rownames(rps.tgt)[grep("B",rownames(rps.tgt))]])[1]
      refpts[3,2] = an(rps[rownames(rps)%in%c("Bpa","Bthr","Btri")])[1]
      refpts[4,2] = an(rps[rownames(rps)%in%c("Blim")])[1]
     }    
  
    return(list(timeseries=timeseries,refpts=refpts))
  }
  
  #' stock2ratios()
  #' @param object of class *FLStockR*  
  #' @return FLStockR with ratios F/Ftgt and B/Btgt
  #' @export
  stock2ratios <- function(object){
       
      stks= TRUE
      if(class(object)=="FLStockR"){
        object= FLStocks(stk=object)
        stks=FALSE
      }
      
      out =FLStocks(lapply(object,function(x){
      
    #x = object[[2]]
      
      B = as.FLQuant(unitSums(ssb(x))/x@refpts[[2]])
      H = unitMeans(fbar(x))/x@refpts[[1]]
      R = unitSums(rec(x))
      C = unitSums(computeCatch(x))
      dimnames(B)$age = "1"
      dimnames(C)$age = "1"
      dimnames(H)$age = "1"
      dimnames(R)$age = "1"
      year = an(dimnames(x)$year)
      iters = an(dimnames(x)$iter)
      
       stk = FLStockR(stock.n=FLQuant(R, dimnames=list(age="1", year = (year),iter=iters)),
      catch.n = C,
      landings.n = C,
      discards.n = FLQuant(0, dimnames=list(age="1", year = (year),iter=iters)),
      stock.wt=FLQuant(1, dimnames=list(age="1", year = (year),iter=iters)),
      landings.wt=FLQuant(1, dimnames=list(age="1", year = year,iter=iters)),
      discards.wt=FLQuant(1, dimnames=list(age="1", year = year,iter=iters)),
      catch.wt=FLQuant(1, dimnames=list(age="1", year = year,iter=iters)),
      mat = as.FLQuant(data.frame(age=1,year=year,unit="unique",
                                  season="all",area="unique",iter=iters,data=an(B/R))),
      
      #mat=B/R,
      m=FLQuant(0.0001, dimnames=list(age="1", year = year)),
      harvest = H,
      m.spwn = FLQuant(0, dimnames=list(age="1", year = year)),
      harvest.spwn = FLQuant(0.0, dimnames=list(age="1", year = year))
    )
    units(stk) = standardUnits(stk)
    stk@catch = unitSums(computeCatch(stk))
    stk@landings = unitSums(computeLandings(stk))
    stk@discards = unitSums(computeStock(stk))
    stk@stock = unitSums(ssb(x))
    br = c("Bthr","Blim","Bpa")
    stk@refpts = FLPar(Ftgt=1,Btgt=1)
    if(any(rownames(x@refpts)%in%br)){
      stk@refpts = rbind(stk@refpts,x@refpts[rownames(x@refpts)%in%br]/x@refpts[[2]])
      
    } 
    
    
    row.names(stk@refpts)[1:2] = row.names(x@refpts)[1:2] 
    stk@desc = x@desc
    return(stk)
    }))
    if(!stks){
      out = out[[1]]
    }
    
    return(out)
  }


  #' Function to summarise Type 1 GFCM fishing opportunities
  #' @param stk *FLStocks* with list of *FLStockR* objects 
  #' @param uncertainty *FLStocks* with list of *FLStockR* objects and iters 
  #' @param eval.yrs evaluation years of forecast 
  #' @param refyr change in percentage biomass to refyr
  #' @return data.frame
  #' @export
  fwd2type1<- function(stock,uncertainty,eval.yrs=NULL,dB=NULL,refyr=NULL){
    
    object= stock
    if(!class(object)=="FLStocks"){
      object = FLStocks(forecast=stock)
    }
    
    if(is.null(refyr)){
      refyr = an(range(object[[1]])["maxyear"])-1
    }
    
    if(is.null(eval.yrs)){
      eval.yrs = an(range(object[[1]])["maxyear"])
    }
    
    df = do.call(rbind,Map(function(x,y,z){
      stk = window(x,end=max(eval.yrs))
      flqs = FLQuants(
        Cy = round(catch(x)[,ac(eval.yrs)],1),
        Fy = round(fbar(x)[,ac(eval.yrs)],3),
        By = round(ssb(x)[,ac(eval.yrs)],1),
        dB = round(100*((ssb(x)[,ac(eval.yrs)])/(ssb(x)[,ac(eval.yrs-1)])-1),2),       
        PBlim = apply((round(ssb(z)[,ac(eval.yrs)],1)<Blim)*100,2,mean)
      )
      out = as.data.frame(flqs)
      data.frame(scenario=y, t(as.matrix(out$data)))
    },x=object,y=names(object),z=uncertainty)
    )
    names(df) = c("Basis",paste0("C_",eval.yrs),
                  paste0("F_",eval.yrs),
                  paste0("SSB_",eval.yrs),
                  "SSB_change","P(SSB<Blim)")
    
    
    
    rownames(df) = 1:nrow(df)
    return(df)  
  } # End of function
  
  
  
