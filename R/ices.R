

#' Function to summarise forecast results
#' @param object *FLStocks* with list of *FLStockR* objects 
#' @param eval.yrs evaluation years of forecast 
#' @param rel if TRUE ratios B/Btgt and F/Ftgt are shown
#' @param dB computes change in percentage biomass to refyr
#' @return data.frame
#' @export
fwd2ices <- function(object,eval.yrs=NULL, rel=NULL,dB=NULL,refyr=NULL){

  
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
        flqs = FLQuants(c(flqs,FLQuants(Btgt=round(ssb(x)[,ac(eval.yrs)]/an(object[[1]]@refpts[2]),3))))
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
    
    if(any(c("Bmsy","Btgt")%in%rownames(object[[1]]@refpts))){
      nam = c(nam, paste0("B",eval.yrs,"/",names(object[[1]])[2]))
    }
    
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




#' ss2ices()
#' @param mvln output of ssmvln() 
#' @param output choice c("iters","mle")[1]
#' @param Fmsy if specified the ratio F/Fmsy is calculated (required for ensembles)
#' @param Btrigger if specified the ratio SSB/Btrigger is calculated (required for ensembles)
#' @param quantiles default is 95CIs as c(0.025,0.975)
#' @return ICES list with $timeseris and $refpts
#' @export
ss2ices <- function(mvln,quantiles = c(0.05,0.95)){
  
    kb = mvln$kb
    mle = mvln$mle
  
    quants=c("SSB","F","Catch","stock","harvest","Recr")
    
    mu = aggregate(cbind(stock,harvest,SSB,F,Catch,Recr)~year+run,kb,
                   quantile,c(0.5,quantiles))
    
    timeseries =  data.frame(year=mle$year,
                             Rec_low=mu[,quants[6]][,2],
                             Rec=mle$Recr,
                             Rec_high=mu[,quants[6]][,3],
                             SSB_low=mu[,quants[1]][,2],
                             SSB=mle$SSB,
                             SSB_high=mu[,quants[1]][,3],
                             
                             
                             Landings=mle$Landings,
                             Discards=mle$Discards,
                             
                             F_low=mu[,quants[2]][,2],
                             F=mle$F,
                             F_upp=mu[,quants[2]][,3]
                             
                        
    )
    
    timeseries[,-1] =     round(timeseries[,-1],3)  
    if("forecast"%in%mle$type){
      timeseries[mle$type=="forecast",-c(1:7)] = NA 
    }
  
  return(timeseries)
}

#' applyAR()
#' Function to apply ICES Advice rule for F < MSY Btrigger
#' @param b biomass
#' @param fmsy Fmsy target
#' @param btrigger if b below btrigger F is reduced
#' @param fmin  F below bmin
#' @param bmin  Option for de facto fishing closure (e.g. bmin = blim)
#' @return F advice 
#' @export
applyAR <- function(b, btrigger, fmsy, bmin=0,fmin=0.001){
  
  # BELOW lim
  out <- c(ifelse(b <= bmin, fmin,ifelse(b< btrigger,
                         (b - bmin) * ((fmsy - fmin) / 
                                          (btrigger - bmin)) + fmin,fmsy)))
 return(out)
}
  

#' fwdF4B ()
#' Bisection Function to search F for a given biomass
#' @param stock FLStock object 
#' @param sr stock recruitment relationship
#' @param btgt target biomass of the search
#' @param nfy number of forecast (default 3)
#' @param niy number of intermediate years (default 1)
#' @param ival intermediate year value
#' @param imet intermediate year metric ("TAC","F")
#' @param ftune tuning limits for F search c(0.1,2) 
#' @param tol precision tolerance
#' @param verbose 
#' @return  FLStock
#' @export

fwdF4B =  function(stock,sr,btgt,nfy=3,niy=1,ival=NULL,imet="TAC",ftune=c(0,2),tol=0.0001,verbose=TRUE){
  rpts = NULL
  if(class(stock)=="FLStockR"){
    rpts = stock@refpts
    stock = as(stock,"FLStock")
  }
  
  fyrs = (dims(stock)$maxyear+1):(dims(stock)$maxyear+nfy) 
  nfy = length(fyrs)
  stkf = stf(stock,nfy)
  iyrs = (dims(stock)$maxyear+1):(dims(stock)$maxyear+niy)
  if(niy>0){
    val = ifelse(length(ival)==niy,ival,rep(ival[1],niy))
    ictrl <- fwdControl(data.frame(year = iyrs,
                                   quant = ifelse(imet=="TAC","catch","f"),
                                   value = val))
    stkf =fwd(stkf,sr=sr,control = ictrl)
  }
  
  
  statistic <- list(FBtgt=list(~yearMeans((SB/SBtgt) < 1), name="F4B",
                              desc="Search Btgt"))
  out <- bisect(stkf, sr=sr, metrics=list(SB=ssb), 
                    refpts=FLPar(SBtgt=btgt), statistic=statistic, years=fyrs[-1],pyears=max(fyrs), tune=list(fbar=ftune), prob=0.5,tol=0.0001,verbose = verbose)
  if(!is.null(rpts)){
    out = FLStockR(out)
    out@refpts = rpts
  }
  return(out)
}

#' Function to summarise ICES fishing opportunities
#' @param stk *FLStocks* with list of *FLStockR* objects 
#' @param uncertainty *FLStocks* with list of *FLStockR* objects and iters 
#' @param eval.yrs evaluation years of forecast 
#' @param refyr change in percentage biomass to refyr
#' @return data.frame
#' @export
fwd2ices<- function(stock,uncertainty,eval.yrs=NULL,dB=NULL,refyr=NULL){
  
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
        Cy = round(catch(x)[,ac(eval.yrs-1)],1),
        Fy = round(fbar(x)[,ac(eval.yrs-1)],3),
        By = round(ssb(x)[,ac(eval.yrs)],1),
        dB = 100*(round(ssb(x)[,ac(eval.yrs)],1)/round(ssb(x)[,ac(eval.yrs-1)],1)-1),       
        PBlim = apply((round(ssb(z)[,ac(eval.yrs)],1)<Blim)*100,2,mean)
        )
      out = as.data.frame(flqs)
      data.frame(scenario=y, t(as.matrix(out$data)))
    },x=object,y=names(object),z=uncertainty)
    )
    names(df) = c("Basis",paste0("C_",eval.yrs-1),
                  paste0("F_",eval.yrs-1),
                  paste0("SSB_",eval.yrs),
                  "SSB_change","P(SSB<Blim)")
  

  
  rownames(df) = 1:nrow(df)
  return(df)  
} # End of function


  

