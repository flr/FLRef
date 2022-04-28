#{{{
#' Fsim()
#
#' Simulates stochastic stock dynamics under under constant Fbrp
#'
#' @param brp output object from computeFbrp() of class FLBRP 
#' @param sigmaR lognormal recruitment standard deviation
#' @param rho AR1 recruitment autocorrelation coefficient 
#' @param nyears number of simulation years
#' @param iters number simulation iterations
#' @param yrs.eval last years to be used evaluation period, default nyears/2
#' @param verbose cat comments
#' @return list of  FLPar, FLStock and FLBRP objects
#' @export

Fsim <- function(brp,sigmaR=0.5,rho=0.,nyears=100,iters=250,yrs.eval=NULL,verbose=TRUE){
  fbar(brp)[] = 0.01
  stock=as(brp,"FLStock")
  ref = rownames(refpts(brp))[grep("F",rownames(refpts(brp)))][1]
  if(nyears<100){
    stock = window(stock,end=nyears+1)
  }
  if(nyears>100){
    stock = stf(stock,nyears-100+1)
  }
  if(is.null(yrs.eval)) yrs.eval=ceiling(nyears/2)
  sr = as(brp,"FLSR")
  stock = propagate(stock,iters)
  stock.n(stock)[1:(dim(stock)[1]-1),1] = stock.n(stock)[1:(dim(stock)[1]-1),1]*rlnorm(iters*(dim(stock)[1]-1),0-0.5*sigmaR^2,sigmaR)
  Fbrp = an(refpts(brp)[ref,"harvest"])
  devs = ar1rlnorm(rho=rho, years=1:(nyears+1), iters=iters,meanlog= 0, sdlog=sigmaR)
  
  
  run = ffwd(stock, sr=sr,
           fbar=FLQuant(Fbrp, dimnames=list(year=2:(nyears+1))),
           deviances=devs)
  run=window(run,start=2,end=(nyears+1))
  #range(run)[4:5]=c(1,nyears)
  
  run@landings = computeLandings(run)
  run@discards = computeDiscards(run)
  run@catch = computeCatch(run)
  
  statistic <- list(FP05=list(~apply(iterMeans((SB/SBlim) < 1), c(1, 3:6), max),name="Prisk", desc="ICES Prisk"))
  pyrs = ((nyears-yrs.eval+1):nyears)
  
  Prisk= mean(mse::performance(run, metrics=list(SB=ssb), statistics=statistic, refpts=FLPar(SBlim=an(refpts(brp)["Blim","ssb"])), years=pyrs)$data)
  
  if(verbose) cat(paste0("Risk3 for ",ref," is P(SSB<Blim) = ",Prisk*100,"%"),"\n")
  
  out = list()
  out$params= FLPar(Fbrp=median(fbar(run)[,pyrs]),
                      MMY = median(landings(run)[,pyrs]),
                      Btrg= median(ssb(run)[,pyrs]),
                      Prisk=Prisk,styr=pyrs[1],endyr=tail(pyrs,1))
  
  rownames(out$params)[1] = paste(ref)
  out$brp = brp
  out$sr =sr 
  out$devs =devs 
  out$stock = run
  
  
  return(out)
  
}  
#}}}

#{{{
# computeFp05()
#
#' Calculates the Fbar value giving a maximum probability of ssb being below Blim of 5 percent
#'
#' @param object output from Fsim()
#' @param range range of Fbar value to be evaluated
#' @param iters Number of iterations, cannot exceed input object
#' @param verbose Should progress be shown, TRUE.
#' @return Fp.05
#' @export

computeFp05 <- function(object,iters="missing",range="missing",tol=0.001,maxit=15,verbose=TRUE){
stock = object$stock 
years = (dims(stock)$minyear:dims(stock)$maxyear)[-1]
pyrs = an(object$params["styr"]) :an(object$params["endyr"])
Fbrp = an(object$params[1])
Flim = an(refpts(object$brp)["Blim","harvest"])
Prisk = an(object$params["Prisk"])

if(missing(range)){
if(Prisk<0.05){ 
  range = c(0.8*Fbrp,1.1*Flim)} else {
  range = c(0.1*Fbrp,1.1*Fbrp)
  }
  }  

if(!missing(iters)){
if(dims(stock)$iter < iters) stock = iter(stock,1:iters)
}
statistic <- list(FP05=list(~apply(iterMeans((SB/SBlim) < 1), c(1, 3:6), max),
                            name="P.05", desc="ICES P.05"))
run <- mse::bisect(stock, sr=object$sr, refpts=FLPar(SBlim=refpts(object$brp)["Blim","ssb"]), deviances=object$devs,
              metrics=list(SB=ssb), statistic=statistic, years=years, pyears=pyrs, 
              tune=list(fbar=range), prob=0.05, tol=tol, verbose=verbose)


Prisk= mean(mse::performance(run, metrics=list(SB=ssb), statistics=statistic, refpts=FLPar(SBlim=an(refpts(brp)["Blim","ssb"])), years=pyrs)$data)


out = list()
out$params= FLPar(Fp05=median(fbar(run)[,ac(pyrs)]),
                  Yp05 = median(landings(run)[,ac(pyrs)]),
                  Btrg= median(ssb(run)[,ac(pyrs)]),
                  Prisk=Prisk,styr=pyrs[1],endyr=tail(pyrs,1))

out$brp = brp
out$sr =sr 
out$devs =devs 
out$stock = run
return(out)  

fp05 <- mean(fbar(run)[,ac(100)])


if(verbose) cat(paste0("Fp.05 = ",round(fp05,3)," is ",ifelse(fp05<Fbrp,"smaller","larger")," than Fbrp = ",round(Fbrp,3)),"\n")
return(out)
}
#}}}


# opt.bisect {{{

#' Bisection approach to optimise x for maximising y  
#'
#'
#' The plain bisection algorithm (Burden & Douglas, 1985) is employed here to
#' find the value of a given forecast target quantity (e.g. `fbar`) for which
#' a selected value of a performance statistic is obtained over a chosen period.
#' @references {Burden, Richard L.; Faires, J. Douglas (1985), "2.1 The Bisection Algorithm", Numerical Analysis (3rd ed.), PWS Publishers, ISBN 0-87150-857-5}
#' @param stock object class FLStock
#' @param sr object class FLSR
#' @param metrics FLQuant of FLStock to be defined  
#' @param statistic 
#' @param years years to be evaluated
#' @param tune range for input x
#' @param tol tolerance level
#' @param maxit number of optimisation steps
#' @param log if TRUE, optimise on log-scale  
#' @author Credits to Iago Mosqueira
#' @examples
#' data(ple4)
#' stock <- propagate(stf(ple4, end=2118), 200)
#' srr <- predictModel(model=rec ~ ifelse(ssb <= b, a * ssb, a * b), params=FLPar(a=1.29, b=1.35e+06))
#' # GENERATE SRR deviances
#' devs <- ar1rlnorm(rho=0.4, 2018:2118, iters=200, meanlog=0, sdlog=0.5)
#' # DEFINE MMY statistic
#' statistic <- list(MMY=list(~apply(L,1,median), name="MMY",
#'   desc="ICES Maximum Median Yield"))
#' # CALL bisect over 100 years, Fmmy calculated over last 50.
#' fmmy <- opt.bisect(stock, sr=srr, deviances=devs, metrics=list(L=landings), 
#' statistic=statistic, years=2018:2118,
#' pyears=2069:2118, tune=list(fbar=c(0.01, 0.2)))
#' # fmmy
#' mean(fbar(fmmy)[,ac(2069:2118)])  

opt.bisect <- function(stock, sr, deviances=rec(stock) %=% 1, metrics,
                       statistic, years, pyears=years, tune, tol=0.001, maxit=15,log=TRUE, verbose=TRUE) {
  
  # CHOOSE method
  
  if(names(tune)[1] %in% c("f", "fbar"))
    foo <- ffwd
  else
    foo <- fwd
  
  # --- RUN at min
  
  cmin <- fwdControl(year=years, quant=names(tune)[1], value=unlist(tune)[1])
  
  # PRINT at top
  if(verbose)
    cat(paste0("[1] ", names(tune), ": ", unlist(tune)[1]))
  rmin <- foo(stock, sr=sr, control=cmin, deviances=deviances)
  rmin@landings = computeLandings(rmin) # TODO add to foo
  rmin@discards = computeDiscards(rmin) # TODO add to foo
  
  pmin <- performance(rmin, metrics=metrics, 
                      statistics=statistic, years=pyears)
  
  obmin <- mean(pmin$data, na.rm=TRUE)
  if(log) obmin = log(obmin) 
  
  if(verbose)
    cat(" Min range value", mean(pmin$data, na.rm=TRUE),"\n")
  
  
  
  # --- RUN at max
  
  cmax <- fwdControl(year=years, quant=names(tune)[1], value=unlist(tune)[2])
  
  # PRINT at top
  if(verbose)
    cat(paste0("[2] ", names(tune), ": ", unlist(tune)[2]))
  
  rmax <- foo(stock, sr=sr, control=cmax, deviances=deviances)
  rmax@landings = computeLandings(rmax)
  rmax@discards = computeDiscards(rmax)
  
  
  pmax <- performance(rmax, metrics=metrics, 
                      statistics=statistic, years=pyears)
  
  obmax <- mean(pmax$data, na.rm=TRUE)
  if(log) obmax = log(obmax) 
  
  if(verbose)
    cat(" Max range value:", mean(pmax$data, na.rm=TRUE),"\n")
  
  # CHECK cmax result
  
  # --- LOOP bisecting
  maxi = max(obmin,obmax)
  
  count <- 0
  while(count <= maxit) {
    
    # RUN at mid
    cmid <- control
    cmid <- fwdControl(year=years, quant=names(tune)[1],
                       value=(cmin$value + cmax$value) / 2)
    
    # PRINT at mid
    if(verbose)
      cat(paste0("[", count + 3, "] ", names(tune), ": ", cmid$value[1]))
    
    rmid <- foo(stock, sr=sr, control=cmid, deviances=deviances)
    rmid@landings = computeLandings(rmid)
    rmid@discards = computeDiscards(rmid)
    
    
    pmid <- performance(rmid, metrics=metrics, 
                        statistics=statistic, years=pyears)
    obmid <- mean(pmid$data, na.rm=TRUE)
    if(log) obmid=log(obmid)
    
    if(verbose)
      cat(" new value:", mean(pmid$data, na.rm=TRUE), " - diff: ", obmid-maxi, "\n")
    
    # CHECK and RETURN cmid result
    if(isTRUE(all.equal(obmid-maxi, 0, tolerance=tol))) {
      return(rmid)
    }
    
    
    # TEST LEFT
    if(obmin > obmax) {
      
      # SET max as new mid
      cmax <- cmid
      obmax <- obmid
      if(isTRUE(all.equal(cmin$value[1], cmid$value[1], tolerance=tol))) {
        return(rmid)
      }
    } else {
      
      # SET min as new mid
      cmin <- cmid
      obmin <- obmid
      if(isTRUE(all.equal(cmid$value[1], cmax$value[1], tolerance=tol))) {
        return(rmid)
      }
    }
    maxi = max(obmin,obmax)
    
    
    count <- count + 1
  }
  
  warning("Solution not found within 'maxit', check 'range', 'maxit' or 'tol'.")
  
  return(rmid)
  
} # }}}

#{{{
#' computeFmmy()
#
#' Uses opt.bisect to derive the F at Maximum Median Yield from stochastic simulations  
#'
#' @param brp output object from computeFbrp() of class FLBRP 
#' @param sigmaR lognormal recruitment standard deviation
#' @param rho AR1 recruitment autocorrelation coefficient 
#' @param nyears number of simulation years
#' @param iters number simulation iterations
#' @param yrs.eval last years to be used evaluation period, default nyears/2
#' @param range range of Fbar value to be evaluated
#' @param tol tolerance
#' @param maxit number of steps
#' @param verbose cat comments
#' @return list of  FLPar, FLStock and FLBRP objects
#' @export
#' @examples 
#' data(ple4)
#' srr = srrTMB(as.FLSR(ple4,model=segreg),spr0=spr0y(ple4),plim=0.15)
#' brp=computeFbrp(ple4,st,proxy = c("msy","sprx"),x=35,blim=params(st)[[2]])
#' ploteq(brp)
#' # DEFINE MMY statistic
#' statistic <- list(MMY=list(~apply(L,1,median), name="MMY",
#'   desc="ICES Maximum Median Yield"))
#' # CALL bisect over 100 years, Fmmy calculated over last 50.
#' fmmy <- opt.bisect(stock, sr=srr, deviances=devs, metrics=list(L=landings), 
#' statistic=statistic, years=2018:2118,
#' pyears=2069:2118, tune=list(fbar=c(0.01, 0.2)))
#' # fmmy
#' mean(fbar(fmmy)[,ac(2069:2118)]) 
#' 
#' 
computeFmmy <- function(brp,sigmaR=0.5,rho=0.0,nyears=100,iters=250,yrs.eval=NULL,range="missing",tol=0.01,maxit=15,verbose=TRUE){
  fbar(brp)[] = 0.01
  stock=as(brp,"FLStock")
  sr = as(brp,"FLSR")
  fbrp = computeFbrp(stock,sr,proxy="msy",blim=an(Fbrp(brp)["Blim"]),verbose=FALSE)
  
  ref = rownames(refpts(brp))[grep("F",rownames(refpts(brp)))][1]
  if(nyears<100){
    stock = window(stock,end=nyears+1)
  }
  if(nyears>100){
    stock = stf(stock,nyears-100+1)
  }
  if(is.null(yrs.eval)) yrs.eval=ceiling(nyears/2)
  pyrs = ((nyears-yrs.eval+1):nyears)
  stock = window(stock,start=1,end=(nyears))
  
  stock = propagate(stock,iters)
  stock.n(stock)[1:(dim(stock)[1]-1),1] = stock.n(stock)[1:(dim(stock)[1]-1),1]*rlnorm(iters*(dim(stock)[1]-1),0-0.5*sigmaR^2,sigmaR)
  Fbrp = an(refpts(brp)[ref,"harvest"])
  devs = ar1rlnorm(rho=rho, years=2:(nyears), iters=iters,meanlog= 0, sdlog=sigmaR)
  statistics <- list(MMY=list(~apply(L,1,median), name="MMY",desc="ICES Maximum Median Yield"))
  if(missing(range)){
    range = an(c(0.2*Fbrp(fbrp)[[1]],Fbrp(fbrp)[[1]]*1.3))
  }
  run  = opt.bisect(stock, sr, deviances=devs, metrics=list(L=landings),
                    statistic=statistics, years=2:(nyears), pyears=pyrs, tune=list(fbar=range), tol=tol, maxit=maxit, verbose=TRUE,log=TRUE) 
  
  run= window(res,start=2,end=(nyears))
  
  if(verbose) cat(paste0("Fmmy = ",round(median(fbar(run)[,ac(pyrs)]),3)," vs Fmsy = ",round(Fbrp(fbrp)[[1]],3)),"\n")
  
  statistic <- list(FP05=list(~apply(iterMeans((SB/SBlim) < 1), c(1, 3:6), max),name="Prisk", desc="ICES Prisk"))
  
  Prisk= mean(mse::performance(run, metrics=list(SB=ssb), statistics=statistic, refpts=FLPar(SBlim=an(refpts(brp)["Blim","ssb"])), years=pyrs)$data)
  
  if(verbose) cat(paste0("Risk3 for Fmmy = ",round(median(fbar(run)[,ac(pyrs)]),3)," is P(SSB<Blim) = ",Prisk*100,"%"),"\n")
  
  
  out = list()
  out$params= FLPar(Fmmy=median(fbar(run)[,ac(pyrs)]),
                    MMY = median(landings(run)[,ac(pyrs)]),
                    Btrg= median(ssb(run)[,ac(pyrs)]),
                    Prisk=Prisk,styr=pyrs[1],endyr=tail(pyrs,1))
  
  out$brp = brp
  out$sr =sr 
  out$devs =devs 
  out$stock = run
  return(out)  
  
}  
#}}}

#' getF()
#' 
#' Helper functio to extract F from various FLRef output
#' @param x output object from computeFbrp() of class FLBRP 
#' @export 
getF <- function(x){
  if(class(x)=="FLPar") y = an(x[[1]])
  if(class(x)=="list") y = an(x$params[[1]])
  if(class(x)=="FLBRP") y = an(Fbrp(x)[[1]])
  return(y)
}




