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

Fsim <- function(brp,sigmaR=0.5,rho=0.,nyears=100,iters=1000,yrs.eval=NULL,verbose=TRUE){
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
  
  Prisk= mse::performance(run, metrics=list(SB=ssb), statistics=statistic, refpts=FLPar(SBlim=an(refpts(brp)["Blim","ssb"])), years=list(pyrs))  
  
  if(verbose) cat(paste0("Risk3 for ",ref," is P(SSB<Blim) = ",Prisk$data*100,"%"),"\n")
  
  out = list()
  out$params= FLPar(Fbrp=median(fbar(run)[,pyrs]),
                      MMY = median(landings(run)[,pyrs]),
                      Btrg= median(ssb(run)[,pyrs]),
                      Prisk=Prisk$data,styr=pyrs[1],endyr=pyrs[2])
  
  rownames(out$params)[1] = paste(ref)
  out$brp = brp
  out$sr =sr 
  out$devs =devs 
  out$stock = run
  
  
  return(out)
  
}  
#}}}

#{{{
# Fp05()
#
#' Calculates the Fbar value giving a maximum probability of ssb being below Blim of 5 percent
#'
#' @param object output from Fsim()
#' @param range range of Fbar value to be evaluated
#' @param iters Number of iterations, cannot exceed input object
#' @param verbose Should progress be shown, TRUE.
#' @return Fp.05
#' @export

Fp05 <- function(object,iters="missing",range="missing",verbose=TRUE){
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
}}  

if(!missing(iters)){
if(dims(stock)$iter < iters) stock = iter(stock,1:iters)
}
statistic <- list(FP05=list(~apply(iterMeans((SB/SBlim) < 1), c(1, 3:6), max),
                            name="P.05", desc="ICES P.05"))
res <- mse::bisect(stock, sr=object$sr, refpts=FLPar(SBlim=refpts(object$brp)["Blim","ssb"]), deviances=object$devs,
              metrics=list(SB=ssb), statistic=statistic, years=years, pyears=pyrs, 
              tune=list(fbar=range), prob=0.05, tol=0.01, verbose=verbose)

fp05 <- mean(fbar(res)[,100])


if(verbose) cat(paste0("Fp.05 = ",round(fp05,3)," is ",ifelse(fp05<Fbrp,"smaller","larger")," than Fbrp = ",round(Fbrp,3)),"\n")
return(fp05)
}
#}}}

