
#{{{
#' Fsim()
#
#' Simulates stochastic stock dynamics under under constant Fbrp
#'
#' @param brp output object from computeFbrp() of class `FLBRP` 
#' @param sigmaR lognormal recruitment standard deviation
#' @param rho AR1 recruitment autocorrelation coefficient 
#' @param nyears number of simulation years
#' @param iters number simulation iterations
#' @param yrs.eval last years to be used evaluation period, default nyears/2
#' @param verbose
#' @return stock object `FLStock` with iterations 
#' @export

Fsim <- function(brp,sigmaR=0.5,rho=0.,nyears=100,iters=1000,yrs.eval=NULL,verbose=TRUE){
  stock=as(brp,"FLStock")
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
  Fbrp = an(refpts(brp)["Fbrp","harvest"])
  run = FLasher::ffwd(stock, sr=sr,
           fbar=FLQuant(Fbrp, dimnames=list(year=2:(nyears+1))),
           deviances=ar1rlnorm(rho=rho, years=1:(nyears+1), iters=iters,meanlog= 0, sdlog=sigmaR))
  run=window(run,start=2,end=(nyears+1))
  range(run)[4:5]=c(1,nyears)
  
  statistic <- list(FP05=list(~apply(iterMeans((SB/SBlim) < 1), c(1, 3:6), max),name="Prisk", desc="ICES Prisk"))
  pyrs = (nyears-yrs.eval+1):nyears
  
  Prisk= mse::performance(run, metrics=list(SB=ssb), statistics=statistic, refpts=FLPar(SBlim=an(refpts(brp)["Blim","ssb"])), years=list(pyrs))  
  
  if(verbose) cat(paste0("Risk3 P(SSB<Blim) = ",Prisk$data*100,"%"))
  
  
  return(run)
  
}  
#}}}



