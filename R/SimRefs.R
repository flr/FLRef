# F sim {{{

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
#' @examples
#' data(ple4)
#' hs = srrTMB(as.FLSR(ple4,model=segreg),spr0=spr0y(ple4),lplim=0.05,uplim=0.25)
#' blim = params(hs)[[2]]
#' brp = computeFbrp(ple4,hs,proxy=c("sprx","f0.1","msy"),x=40,blim=blim)
#' ploteq(brp)
#' fsim = Fsim(brp,sigmaR=0.7,rho=0.3)
#' plotFsim(fsim)
#' plotFsim(fsim,panels=2)

Fsim <- function(brp,Ftgt=NULL,sigmaR=0.5,rho=0.,nyears=100,iters=250,yrs.eval=NULL,verbose=TRUE){
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
  if(is.null(Ftgt)){
  Fbrp = an(refpts(brp)[ref,"harvest"])
  } else {
  Fbrp = Ftgt  
  }
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
                      Btgt= median(ssb(run)[,pyrs]),
                      Prisk=Prisk,styr=pyrs[1],endyr=tail(pyrs,1))
  
  rownames(out$params)[1] = paste(ref)
  out$brp = brp
  out$sr =sr 
  out$devs =devs 
  out$stock = run
  
  
  return(out)
  
}  
# }}}

# Fp05 {{{

# Fp05()
#
#' Calculates the Fbar value giving a maximum probability of ssb being below Blim of 5 percent
#'
#' @param object output from Fsim()
#' @param range range of Fbar value to be evaluated
#' @param iters Number of iterations, cannot exceed input object
#' @param verbose Should progress be shown, TRUE.
#' @return list
#' @export
#' @examples 
#' data(ple4)
#' bh = srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=spr0y(ple4))
#' brp = computeFbrp(ple4,bh,proxy="bx",x=35,blim=0.2) # set Blim higher
#' fsim = Fsim(brp,sigmaR=0.7,rho=0.3,iters=500)
#' plotFsim(fsim)
#' fp.05 = Fp05(fsim)
#' plotFsim(fp.05,panels=c(2,4)) # black line is Fp0.05
#' getF(fp.05) 

Fp05 <- function(object,iters="missing",range="missing",tol=0.001,maxit=20,verbose=TRUE){
stock = object$stock 
years = (dims(stock)$minyear:dims(stock)$maxyear)[-1]
pyrs = an(object$params["styr"]) :an(object$params["endyr"])
Fbrp = an(object$params[1])
Flim = an(refpts(object$brp)["Blim","harvest"])
Prisk = an(object$params["Prisk"])
brp =  object$brp
sr = object$sr 


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
run <- bisect(stock, sr=object$sr, refpts=FLPar(SBlim=refpts(object$brp)["Blim","ssb"]), deviances=object$devs,
              metrics=list(SB=ssb), statistic=statistic, years=years, pyears=pyrs, 
              tune=list(fbar=range), prob=0.05, tol=tol, verbose=verbose)


Prisk= mean(mse::performance(run, metrics=list(SB=ssb), statistics=statistic, refpts=FLPar(SBlim=an(refpts(object$brp)["Blim","ssb"])), years=pyrs)$data)


out = list()
out$params= FLPar(Fp05=median(fbar(run)[,ac(pyrs)]),
                  Yp05 = median(landings(run)[,ac(pyrs)]),
                  Btgt= median(ssb(run)[,ac(pyrs)]),
                  Prisk=Prisk,styr=pyrs[1],endyr=tail(pyrs,1))

out$brp = object$brp
out$sr =object$sr 
out$devs =object$devs
out$stock = run


fp05 <- mean(fbar(run)[,ac(100)])
if(verbose) cat(paste0("Fp.05 = ",round(fp05,3)," is ",ifelse(fp05<Fbrp,"smaller","larger")," than Fbrp = ",round(Fbrp,3)),"\n")

return(out)
}
#}}}

# ref.bisect {{{

#' Bisection approach to find target F for B (sex-structured)
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
#' @param target target with default 1 for ratios
#' @param tol tolerance level
#' @param maxit number of optimisation steps
#' @param log if TRUE, optimise on log-scale  
#' @export
#' @author Credits to Iago Mosqueira

ref.bisect = function (stock, sr, deviances = rec(stock) %=% 1, metrics, refpts, 
                       statistic, years, pyears = years, tune, target=1, tol = 0.01, 
                       maxit = 15, verbose = TRUE) 
{
  if (names(tune)[1] %in% c("f", "fbar")) 
    foo <- ffwd
  else foo <- fwd
  cmin <- fwdControl(year = years, quant = names(tune)[1], 
                     value = unlist(tune)[1])
  if (verbose) 
    cat(paste0("[1] ", names(tune), ": ", unlist(tune)[1]))
  rmin <- foo(stock, sr = sr, control = cmin, deviances = deviances)
  pmin <- performance(rmin, metrics = metrics, statistics = statistic, 
                      refpts = refpts, years = pyears)
  obmin <- (target-mean(pmin$data, na.rm = TRUE))
  if (verbose) 
    cat(" - target:", mean(pmin$data, na.rm = TRUE), " - diff: ", 
        obmin, "\n")
  if (isTRUE(all.equal(obmin, 0, tolerance = tol))) 
    return(rmin)
  cmax <- fwdControl(year = years, quant = names(tune)[1], 
                     value = unlist(tune)[2])
  if (verbose) 
    cat(paste0("[2] ", names(tune), ": ", unlist(tune)[2]))
  rmax <- foo(stock, sr = sr, control = cmax, deviances = deviances)
  pmax <- performance(simplify(rmax,weighted=TRUE), metrics = metrics, statistics = statistic, 
                      refpts = refpts, probs = NULL, years = pyears)
  obmax <- (target -  mean(pmax$data, na.rm = TRUE))
  if (verbose) 
    cat(" - target:", mean(pmax$data, na.rm = TRUE), " - diff: ", 
        obmax, "\n")
  if (isTRUE(all.equal(obmax, 0, tolerance = tol))) 
    return(rmax)
  if ((obmin * obmax) > 0) {
    warning("Range of hcr param(s) cannot achieve requested tuning objective probability")
    return(list(min = rmin, max = rmax))
  }
  count <- 0
  while (count <= maxit) {
    cmid <- control
    cmid <- fwdControl(year = years, quant = names(tune)[1], 
                       value = (cmin$value + cmax$value)/2)
    if (verbose) 
      cat(paste0("[", count + 3, "] ", names(tune), ": ", 
                 cmid$value[1]))
    rmid <- foo(stock, sr = sr, control = cmid, deviances = deviances)
    pmid <- performance(rmid, metrics = metrics, statistics = statistic, 
                        refpts = refpts, probs = NULL, years = pyears)
    obmid <- target-mean(pmid$data, na.rm = TRUE) 
    if (verbose) 
      cat(" - target:", mean(pmid$data, na.rm = TRUE), " - diff: ", 
          obmid, "\n")
    if (isTRUE(all.equal(obmid, 0, tolerance = tol))) {
      return(rmid)
    }
    if ((obmin * obmid) < 0) {
      cmax <- cmid
      obmax <- obmid
      if (isTRUE(all.equal(cmin$value[1], cmid$value[1], 
                           tolerance = tol))) {
        return(rmid)
      }
    }
    else {
      cmin <- cmid
      obmin <- obmid
      if (isTRUE(all.equal(cmid$value[1], cmax$value[1], 
                           tolerance = tol))) {
        return(rmid)
      }
    }
    count <- count + 1
  }
  warning("Solution not found within 'maxit', check 'range', 'maxit' or 'tol'.")
  return(rmid)
}
# }}}

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
#' @export
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

# Fmmy {{{

#' Fmmy()
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
#' bh = srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=spr0y(ple4))
#' brp = computeFbrp(ple4,bh,proxy=c("bx","msy"),x=35,blim=0.1)
#' fmmy = Fmmy(brp,sigmaR=0.7,rho=0.3)
#' getF(fmmy) # FMMY value 
#' plotFsim(fmmy)
#' brpfmmy = computeFbrp(ple4,bh,proxy=getF(fmmy),blim=0.1)
#' fsim = Fsim(brpfmmy,sigmaR=0.7,rho=0.3)
#' plotFsim(fsim)

Fmmy <- function(brp,sigmaR=0.5,rho=0.0,nyears=100,iters=250,yrs.eval=NULL,range="missing",tol=0.001,maxit=15,verbose=TRUE){
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
  devs = rlnormar1 (rho=rho, year=2:(nyears), n=iters,meanlog= 0, sdlog=sigmaR)
  statistics <- list(MMY=list(~apply(L,1,median), name="MMY",desc="ICES Maximum Median Yield"))
  if(missing(range)){
    range = an(c(0.2*Fbrp(fbrp)[[1]],Fbrp(fbrp)[[1]]*1.3))
  }
  run  = opt.bisect(stock, sr, deviances=devs, metrics=list(L=landings),
                    statistic=statistics, years=2:(nyears), pyears=pyrs, tune=list(fbar=range), tol=tol, maxit=maxit, verbose=TRUE,log=TRUE) 
  
  run= window(run,start=2,end=(nyears))
  
  if(verbose) cat(paste0("Fmmy = ",round(median(fbar(run)[,ac(pyrs)]),3)," vs Fmsy = ",round(Fbrp(fbrp)[[1]],3)),"\n")
  
  statistic <- list(FP05=list(~apply(iterMeans((SB/SBlim) < 1), c(1, 3:6), max),name="Prisk", desc="ICES Prisk"))
  
  Prisk= mean(mse::performance(run, metrics=list(SB=ssb), statistics=statistic, refpts=FLPar(SBlim=an(refpts(brp)["Blim","ssb"])), years=pyrs)$data)
  
  if(verbose) cat(paste0("Risk3 for Fmmy = ",round(median(fbar(run)[,ac(pyrs)]),3)," is P(SSB<Blim) = ",Prisk*100,"%"),"\n")
  
  
  out = list()
  out$params= FLPar(Fmmy=median(fbar(run)[,ac(pyrs)]),
                    MMY = median(landings(run)[,ac(pyrs)]),
                    Btgt= median(ssb(run)[,ac(pyrs)]),
                    Prisk=Prisk,styr=pyrs[1],endyr=tail(pyrs,1))
  
  out$brp = brp
  out$sr =sr 
  out$devs =devs 
  out$stock = run
  return(out)  
  
}  
# }}}

# getF {{{

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
# }}}

# bisect {{{

#' Bisection search for a forecast target matching a performance statistic
#'
#' Uses a bisection algorithm to find the value of a forecast control quantity
#' that gives a specified target value for a performance statistic. Typical use
#' cases include finding the constant \code{fbar} that gives a chosen risk level,
#' such as \eqn{P(SSB < B_{lim}) = 0.05}, over a projection period.
#'
#' The function evaluates the forecast at the lower and upper bounds supplied in
#' \code{tune}. If the requested probability lies between the two outcomes, the
#' interval is repeatedly bisected until the performance statistic is within
#' \code{tol} of \code{prob}, or until \code{maxit} iterations are reached.
#'
#' @param stock An \code{FLStock} object used as the starting point for the
#'   forecast. The object can be propagated over iterations before calling
#'   \code{bisect()}.
#' @param sr A stock-recruitment object or model passed to \code{fwd()} or
#'   \code{ffwd()} through the \code{sr} argument.
#' @param deviances Recruitment deviances used in the stochastic forecast.
#'   Defaults to \code{rec(stock) \%=\% 1}, i.e. deterministic recruitment
#'   multipliers of one. For stochastic projections this is usually an
#'   \code{FLQuant} of lognormal or autocorrelated recruitment deviances.
#' @param metrics A named list of metric functions passed to
#'   \code{performance()}, for example \code{list(SB = ssb)}.
#' @param refpts Reference points passed to \code{performance()}, typically an
#'   \code{FLPar} object such as \code{FLPar(SBlim = 150000)}.
#' @param statistic A named list defining the performance statistic to evaluate,
#'   passed to \code{performance()}. For example, a risk statistic based on
#'   \code{SB / SBlim < 1}.
#' @param years Projection years over which the forecast control is applied.
#' @param pyears Years over which the performance statistic is evaluated.
#'   Defaults to \code{years}. This is often a subset of the projection period,
#'   for example the final 50 years of a long stochastic forecast.
#' @param tune A named numeric vector or list giving the forecast quantity to
#'   tune and its lower and upper bounds. For example,
#'   \code{list(fbar = c(0.1, 1.0))}. If the tuned quantity is named
#'   \code{"f"} or \code{"fbar"}, \code{ffwd()} is used; otherwise
#'   \code{fwd()} is used.
#' @param prob Numeric. Target value of the performance statistic. For example,
#'   \code{0.05} for a 5 percent risk threshold.
#' @param tol Numeric. Absolute tolerance used to decide whether the target
#'   statistic has been reached. Default is \code{0.01}.
#' @param maxit Integer. Maximum number of bisection iterations. Default is
#'   \code{15}.
#' @param verbose Logical. If \code{TRUE}, print the value of the tuned control,
#'   the resulting statistic, and the difference from the target at each
#'   bisection step.
#'
#' @return If a solution is found, an \code{FLStock} object returned by
#'   \code{fwd()} or \code{ffwd()} for the tuned control value. If the requested
#'   target is already achieved at the lower or upper bound, the corresponding
#'   forecast object is returned. If the supplied range does not bracket the
#'   target statistic, a list with elements \code{min} and \code{max} is
#'   returned, containing the forecasts at the lower and upper bounds.
#'
#' @details
#' The bisection search requires that the target objective is bracketed by the
#' two values supplied in \code{tune}. In practice, this means that:
#'
#' \deqn{
#' \left[g(x_{min}) - p\right]
#' \left[g(x_{max}) - p\right] < 0
#' }
#'
#' where \eqn{g(x)} is the performance statistic produced by forecasting with
#' control value \eqn{x}, and \eqn{p} is the target value supplied by
#' \code{prob}. If both bounds give performance statistics on the same side of
#' the target, the function cannot identify a unique bracketed solution.
#'
#' The function is most useful for risk-based advice calculations, such as
#' finding the fishing mortality that gives a pre-specified probability of
#' falling below \code{Blim} over a future evaluation period.
#'
#' @references
#' Burden, R. L. and Faires, J. D. (1985). Numerical Analysis. 3rd edition.
#' PWS Publishers. Section 2.1: The Bisection Algorithm.
#'
#' @examples
#' \dontrun{
#' data(ple4)
#'
#' # Extend and propagate stock
#' stock <- propagate(stf(ple4, end = 2118), 100)
#'
#' # Define a Beverton-Holt stock-recruitment relationship
#' srr <- predictModel(
#'   model = rec ~ a * ssb * exp(-b * ssb),
#'   params = FLPar(a = 5.20, b = 1.65e-6)
#' )
#'
#' # Generate autocorrelated recruitment deviances
#' devs <- ar1rlnorm(
#'   rho = 0.4,
#'   years = 2018:2118,
#'   iters = 100,
#'   meanlog = 0,
#'   sdlog = 0.5
#' )
#'
#' # Define an ICES-style risk statistic
#' statistic <- list(
#'   FP05 = list(
#'     ~ yearMeans((SB / SBlim) < 1),
#'     name = "P.05",
#'     desc = "ICES P.05"
#'   )
#' )
#'
#' # Find fbar that gives P(SSB < Blim) = 0.05
#' fp05fwd <- bisect(
#'   stock = stock,
#'   sr = srr,
#'   deviances = devs,
#'   metrics = list(SB = ssb),
#'   refpts = FLPar(SBlim = 150000),
#'   statistic = statistic,
#'   years = 2018:2118,
#'   pyears = 2069:2118,
#'   tune = list(fbar = c(0.1, 1.0)),
#'   prob = 0.05
#' )
#' }
#'
#' @export

bisect <- function(stock, sr, deviances=rec(stock) %=% 1, metrics, refpts,
  statistic, years, pyears=years, tune, prob, tol=0.01, maxit=15, verbose=TRUE) {

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
  
  pmin <- performance(rmin, metrics=metrics, 
    statistics=statistic, refpts=refpts, probs=NULL, years=pyears)

  obmin <- mean(pmin$data, na.rm=TRUE) - prob

  if(verbose)
    cat(" - prob:", mean(pmin$data, na.rm=TRUE), " - diff: ", obmin, "\n")
  
  # CHECK cmin result
  if(isTRUE(all.equal(obmin, 0, tolerance=tol)))
    return(rmin)
  
  # --- RUN at max

  cmax <- fwdControl(year=years, quant=names(tune)[1], value=unlist(tune)[2])

  # PRINT at top
  if(verbose)
    cat(paste0("[2] ", names(tune), ": ", unlist(tune)[2]))

  rmax <- foo(stock, sr=sr, control=cmax, deviances=deviances)
  
  pmax <- performance(rmax, metrics=metrics, 
    statistics=statistic, refpts=refpts, probs=NULL, years=pyears)
  obmax <- mean(pmax$data, na.rm=TRUE) - prob

  if(verbose)
    cat(" - prob:", mean(pmax$data, na.rm=TRUE), " - diff: ", obmax, "\n")
  
  # CHECK cmax result
  if(isTRUE(all.equal(obmax, 0, tolerance=tol)))
    return(rmax)
 
  # --- CHECK range includes 0
  if((obmin * obmax) > 0) {
    warning("Range of hcr param(s) cannot achieve requested tuning objective probability")
    return(list(min=rmin, max=rmax))
  } 

  # --- LOOP bisecting

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

    pmid <- performance(rmid, metrics=metrics, 
      statistics=statistic, refpts=refpts, probs=NULL, years=pyears)
    obmid <- mean(pmid$data, na.rm=TRUE) - prob

    if(verbose)
      cat(" - prob:", mean(pmid$data, na.rm=TRUE), " - diff: ", obmid, "\n")

    # CHECK and RETURN cmid result
    if(isTRUE(all.equal(obmid, 0, tolerance=tol))) {
      return(rmid)
    }

    # TEST LEFT
    if((obmin * obmid) < 0) {

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
    count <- count + 1
  }

  warning("Solution not found within 'maxit', check 'range', 'maxit' or 'tol'.")

  return(rmid)

} # }}}


