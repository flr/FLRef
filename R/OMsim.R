
# rffwd {{{

#' rffwd() Project forward an FLStock with evolutionary Fbar
#'
#' @param object An *FLStock*
#' @param sr A stock-recruit relationship, *FLSR* or *predictModel*.
#' @param fbar Yearly target for average fishing mortality, *FLQuant*.
#' @param control Yearly target for average fishing mortality, *FLPar*.
#' @param deviances Deviances for the strock-recruit relationsip, *FLQuant*.
#'
#' @return The projected *FLStock* object.
#' @export
#' @examples 
#' data(ple4)
#' sr <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=mean(spr0y(ple4)))
#' brp = computeFbrp(ple4,sr,proxy="msy") 
#' fbar(brp) = FLQuant(rep(0.01,70))
#' stk = as(brp,"FLStock")
#' units(stk) = standardUnits(stk)
#' its = 100
#' stk <- FLStockR(propagate(stk, its))
#' stk@refpts= Fbrp(brp)
#' b0=an(Fbrp(brp)["B0"])
#' control = FLPar(Feq=0.15,Frate=0.1,Fsigma=0.15,SB0=b0,minyear=2,maxyear=70,its=its)
#' run <- rffwd(stk, sr=sr,control=control,deviances=ar1rlnorm(0.3, 1:70, its, 0, 0.6))
#' plotAdvice(run)

rffwd <- function(object, sr, fbar=control, control=fbar, deviances="missing") {
  
  # DIMS
  dm <- dim(object)
  
  # EXTRACT slots
  sn <- stock.n(object)
  sm <- m(object)
  sf <- harvest(object)
  se <- catch.sel(object)
  
  # DEVIANCES
  if(missing(deviances)) {
    deviances <- rec(object) %=% 1
  }
  
  # HANDLE fwdControl
  if(is(fbar, "fwdControl")) {
    # TODO CHECK single target per year & no max/min
    # TODO CHECK target is fbar/f
    fbar <- sf[1, fbar$year] %=% fbar$value
  }
  
  if(is(fbar,"FLPar")){
    control = fbar
    fbar = FLQuant(an(control["Feq"]),
      dimnames=list(year=control["minyear"]:control["maxyear"]))
    fbar= propagate(fbar,control["its"])
  }
  
  # SET years
  yrs <- match(dimnames(fbar)$year, dimnames(object)$year)
  
  # COMPUTE harvest
  fages <- range(object, c("minfbar", "maxfbar"))
  sf[, yrs] <- (se[, yrs] %/%
    quantMeans(se[seq(fages[1], fages[2]), yrs])) %*% fbar
  
  # COMPUTE TEP
  sw <- stock.wt(object)
  ma <- mat(object)
  ms <- m.spwn(object)
  fs <- harvest.spwn(object)
  ep <- exp(-(sf * fs) - (sm * ms)) * sw * ma
  
  # LOOP over years
  for (i in yrs - 1) {
    # rec * deviances
    sn[1, i + 1] <- eval(sr@model[[3]],   
      c(as(sr@params, 'list'), list(ssb=c(colSums(sn[, i] * ep[, i]))))) *
      c(deviances[, i + 1])
    # n
    sn[-1, i + 1] <- sn[-dm[1], i] * exp(-sf[-dm[1], i] - sm[-dm[1], i])
    # pg
    sn[dm[1], i + 1] <- sn[dm[1], i + 1] +
      sn[dm[1], i] * exp(-sf[dm[1], i] - sm[dm[1], i])
    
    # Endogenous F 
    if(is(control,"FLPar")){
      fbar[,i] =  quantMeans(sf[seq(fages[1], fages[2]), i]) *
        ((c(colSums(sn[, i] * ep[, i]))/an(control["Feq"] *
        control["SB0"]))^an(control["Frate"])) *
        exp(rnorm(an(control["its"]), mean=-control["Fsigma"]^2/2,
          sd=control["Fsigma"]))
      sf[, i+1] <- (se[,i+1] %/% quantMeans(se[seq(fages[1], fages[2]), i+1])) %*%
        fbar[,i]
      ep[,i+1] <- exp(-(sf[,i+1] * fs[,i+1]) - (sm[,i+1] * ms[,i+1])) * sw[,i+1] *
        ma[,i+1]
    }   
  }
  
  # UPDATE stock.n & harvest
  stock.n(object) <- sn
  harvest(object) <- sf
  
  # UPDATE stock, ...
  stock(object) <- computeStock(object)
  
  # catch.n
  catch.n(object)[,-1] <- (sn * sf / (sm + sf) * (1 - exp(-sf - sm)))[,-1]
  
  # landings.n & discards.n
  landings.n(object)[is.na(landings.n(object))] <- 0
  discards.n(object)[is.na(discards.n(object))] <- 0
  
  landings.n(object) <- catch.n(object) * (landings.n(object) / 
    (discards.n(object) + landings.n(object)))
  
  discards.n(object) <- catch.n(object) - landings.n(object)
  
  # landings & discards
  landings(object) <- computeLandings(object)
  discards(object) <- computeDiscards(object)
  
  # catch.wt
  catch.wt(object) <- (landings.wt(object) * landings.n(object) + 
    discards.wt(object) * discards.n(object)) / catch.n(object)
  
  # catch
  catch(object) <- quantSums(catch.n(object) * catch.wt(object))
  
  return(object)
}
# }}}

# newselex() {{{
#
#' generates flexible 5-paramater selex curves 
#'
#' @param object FLQuant from catch.sel() or sel.pattern()
#' @param selexpars Selectivity Parameters selexpars S50, S95, Smax, Dcv, Dmin 
#' \itemize{
#'   \item S50:  age at 50% selectivity 
#'   \item S95:  age at 50% selectivity
#'   \item Smax: age at peak of selectivity before descending limb 
#'   \item Dcv: CV demeterming the steepness of the descending half-normal slope 
#'   \item Dmin: determines the minimum retention of oldest fishes
#' }    
#' @return FLquant with selectivity pattern
#' @export 
#' @examples 
#' data(ple4)
#' sel = newselex(catch.sel(ple4),FLPar(S50=2,S95=3,Smax=4.5,Dcv=0.6,Dmin=0.3))
#' ggplot(sel)+geom_line(aes(age,data))+ylab("Selectivity")+xlab("Age")
#' # Simulate
#' harvest(ple4)[] = sel
#' sr <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=mean(spr0y(ple4)))
#' brp = computeFbrp(ple4,sr,proxy="msy") 
#' fbar(brp) = FLQuant(rep(0.01,70))
#' stk = as(brp,"FLStock")
#' units(stk) = standardUnits(stk)
#' its = 100
#' stk <- FLStockR(propagate(stk, its))
#' stk@refpts= Fbrp(brp)
#' b0=an(Fbrp(brp)["B0"])
#' control = FLPar(Feq=0.15,Frate=0.1,Fsigma=0.15,SB0=b0,minyear=2,maxyear=70,its=its)
#' run <- rffwd(stk, sr=sr,control=control,deviances=ar1rlnorm(0.3, 1:70, its, 0, 0.6))
#' plotAdvice(run)

newselex<- function(object,selexpars){
  age= dims(object)$min:dims(object)$max
  pars = selexpars
  if(length(selexpars)<3) pars= rbind(sp,FLPar(Smax=max(age)+1,Dcv=0.1,Dmin=0))
  S50 = pars[[1]]
  S95 = pars[[2]]
  Smax =pars[[3]]
  Dcv =pars[[4]]
  Dmin =pars[[5]]
  psel_a = 1/(1+exp(-log(19)*(age-S50)/(S95-S50)))
  psel_b = dnorm(age,Smax,Dcv*Smax)/max(dnorm(age,Smax,Dcv*Smax))
  psel_c = 1+(Dmin-1)*(psel_b-1)/-1
  psel = ifelse(age>=max(Smax),psel_c,psel_a)
  psel = psel/max(psel)
  res = object
  res[] = psel
  return(res)
}
# }}}

# bioidx.sim() {{{
#
#' generates FLIndexBiomass with random observation error from an FLStock
#'
#' @param object FLStock
#' @param sel FLQuant with selectivity.pattern 
#' @param sigma observation error for log(index) 
#' @param q catchability coefficient for scaling
#' @return FLIndexBiomass 
#' @export 
#' @examples 
#' data(ple4)
#' sel = newselex(catch.sel(ple4),FLPar(S50=1.5,S95=2.1,Smax=4.5,Dcv=1,Dmin=0.1))
#' ggplot(sel)+geom_line(aes(age,data))+ylab("Selectivity")+xlab("Age")
#' object = propagate(ple4,10)
#' sel = newselex(catch.sel(object),FLPar(S50=2.5,S95=3.2,Smax=3.5,Dcv=0.6,Dmin=0.2))
#' idx = bioidx.sim(object,sel=sel,q=0.0001)
#' # Checks
#' ggplot(idx@sel.pattern)+geom_line(aes(age,data))+ylab("Selectivity")+xlab("Age")
#' ggplot(idx@index)+geom_line(aes(year,data,col=ac(iter)))+theme(legend.position = "none")+ylab("Index")

bioidx.sim <- function(object,sel=catch.sel(object),sigma=0.2,q=0.001,rho=0){
  
  sel = sel%/%apply(sel,2:6,max)
  if(dims(sel)$unit<dims(object)$unit){
    sel <- expand(sel, unit = c("F", "M"))
  }
  selinp =  catch.sel(object)
  years = ac(dimnames(object)$year)
  selinp[] = sel
  idx = survey(object,sel=selinp,biomass=TRUE)
  qdevs =  FLQuant(rlnormar1(rho=rho,years=years,n= dims(object)$iter,meanlog=0,sdlog=sigma),quant="age")
    
  idx@index.q = q*qdevs

  idx@index.var[] = sigma
  idx@index = idx@index.q*idx@index
  
  units(idx)[] = "NA"
  return(idx)
}
# }}}

# idx.sim {{{

#' generates FLIndex with lognormal annual and multinomial age composition observation error 
#' @param object FLStock
#' @param sel FLQuant with selectivity.pattern 
#' @param ess effective sample size for age composition sample
#' @param sigma annual observation error for log(q)
#' @param ages define age range 
#' @param years define year range
#' @param q catchability coefficient for scaling
#' @return FLIndex
#' @export 
#' @examples 
#' data(ple4)
#' sel = newselex(catch.sel(ple4),FLPar(S50=1.5,S95=2.1,Smax=4.5,Dcv=1,Dmin=0.1))
#' ggplot(sel)+geom_line(aes(age,data))+ylab("Selectivity")+xlab("Age")
#' object = propagate(ple4,10)
#' idx = idx.sim(object,sel=sel,ess=200,sigma=0.2,q=0.01,years=1994:2017)
#' # Checks
#' ggplot(idx@sel.pattern)+geom_line(aes(age,data))+ylab("Selectivity")+xlab("Age")
#' ggplot(idx@index)+geom_line(aes(year,data,col=ac(iter)))+facet_wrap(~age,scales="free_y")+
#' theme(legend.position = "none")+ylab("Index")

idx.sim <- function(object,sel=catch.sel(object),ages=NULL,years=NULL,ess=200,sigma=0.2,q=0.01){
  if(is.null(ages)){
    ages = an(dimnames(object)$age)
  }
  if(is.null(years)){
    years = an(dimnames(object)$year)
  }
  if(dims(sel)$unit<dims(object)$unit){
    sel <- expand(sel, unit = c("F", "M"))
  }
  selinp =  catch.sel(object)
  selinp[] = sel
  sel = trim( selinp,age=ages,year=years)
  object = trim(object,age=ages,year=years)
  idx = index = trim(survey(object,ages=ac(ages),sel=sel,biomass=F))
  devs = rlnorm(dims(object)$iter*dims(object)$year,log(q),sigma)
  for(i in seq(ages)){
    idx@index.q[i,] = devs 
  }
  res = idx@index
  # Sample age-comp from multinomial
  for(i in seq(dims(object)$iter)){
    for(y in seq(dims(object)$year)){
      prob =  c(res[,y,,,,i]%/%apply(res[,y,,,,i],2,sum))
      res[, y, , , , i] <- apply(rmultinom(ess, 1, prob = prob),1,sum)
    }
  }
  
  idx@index.var[] = sigma
  fac = apply(index@index,2:6,sum)/apply(res,2:6,sum)
  res = res%*%fac
  units(idx@index.q) = "1"
  idx@index = idx@index.q*res
  return(idx)
}
# }}}

# pgquant {{{
#
#' sets plus group on FLQuant
#' @param object FLQuant
#' @param pg 
#' @return FLQuant
#' @export 

pgquant <- function(object,pg){
  ages = an(dimnames(object)$age)
  age = ages[ages<=pg]
  plus = ac(ages[ages>=pg])
  res = trim(object,age=age)  
  res[ac(pg),] = quantSums(object[plus,]) 
  return(res)
}
# }}}



# ca.sim() {{{

# Adds an optional `u`
# argument for CRN-safe sampling; existing call sites that don't pass `u`
# are unaffected (falls back to the original stats::rmultinom() behaviour).

#' Generate catch.n with lognormal annual and multinomial age composition
#' observation error
#'
#' @param object `FLQuant` of true numbers/catch-at-age.
#' @param ess Effective sample size for the age-composition sample.
#' @param u Optional. A `[ess, year, iter]` array of pre-generated
#'   `Uniform(0,1)` deviates (e.g. from [rUnif_lfd()]) for CRN-safe
#'   sampling via [rmultinom_crn()]. If `NULL` (default), sampling falls
#'   back to `stats::rmultinom()` against the live RNG, i.e. the original
#'   `ca.sim()` behaviour is unchanged.
#' @param what One of `"catch"`, `"landings"`, `"discards"` (currently
#'   unused inside the function body, kept for interface compatibility
#'   with the original `ca.sim()`).
#' @param rescale Logical. If `TRUE`, rescale the sampled age composition
#'   back to the abundance scale of `object` (sampled proportions x true
#'   total), rather than returning raw counts on the `ess` scale.
#'
#' @return `FLQuant`, same dimensions as `object`, with sampled
#'   catch-at-age (or rescaled proportions x true total, if
#'   `rescale=TRUE`).
#'
#' @examples
#' \dontrun{
#' data(ple4)
#' object <- propagate(catch.n(ple4), 10)
#'
#' ## original behaviour, unchanged
#' ca1 <- ca.sim(object, ess = 200)
#'
#' ## CRN-safe: same u -> same sample, even if called twice
#' u <- array(runif(200 * dim(object)[2] * 10),
#'            dim = c(200, dim(object)[2], 10))
#' ca2 <- ca.sim(object, ess = 200, u = u)
#' ca3 <- ca.sim(object, ess = 200, u = u)
#' stopifnot(identical(ca2, ca3))
#' }
#'
#' @export
ca.sim <- function(object, ess = 200, u = NULL,
                   what = c("catch", "landings", "discards")[1],
                   rescale = FALSE) {
  
  res <- object
  ref <- res
  
  for (i in seq(dims(object)$iter)) {
    for (y in seq(dims(object)$year)) {
      
      prob <- c(res[, y, , , , i] %/% apply(res[, y, , , , i], 2, sum))
      
      if (is.null(u)) {
        res[, y, , , , i] <- apply(rmultinom(ess, 1, prob = prob), 1, sum)
      } else {
        res[, y, , , , i] <- rmultinom_crn(u[, y, i], prob)
      }
    }
  }
  
  if (rescale) {
    res <- quantSums(ref) %*% (res %/% quantSums(res))
  }
  
  res
}

# }}}

# iALK() {{{
#
#' inverse ALK function with lmin added to FLCore::invALK 
#' @param params growth parameter, default FLPar(linf,k,t0)
#' @param model growth model, only option currently vonbert
#' @param age age vector
#' @param cv of length-at-age
#' @param lmax maximum upper length specified lmax*linf
#' @param max maximum size value
#' @param lmin minimum length
#' @param reflen evokes fixed sd for L_a at sd = cv*reflen
#' @param bin length bin size, dafault 1
#' @param timing t0 assumed 1st January, default seq(0,11/12,1/12), but can be single event 0.5
#' @param unit default is "cm"
#' @return FLPar age-length matrix
#' @export 

iALK <- function(params, model=vonbert, age, cv=0.1,lmin=5, lmax=1.2, bin=1,
                   max=ceiling(linf * lmax), reflen=NULL) {
  
  linf <- c(params['linf'])
  
  # FOR each age
  bins <- seq(lmin, max, bin)
  
  # METHOD
  if(isS4(model))
    len <- do.call(model, list(age=age,params=params))
  else {
    lparams <- as(FLPar(params), "list")
    len <- do.call(model, c(list(age=age),
                            lparams[names(lparams) %in% names(formals(model))]))
  }
  
  if(is.null(reflen)) {
    sd <- abs(len * cv)
  } else {
    sd <- reflen * cv
  }
  
  probs <- Map(function(x, y) {
    p <- c(pnorm(1, x, y),
           dnorm(bins[-c(1, length(bins))], x, y),
           pnorm(bins[length(bins)], x, y, lower.tail=FALSE))
    return(p / sum(p))
  }, x=len, y=sd)
  
  res <- do.call(rbind, probs)
  
  alk <- FLPar(array(res, dim=c(length(age), length(bins), 1)),
               dimnames=list(age=age, len=bins, iter=1), units="")
  
  return(alk)
} 
# }}}

# ALK() {{{
#
#' ALK function
#' @param N_a numbers at age sample for single event
#' @param iALK from iALK() outout
#' @return FLPar of ALK
#' @export 
ALK <- function(N_a,iALK){
  alk = iALK
  alk[] = N_a
  alk = alk*iALK
  alksum =apply(alk,2,sum)
  for(i in 1:dim(alk)[1]){
    alk[i,]=alk[i,]/an(alksum)
  }   
  return(alk)
}
# }}}

# alk.sample() {{{
#
#' generates annual ALK sample with length stratified sampling
#' @param lfds length frequency *FLQuant*
#' @param alks annual ALK proportions at age output form ALKs() *FLPars*
#' @param nbin number of samples per length bin
#' @param n.sample sample size of lfd 
#' @return FLPars of sampled ALK
#' @export 

alk.sample <- function(lfds,alks,nbin = 20,n.sample=1){
  res = alks
  if(n.sample>1){
    lfdn = lfds%/%apply(lfds,2:6,sum)*n.sample
  } else {
    lfdn = lfds
  }
  for(i in seq(dims(lfds)$iter)){
    for(y in seq(dims(lfds)$year)){
      nL = pmin(lfdn[,y],nbin )  # check len samples
      for(l in seq(dim(alks[[1]])[2])){
        if(nL[l]>1){  
          res[[y]][,l][] = apply(rmultinom(nL[l], 1, prob = c(alks[[y]][,l])),1,sum)
        } else {
          res[[y]][,l][] = 0}
        
      }
    }
  }
  res
}  
# }}}

# ALKs() {{{
#
#' annual ALK function
#' @param object FLQuant with numbers at age
#' @param iALK from iALK() outout
#' @return FLPars of ALK
#' @export 
ALKs <- function(object,iALK){
  it = dim(object)[6]
  nyr= dim(object)[2]
  year = (dimnames(object)$year)
  alks = FLPars(lapply(year,function(x){
  alk = propagate(iALK,it)
  alk[] = object[,ac(x)]
  for(i in seq(it)){
  iter(alk,i) = iter(alk,i) *iALK
  iter(alk,i) =  sweep(iter(alk,i),2,apply(iter(alk,i),2,sum),"/")
  }
  
  return(alk)
  }))
  names(alks) = ac(year)
 return(alks)
}
# }}}

# applyALK() {{{
#
#' applyALK function to length to age
#' @param lfd *FLQuant* with numbers at length
#' @param alks *FLPars* annual ALKs
#' @return FLQuant for numbers at age
#' @export 

applyALK <- function(lfds,alks){
  yr = dimnames(lfds)$year
  year = an(yr)
  if(class(alks)=="FLPar"){
   alks = FLPars(lapply(yr,function(x){
      alks
    }))
   names(alks) = yr
  }
  age = ac(dimnames(alks[[1]])$age)
  
  if(any(dimnames(lfds)$year !=  names(alks)))
     stop("ALKs list must have the same years as length data")
  its = dims(lfds)$iter
  
  res <- FLQuant(NA, units="1000",
                 dimnames=list(age=age, year=year,iter=1:its))
  
  for(i in seq(its)){
  for(y in 1:length(year)){
    mf = (model.frame(iter(alks[[y]],i))[,seq(dim(alks[[y]])[1])])
    mf = mf/pmax(apply(mf,1,sum),0.01)
    iter(res,i)[,y] = apply(t(mf[,seq(dim(alks[[y]])[1])])%*%c(iter(lfds[,y],i)),1,sum)
  }
  }
  return(res)
}

#{{{
# condition
#' Fold length selectivity into an inverse ALK
#'
#' Shared helper: re-weights an inverse age-length key (age x len, rows
#' sum to 1) by a length-selectivity vector on the same length grid, and
#' renormalises each age row back to sum to 1. Ages with zero selected
#' probability get an all-zero row rather than a divide-by-zero. Used by
#' both `build_gear()` and the gear-aware `len.sim()` so the two stay
#' consistent rather than each carrying its own copy of this logic.
#'
#' @param ialk_mat Numeric matrix, age (rows) x len (cols), rows sum to 1.
#' @param sel_len_vec Numeric vector of length-selectivity, named by
#'   length. Reordered internally to match `colnames(ialk_mat)`.
#'
#' @return Numeric matrix, same dimensions as `ialk_mat`, rows sum to 1
#'   (or all-zero where no length was selected for that age).
#'
#' @export
condition_alk <- function(ialk_mat, sel_len_vec) {
  
  sel_len_vec <- sel_len_vec[colnames(ialk_mat)]
  
  weighted <- sweep(ialk_mat, 2, sel_len_vec, "*")
  row_tot <- rowSums(weighted)
  cond <- weighted
  valid <- row_tot > 0
  cond[valid, ] <- weighted[valid, , drop = FALSE] / row_tot[valid]
  cond[!valid, ] <- 0
  cond
}
#}}}

#' Generate survey (pulse) or continuous length-frequency samples
#'
#' Samples a length-frequency distribution from `N_a` by projecting it
#' through the biological inverse age-length key at one or more
#' within-year timing events (`timing`), optionally re-weighted by a
#' gear's length selectivity.
#'
#' Unlike [lfd.sim()]'s two-stage fishery-dependent design (a noisy
#' catch-at-age sample feeding a length draw), `len.sim()` samples length
#' directly from `N_a` at each timing event and sums across events -- the
#' single-stage design appropriate for a survey or continuous sampling
#' programme, where `N_a` is already the relevant numbers/catch-at-age
#' for that observation process rather than something to be re-sampled
#' for age first.
#'
#' Because growth (`t0 + timing[t]`) shifts within a year, the ALK -- and
#' therefore any gear-selectivity weighting applied to it via
#' [condition_alk()] -- is rebuilt fresh at each timing event, rather
#' than reusing a single pre-built `gear$condALK` as `lfd.sim()` does.
#'
#' @param N_a Numbers (or catch)-at-age `FLQuant`.
#' @param params Growth parameters, `FLPar(linf, k, t0)`.
#' @param model Growth model. Default `vonbert`.
#' @param ess Effective sample size, split evenly across `timing` events
#'   (`round(ess / length(timing))` per event, as in the original).
#'   Default `250`. If using `gear`, consider passing `gear$ess_len`
#'   explicitly.
#' @param timing Within-year timing of sampling events (as fractions of a
#'   year added to `t0`). Default `seq(0, 11/12, 1/12)` (roughly monthly,
#'   continuous sampling); pass a single value (e.g. `0.5`) for a single
#'   survey pulse.
#' @param unit Length unit label. Default `"cm"`.
#' @param scale Logical. If `TRUE` (default), rescale the output so its
#'   per-year/iter total matches `N_a`'s true total.
#' @param reflen,bin,cv,lmin,lmax Passed to `iALK()`, as in the original.
#' @param gear Optional gear object as returned by [build_gear()]. If
#'   supplied, its `sel_len` is folded into the ALK at every timing event
#'   via [condition_alk()]. If `NULL` (default), no selectivity weighting
#'   is applied -- reproduces the original `len.sim()` behaviour exactly.
#' @param u Optional list of length `length(timing)`, each element a
#'   `[ess_t, year, iter]` array of pre-generated `Uniform(0,1)` deviates
#'   (`ess_t = round(ess / length(timing))`) for CRN-safe sampling via
#'   [rmultinom_crn()] -- see [rUnif_len()]. If `NULL` (default), falls
#'   back to `stats::rmultinom()` against the live RNG, i.e. the
#'   original behaviour.
#'
#' @return `FLQuant` with a `len` dimension, the sampled length-frequency
#'   by year and iter.
#'
#' @examples
#' \dontrun{
#' ## original behaviour, unchanged: no gear, no CRN
#' lfd_survey <- len.sim(stock.n(om)[, "2024"], params = lhpar, ess = 300,
#'                        timing = 0.5)
#'
#' ## gear-aware, CRN-safe
#' u_len <- rUnif_len(trawl, timing = seq(0, 11/12, 1/12),
#'                     years = an(dimnames(N_a)$year), nsim = dims(N_a)$iter,
#'                     seed = 456)
#' lfd_survey_g <- len.sim(N_a, params = lhpar, ess = trawl$ess_len,
#'                          timing = seq(0, 11/12, 1/12),
#'                          gear = trawl, u = u_len)
#' }
#'
#' @export
len.sim <- function(N_a, params, model = vonbert, ess = 250,
                    timing = seq(0, 11/12, 1/12), unit = "cm", scale = TRUE,
                    reflen = NULL, bin = 1, cv = 0.1, lmin = 5, lmax = 1.2,
                    gear = NULL, u = NULL) {
  
  gp <- c(params)
  age <- an(dimnames(N_a)$age)
  its <- dims(N_a)$iter
  yrs <- dimnames(N_a)$year
  ess_t <- round(ess / length(timing))
  
  out <- NULL   # built on the first timing event, once we know the len grid
  
  for (t in seq_along(timing)) {
    
    ialk <- iALK(
      params = c(linf = gp[["linf"]], k = gp[["k"]], t0 = gp[["t0"]] + timing[t]),
      model = model, age = age, lmax = lmax, reflen = reflen, bin = bin, lmin = lmin
    )
    ialk_mat <- c(ialk)
    dim(ialk_mat) <- dim(ialk)[1:2]
    dimnames(ialk_mat) <- list(age = dimnames(ialk)$age, len = dimnames(ialk)$len)
    
    if (is.null(gear)) {
      cond <- ialk_mat                      # no selectivity: original behaviour
    } else {
      sel_len_vec <- as.numeric(gear$sel_len)
      names(sel_len_vec) <- dimnames(gear$sel_len)$len
      cond <- condition_alk(ialk_mat, sel_len_vec)
    }
    
    if (is.null(out)) {
      out <- FLQuant(
        0,
        dimnames = list(
          len = colnames(cond), year = yrs, unit = "unique",
          season = "all", area = "unique", iter = seq_len(its)
        )
      )
    }
    
    for (y in seq_along(yrs)) {
      for (i in seq_len(its)) {
        
        n_a <- as.numeric(N_a[, y, , , , i])
        if (sum(n_a) == 0) next
        
        p_len <- as.numeric((n_a / sum(n_a)) %*% cond)
        
        if (is.null(u)) {
          len_n <- apply(rmultinom(ess_t, 1, prob = p_len), 1, sum)
        } else {
          len_n <- rmultinom_crn(u[[t]][, y, i], p_len)
        }
        
        out[, y, , , , i] <- out[, y, , , , i] + len_n
      }
    }
  }
  
  units(out) <- unit
  
  if (scale) {
    fac <- apply(N_a, 2:6, sum) / apply(out, 2:6, sum)
    fac[!is.finite(fac)] <- 0
    out <- out %*% fac
  }
  
  out
}

#' Pre-generate CRN uniform-deviate streams for len.sim()
#'
#' Generates one `Uniform(0,1)` array per within-year timing event, each
#' shaped `[ess_t, year, iter]`, for use with [len.sim()]'s `u` argument
#' via [rmultinom_crn()]. Analogous to [rUnif_lfd()], but shaped for
#' `len.sim()`'s timing-event loop rather than the two-stage age/length
#' design.
#'
#' @param gear A gear object as returned by [build_gear()]; used for
#'   `ess_len` (per [len.sim()]'s `ess` argument, split across timing
#'   events the same way `len.sim()` does internally).
#' @param timing Within-year timing vector, as passed to [len.sim()].
#' @param years Integer vector of years to generate deviates for.
#' @param nsim Integer. Number of OM iterations.
#' @param seed Optional integer seed. If `NULL` (default), the current
#'   RNG state is used as-is.
#'
#' @return A list of length `length(timing)`, each element an array of
#'   dimension `[round(gear$ess_len / length(timing)), length(years), nsim]`.
#'
#' @examples
#' \dontrun{
#' u_len <- rUnif_len(trawl, timing = seq(0, 11/12, 1/12),
#'                     years = 2025:2054, nsim = 100, seed = 456)
#' length(u_len)          # one array per timing event
#' dim(u_len[[1]])        # ess_t x years x nsim
#' }
#'
#' @export
rUnif_len <- function(gear, timing, years, nsim, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ess_t <- round(gear$ess_len / length(timing))
  yrs <- as.character(years)
  
  lapply(seq_along(timing), function(t) {
    array(
      stats::runif(ess_t * length(yrs) * nsim),
      dim = c(ess_t, length(yrs), nsim),
      dimnames = list(draw = seq_len(ess_t), year = yrs, iter = seq_len(nsim))
    )
  })
}

# {{{
# lfd.sim
#
# Depends on ca.sim() (ca_sim_crn.R) and rmultinom_crn() (rmultinom_crn.R)
# already being loaded.

#' Two-stage gear-selective length-frequency sample from an operating model
#'
#' Generates a length-frequency sample for one gear in two stages:
#' \enumerate{
#'   \item a low-effective-sample-size catch-at-age sample is drawn from
#'     the gear's true catch-at-age (via [ca.sim()]);
#'   \item that sampled age composition is projected through the gear's
#'     conditional inverse age-length key,
#'     \eqn{P(l \mid a, \text{caught by gear})} (built once by
#'     [build_gear()], so selectivity is already folded in), and a larger
#'     length sample is drawn from the resulting expected length
#'     distribution.
#' }
#' This mimics a sampling design where a small, possibly more costly,
#' age-reading sample (e.g. observer-collected otoliths) informs the age
#' structure applied to a larger, cheaper length sample (e.g. market or
#' landings-site length measurements) -- rather than treating the two as
#' independent draws from the true population.
#'
#' @param object `FLQuant` of the gear's true catch (or numbers)-at-age,
#'   e.g. one element of `catch_n_gear` from the operating model.
#' @param gear A gear object as returned by [build_gear()]; must contain
#'   `condALK` (age x len matrix, rows sum to 1) and, unless overridden
#'   below, `ess_age`/`ess_len`.
#' @param ess_age,ess_len Integer. Effective sample sizes for the two
#'   stages. Default to `gear$ess_age`/`gear$ess_len`.
#' @param u_age Optional `[ess_age, year, iter]` array of pre-generated
#'   `Uniform(0,1)` deviates for the stage-1 draw (see [rUnif_lfd()]). If
#'   `NULL` (default), stage 1 falls back to `ca.sim()`'s own
#'   `stats::rmultinom()` behaviour (no CRN).
#' @param u_len Optional `[ess_len, year, iter]` array of pre-generated
#'   `Uniform(0,1)` deviates for the stage-2 draw. If `NULL` (default),
#'   stage 2 falls back to a plain `stats::rmultinom()` draw (no CRN).
#' @param scale Logical. If `TRUE` (default), rescale the returned length
#'   sample so its per-year/iter total matches the true total in `object`,
#'   rather than returning raw counts on the `ess_len` scale.
#'
#' @return `FLQuant` with a `len` dimension (matching `gear$len_bins`),
#'   the sampled length frequencies by year and iter.
#'
#' @examples
#' \dontrun{
#' ## trawl <- build_gear(...) from build_gear.R examples
#' ## catch_n_gear$Trawl : FLQuant, true catch-at-age for Trawl, from the OM
#'
#' lfd_trawl <- lfd.sim(catch_n_gear$Trawl, trawl)
#'
#' ## sanity check: shape should show the dome, not just decline with length
#' plot(lfd_trawl)
#'
#' ## CRN-safe version, reproducible across calls
#' devs <- rUnif_lfd(om_gears, years = an(dimnames(catch_n_gear$Trawl)$year),
#'                    nsim = dim(catch_n_gear$Trawl)[6], seed = 456)
#' lfd_trawl_a <- lfd.sim(catch_n_gear$Trawl, trawl,
#'                         u_age = devs$age$Trawl, u_len = devs$len$Trawl)
#' lfd_trawl_b <- lfd.sim(catch_n_gear$Trawl, trawl,
#'                         u_age = devs$age$Trawl, u_len = devs$len$Trawl)
#' stopifnot(identical(lfd_trawl_a, lfd_trawl_b))
#' }
#'
#' @export
lfd.sim <- function(object, gear,
                    ess_age = gear$ess_age, ess_len = gear$ess_len,
                    u_age = NULL, u_len = NULL, scale = TRUE) {
  
  age <- an(dimnames(object)$age)
  years <- dimnames(object)$year
  its <- dims(object)$iter
  
  condALK <- gear$condALK
  if (!identical(dimnames(condALK)$age, as.character(age))) {
    stop(
      "'gear$condALK' age dimension does not match 'object' age dimension. ",
      "Check 'object' and 'gear' were built for the same age range."
    )
  }
  
  len_lower <- colnames(condALK)
  out <- FLQuant(
    NA_real_,
    dimnames = list(
      len = len_lower, year = years, unit = "unique",
      season = "all", area = "unique", iter = seq_len(its)
    )
  )
  
  ## --- stage 1: low-ESS catch-at-age sample ------------------------------
  age_n_flq <- ca.sim(object, ess = ess_age, u = u_age)
  
  ## --- stage 2: larger length sample via the conditional ALK -------------
  for (y in seq_along(years)) {
    for (i in seq_len(its)) {
      
      age_n <- as.numeric(age_n_flq[, y, , , , i])
      
      if (sum(age_n) == 0) {
        ## degenerate stage-1 draw (e.g. ~zero catch that year/iter):
        ## fall back on the true age proportions so stage 2 isn't handed
        ## an all-zero mixture
        age_n <- as.numeric(object[, y, , , , i])
      }
      
      age_p <- age_n / sum(age_n)
      p_len <- as.numeric(age_p %*% condALK)
      
      if (is.null(u_len)) {
        len_n <- apply(rmultinom(ess_len, 1, prob = p_len), 1, sum)
      } else {
        len_n <- rmultinom_crn(u_len[, y, i], p_len)
      }
      
      out[, y, , , , i] <- len_n
    }
  }
  
  units(out) <- "cm"
  
  if (scale) {
    fac <- apply(object, 2:6, sum) / apply(out, 2:6, sum)
    fac[!is.finite(fac)] <- 0
    out <- out %*% fac
  }
  
  out
}

# {{{
# rUnif_lfd 
#
#' Pre-generate CRN uniform-deviate streams for a two-stage length OEM
#'
#' Generates one `Uniform(0,1)` array per gear per sampling stage (age,
#' length), each shaped `[draw, year, iter]`, for use with
#' [rmultinom_crn()]. Intended to be called once per MSE run/scenario set
#' and reused everywhere the length OEM draws a sample, so the same
#' underlying randomness is shared across MP variants (CRN).
#'
#' @param om_gears Named list of gear objects as returned by
#'   [build_gear()]; must each contain `ess_age` and `ess_len`.
#' @param years Integer vector of years to generate deviates for
#'   (typically the full projection horizon).
#' @param nsim Integer. Number of OM iterations.
#' @param seed Optional integer seed, set via `set.seed()` before
#'   generation for reproducibility across sessions. If `NULL` (default),
#'   the current RNG state is used as-is.
#'
#' @return A list with two elements, `age` and `len`, each a named list
#'   (one array per gear) of dimension `[ess, length(years), nsim]`.
#'
#' @examples
#' \dontrun{
#' om_gears <- list(Trawl = trawl, Gillnet = gillnet)   # from build_gear()
#' devs <- rUnif_lfd(om_gears, years = 2025:2054, nsim = 100, seed = 456)
#'
#' dim(devs$age$Trawl)     # 50 x 30 x 100  (ess_age x years x nsim)
#' dim(devs$len$Gillnet)   # 300 x 30 x 100 (ess_len x years x nsim)
#' }
#'
#' @export
rUnif_lfd <- function(om_gears, years, nsim, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  gears <- names(om_gears)
  yrs <- as.character(years)
  
  mk <- function(field) {
    stats::setNames(lapply(gears, function(g) {
      n <- om_gears[[g]][[field]]
      array(
        stats::runif(n * length(yrs) * nsim),
        dim = c(n, length(yrs), nsim),
        dimnames = list(draw = seq_len(n), year = yrs, iter = seq_len(nsim))
      )
    }), gears)
  }
  
  list(age = mk("ess_age"), len = mk("ess_len"))
}
#}}}



#' updsr()
#' 
#' updates sr in brp after changing biology 
#' @param object An *FLBRP*
#' @param s assumed steepness s
#' @param v input option new SB0
#' @return FLBRP
#' @export
#' @examples
#' data(ple4)
#' sr <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=mean(spr0y(ple4)))
#' brp = FLBRP(ple4,sr)
#' s = sr@SV[[1]]
#' params(brp)
#' # change
#' m(brp) = Mlorenzen(stock.wt(brp),Mref=0.15)
#' brpupd =updsr(brp,s)
#' params(brp)

updsr <- function(object,s=0.7,v=NULL){
  if(is.null(v))
  v = refpts(object)["virgin","ssb"]
  sr=SRModelName(model(object))
  par=FLPar(s=s,v = v)
  params(object)=FLCore::ab(par[c("s","v")],sr,spr0=spr0(object))[c("a","b")]
  brp(object)
  return(object)
}  
# }}}

# fudc() {{{

#' generates an up-down-constant F-pattern 
#' @param object An *FLStock*
#' @param fref reference denominator for fbar 
#' @param fhi factor for high F as fhi = fbar/fref
#' @param flo factor for low F as flo = fbar/fref
#' @param sigmaF variation on fbar
#' @param breaks relative location of directional change
#' @return FLQuant
#' @export
#' @examples
#' data(ple4)
#' sr <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=mean(spr0y(ple4)))
#' brp = computeFbrp(ple4,sr,proxy="msy")
#' fmsy = Fbrp(brp)["Fmsy"]
#' stki = propagate(ple4,100)
#' fy = fudc(ple4,fhi=2,flo=0.9,fref=fmsy,sigmaF=0)
#' fyi = fudc(stki,fhi=2,flo=0.9,fref=fmsy,sigmaF=0.2)
#' plot(fy,fyi)+ylab("F")
#' #Forcasting
#' om <- FLStockR(ffwd(stki,sr,fbar=fyi))
#' om@refpts = Fbrp(brp)
#' plotAdvice(window(om,start=1960))

fudc = function(object,fref=0.2,fhi=2.5,flo=0.8,sigmaF=0.2,breaks=c(0.5,0.75)){
  fref = c(fref)
  f0 = median(c(fbar(object)[,1]))
  f=(fbar(object))
  x = an(dimnames(object)$year)
  steps = length(x)
  f[] = f0
  y0 = x[1]
  y1 = x[floor(steps*breaks[1])]
  y2 = x[ceiling(steps*breaks[2])]
  f[,x > y0 & x <= y1] = ((fhi*fref-f0)/(y1-y0))*(x[x > y0 & x <= y1] - y0) +f0
  f[,x > y1 & x <= y2] = (-(fhi*fref-flo*fref)/(y2-y1))*(x[x > y1 & x <= y2]-y1) +fhi*fref
  f[,x > y2] = fref*flo
  flq=f*rlnorm(f,0,sigmaF)
  units(flq) ="f"
  return(flq[,-1])
}

# }}}

# schaefer.sim() {{{

#' schaefer.sim()
#' 
#' generates a Schafer surplus production model with process and observation error
#' @param k carrying capacity
#' @param r intrinsic rate of population increase
#' @param q catchability coefficient 
#' @param pe process error 
#' @param oe process error 
#' @param bk initial fraction of b/k
#' @param years time horizon 
#' @param f0 factor for initial year as f0 = f/fmsy
#' @param fhi factor for high F as fhi = f/fmsy
#' @param flo factor for low F as flo = fbar/fmsy
#' @param sigmaF variation on f trajectory
#' @param iters number of iterations
#' @param rel if TRUE metrics B/Bmsy and F/Fmsy are produced
#' @return FLQuants
#' @export
#' @examples
#' stk = schaefer.sim(iters=100,q=0.5) 
#' plotAdvice(stk)
#' plot(FLIndex(index=iter(stk@stock,1))) # index

schaefer.sim <- function(k=10000,r=0.3,q=0.5,pe=0.1,oe=0.2,bk=0.9,
                         years=1980:2022,f0=0.2,fhi=2.2,flo=0.8,sigmaF=0.15,iters=1,
                         blim=0.3,bthr=0.5,rel=FALSE){
  fmsy= fref= r/2
  bmsy = k/2
  f0 = f0*fmsy
  f = propagate(FLQuant(f0,dimnames=list(year= years,age=1),units="f"),iters)
  x = an(dimnames(f)$year)
  steps = length(x)
  y0 = x[1]
  y1 = x[floor(steps/2)]
  y2 = x[ceiling(3*steps/4)]
  f[,x > y0 & x <= y1] = ((fhi*fref-f0)/(y1-y0))*(x[x > y0 & x <= y1] - y0) +f0
  f[,x > y1 & x <= y2] = (-(fhi*fref-flo*fref)/(y2-y1))*(x[x > y1 & x <= y2]-y1) +fhi*fref
  f[,x > y2] = fref*flo
  f=f*rlnorm(f,0,sigmaF)
  units(f) ="f"
  b = propagate(FLQuant(k*bk,dimnames=list(year= c(years,max(years)+1),age=1),units="t"),iters)
  pdevs = rlnorm(b,0,pe)
  b[,1] = b[,1]*pdevs[,1] 
  for(y in 2:(steps+1)){
    b[,y] = (b[,y-1]+r*b[,y-1]*(1-b[,y-1]/k)-b[,y-1]*f[,y-1])*pdevs[,y-1] 
  }
  B = b[,ac(years)]
  C = b[,ac(years)]*f 
  units(C) = "t"
  H = f
  if(rel){
    B = B/bmsy
    H = H/fmsy
  }
  
  df = as.data.frame(B)
  year = unique(df$year)
  N = as.FLQuant(data.frame(age=1,year=df$year,unit="unique",
                            season="all",area="unique",iter=df$iter,data=1))
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
    Fmsy = fmsy,
    Bmsy = bmsy,
    MSY = r*k/4,
    Blim= bmsy*blim,
    Bthr= bmsy*bthr,
    B0 = k,
  )
  
  # index
  obsdev = rlnorm(b[,ac(years)],0,oe)
  index = b[,ac(years)]*q*obsdev
  
  stk@stock = index
  
  if(rel){
    stk@refpts[1:2] =1 
    stk@refpts["Blim"] = blim
    stk@refpts["Bthr"] = bthr
    stk@refpts["B0"] = k/bmsy
  }
  stk@desc = "spm"
  
  
  
  
  return(stk)
}       
# }}}


# updCatch.n {{{

#' computes catch.n for a given harvest 
#' @param object An *FLStock*
#' @return FLQuant
#' @export

updCatch.n <- function(stk) {
  z <- harvest(stk) + m(stk)
  cn <- stock.n(stk) * harvest(stk) / z * (1 - exp(-z))
  cn[z == 0] <- 0
  catch.n(stk) <- cn
  stk
}



# {{{
# SOP corrections 
#
#' scales catch-at-age to total catch with error (optional)
#' @param object FLQuant catch.n, discard.n, landings.n
#' @param stock FLStock
#' @param sigma observation error
#' @param what type c("catch", "landings", "discards")
#' @return FLQuant
sops <- function(object,stock,sigma=0.1,what=c("catch","landings","discards")[1]){
  dmo = dimnames(object) 
  stock = stock[dmo$age,dmo$year]
  if(what=="catch") out = (catch(stock)/quantSums(object*catch.wt(stock)))%*%object
  if(what=="landings") out = (landings(stock)/quantSums(object*landings.wt(stock)))%*%object
  if(what=="landings") out = (discards(stock)/quantSums(object*discards.wt(stock)))%*%object
  
  devs = rlnorm(dims(object)$iter*dims(object)$year,0,sigma)
  for(i in seq(dims(object)$age)){
    out[i,] = out[i,]*devs 
  }
  out = out*devs
  #out[out==0] = NA
  return(out)
}
