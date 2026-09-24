
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
  return(flq)
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
