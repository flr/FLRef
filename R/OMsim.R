
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
    fbar = FLQuant(an(control["Feq"]), dimnames=list(year=control["minyear"]:control["maxyear"]))
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
        ((c(colSums(sn[, i] * ep[, i]))/an(control["Feq"]*control["SB0"]))^an(control["Frate"]))*
        exp(rnorm(an(control["its"]), mean=-control["Fsigma"]^2/2, sd=control["Fsigma"]))
      sf[, i+1] <- (se[,i+1]%/%quantMeans(se[seq(fages[1], fages[2]), i+1])) %*% fbar[,i]
      ep[,i+1] <- exp(-(sf[,i+1] * fs[,i+1]) - (sm[,i+1] * ms[,i+1])) * sw[,i+1] * ma[,i+1]
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
  
  # catch.wt
  catch.wt(object) <- (landings.wt(object) * landings.n(object) + 
                         discards.wt(object) * discards.n(object)) / catch.n(object)
  
  # catch
  catch(object) <- quantSums(catch.n(object) * catch.wt(object))
  
  return(object)
}
# }}}


# {{{
# newselex()
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



# {{{
# bioidx.sim()
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

bioidx.sim <- function(object,sel=catch.sel(object),sigma=0.2,q=0.001){
  idx = survey(object,sel=sel,biomass=T)
  
  idx@index.q[] =  rlnorm(dims(object)$iter*dims(object)$year,log(q),sigma)
  idx@index.var[] = sigma
  idx@index = idx@index.q*idx@index
  return(idx)
}



#' asem2spm()
#' @param object An *FLBRP*
#' @param quant choose between vb and ssb 
#' @param spcurve if TRUE a data.frame is added
#' @param rel if TRUE ratios are produced for spcurve
#' @return prior means for r and m *FLPar*
#' @export
#' @examples
#' data(ple4)
#' sr <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=mean(spr0y(ple4)))
#' brp = FLBRP(ple4,sr)
#' asem2spm(brp)[1:4]
#' plotpf(brp)
#' plotpf(brp,rel=T)

asem2spm <- function(object,quant=c("vb","ssb"),fmsy=NULL,rel=FALSE,spcurve=FALSE){
  quant = quant[1]
  pbrp = brp=object
  
  if(is.null(fmsy)){
    fmsy = refpts(brp)["msy","harvest"] 
  }
  fbar(brp) = FLQuant(c(0.,fmsy))
  MSY = an(landings(brp)[,2])
  if(quant=="vb") bio = vb(brp)
  if(quant=="ssb") bio = ssb(brp)
 
  Fmsy=an(MSY/bio[,2])
  BmsyK =bio[,2]/bio[,1]
  mi = seq(0.0001,5,0.0001)
  m = mi[abs(mi^(-1/(mi-1))-an(BmsyK))==min(abs(mi^(-1/(mi-1))-an(BmsyK)))][1] 
  r = Fmsy*(m-1)/(1-1/m) 
  pars= FLPar(r=r,m=m,BmsyK =BmsyK,Fmsy=Fmsy,MSY=MSY ,Bmsy=bio[,2],B0=bio[,1])  


   if(spcurve){
    if(quant=="vb")
          out = data.frame(Biomass=c(vb(pbrp)),SP=c(landings(pbrp)),quant=quant)
    if(quant=="ssb") 
          out =data.frame(Biomass=c(ssb(pbrp)),SP=c(landings(pbrp)),quant=quant)
   
     if(rel){
       out$Biomass = out$Biomass/an(pars["B0"])  
       out$SP = out$SP/an(pars["MSY"])  
     }
  }
  
  if(!spcurve)
    return(pars)
  
  if(spcurve)
    return(list(pars=pars,curve=out))
  
} #}}}



#' plotpf()
#' 
#' plots production functions 
#' @param object An *FLBRP*
#' @param quant choose between vb and ssb or both
#' @param fmsy default if Fmsy
#' @param rel if TRUE ratios are produced for spcurve
#' @return ggplot
#' @export
#' @examples
#' data(ple4)
#' sr <- srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=mean(spr0y(ple4)))
#' brp = FLBRP(ple4,sr)
#' asem2spm(brp)[1:4]
#' plotpf(brp)
#' plotpf(brp,rel=T)
plotpf <- function(object,quant=c("vb","ssb"),fmsy=NULL,rel=FALSE){

pf = function(asem,quant,rel=FALSE){
  quant = quant[1]
  bio = asem$curve$Biomass
  r= c(asem$pars["r"])
  m= c(asem$pars["m"])
  k= ifelse(rel,1,c(asem$pars["B0"]))
  SP = bio*r/(m-1)*(1-(bio/k)^(m-1))
  if(rel) SP= SP/max(SP,na.rm = T)
  return(data.frame(Biomass=bio,SP=SP,quant=quant,method="spm"))
}

asem = asem2spm(object,quant=quant[1],fmsy=fmsy,spcurve = TRUE,rel=rel)

y = an(asem$pars["MSY"])
x = an(asem$pars["Bmsy"])

if(rel){
 x = an(asem$pars["BmsyK"])
 y =1
}
xy = data.frame(x=x,y=y,method="spm",quant=quant[1])
df = rbind(data.frame(asem$curve,method="asem"),pf(asem,quant[1],rel=rel))

if(length(quant)>1){
  asem = asem2spm(object,quant=quant[2],fmsy=fmsy,spcurve = TRUE,rel=rel)
  y = an(asem$pars["MSY"])
  x = an(asem$pars["Bmsy"])
  if(rel){
    y =1
    x = an(asem$pars["BmsyK"])
  }
  xy = rbind(xy,data.frame(x=x,y=y,method="spm",quant=quant[2]))
  df = rbind(df,data.frame(asem$curve,method="asem"),pf(asem,quant[2],rel=rel)) 
  
}

p = ggplot(df,aes(Biomass,SP,linetype=method,col=quant))+
    geom_line()+
    geom_segment(data= xy,aes(x =x,y=0,xend = x, yend = y),linetype=2)+
    theme_bw()+theme(legend.title = element_blank())+
    ylab("Surplus Production")+
    xlab(ifelse(!rel,"Biomas",expression(B/B[0])))+
    scale_x_continuous(expand = expansion(mult = c(0, .05)),labels=human_numbers, limits=c(0, NA))+
   scale_y_continuous(expand = expansion(mult = c(0, .1)))
   return(p)
    
}


#' updsr()
#' 
#' updates sr in brp after changing biology 
#' @param object An *FLBRP*
#' @param s assumed steepness s
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

updsr <- function(object,s=0.7){
  sr=SRModelName(model(object))
  par=FLPar(s=s,v = refpts(object)["virgin","ssb"])
  params(object)=FLCore::ab(par[c("s","v")],sr,spr0=spr0(eql))[c("a","b")]
  return(object)
}  
#}}}



#' fudc()
#' 
#' generates an up-down-constant F-pattern 
#' @param object An *FLStock*
#' @param fref reference denominator for fbar 
#' @param fhi factor for high F as fhi = fbar/fref
#' @param flo factor for low F as flo = fbar/fref
#' @param sigmaF variation on fbar
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

fudc = function(object,fref=0.2,fhi=2.5,flo=0.8,sigmaF=0.2){
  fref = c(fref)
  f0 = median(c(fbar(object)[,1]))
  f=(fbar(object))
  x = an(dimnames(object)$year)
  steps = length(x)
  f[] = f0
  y0 = x[1]
  y1 = x[floor(steps/2)]
  y2 = x[ceiling(3*steps/4)]
  f[,x > y0 & x <= y1] = ((fhi*fref-f0)/(y1-y0))*(x[x > y0 & x <= y1] - y0) +f0
  f[,x > y1 & x <= y2] = (-(fhi*fref-flo*fref)/(y2-y1))*(x[x > y1 & x <= y2]-y1) +fhi*fref
  f[,x > y2] = fref*flo
  flq=f*rlnorm(f,0,sigmaF)
  units(flq) ="f"
  return(flq[,-1])
}

