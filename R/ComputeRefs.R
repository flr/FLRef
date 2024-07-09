
#{{{
#' computeFbrp()
#
#' Computes biological reference points corresponding to the proxy Fbrp
#'
#' @param stock object of class FLStock 
#' @param sr stock recruitment model of class FLSR
#' @param proxy choice of Fmsy proxies (combinations permitted) 
#' \itemize{
#'   \item "sprx"  spawning potential ratio spr/spr0 with basis x 
#'   \item "bx" SSB as fraction xSSB0
#'   \item "f0.1" 10% slope of yield-per-recruit curve
#'   \item "fe40" Patterns estimator for Fmsy
#'   \item "msy"  maximum surplus production (not defined for segreg)
#'   \item  numeric user value 
#' }
#' @param x basis in percent for sprx and bx, e.g. 40 for spr40
#' @param blim values < 1 are taken as fraction to B0 and blim > 1 as absolute values unless specified otherwise
#' @param type type of blim input, values < 1 are  
#' \itemize{
#'   \item "b0" fraction to B0  
#'   \item "btgt" fraction to Btarget (first occurring in proxy)    
#'   \item "value" absolute value
#' }
#' @param btri Btrigger can specified as absolute value
#' @param bpa Bpa can specified as absolute value
#' @param bthresh Bthresh (GFCM) interchangeable use with Bpa
#' @param fmax maximum Flim = max(Flim,fmax*Fbrp)
#' @param verbose   
#' @return brp object of class FLBRP with computed Fbrp reference points   
#' @export
#' @examples
#' data(ple4)
#' srr = srrTMB(as.FLSR(ple4,model=rickerSV),spr0=spr0y(ple4))
#' brp = computeFbrp(stock=ple4,sr=srr,proxy=c("sprx","f0.1"),blim=0.1,type="b0")
#' ploteq(brp,obs=TRUE,refpts="msy")

computeFbrp <- function(stock,sr='missing',proxy=NULL,x=NULL,blim=0.1,type=c("b0","btgt","value"),btri="missing",bpa="missing",bthresh="missing",verbose=T,fmax=10, ...){
  if(type[1]=="btrg") type="btgt"
  
  if(class(stock)=="FLStockR") stock = as(stock,"FLStock")
      
  
  if(missing(sr)){
    if(verbose)cat(paste0("Computing Per-Recruit Quantities"),"\n")
    sr = FLSR(params=FLPar(1, params='a'),model=formula(rec~a))
    
  }  
    srmod = SRModelName(model(sr))
 
    if(srmod%in%c("bevholtSV","rickerSV")){
        if(is.null(proxy)) proxy=c("bx","f0.1","fe40","msy")
        if(is.null(x)) x=35
    } else {
      if(is.null(proxy)) proxy=c("sprx","f0.1","fe40","msy")
      if(is.null(x)) x=40
    }
    
  
  
 
  type = type[1]
  brp = brp(FLBRP(stock,sr, ...))
  
  if(is.numeric(proxy)){
    Fbrp = FLPar(Fref=proxy)
    if(verbose)cat(paste0("Computing Fref = ",proxy," with Btgt = Bref"),"\n")
  } else {
  
  Fbrp = FLPar(none=NA)
  if(c("sprx")%in%proxy){
    pr = brp(FLBRP(stock)) # per-recruit
    Fbrp = rbind(Fbrp,FLPar(Fspr = refpts(pr+FLPar(Btgt=refpts(pr)["virgin","ssb"]*x*0.01))["Btgt","harvest"]))
    if(verbose)cat(paste0("Computing Fspr",x," with Btgt = Bspr",x),"\n")
  }
  if(c("bx")%in%proxy){
    Fbrp = rbind(Fbrp,FLPar(Fb = refpts(brp+FLPar(Btgt=refpts(brp)["virgin","ssb"]*x*0.01))["Btgt","harvest"]))
    if(verbose)cat(paste0("Computing Fsb",x," with Btgt = Bsb",x),"\n")
  }  
  if(c("f0.1")%in%proxy){
    Fbrp = rbind(Fbrp,FLPar(F0.1=refpts(brp)["f0.1","harvest"]))
    if(verbose)cat(paste0("Computing F0.1 with Btgt = B[F.01]"),"\n")
  }
  
  if(c("fe40")%in%proxy){
    Fbrp = rbind(Fbrp,FLPar(Fe40=Fe40(stock)))
    if(verbose)cat(paste0("Computing Fe40 with Btgt = B[Fe40]"),"\n")
  }
  
  if(c("msy")%in%proxy){
    if(SRModelName(model(sr))%in%c("mean")){
      if(verbose)cat(paste0("Warning: Fmsy = Fmax is not unambiguously defined for ",SRModelName(model(sr))," SR model","\n"))
      Fbrp = rbind(Fbrp,FLPar(Fmsy=refpts(brp)["msy","harvest"]))
    } else if(SRModelName(model(sr))%in%c("segregA1","segreg")){  
      flim = refpts(brp+FLPar(Blim=sr@params[[2]]))["Blim","harvest"]
      Fbrp = rbind(Fbrp,FLPar(Fmsy=min(flim,refpts(brp)["fmax","harvest"])))
      if(verbose)cat(paste0("Warning: Fmsy = min(Flim,Fmax) is not unambiguously defined for ",SRModelName(model(sr))," SR model","\n"))
    } else {
    #add warning for segreg
    
    Fbrp = rbind(Fbrp,FLPar(Fmsy=refpts(brp)["msy","harvest"]))
      if(verbose)cat(paste0("Computing Fmsy with Btgt = Bmsy"),"\n")
    }
  }
  
  ord = c("Fspr","Fb","F0.1","Fe40","Fmsy")[match(proxy,c("sprx","bx","f0.1","fe40","msy"))]
  Fbrp = Fbrp[-1]
  Fbrp = Fbrp[match(rownames(Fbrp),ord)]
  }
  
  B0 = an(refpts(brp)["virgin","ssb"])
  
  brpf =  brp+Fbrp
  if(blim<=1 & type!="value"){
    if(type%in%c("b0")){
      Blim =an(refpts(brpf)["virgin","ssb"])*blim
      if(srmod=="segreg"){
        if(verbose & Blim<sr@params[[2]]) cat("\n",paste0(" Upward adjusted Blim to breakpoint  = ",round(sr@params[[2]],0)," of segreg"),"\n")
        if(verbose & Blim>=sr@params[[2]]) cat("\n",paste0(" Blim = ",blim,"B0"),"\n")
        Blim = max(Blim,sr@params[[2]])
        } else {
      if(verbose)cat("\n",paste0(" Blim = ",blim,"B0"),"\n")
        }
        }
    if(type%in%c("btgt")){
      Blim =an(refpts(brpf)[rownames(Fbrp)[1],"ssb"])*blim
      if(srmod=="segreg"){
        if(verbose & Blim<sr@params[[2]]) cat("\n",paste0(" Upward adjusted Blim to breakpoint  = ",round(sr@params[[2]],0)," of segreg"),"\n")
        if(verbose & Blim>=sr@params[[2]]) cat("\n",paste0("Blim = ",blim," with Btgt corresponding to ", rownames(Fbrp)[1]),"\n")
        Blim = max(Blim,sr@params[[2]])
      } else {
        if(verbose)cat("\n",paste0("Blim = ",blim," with Btgt corresponding to ", rownames(Fbrp)[1]),"\n")
      }

    }
  } else {
    Blim = blim
    if(srmod=="segreg"){
      if(verbose & Blim<sr@params[[2]]) cat("\n",paste0(" Upward adjusted Blim to breakpoint  = ",round(sr@params[[2]],0)," of segreg"),"\n")
      if(verbose & Blim>=sr@params[[2]]) cat("\n",paste0("Blim as input value Blim = ",blim),"\n")
      Blim = max(Blim,sr@params[[2]])
    } else {
      if(verbose)cat("\n",paste0("Blim as input value Blim = ",blim),"\n")
    }
    
    Blim = blim
  }
  

  
  # rename
  rownames(Fbrp)[which(rownames(Fbrp)%in%c("Fspr","Fb"))] = paste0(rownames(Fbrp)[which(rownames(Fbrp)%in%c("Fspr","Fb"))],x) 
  fref  = rownames(Fbrp)
  refs = rbind(FLPar(Blim=Blim),Fbrp)
  
  
  if(!missing(btri)){
    refs= rbind(refs,FLPar(Btri=btri))
  }
  
  if(!missing(bthresh)){
  refs= rbind(refs,FLPar(Bthr=bthresh))
  }
  
  if(!missing(bpa)){
    refs= rbind(refs,FLPar(Bpa=bpa))
  }
  
  refs=rbind(refs,FLPar(B0=0.99*B0))
  
  # do check
  check = brp+refs
  if(refpts(check)["Blim","harvest"]>fmax*refpts(check)[paste(fref[1]),"harvest"]){
    refs["Blim"] = refpts(brp+FLPar(Flim=fmax*refpts(check)[paste(fref[1]),"harvest"]))["Flim","ssb"] 
  }
    
  brp = brp+refs
  refpts(brp)["B0"] = refpts(brp)["virgin"]
  #if(fmax*refpts(check)["Fbrp","harvest"]>max(fbar(brp))){
    fl = min(c(fmax*refpts(check)[paste(fref),"harvest"],refpts(check)["Blim","harvest"]*1.5))
    fl = max(fbar(stock),fl,na.rm=TRUE)
    fbar(brp) = seq(0,fl,fl/100)
  #} 
  
  return(brp)
}
#}}}

# {{{

#' computeFbrps()
#
#' Computes biological reference points corresponding to the proxy Fbrp
#'
#' @param stock object of class FLStock 
#' @param sr stock recruitment model of class FLSR
#' @param proxies choice of Fmsy proxies
#' \itemize{
#'   \item "all"  both sprx and bx 
#'   \item "sprx"  spawning potential ratio spr/spr0 with basis x 
#'   \item "bx" SSB as fraction xSSB0
#' }    
#' @param fmsy if TRUE, Fmsy is computed (not suggest for segreg or geomean sr) 
#' @param f0.1 if TRUE, F0.1 is computed
#' @param fmax maximum Flim = minfmax*Fbrp)
#' @param verbose   
#' @return brp object of class FLBRP with computed Fbrp reference points   
#' @export

computeFbrps <- function(stock,sr="missing",proxy=c("sprx","bx","all"),fmsy=FALSE,f0.1=TRUE,fmax=5,verbose=T, ...){
  
    
    # use geomean sr if sr = NULL (only growth overfishing)
    if(missing(sr)){
      if(verbose)cat(paste0("Computing Per-Recruit Quantities"),"\n")
      sr = FLSR(params=FLPar(1, params='a'),model=formula(rec~a))
      
    }  
    
  proxy=proxy[1] 
 
  brp = brp(FLBRP(stock,sr, ...))
  if(proxy%in%c("sprx","all")){
    pr = brp(FLBRP(stock)) # per-recruit
    Fsprs = FLPar(
      Fspr35 = an(refpts(pr+FLPar(Btgt=refpts(pr)["virgin","ssb"]*35*0.01))["Btgt","harvest"]),
      Fspr40 = an(refpts(pr+FLPar(Btgt=refpts(pr)["virgin","ssb"]*40*0.01))["Btgt","harvest"]),
      Fspr45 = an(refpts(pr+FLPar(Btgt=refpts(pr)["virgin","ssb"]*45*0.01))["Btgt","harvest"]),
      Fspr50 = an(refpts(pr+FLPar(Btgt=refpts(pr)["virgin","ssb"]*50*0.01))["Btgt","harvest"])
    )  
      if(verbose)cat(paste0("Computing Fspr% 35-50 with Btgt = Bspr"),"\n")
  }
  
  
  
  if(proxy%in%c("bx","all")){
    Bx = FLPar(
          B30=refpts(brp)["virgin","ssb"]*30*0.01,
          B35=refpts(brp)["virgin","ssb"]*35*0.01,
          B40=refpts(brp)["virgin","ssb"]*40*0.01,
          B45=refpts(brp)["virgin","ssb"]*45*0.01)
    fbrps = refpts(brp+Bx)[8:11,"harvest"] 
    FBs = FLPar(Fb30= fbrps[1,],Fb35= fbrps[2,],Fb40= fbrps[3,],Fb45= fbrps[4,] )
    
    if(verbose)cat("\n",paste0("Computing Fsb% 30-45 with Btgt = Bsb"),"\n")
  }  
  if(proxy=="all") Fbrps=rbind(Fsprs,FBs) 
  if(proxy=="sprx") Fbrps=Fsprs 
  if(proxy=="bx") Fbrps=FBs 
  if(fmsy) Fbrps = rbind(Fbrps,FLPar(Fmsy=refpts(brp)["msy","harvest"]))
  if(f0.1) Fbrps = rbind(Fbrps,FLPar(F0.1=refpts(brp)["f0.1","harvest"]))
  
  fref = rownames(Fbrps)
  
  Fbrps = rbind(Fbrps, FLPar(B0 = 0.99*an(refpts(brp)["virgin","ssb"])))
  
  
  brp =  brp(brp+Fbrps)
  # Fix
  refpts(brp)["B0"] =   refpts(brp)["virgin"] 
  
  
  fl = min(c(fmax*refpts(brp)[paste(fref),"harvest"]),na.rm=TRUE)
  fl = max(fbar(stock),fl,na.rm=TRUE)
  fbar(brp) = seq(0,fl,fl/100)
  
  
  return(brp)
}
#}}}



#{{{
#' Fbrp()
#
#' Extract Fbrp based reference points from output of computeFbrp
#' @param brp input of class FLBRP from ComputeFbrp    
#' @return FLPar object with computed Fbrp reference points   
#' @export

Fbrp <- function(brp){
  rpt = refpts(brp)
  ref = rownames(rpt)[grep("F",rownames(rpt))][1]
  out = FLPar(Fbrp = rpt[ref,"harvest"],
        Btgt=rpt[ref,"ssb"],
        Blim=rpt["Blim","ssb"],
        Flim=rpt["Blim","harvest"],
        Yeq = rpt[ref,"yield"],
        B0 = rpt["virgin","ssb"], 
        R0= rpt["virgin","rec"]
        )
  rownames(out)[1] = ref
  if("Fp.05"%in%rownames(rpt)){
    out = rbind(out,FLPar(Fp.05=rpt["Fp.05","harvest"]))
  }
  if("Btri"%in%rownames(rpt)){
    out = rbind(out,FLPar(Btri=rpt["Btri","ssb"]))
  }
 
  if("Bpa"%in%rownames(rpt)){
    out = rbind(out,FLPar(Bpa=rpt["Bpa","ssb"],Fpa=rpt["Bpa","harvest"]))
  }
  
  
  if("Bthr"%in%rownames(rpt)){
    out = rbind(out,FLPar(Bthr=rpt["Bthr","ssb"],Fthr=rpt["Bthr","harvest"]))
  }
  
  return(out)
}
#}}}


#{{{
#' Fe40()
#
#' Patterson estimator for Fmsy
#' @param stock input of class FLStock     
#' @param nyears number of years to average      
#' @return value  
#' @export
Fe40 = function(stock,nyears=3){
  fbar.range = range(stock)[c("minfbar")]:range(stock)[c("maxfbar")]
  Mbar = apply(m(stock)[ac(fbar.range),],2:6,mean)
  return(mean(tail(0.4/0.6*Mbar,nyears)))
  }
#}}}


#{{{
#' ABItgt()
#
#' Computes ABI for target F, e.g. ABImsy (Griffith et al. 2023)
#'
#' @param stock object of class FLStock 
#' @param ftgt target F at equilibrium, e.g. Fmsy
#' @param thresh quantile ageref treshold, default 0.9
#' @return *FLQuant* 
#' @export
#' @examples
#' data(ple4)
#' ABImsy = ABItgt(ple4,ftgt=0.22,thresh=0.9)
#' plot(ABImsy)+ylim(0,2)+
#'  geom_hline(yintercept = c(0.8,1),col=c(2,1),linetype=c(2,1))+ylab(expression(ABI[MSY]))

ABItgt <- function(stock,ftgt=0.2,thresh=0.9, ...){
  eqstk = brp(FLBRP(stock, ...))
  fbar(eqstk)[,1][] = 0.00001 # compute for eq Fmsy
  fbar(eqstk)[,1:101][] = ftgt # compute for equilibrium Fmsy
  eqstk = brp(eqstk)
  eqstk = window(as(eqstk, "FLStock"),start=2,end=2) # year1 F=0, year2=Fmsy
  eqstk@name = stock@name # name stk
  n_a = stock.n(eqstk)[-1,] # remove first age
  ages = dims(n_a)$min:dims(n_a)$max
  cums = apply(n_a,2:6,cumsum)
  n_thresh = sum(n_a*thresh)
  aref = min(ages[which((n_thresh-cums)^2==min((n_thresh-cums)^2))]+1,range(eqstk)["plusgroup"]-1)
  rp = sum(n_a[ac(aref:max(ages)),])/sum(n_a) # ref proportion
  flq= quantSums(stock.n(stock)[ac(aref:range(stock)[2]),])/quantSums(stock.n(stock)[-1,])/rp
  
  return(flq)
}
#}}}

