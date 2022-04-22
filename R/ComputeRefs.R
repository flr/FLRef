
#{{{
#' computeFbrp()
#
#' Computes biological reference points corresponding to the proxy Fbrp
#'
#' @param stock object of class FLStock 
#' @param srr stock recruitment model of class FLSR
#' @param proxies choice of Fmsy proxies (combinations permitted) 
#' \itemize{
#'   \item "sprx"  spawning potential ratio spr/spr0 with basis x 
#'   \item "bx" SSB as fraction xSSB0
#'   \item "f0.1" 10% slope of yield-per-recruit curve
#'   \item "msy"  maximum surplus production (not defined for segreg)
#' }    
#' @param x basis in percent for sprx and bx 
#' @param blim values < 1 are taken as fraction to B0 and blim > 1 as absolute values unless specified otherwise
#' @param type type of blim input, values < 1 are  
#' \itemize{
#'   \item "b0" fraction to B0  
#'   \item "btrg" fraction to Btarget (first occurring in proxy)    
#'   \item "value" absolute value
#' }
#' @param btri Btrigger can specified as absolute value
#' @param bpa Bpa can specified as absolute value
#' @param fmax maximum Flim = max(Flim,fmax*Fbrp)
#' @param verbose   
#' @return brp object of class FLBRP with computed Fbrp reference points   
#' @export

computeFbrp <- function(stock,srr=NULL,proxy=c("sprx","bx","f0.1","msy"),x=40,blim=0.1,type=c("b0","btrg","value"),btri="missing",bpa="missing",verbose=T,fmax=10){
   
  # use geomean sr if sr = NULL (only growth overfishing)
  if(is.null(srr)){
    srr = fmle(as.FLSR(stock,model=geomean),method="BFGS")
    if(verbose)cat(paste0("Computing geomean from S-R data in the absense of a specified SRR"),"\n")
  }  
  
  Fbrp = FLPar(none=NA)
  
  
  type = type[1]
  brp = brp(FLBRP(stock,srr))
  if(c("sprx")%in%proxy){
    pr = brp(FLBRP(stock)) # per-recruit
    Fbrp = rbind(Fbrp,FLPar(Fspr = refpts(pr+FLPar(Btrg=refpts(pr)["virgin","ssb"]*x*0.01))["Btrg","harvest"]))
    if(verbose)cat(paste0("Computing Fspr",x," with Btrg = Bspr",x),"\n")
  }
  if(c("bx")%in%proxy){
    Fbrp = rbind(Fbrp,FLPar(Fb = refpts(brp+FLPar(Btrg=refpts(brp)["virgin","ssb"]*x*0.01))["Btrg","harvest"]))
    if(verbose)cat(paste0("Computing Fsb",x," with Btrg = Bsb",x),"\n")
  }  
  if(c("f0.1")%in%proxy){
    Fbrp = rbind(Fbrp,FLPar(F0.1=refpts(brp)["f0.1","harvest"]))
    if(verbose)cat(paste0("Computing F0.1 with Btrg = B[F.01]"),"\n")
  }
  
  if(c("msy")%in%proxy){
    if(SRModelName(model(srr))%in%c("segregA1","segreg","mean")){
       warning(paste0("MSY is not unambiguously defined for ",SRModelName(model(srr))," SR model","\n"))
    } else {
    #add warning for segreg
    Fbrp = rbind(Fbrp,FLPar(Fmsy=refpts(brp)["msy","harvest"]))
      if(verbose)cat(paste0("Computing Fmsy with Btrg = Bmsy"),"\n")
    }
  }
  
  ord = c("Fspr","Fb","F0.1","Fmsy")[match(proxy,c("sprx","bx","f0.1","msy"))]
  
  Fbrp = Fbrp[-1]
  
  Fbrp = Fbrp[match(rownames(Fbrp),ord)]
  
  B0 = an(refpts(brp)["virgin","ssb"])
  
  brpf =  brp+Fbrp
  if(blim<=1 & type!="value"){
    if(type%in%c("b0")){
      Blim =an(refpts(brpf)["virgin","ssb"])*blim
      if(verbose)cat("\n",paste0(" Blim = ",blim,"B0"),"\n")
    }
    if(type%in%c("btrg")){
      Blim =an(refpts(brpf)[rownames(Fbrp)[1],"ssb"])*blim
      if(verbose)cat("\n",paste0("Blim = ",blim," with Btrg corresponding to ", rownames(Fbrp)[1]),"\n")
    }
  } else {
    if(verbose)cat("\n",paste0("Blim as input value Blim = ",blim),"\n")
    Blim = blim
  }
  
  #subs = c("spr","b","0.1","msy")[which(c("sprx","bx","f0.1","msy")%in%proxy)] 
  #xi = which(c("sprx","bx")%in%proxy)
  
  
  # rename
  rownames(Fbrp)[which(rownames(Fbrp)%in%c("Fspr","Fb"))] = paste0(rownames(Fbrp)[which(rownames(Fbrp)%in%c("Fspr","Fb"))],x) 
  fref  = rownames(Fbrp)
  refs = rbind(Fbrp,FLPar(Blim=Blim,B0=B0))
  
  
  if(!missing(btri)){
    refs= rbind(refs,FLPar(Btri=btri))
  }
  
  if(!missing(bpa)){
    refs= rbind(refs,FLPar(Bpa=bpa))
  }
  
  
  # do check
  check = brp+refs
  if(refpts(check)["Blim","harvest"]>fmax*refpts(check)[paste(fref[1]),"harvest"]){
    refs["Blim"] = refpts(brp+FLPar(Flim=fmax*refpts(check)[paste(fref[1]),"harvest"]))["Flim","ssb"] 
  }
    
  brp = brp+refs
  #if(fmax*refpts(check)["Fbrp","harvest"]>max(fbar(brp))){
    fl = min(c(fmax*refpts(check)[paste(fref),"harvest"],refpts(check)["Blim","harvest"]*1.5))
    fbar(brp) = seq(0,fl,fl/100)
  #} 
  
  return(brp)
}
#}}}

#{{{
#' computeFbrps()
#
#' Computes biological reference points corresponding to the proxy Fbrp
#'
#' @param stock object of class FLStock 
#' @param srr stock recruitment model of class FLSR
#' @param proxies choice of Fmsy proxies
#' \itemize{
#'   \item "all"  both sprx and bx 
#'   \item "sprx"  spawning potential ratio spr/spr0 with basis x 
#'   \item "bx" SSB as fraction xSSB0
#' }    
#' @param fmsy if TRUE, Fmsy is computed (not suggest for segreg or geomean srr) 
#' @param f0.1 if TRUE, F0.1 is computed
#' @param verbose   
#' @return brp object of class FLBRP with computed Fbrp reference points   
#' @export

computeFbrps <- function(stock,srr=NULL,proxy=c("all","sprx","bx"),fmsy=FALSE,f0.1=TRUE,verbose=T){
  
  # use geomean sr if sr = NULL (only growth overfishing)
  if(is.null(srr)){
    srr = fmle(as.FLSR(stock,model=geomean),method="BFGS")
    if(verbose)cat(paste0("Computing geomean from S-R data in the absense of a specified SRR","\n"))
  }  
  proxy=proxy[1] 
 
  brp = brp(FLBRP(stock,srr))
  if(proxy%in%c("sprx","all")){
    pr = brp(FLBRP(stock)) # per-recruit
    Fsprs = FLPar(
      Fspr35 = an(refpts(pr+FLPar(Btrg=refpts(pr)["virgin","ssb"]*35*0.01))["Btrg","harvest"]),
      Fspr40 = an(refpts(pr+FLPar(Btrg=refpts(pr)["virgin","ssb"]*40*0.01))["Btrg","harvest"]),
      Fspr45 = an(refpts(pr+FLPar(Btrg=refpts(pr)["virgin","ssb"]*45*0.01))["Btrg","harvest"]),
      Fspr50 = an(refpts(pr+FLPar(Btrg=refpts(pr)["virgin","ssb"]*50*0.01))["Btrg","harvest"])
    )  
      if(verbose)cat(paste0("Computing Fspr% 35-50 with Btrg = Bspr"),"\n")
  }
  
  
  if(proxy%in%c("bx","all")){
    Bx = FLPar(
          B30=refpts(brp)["virgin","ssb"]*30*0.01,
          B35=refpts(brp)["virgin","ssb"]*35*0.01,
          B40=refpts(brp)["virgin","ssb"]*40*0.01,
          B45=refpts(brp)["virgin","ssb"]*45*0.01)
    fbrps = refpts(brp+Bx)[8:11,"harvest"] 
    FBs = FLPar(Fb30= fbrps[1,],Fb35= fbrps[2,],Fb40= fbrps[3,],Fb45= fbrps[4,] )
    
    if(verbose)cat("\n",paste0("Computing Fsb% 30-45 with Btrg = Bsb"),"\n")
  }  
  if(proxy=="all") Fbrps=rbind(Fsprs,FBs) 
  if(proxy=="sprx") Fbrps=Fsprs 
  if(proxy=="bx") Fbrps=FBs 
  if(fmsy) Fbrps = rbind(Fbrps,FLPar(Fmsy=refpts(brp)["msy","harvest"]))
  if(f0.1) Fbrps = rbind(Fbrps,FLPar(F0.1=refpts(brp)["f0.1","harvest"]))
  
  
  Fbrps = rbind(Fbrps, FLPar(B0 = an(refpts(brp)["virgin","ssb"])))
  
  brp =  brp(brp+Fbrps)
  
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
        Btrg=rpt[ref,"ssb"],
        Blim=rpt["Blim","ssb"],
        Flim=rpt["Blim","harvest"],
        Yeq = rpt[ref,"yield"],
        B0 = rpt["virgin","ssb"] 
        )
  rownames(out)[1] = ref
  if("Fp.05"%in%rownames(rpt)){
    out = rbind(out,FLPar(Fp.05=rpt["Fp.05","harvest"]))
  }
  if("Btri"%in%rownames(rpt)){
    out = rbind(out,FLPar(Btri=rpt["Btri","ssb"]))
  }
 
  
  
  return(out)
}
#}}}

