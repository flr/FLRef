
#{{{
#' computeFbrp()
#
#' Computes biological reference points corresponding to the proxy Fbrp
#'
#' @param stock object of class FLStock 
#' @param sr stock recruitment model of class FLSR
#' @param proxies choice of plots 1:4
#' \itemize{
#'   \item sprx  spawning potential ratio spr/spr0 with basis x 
#'   \item bx SSB as fraction xSSB0
#'   \item f0.1 10% slope of yield-per-recruit curve
#'   \item msy  maximum surplus production (not defined for segreg)
#' }    
#' @param x basis in percent for sprx and bx 
#' @param blim values < 1 are taken as fraction to B0 and blim > 1 as absolute values unless specified otherwise
#' @param type type of blim input, values < 1 are  
#' \itemize{
#'   \item b0  fraction to B0  
#'   \item btrg fraction to Btarget    
#'   \item value absolute value
#' }
#' @param btri Btrigger can specified as ratio to Btrg
#' @param verbose   
#' @return brp object of class FLBRP with computed Fbrp reference points   
#' @export

computeFbrp <- function(stock,sr=NULL,proxy=c("sprx","bx","f0.1","msy"),x=40,blim=0.1,type=c("b0","btrg","value"),btri="missing",verbose=T){
 
  # use geomean sr if sr = NULL (only growth overfishing)
  if(is.null(sr)){
    sr = fmle(as.FLSR(stock,model=geomean),method="BFGS")
    if(verbose)cat(paste0("Computing geomean from S-R data in the absense of a specified SRR","\n"))
  }  
  proxy=proxy[1] 
  type = type[1]
  brp = brp(FLBRP(stock,sr))
  if(proxy%in%c("sprx")){
    pr = brp(FLBRP(stock)) # per-recruit
    Fbrp = an(refpts(pr+FLPar(Btrg=refpts(pr)["virgin","ssb"]*x*0.01))["Btrg","harvest"])
    if(verbose)cat(paste0("Computing Fspr",x," with Btrg = Bspr",x," and"))
  }
  if(proxy%in%c("bx")){
    Fbrp = an(refpts(brp+FLPar(Btrg=refpts(brp)["virgin","ssb"]*x*0.01))["Btrg","harvest"])
    if(verbose)cat(paste0("Computing Fsb",x," with Btrg = Bsb",x, " and"))
  }  
  if(proxy%in%c("f0.1")){
    Fbrp = an(refpts(brp)["f0.1","harvest"])
    if(verbose)cat(paste0("Computing F0.1 with Btrg = B[F.01]"," and"))
  }
  
    
  if(proxy%in%c("msy")){
    if(SRModelName(model(sr))%in%c("segregA1","segreg","geomean")) 
       stop(paste0("MSY is not unambiguously defined for ",SRModelName(model(sr))," SR model","\n"))
    
    #add warning for segreg
    Fbrp = an(refpts(brp)["msy","harvest"])
    if(verbose)cat(paste0("Computing Fmsy with Btrg = Bmsy"," and"))
  }
  
  B0 = an(refpts(brp)["virgin","ssb"])
  
  brpf =  brp+FLPar(Fbrp=Fbrp)
  if(blim<=1 & type!="value"){
  if(type%in%c("b0")){
    Blim =an(refpts(brpf)["virgin","ssb"])*blim
    if(verbose)cat(paste0(" Blim = ",blim,"B0"),"\n")
  }
  if(type%in%c("btrg")){
    Blim =an(refpts(brpf)["Fbrp","ssb"])*blim
    if(verbose)cat(paste0("Blim = ",blim,"Btrg"),"\n")
  }
  } else {
    if(verbose)cat(paste0("Blim as input value Blim = ",blim),"\n")
    Blim = blim
  }
  
  refs = FLPar(Fbrp=Fbrp,Blim=Blim,B0=B0)
  if(!missing(btri)){
    refs= rbind(refs,FLPar(Btri=btri*an(refpts(brpf)["Fbrp","ssb"])))
  }
  
  brp = brp(brp+refs)
  
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
  out = FLPar(Fbrp = rpt["Fbrp","harvest"],
        Btrg=rpt["Fbrp","ssb"],
        Blim=rpt["Blim","ssb"],
        Flim=rpt["Blim","harvest"],
        Yeq = rpt["Fbrp","yield"],
        B0 = rpt["virgin","ssb"] 
        )
  if("Fp.05"%in%rownames(rpt)){
    out = rbind(out,FLPar(Fp.05=rpt["Fp.05","harvest"]))
  }
  if("Btri"%in%rownames(rpt)){
    out = rbind(out,FLPar(Fp.05=rpt["Btri","ssb"]))
  }
  
  return(out)
}
#}}}

