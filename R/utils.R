

# rGclass {{{ 
#' Function to characterize Productivity and refpts based on r and Generation
#'
#' @param r value of the intrinsic rate of population increase
#' @param gt generation time G
#'
#' @return list with Productivity category and suggest Fbrps
#' @export
rGclass = function(r=NULL,gt=NULL){
  rg = data.frame(Low=c(0.00001,0.15),Medium=c(0.150001,0.5),High=c(0.500001,1))
  gg = data.frame(Low=c(50,10.001),Medium=c(10,5),High=c(5,0))
  
  Fspr = c(40,40,50)
  Fsb = c(35,35,40)
  
  selr = selg = 100
  if(is.null(r)==FALSE){
    mur = apply(rg,2,min)
    selr = max(which(r>mur))
  }
  if(is.null(gt)==FALSE){
    mug = apply(gg,2,min)  
    selg = min(which(gt>mug))
  }
  sel = min(selg,selr)
  
  category = names(rg)[sel]  
  return(list(class=category,Fspr=Fspr[sel],Fsb=Fsb[sel]))
}

# #{{{ color options
#' r4sscol
#' @param n number of colors
#' @param alpha transluscency 
#' @return vector of color codes
#' @export
rc4 <- function(n,alpha=1){
  # a subset of rich.colors by Arni Magnusson from the gregmisc package
  # a.k.a. rich.colors.short, but put directly in this function
  # to try to diagnose problem with transparency on one computer
  x <- seq(0, 1, length = n)
  r <- 1/(1 + exp(20 - 35 * x))
  g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
  b <- dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))
  rgb.m <- matrix(c(r, g, b), ncol = 3)
  rich.vector <- apply(rgb.m, 1, function(v) rgb(v[1], v[2], v[3], alpha=alpha))
  
  return(rich.vector)
}

#' ss3col
#' @param n number of colors
#' @param alpha transluscency 
#' @return vector of color codes
#' @export
ss3col <- function(n,alpha=1){
  if(n>3) col <- rc4(n+1)[-1]
  if(n<3)  col <- rc4(n)
  if(n==3) col <- c("blue","green","red")
  if(alpha<1){
    # new approach thanks to Trevor Branch
    cols <- adjustcolor(col, alpha.f=alpha)
  } else {
    cols=col
  }
  return(cols)
}


#' huecol
#' @param n number of colors
#' @param alpha transluscency 
huecol <- function(n,alpha=1) {
  hues = seq(15, 375, length = n + 1)
  adjustcolor(hcl(h = hues, l = 65, c = 100)[1:n],alpha.f=alpha)
}

# }}}# color options


#' stockMedians
#' @param stock class FLStock or FLStockR
#' @return medians of all FLstock FlQuants 
#' @export
stockMedians <- function(stock){
  stk = stock
  stk@catch.n = iterMedians(stk@catch.n)  
  stk@stock.n = iterMedians(stk@stock.n)  
  stock@landings.n =iterMedians(stk@catch.n)
  stock@discards.n =iterMedians(stk@discards.n)
  stk@catch.wt = iterMedians(stk@catch.wt)  
  stk@stock.wt = iterMedians(stk@stock.wt)  
  stock@landings.wt =iterMedians(stk@landings.wt)
  stock@discards.wt =iterMedians(stk@discards.wt)
  stock@m =iterMedians(stk@m)
  stock@mat =iterMedians(stk@mat)
  stock@harvest = iterMedians(stk@harvest)
  stock@harvest.spwn = iterMedians(stk@harvest.spwn)
  stock@m.spwn = iterMedians(stk@m.spwn)
  if(class(stock)=="FLStockR"){
    stk@refpts = iterMedians(stk@refpts)
  }
  stk@catch = computeCatch(stk)
  stk@discards = computeDiscards(stk)
  stk@landings = computeLandings(stk)
  stk@stock = computeStock(stk)
  return(iter(stk,1))
}


#' ssmvln()
#'
#' function to generate uncertainty for Stock Synthesis 
#'
#' @param ss3rep from r4ss::SS_output
#' @param out choice c("iters","mle")
#' @param Fref  Choice of Fratio c("MSY","Btgt","SPR","F01"), correponding to F_MSY and F_Btgt                                                               
#' @param years single year or vector of years for mvln   
#' @param virgin if FALSE (default) the B0 base for Bratio is SSB_unfished
#' @param mc number of monte-carlo simulations   
#' @param weight weighting option for model ensembles weight*mc 
#' @param run qualifier for model run
#' @param plot option to show plot
#' @param ymax ylim maximum
#' @param xmax xlim maximum
#' @param addprj include forecast years
#' @param legendcex=1 Allows to adjust legend cex
#' @param verbose Report progress to R GUI?
#' @param seed retains interannual correlation structure like MCMC 
#' @return output list of quant posteriors and mle's
#' @author Henning Winker (GFCM)
#' @export

ssmvln = function(ss3rep,Fref = NULL,years=NULL,virgin=FALSE,mc=1000,weight=1,run="MVLN",
                       addprj=FALSE,ymax=NULL,xmax=NULL,legendcex=1,verbose=TRUE,seed=123){
  
  
  status=c('Bratio','F')
  quants =c("SSB","Recr")
  mc = round(weight*mc,0)
  hat = ss3rep$derived_quants
  cv = cv1 = ss3rep$CoVar
  
  
  
  if(is.null(cv)) stop("CoVar from Hessian required")
  # Get years
  allyrs = unique(as.numeric(gsub(paste0(status[1],"_"),"",hat$Label[grep(paste0(status[1],"_"), hat$Label)])))[-1]
  allyrs = allyrs[!is.na(allyrs)] 
  
  if(is.null(years) & addprj==TRUE) yrs = allyrs   
  if(is.null(years) & addprj==FALSE) yrs = allyrs[allyrs<=ss3rep$endyr]
  if(is.null(years)==FALSE) yrs = years[years%in%allyrs==TRUE]
  estimate = ifelse(yrs<=ss3rep$endyr,"fit","forecast")
  
  # brp checks for starter file setting
  refyr = max(yrs)
  endyr = ss3rep$endyr
  bt = hat[hat$Label==paste0("SSB_",endyr),2]
  ft = hat[hat$Label==paste0("F_",endyr),2]
  # virgin
  bv =hat[hat$Label%in%c("SSB_virgin","SSB_Virgin"),2]
  # unfished
  b0 =hat[hat$Label%in%c("SSB_unfished","SSB_Unfished"),2]
  fmsy = hat[hat$Label%in%c("Fstd_MSY","annF_MSY"),2]
  fspr = hat[hat$Label%in%c("Fstd_SPR","annF_SPR"),2]
  bmsy = hat[hat$Label==paste0("SSB_MSY"),2]
  
  bf01= f01=option.btgt = FALSE
  
  if("SSB_Btgt"%in%hat$Label){
    btgt = hat[hat$Label==paste0("SSB_Btgt"),2]
    if(!is.null(Fref)) if(Fref=="F01") stop("F01 not defined: choose Fref = MSY or Btgt ")
  }
  if("SSB_F01"%in%hat$Label){
    f01=TRUE
    if(is.null(Fref)) Fref = "Btgt"
    btgt = hat[hat$Label==paste0("SSB_F01"),2]*b0 #  Why SSB_F01 ratio???
    if(round(hat[hat$Label%in%c("annF_F01"),2],2)==round(fmsy,2)){
      option.btgt=TRUE 
      if(!is.null(Fref))if(Fref=="MSY") stop("FMSY not defined: choose Fref = F01 ")
      bf01= TRUE}
  }
  
  
  bratio = hat[hat$Label==paste0("Bratio_",endyr),2]
  bb.check = c(bt/bv,bt/bmsy,bt/btgt)
  
  option.btgt = FALSE
  if(btgt==bmsy){
    option.btgt = TRUE
    if(is.null(Fref)) Fref = "Btgt"
    if(Fref=="MSY") stop("FMSY not defined: choose Fref = Btgt ")
    
  }
  if(fmsy==fspr){
    if(is.null(Fref)) Fref = "SPR"
    if(Fref=="MSY") stop("FMSY not defined: choose Fref = SPR")
  }
  
  # bratio definition
  bb = max(which(abs(bratio-bb.check)==min(abs(bratio-bb.check))))   
  if(bf01) bb=4
  
  bbasis  = c("SSB/SSB0","SSB/SSBMSY","SSB/SSBtgt","SSB/SSBF01")[bb]
  fbasis = strsplit(ss3rep$F_report_basis,";")[[1]][1]
  if(is.na(ss3rep$btarg)) ss3rep$btarg=0
  gettrg = ifelse(ss3rep$btarg>0,ss3rep$btarg,round(btgt/b0,2))
  if(fbasis%in%c("_abs_F","(F)/(Fmsy)",paste0("(F)/(F_at_B",ss3rep$btarg*100,"%)"),paste0("(F)/(F",ss3rep$btarg*100,"%SPR)"))){
    fb = which(c("_abs_F","(F)/(Fmsy)",paste0("(F)/(F_at_B",ss3rep$btarg*100,"%)"),
                 paste0("(F)/(F",ss3rep$btarg*100,"%SPR)"))%in%fbasis)
  } else { stop("F_report_basis is not defined, please rerun Stock Synthesis with recommended starter.ss option for F_report_basis: 1")}
  
  if(is.null(Fref) & fb%in%c(1,2)) Fref = "MSY"
  if(is.null(Fref) & fb%in%c(3)) Fref = "Btgt"
  if(is.null(Fref) & fb%in%c(4)) Fref = "SPR"
  if(Fref=="F01") Fref="Btgt" # hack to avoid redundancy later
  
  if(verbose) cat("\n","starter.sso with Bratio:",bbasis,"and F:",fbasis,"\n","\n")
  
  bref  = gettrg
  
  if(fb==4 & Fref[1] %in% c("Btgt","MSY")) stop("Fref = ",Fref[1]," option conflicts with ",fbasis," in starter.sso, please choose Fref = SPR")
  if(fb==2 & Fref[1] %in% c("Btgt","SPR")) stop("Fref = ",Fref[1]," option conflicts with ",fbasis," in starter.sso, please choose Fref = MSY")
  if(fb==3 & Fref[1] %in% c("Btgt","MSY")) stop("Fref = ",Fref[1]," option conflicts with ",fbasis,", in starter.sso, please choose Fref = Btgt")
  if(fb%in%c(1,2) &  Fref[1] =="MSY") Fquant = "MSY"
  if(fb%in%c(1,3) & Fref[1] =="Btgt") Fquant = "Btgt"
  if(fb%in%c(1,4) & Fref[1] =="SPR") Fquant = "SPR"
  if(Fquant == "Btgt" & f01) Fquant = "F01"
  
  # check ss3 version
  if("Fstd_MSY"%in%hat$Label){Fname = "Fstd_"} else {Fname="annF_"}
  cv <- cv[cv$label.i %in% paste0(status,"_",yrs),]
  cv1 = cv1[cv1$label.i%in%paste0(Fname,Fquant) & cv1$label.j%in%paste0(status,"_",yrs),]
  fref = hat[hat$Label==paste0(Fname,Fquant),]
  cv$label.j[cv$label.j=="_"] <- cv$label.i[cv$label.j=="_"]
  
  if(is.null(hat$Label)){ylabel = hat$LABEL} else {ylabel=hat$Label}
  kb=mle = NULL
  for(yi in 1:length(yrs)){ 
    set.seed(seed)
    yr = yrs[yi]
    x <- cv[cv$label.j %in% paste0(status[2],"_",c(yr-1,yr,yr+1)) & cv$label.i %in% paste0(status[1],"_",c(yr-1,yr,yr+1)),]
    x1 = cv1[cv1$label.j %in% paste0(status[1],"_",c(yr-1,yr,yr+1)),] 
    x2 = cv1[cv1$label.j %in% paste0(status[2],"_",c(yr-1,yr,yr+1)),] 
    y = hat[ylabel %in% paste0(status,"_",yr),] # old version Label not LABEL
    y$Value[1] = ifelse(y$Value[1]==0,0.001,y$Value[1])
    varF = log(1+(y$StdDev[1]/y$Value[1])^2) # variance log(F/Fmsy)  
    varB = log(1+(y$StdDev[2]/y$Value[2])^2) # variance log(SSB/SSBmsy)  
    varFref = log(1+(fref$StdDev[1]/fref$Value)^2) # variance log(F/Fmsy)  
    cov = log(1+mean(x$corr)*sqrt(varF*varB)) # covxy
    cov1 = log(1+mean(x1$corr)*sqrt(varB*varFref)) # covxy
    cov2 = log(1+mean(x2$corr)*sqrt(varF*varFref)) # covxy
    # MVN means of SSB/SBBmsy, Fvalue and Fref (Ftgt or Fmsy) 
    mvnmu = log(c(y$Value[2],y$Value[1],fref$Value)) # Assume order F_ then Bratio_ 
    # Create MVN-cov-matrix
    mvncov = matrix(NA,ncol=3,nrow=3)
    diag(mvncov) = c(varB,varF,varFref)
    mvncov[1,2] = mvncov[2,1] = cov 
    mvncov[2,3] = mvncov[3,2] = cov1 
    mvncov[1,3] = mvncov[3,1] = cov2 
    kb.temp = data.frame(year=yr,run=run,type=estimate[yi],iter=1:mc,exp(mvtnorm::rmvnorm(mc ,mean = mvnmu,sigma = mvncov,method=c( "svd")))) # random  MVN generator
    colnames(kb.temp) = c("year","run","type","iter","stock","harvest","F")
    if(length(quants)>0){
      quant=NULL
      for(qi in 1:length(quants)){
        qy = hat[ylabel %in% paste0(quants[qi],"_",yr),]
        qsd = sqrt(log(1+(qy$StdDev[1]/qy$Value[1])^2))
        quant = cbind(quant,rlnorm(mc,log(qy$Value[1])-0.5*qsd*qsd,qsd))
      }     
      colnames(quant) = quants    
      kb.temp = cbind(kb.temp,quant)
    }
    kb = rbind(kb,cbind(kb.temp))
    mle = rbind(mle,data.frame(year=yr,run=run,type=estimate[yi],stock=y$Value[2],harvest=y$Value[1],F=fref$Value[1])) 
  }
  # add mle quants
  qmles = NULL
  for(qi in 1:length(quants)){
    qmles = cbind(qmles, hat[ylabel %in% paste0(quants[qi],"_",yrs),]$Value)
  } 
  colnames(qmles) = quants    
  mle = cbind(mle,qmles)
  
  mle = mle[,c(1:5,7,6,8)]
  kb = kb[,c(1:6,8,7,9)]
  
  # virgin or unfished?
  if(!virgin){
    if(bb%in%c(1,3)){
      if(bb==1 | bb==3 & !option.btgt){
        kb[,"stock"] = kb[,"stock"]*(bv/b0)
        mle[,"stock"] = mle[,"stock"]*(bv/b0)
      }}}
  if(virgin){ # reverse correction
    if(bb==3 & option.btgt){
      kb[,"stock"] = kb[,"stock"]*(b0/bv)
      mle[,"stock"] = mle[,"stock"]*(b0/bv)
    }}
  
  
  # Take ratios
  if(bb==1){
    kb[,"stock"] = kb[,"stock"]/bref  
    mle[,"stock"] = mle[,"stock"]/bref
  }
  
  if(fb> 1){
    kb[,"F"] = kb[,"F"]*kb[,"harvest"] 
    mle[,"F"] = mle[,"F"]*mle[,"harvest"] 
    
  } else {
    fi = kb[,"harvest"]
    fm = mle[,"harvest"]
    kb[,"harvest"] = kb[,"harvest"]/kb[,"F"]
    kb[,"F"] = fi
    mle[,"harvest"] = mle[,"harvest"]/mle[,"F"]
    mle[,"F"] = fm
  }    
  
  
  # Add catch
  C_obs = aggregate(Exp~Yr,ss3rep$catch,sum)
  #colnames(C_obs) = c("Yr","Obs")
  Cobs = C_obs[C_obs$Yr%in%yrs,]
  foreyrs = unique(as.numeric(gsub(paste0("ForeCatch_"),"",hat$Label[grep(paste0("ForeCatch_"), hat$Label)])))
  Cfore = data.frame(Yr=foreyrs,Obs=hat$Value[hat$Label%in%paste0("ForeCatch_",foreyrs)] )
  names(Cfore) = names(Cobs)
  Catch = rbind(Cobs,Cfore)
  Catch = Catch[Catch$Yr%in%yrs,]
  kb$Catch = rep(Catch[,2],each=max(kb$iter))
  mle$Catch = Catch[,2]
  trg =round(bref*100,0)
  spr = round(ss3rep$sprtarg*100,0)
  xlab = c(bquote("SSB/SSB"[.(trg)]),expression(SSB/SSB[MSY]),bquote("SSB/SSB"[.(trg)]),expression(SSB/SSB[F0.1]))[bb] 
  ylab = c(expression(F/F[MSY]),
           bquote("F/F"[SB~.(trg)]),
           bquote("F/F"[SPR~.(spr)]),
           expression(F/F[0.1])
  )[which(c("MSY","Btgt","SPR","F01")%in%Fquant)] 
  
  labs = ifelse(quants=="Recr","Recruits",quants)
  refB = c(paste0("B",trg),"Bmsy",paste0("B",trg),"BF0.1")[bb] 
  refF = c("Fmsy",
           paste0("Fb",trg),
           paste0("F",spr),
           "F0.1")[which(c("MSY","Btgt","SPR","F01")%in%Fquant)] 
  refpts = data.frame(RefPoint=c("Ftgt","Btgt","MSY","B0","R0"),value=c((mle$F/mle$harvest)[1],
                      (mle$SSB/mle$stock)[1],
                      MSY=hat$Value[hat$Label=="Dead_Catch_MSY"],
                      B0=b0,
                      MSY=hat$Value[hat$Label=="Recr_unfished"]))
  
  return(list(kb=kb,mle=mle,refpts=refpts, quants=c("stock","harvest","SSB","F","Recr","Catch"),
              labels=c(xlab,ylab,labs[1],"F",labs[2],"Catch"),Btgtref = bref))
} # End 

# {{{
#' blag()
#'
#' function to assign B[y+1] to B[y]. Warning correlation structure of B[y+1] and F[y] is meaningless
#'
#' @param mvn 
#' @return output list of quant posteriors and mle's
#' @author Henning Winker (GFCM)
#' @export

blag <- function(mvn,verbose=TRUE){
  d1 = mvn$mle
  d2 = mvn$kb
  endyr = max(mvn$mle$year)  
  styr = min(mvn$mle$year)
  d1B = d1[d1$year!=styr,]
  d2B = d2[d2$year!=styr,]
  d1 = d1[d1$year!=endyr,] 
  d2 = d2[d2$year!=endyr,] 
  d1[,c("stock","SSB")] = d1B[,c("stock","SSB")] 
  d2[,c("stock","SSB")] = d2B[,c("stock","SSB")] 
  out = mvn
  out$mle =d1
  out$kb=d2
  return(out)
}
# }}}

# {{{
#' Mlorenzen
#' 
#' computes Lorenzen M with scaling option
#' @param object weight-at-age of class *FLQuant* 
#' @param Mref reference M for scaling  
#' @param Aref reference Age for scaling                                       
#' @return FLQuant m()
#' @export
#' @examples 
#' data(ple4)
#' Ml = Mlorenzen(stock.wt(ple4))
#' # Scale
#' Ms = Mlorenzen(stock.wt(ple4),Mref=0.2,Aref=2)
#' flqs = FLQuants(Lorenzen=Ml,Scaled=Ms)

Mlorenzen = function(object,Mref="missing",Aref=2){
  Ml = 3*((1000*object)^(-0.288))
  if(missing("Mref"))
    out = Ml
  if(!missing("Mref"))
    out = Ml*an(Mref/Ml[ac(Aref),])
  
  units(out) =  "m"
  return(out)
} 
#}}}



