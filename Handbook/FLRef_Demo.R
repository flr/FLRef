## ----wrap-hook, echo = FALSE--------------------------------------
library(knitr)
hook_output = knit_hooks$get('output')
knit_hooks$set(output = function(x, options) {
  # this hook is used only when the linewidth option is not NULL
  if (!is.null(n <- options$linewidth)) {
    x = knitr:::split_lines(x)
    # any lines wider than n should be wrapped
    if (any(nchar(x) > n)) x = strwrap(x, width = n)
    x = paste(x, collapse = '\n')
  }
  hook_output(x, options)
})


## ---- echo = FALSE------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "  " ,fig.align = 'center', cache=FALSE,tidy.opts=list(width.cutoff=80), tidy=TRUE)



## ---- eval=FALSE--------------------------------------------------
## installed.packages("devtools")
## 
## installed.packages("ggplot2")
## 
## installed.packages("ggpubr")
## 
## installed.packages("TMB")
## 
## 
## devtools::install_github("flr/FLCore")
## 
## devtools::install_github("flr/FLBRP")
## 
## devtools::install_github("flr/FLasher")
## 
## devtools::install_github("flr/mse")
## 
## devtools::install_github("flr/ggplotFL")
## 
## devtools::install_github("flr/FLSRTMB")
## 
## devtools::install_github("henning-winker/FLRef")
## 
## # only for demo
## install.packages("ggpubr")
## 
## 


## ----message=FALSE------------------------------------------------
library(FLCore)
library(FLBRP)
library(FLasher)
library(FLSRTMB)
library(ggplotFL)
library(FLRef)
library(ggpubr) # For this demo


## ---- eval=TRUE---------------------------------------------------

data(ple4)
stk = ple4

## ----fig1, fig.height=4, fig.cap = "Estimated stock trajectories"----
plot(stk)+theme_bw()+ facet_wrap(~qname,scales="free",ncol=2)


## ---- eval=TRUE---------------------------------------------------

fbrps = computeFbrps(stock=ple4, proxy="sprx",f0.1=TRUE, verbose = FALSE)



## ----fig2, fig.height=4.5, fig.cap = "Estimated per-recruit reference points corresponding to the $F_{BPR}$'s $F_{spr35-50}$ and $F_{0.1}$."----

ploteq(fbrps)



## ----fig3, fig.height=4.5, fig.cap = "Estimated per-recruit reference points corresponding to the $F_{BPR}$, $F_{spr35-50}$ and $F_{0.1}$."----

ploteq(fbrps,refpts = "fmax")



## ---- eval=TRUE---------------------------------------------------

fbrp = computeFbrp(stock=stk, proxy=c("f0.1"),blim=0.25,type="btrg", verbose = FALSE)

Fbrp(fbrp)



## ---- eval=TRUE---------------------------------------------------

fbrp = computeFbrp(stock=stk, proxy=c("sprx","f0.1"),x=40,blim=0.25,type="btrg", verbose = FALSE)

Fbrp(fbrp)



## ----fig4, fig.height=4.5, fig.cap = "Estimated per-recruit reference points corresponding to the"----

ploteq(fbrp,refpts = "fmax")



## ----fig5, fig.height=4.5, fig.cap = "Stock advice plot showing the modelled quantities from a per-recruit perspestive relative to per-recruit based reference points"----

plotAdvice(stk,fbrp)



## ---- eval=TRUE---------------------------------------------------

object = as.FLSR(stk,model=geomean)

gm = srrTMB(object)



## ---- eval=TRUE---------------------------------------------------

fbrp = computeFbrp(stock=stk,sr=gm, proxy=c("sprx","f0.1"),x=40,blim=0.25,type="btrg", verbose = FALSE)


## ----fig6, fig.height=4.5, fig.cap = "Estimated  reference points relative to estimates of recruitment, $SSB$, $F$ and landings"----

ploteq(fbrp,refpts="fmax",obs=TRUE)



## ----fig7, fig.height=4.5, fig.cap = "Stock advice plot showing modelled quantities and the corresponding reference points"----

plotAdvice(stk,fbrp)



## ---- eval=TRUE---------------------------------------------------

bh = srrTMB(as.FLSR(stk,model=bevholtSV),spr0=spr0y(stk),verbose = FALSE)

bh@SV


## ---- eval=TRUE---------------------------------------------------

ri = srrTMB(as.FLSR(stk,model=rickerSV),spr0=spr0y(stk),verbose = FALSE)



## ---- eval=TRUE---------------------------------------------------

hs = srrTMB(as.FLSR(stk,model=segreg),spr0=spr0y(stk),lplim=0.05,uplim=0.2)



## ----fig8, fig.height=3.5, fig.cap = "Comparison of the S-R relationship fitted by assuming Beverton-Holt (bh), Ricker (ri) and Hockey-Stick (hs) function. Open circles show the observed S-R pairs with the solid do denoting the final assessment year"----

srs = FLSRs(bh=bh,ri=ri,hs=hs)

plotsrs(srs)



## ---- eval=TRUE---------------------------------------------------

hs@SV[["BlimB0"]]



## ----fig8a, fig.height=3.5, fig.cap = "Annual $SPR_{0_y}$ as function of time-varying $w_{a_y}$ (here), $mat_{a_y}$ and $M_{a_y}$."----

plot(spr0y(stk))+theme_bw()+
  ylab(expression(SPR[0]))+xlab("Year")+
  geom_hline(yintercept = mean(spr0y(stk)),linetype="dashed")



## ----fig8b, fig.height=3.5, fig.cap = "Comparison of the Hockey-Stick specified with (1) $SRP_{5-20}$ and time-varying $SPR_{0,y}$, (2) the same but with the mean of $SPR_{0,y}$ (3) $SRP_{7-20}$ and mean $SPR_{0,y}$."----


hs1 = srs$hs
hs2 = srrTMB(as.FLSR(stk,model=segreg),
          spr0=mean(spr0y(stk)),lplim=0.05,uplim=0.2)
hs3 = srrTMB(as.FLSR(stk,model=segreg),
             spr0=mean(spr0y(stk)),lplim=0.07,uplim=0.2)

plotsrs(FLSRs(plim0.05=hs1,muSPR0=hs2,plim0.07=hs3))



## ---- eval=TRUE---------------------------------------------------

hsblim(hs1)
hsblim(hs2)
hsblim(hs3)




## ---- eval=TRUE---------------------------------------------------

# Extract Blim
blim = c(params(hs)["b"])
# check break-point relative to B0
hsblim(hs)["SRPlim"]
hsblim(hs3)["SRPlim"]



## ----fig9, fig.height=5.5, fig.cap = "Comparison of the S-R relationship fitted by assuming Beverton-Holt (bh), Ricker (ri) and Hockey-Stick (hs) function."----

p1 = plotsrs(srs,path=FALSE)
p2 = plotsrs(srs,path=TRUE)
p3 = plotsrs(srs,b0=TRUE)
p4 = plotsrs(srs,b0=TRUE,rel=TRUE)

ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="right")



## ----fig10, fig.height=5., fig.cap = "Estimated reference points at equilibrium recruitment, $SSB$, $F$ and landings"----
brps = FLBRPs(
  bh = computeFbrp(stk,sr=srs[["bh"]],proxy=c("f0.1","sprx","msy"),x=40,blim=0.25,type="btrg",verbose = FALSE),
  ri = computeFbrp(stk,sr=srs[["ri"]],proxy=c("f0.1","sprx","msy"),x=40,blim=0.25,type="btrg",verbose = FALSE),
  hs = computeFbrp(stk,sr=srs[["hs"]],proxy=c("f0.1","sprx","msy"),x=40,blim=blim,type="value",verbose = FALSE)
)

# plot
ploteq(brps)



## ----fig10a, fig.height=5., fig.cap = "Estimated  reference points relative to estimates of recruitment, $SSB$, $F$ and landings"----

# plot
ploteq(brps,obs=TRUE)



## ----fig11, fig.height=3.5, fig.cap = "Comparison of alternative parameterisation of the Beverton Holt S-R"----

# Fixed steepness
s = c(seq(0.8,0.95,0.05))
fixs = FLSRs(lapply(as.list(s),function(x){
  srrTMB(as.FLSR(stk,model=bevholtSV),spr0=spr0y(stk),s=x,s.est = FALSE)
}))
names(fixs) = paste0("s=",s)
# with prior with mean s=0.85 and s.logitsd = 0.3
s.pr = srrTMB(as.FLSR(stk,model=bevholtSV),spr0=spr0y(stk),s=0.8,s.logitsd = 0.3)
# uncontrained estimate
s.est = srrTMB(as.FLSR(stk,model=bevholtSV),spr0=spr0y(stk),s=0.8)

#combine
bhs = FLSRs(c(s.est=s.est,s.pr=s.pr,fixs))
# add s estimate
names(bhs)[1:2] = c(paste0("s.est(",round(s.est@SV[["s"]],2),")"),
  paste0("s.pr(",round(s.pr@SV[["s"]],2),")"))

plotsrs(bhs)



## ----fig12, fig.height=4.5, fig.cap = "Comparison of equilibriun curves and reference points for alternative parameterisation of the Beverton Holt S-R"----

bh.brps = FLBRPs(lapply(bhs,function(x){
  computeFbrp(stk,x,proxy=c("f0.1","msy"),blim=0.25,type="btrg",verbose=FALSE)
}))

ploteq(bh.brps,obs=TRUE,panels=4)



## ----fig13, fig.height=4.5, fig.cap = "Stock advice plot showing modelled quantities and the corresponding reference points for a Beverton S-R model with estimated $s$"----

plotAdvice(stk,bh.brps[[1]])+
  ggtitle(paste0(stk@name,": BevHolt with s = ",round(bhs[[1]]@SV[["s"]],3)))



## ----fig14, fig.height=4.5, fig.cap = "Stock advice plot showing modelled quantities and the corresponding reference points for a Beverton S-R model with estimated $s$"----

pars= Fbrp(bh.brps[[1]])
pars
p1= plotAR(bh.brps[[1]],obs=stk,kobe=FALSE,bpa=1.4)
p2= plotAR(bh.brps[[1]],obs=stk,kobe=FALSE,
           bpa=1.4,btrigger=0.7)
p3= plotAR(bh.brps[[1]],obs=stk,kobe=TRUE,
           bpa=1.4,btrigger=0.7,bclose=1,fmin=0.01)
p4= plotAR(bh.brps[[1]],obs=stk,kobe=TRUE,
           bpa=1.4,btrigger=0.7,rel=TRUE)

ggarrange(p1, p2, p3, p4, ncol=2, nrow=2)


