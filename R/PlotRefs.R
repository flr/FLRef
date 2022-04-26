
#{{{
#' ploteq()
#
#' Modification of method plot(`FLBRP`) to plot equilibrium output of computeFbrp()  
#'
#' @param brps output object from computeFbrp of class FLBRP  
#' @param refpts Reference points, defaults are computed refpts from computeFbrp()  
#' \itemize{
#'   \item Fbrp  
#'   \item Blim 
#'   \item B0  
#'   \item Btri  
#' }    
#' @param obs Should observations be plotted? Defaults to `FALSE`.
#' @param dashed plots vertical dashed lines to highlight refpts locations
#' @param colours refpts colours, default is designed for computeFbrp() output
#' @param panels plot panel option 1:4 
#' @param ncol number of plot panel columns
#' @return ggplot  
#' @export

ploteq <- function(brps, refpts="missing", obs=FALSE, dashed=TRUE,
                   colours="missing" ,panels=NULL, ncol=2){
            
            if(class(brps)=="FLBRP") brps = FLBRPs(brps)
            defaults = c("virgin","msy","crash","f0.1","fmax","spr.30","mey")
            
            cname= names(brps)
            
            df = Map(function(x,y){
            # EXTRACT metrics
            data.frame(model.frame(metrics(x,list(ssb=ssb, harvest=fbar, rec=rec, yield=landings, profit=profit)),
                              drop=FALSE),cname=y)
            },x=brps,y=as.list(cname))
            
            # refpts
            drps <- as.vector(do.call(c,lapply(brps,function(x){
            dimnames(refpts(x))$refpt
            })))
            
            drps = unique(drps)
            # re-arrange
            drps = drps[c(which(drps%in%defaults),grep("B",drps),grep("F",drps))]
            if("B0"%in%drps) drps = c(drps[!drps=="B0"],"B0") 
            
            
          
            if(missing(refpts)){ 
              refpts = drps[drps%in%defaults==FALSE]
            } else {
              
              excl = defaults[defaults%in%refpts==FALSE]
              refpts = drps[drps%in%excl==FALSE]
            }
          
            
            rps = FLPars(Map(function(x,y){
            drpx =  drps[drps%in%rownames(refpts(x))]
            rp = refpts(x)[drpx %in% refpts,]
            dms <- dimnames(rp)
            rp[!dms$refpt %in% "mey",
                         !dms$quant %in% c("revenue", "cost", "profit")]
              
            },x=brps,y=cname)) 
            
            # estimated?
            rpf <- TRUE
            
            
            # NO economics
            plots <- list(
              P1=c(x="harvest", y="ssb", panel="SSB ~ F", pos=1),
              P2=c(x="ssb", y="rec", panel="Recruitment ~ SSB", pos=2),
              P3=c(x="harvest", y="yield", panel="Yield ~ F", pos=3),
              P4=c(x="ssb", y="yield", panel="Yield ~ SSB", pos=4))
            
            
            
            # SUBSET panels if not NULL
            if(!is.null(panels))
              plots <- plots[panels]
            
            # APPLY over plots to extract x, y and panel for each element
            
            
            dat <- lapply(plots, function(p) {
              pd = NULL
              for(i in 1:length(df)){
                pd = rbind(pd,data.frame(x=df[[i]][,p['x']], y=df[[i]][,p['y']], iter=df[[i]][,'iter'],
                         panel=p['panel'], pos=p['pos'], row.names=NULL,cname=cname[i]))
              }
              pd
              })
            
            # RBIND into single df
            dat <- do.call(rbind, c(dat, list(make.row.names = FALSE)))
            
            # Limit y to 0 in panels 1:4
            dat[dat$pos %in% 1:4 & dat$y<0, "y"] <- 0
            
            # CREATE facet labels vector
            facl <- setNames(unique(dat$panel), nm=unique(dat$pos))
            
            dat$cname = factor(dat$cname,levels=cname)
            
            # PLOT
            p <- ggplot(dat, aes_(x=~x, y=~y))  + theme_bw()+
              facet_wrap(~pos, scales="free", ncol=ncol, labeller=labeller(pos=facl)) +
              xlab("") + ylab("") +
              scale_x_continuous(expand = expansion(mult = c(0, .05)),labels=human_numbers, limits=c(0, NA))+
              scale_y_continuous(expand = expansion(mult = c(0, .1)))+
              theme(legend.title = element_blank())
            if(length(brps)>1) p <- p + geom_line(aes(color=cname),size=0.6)
            if(length(brps)==1) p <- p + geom_line(size=0.5)
            # PLOT observations
           
            if(obs) {
              
              dfo <- Map(function(x,y){ 
                 data.frame( model.frame(metrics(x,
                                         list(ssb=ssb.obs, harvest=fbar.obs, rec=rec.obs, yield=landings.obs,
                                              profit=profit.obs)), drop=FALSE),cname=y)
              },x=brps,y=as.list(cname))
              
              # APPLY over plots to extract x, y and panel for each element
              dato <- lapply(plots, function(p){
                pdo = NULL
                for(i in 1:length(brps)){
                pdo = rbind(pdo,data.frame(x=dfo[[i]][,p['x']], y=dfo[[i]][,p['y']], iter=dfo[[i]][,'iter'],
                           pos=p['pos'],cname=dfo[[i]]$cname, row.names=NULL))
                }
                pdo
                })
              
              
              
              # REMOVE if NA
              idx <- unlist(lapply(dato, function(x) all(is.na(x$y))))
              
              dato <- do.call(rbind, c(dato[!idx], list(make.row.names = FALSE)))
              p <- p + geom_point(data=dato,pch=21,bg="lightgrey",color=1,cex=1)
            
              }
            
            # PLOT refpts
            if(rpf) {
              rpdat <- lapply(plots, function(p) {
                rpd = NULL
                for(i in 1:length(brps)){
                # CBIND x, + refpt, iter ...
                  rpd= rbind(rpd,cbind(as(rps[[i]][, p['x']], 'data.frame')[, -2],
                      # ... y, panel
                      y=c(rps[[i]][,p['y']]), pos=unname(p['pos']),cname=cname[i]))
                  
                }    
                rpd
              })
              rpdat <- do.call(rbind, c(rpdat, list(make.row.names = FALSE)))
              
              # reorganize
              refpts = unique(c(refpts[grep("Blim",refpts)],refpts[grep("F",refpts)],refpts[-grep("F",refpts)]))
              rpdat$refpt = factor(rpdat$refpt,levels=refpts)
              
              
              # CALCULATE ymin per panel
              rpdat$ymin <- ave(rpdat$y, rpdat$pos, FUN=function(x) pmin(min(x), 0))
              
              # SET shapes and colors
              #if(missing(shapes))
              #  shapes <- rep(21:25,10)[1:length(brps)] # c(rep(21,length(refpts)),rep(22,length(refpts)))
              if(missing(colours))
                colours <- rev(ss3col(length(refpts)))
                
                
                
              # ADD rps points
              p <- p + geom_point(data=rpdat, size=2.5,
                                  aes(x=data, y=y, group=refpt, fill=refpt,shape=cname),alpha=0.6,pch=c(21)) +
                scale_fill_manual(values=colours) 
              
              # ADD refpts labels and text
              if(length(refpts) > 0 & is.character(refpts)){
                rpdat <- rpdat[rpdat$refpt %in% refpts,]
                
                # CALCULATE limits of lines
                rpdat$yend <- rpdat$y * 0.99
                rpdat$ymax <- ave(rpdat$y, rpdat$pos, FUN=max)
                rpdat$ystart <- rpdat$ymin + (rpdat$ymax * 0.05)
                rpdat$ymin <- 0
                rpdat$ystart <- 0
                
                
              if(dashed==TRUE)   
                p<-p+geom_segment(data=rpdat, aes_(x=~data, y=~ystart, xend=~data, yend=~yend),linetype="dashed",color="darkgrey") 
              
              }    
                
              
            }
            return(p)
          }
 # }}}




#{{{
#' plotFsim
#
#' Plots stochastic stock dynamics against refpts for constant Fsim()
#'
#' @param object output object from Fsim() 
#' @param worms option to show individual iterations
#' @param thinning thinning rate of iterations shows, e.g. 10 shows every 10th
#' @param probs determine credibility intervals, default 80th, 90th percentiles   
#' @param plotrefs if TRUE reference points are plotted 
#' @param colour color of CIs
#' @param yrs.eval last years to be used evaluation period, default half nyears
#' @param ncol number of plot panel columns
#' @param label.size size of reference points
#' @return ggplot  
#' @export
 
plotFsim <- function(object,worms=TRUE,thinning = 10,probs=c(0.05,0.2,0.50,0.8,0.95),plotrefs=TRUE,
                     colour="dodgerblue",ncol=2,label.size=3,yrs.eval=NULL){
stock = object$stock
brp = object$brp

nyears=dims(stock)$maxyear-dims(stock)$minyear+1
if(is.null(yrs.eval)) yrs.eval=ceiling(nyears/2)
iters = dims(stock)$iter
if(!worms){  
p <- ggplotFL::plot(stock,metrics=list(Rec=rec,SSB=ssb,Landings=landings,F=fbar))+theme_bw()+xlab("Year")+theme(legend.position = "none")
} else {
p <- ggplotFL::plot(stock,metrics=list(Rec=rec,SSB=ssb,Landings=landings,F=fbar),iter=seq(1,iters,thinning))+scale_color_manual(values=c(grey.colors(length(seq(1,iters,thinning)))))+
  theme_bw()+xlab("Year")+theme(legend.position = "none")
}  

p = p +ggplotFL::geom_flquantiles(fill=colour, probs=probs[c(1,3,5)], alpha=0.4) +
  ggplotFL::geom_flquantiles(fill=colour, probs=probs[c(2,3,4)], alpha=0.5)+
  facet_wrap(~qname, scales="free",ncol=ncol)

ref = rownames(object$params)[1]
Fs = FLPar(Fbrp = refpts(brp)[ref,"harvest"],Flim = refpts(brp)["Blim","harvest"])
rownames(Fs)[1] = ref

if(plotrefs){
  
  ggp =ggplotFL::geom_flpar(data=FLPars(SSB  =FLPar(Btrg=refpts(brp)[ref,"ssb"],Blim=refpts(brp)["Blim","ssb"],B0=refpts(brp)["virgin","ssb"]),
                                   F    =Fs,
                                   Landings    =FLPar(Yeq = refpts(brp)[ref,"yield"]),
                                   Rec=FLPar(R0=refpts(brp)["virgin","rec"])),
                       x=c(0.2*nyears),colour = c(c("darkgreen","red","blue"),c("darkgreen","red"),c("darkgreen","blue")))
  ggp[[2]]$aes_params$size=label.size
  
  
p = p + ggp +
  geom_vline(xintercept = nyears-yrs.eval+0.5,linetype="dashed")
}  

return(p)
}
# }}}

#{{{
#' plotAdvice
#
#' Plots stochastic stock dynamics against refpts for constant Fsim()
#'
#' @param object output object from Fsim() 
#' @param plotrefs if TRUE reference points are plotted 
#' @param colour color of CIs
#' @param ncol number of plot panel columns
#' @param label.size size of refpts labels 
#' @return ggplot  
#' @export

plotAdvice <- function(stock,brp,plotrefs=TRUE,ncol=2,label.size=3){

  
  rp = refpts(brp)
  pr = FALSE
  if(model(brp)=="rec ~ a") pr = TRUE
  
  if(pr){
  b0 = c(rp["B0",c("ssb")])
  r0 = exp(mean(log(rec(stock))))
  ref = rownames(rp)[grep("F",rownames(rp))[1]]
  yref = c(rp[ref,"yield"])
  fref = c(rp[ref,"harvest"])
  fls = FLQuants("Recruitment"=rec(stock)/r0,
                 "Spawning Ratio Potential"=(ssb(stock)/r0)/b0,"Yield per Recruit"=(landings(stock)/r0)/yref,"F"=fbar(stock))
  
  names(fls)[names(fls)=="F"] =  paste0("F(",paste0(range(stock, c("minfbar", "maxfbar")),
         collapse="-"),")")
  
  p = ggplotFL::plot(fls)+ ylim(c(0, NA))+ theme_bw()+xlab("Year")+ facet_wrap(~qname, scales="free",ncol=ncol)  
  
  }
  if(!pr){
    fls = FLQuants("Recruitment"=rec(stock),
                   "SSB"=ssb(stock),"Landings"=landings(stock),"F"=fbar(stock))
    names(fls)[names(fls)=="F"] =  paste0("F(",paste0(range(stock, c("minfbar", "maxfbar")),
                                                      collapse="-"),")")
    p = ggplotFL::plot(fls)+ ylim(c(0, NA))+ theme_bw()+xlab("Year")+ facet_wrap(~qname, scales="free",ncol=ncol)  
    
  }
     
  
  xy =quantile(dims(stock)$minyear:dims(stock)$maxyear,c(0.2,0.45,0.75,0.6,0.3))
  
  
  Fs = FLPar(an(refpts(brp)[rownames(refpts(brp))[grep("F",rownames(refpts(brp)))],"harvest"]),
             params=rownames(refpts(brp))[grep("F",rownames(refpts(brp)))])
  Fs= rbind(Fs,FLPar(Flim=refpts(brp)["Blim","harvest"]))
  nf = length(Fs)-1
  Bs = FLPar(an(refpts(brp)[rownames(refpts(brp))[grep("F",rownames(refpts(brp)))],"ssb"]),
             params=paste0(rownames(refpts(brp))[grep("F",rownames(refpts(brp)))],""))
  Bs= rbind(Bs,FLPar(Flim=refpts(brp)["Blim","ssb"],B0=refpts(brp)["virgin","ssb"]))
  rownames(Bs) = gsub("F","B",rownames(Bs))
  if(any(rownames(Bs)%in%"B0.1")) rownames(Bs)[rownames(Bs)%in%"B0.1"] = "Bf0.1" 
  
  Ys = FLPar(an(refpts(brp)[rownames(refpts(brp))[grep("F",rownames(refpts(brp)))],"yield"]),
             params=paste0(rownames(refpts(brp))[grep("F",rownames(refpts(brp)))],""))
  rownames(Ys) = gsub("F","Y",rownames(Ys))
  if(any(rownames(Bs)%in%"B0.1")) rownames(Ys)[rownames(Ys)%in%"Y0.1"] = "Yf0.1" 
  
  if(pr){
    Bs = Bs/b0
    Ys = Ys/yref    
    
  }
  
  xy = xy[1:nf]
  
  if(plotrefs){
  if(!pr){  
    
    fps =FLPars(SSB  =Bs,
           F    =Fs,
           Catch    =Ys,
           Rec=FLPar(R0=refpts(brp)["virgin","rec"]))
    fps@names = names(fls)[c(2,4,3,1)] 
    ggp = ggplotFL::geom_flpar(data=fps,
                               x=c(c(xy,xy[1],xy[1]),c(xy,xy[1]),c(xy),xy[1]),colour = c(c(rep("darkgreen",nf),"red","blue"),c(rep("darkgreen",nf),"red"),c(rep("darkgreen",nf),"blue")))
    ggp[[2]]$aes_params$size=label.size
    p = p +ggp 
  }
  if(pr){  
  
  rownames(Bs)[grep("B0",rownames(Bs))] = "SPR0"  
  fps = FLPars(SSB  =Bs,
               F    =Fs,
               Catch    =Ys,
               Rec=FLPar(R0=refpts(brp)["virgin","rec"]))    
    
  
  fps@names = names(fls)[c(2,4,3,1)] 
  ggp = ggplotFL::geom_flpar(data=fps,
                              x=c(c(xy,xy[1],xy[1]),c(xy,xy[1]),c(xy),xy[1]),colour = c(c(rep("darkgreen",nf),"red","blue"),c(rep("darkgreen",nf),"red"),c(rep("darkgreen",nf),"blue")))
  ggp[[2]]$aes_params$size=label.size
  p = p +ggp 
 }  
    
    
    }  
  
  return(p)
}
# }}}

#{{{
#' plotWKREF
#
#' Plots the new proposed ICES advice rule
#'
#' @param ftrg Target F = min(Fbrp,Fp0.5) 
#' @param btrg Biomass target corresponding to Fbrp
#' @param blim biomass limit
#' @param btrigger biomass trigger below which F is linearly reduced   
#' @param bthresh biomass threshold beyond which biomass is classified sustainable
#' @param bclose biomass that invokes fishing closure 
#' @param fmin minimum allowable (bycatch) fishing mortality under closure 
#' @param obs obtion to show observation with input class `FLStock`
#' @param kobe add kobe colour-coding
#' @param alpha transparency of shading
#' @param xmax multiplier for upper default xlim
#' @param ymax multiplier for upper default ylim
#' @param xlab option customize xlab
#' @param ylab option customize ylab
#' @param rel option to denote x,y labs as relative B/Btrg and F/Ftrg
#' @param expand option to expand the plot area to border - default TRUE
#' @param labels annotate reference point labels
#' @param critical option to highlight critical zone below blim
#' @return ggplot  
#' @export
#' @examples
#' plotWKREF()
#' # Close fishery at Blim and adjust axis labels to relative
#' plotWKREF(blim=0.2,bclose=0.2,rel=TRUE)
#' # Close fishery at Blim, but allow fmin (e.g. bycatch)
#' plotWKREF(blim=0.2,bclose=0.2,fmin=0.1,rel=TRUE)
#' # Change Btrigger above Btrg
#' plotWKREF(blim=0.2,bclose=0.2,fmin=0.1,btrigger=1.1,rel=TRUE)
#' # Plot stock data
#' data(ple4)
#' plotWKREF(ftrg=0.25,btrg=8e+05,btrigger = 0.9*8e+05, blim=2e5,bclose=3e5,fmin=0.03,obs=ple4)



plotWKREF <- function(ftrg = 1,btrg=1,blim=0.2,btrigger=0.9*btrg,bthresh=0.8*btrg,bclose=0,fmin=0, obs="missing", kobe=TRUE,
                      alpha=1,xmax=1.3,ymax=1.5,ylab="missing",xlab="missing",rel=FALSE,expand=TRUE,labels=TRUE,critical=kobe) {
  #Define axis
  metric="ssb"
  output="fbar"
  
  # SET args
  xlim <- btrg * xmax
  ylim <- ftrg * ymax
  
  # GET observations
  if(!missing(obs)) {
    if(class(obs)=="FLStock")
      obs <- model.frame(metrics(obs, list(met=get(metric), out=get(output))))
    
    if(class(obs)=="FLQuants"){
      obs=obs[1:2]
      names(obs) = c("met","out")
    }
    
    xlim <- max(c(obs$met,btrg*1.5)) * 1.05
    ylim <- max(c(obs$out,ftrg*1.3)) * 1.05
  }
  
  # SET met values
  met <- seq(0, xlim, length=200)
  
  # BELOW lim
  out <- ifelse(met <= bclose, fmin,
                # Revised
                ifelse(met < btrigger,
                       pmax(c(((ftrg-fmin)/(btrigger-bclose))*(met - bclose) +fmin,fmin)),
                       # ABOVE btrigger
                       ftrg))
  # DATA
  dat <- data.frame(met=met, out=out)
  
  if(missing(xlab)){ 
    xlab="SSB"
    if(rel) xlab = expression(B/B[trg])
  }
  if(missing(ylab)){ 
    ylab="Fishing Mortality"
   if(rel) ylab = expression(F/F[trg])
  }
  p <- ggplot(dat, aes(x=met, y=out))+theme_bw()  +
    xlab(xlab) + ylab(ylab)  
  if(expand){
    p <- p + scale_x_continuous(expand = c(0, 0), limits = c(0, xlim)) + 
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim))+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.line = element_blank())
    
  }
  
  # PREDICT for threshold
  ythresh <- ifelse(bthresh < btrigger,
                    pmax(c(((ftrg-fmin)/(btrigger-bclose))*(bthresh - bclose) +fmin,fmin)),
                    ftrg)
  ytarget <- ifelse(btrg < btrigger,
                    pmax(c(((ftrg-fmin)/(btrigger-bclose))*(btrg - bclose) +fmin,fmin)),
                    ftrg)
  if(kobe) {
    ylwmin = ifelse(fmin>0,0,bclose)
    
    if(btrigger>=bthresh){
    yell = geom_polygon(data=data.frame(x=c(ylwmin,ylwmin,bclose, bthresh, bthresh,bthresh ,bclose,ylwmin),
                                          y=c(fmin,0, 0, 0,0,ythresh,fmin,fmin)),
                          aes(x=x, y=y), fill="yellow", alpha=alpha)
    } else {
     yell= geom_polygon(data=data.frame(x=c(ylwmin,ylwmin,bclose, bthresh, bthresh,btrigger ,bclose,ylwmin),
                                          y=c(fmin,  0,    0,      0,       ythresh,ythresh,fmin,fmin)),
                          aes(x=x, y=y), fill="yellow", alpha=alpha)  
      
    }
    # RED
    if(btrigger>=bthresh){
      rdd = geom_polygon(data=data.frame(
        x=c(0, bclose,bthresh, bthresh, bthresh, bthresh, 0, 0),
        y=c(fmin, fmin, ythresh, ftrg,ftrg ,ylim, ylim, fmin)),
        aes(x=x, y=y), fill="red", alpha=alpha)
    } else {
      rdd = geom_polygon(data=data.frame(
        x=c(0, bclose,btrigger, bthresh, bthresh, bthresh, 0, 0),
        y=c(fmin, fmin, ythresh, ftrg,ftrg ,ylim, ylim, fmin)),
        aes(x=x, y=y), fill="red", alpha=alpha)
    }
    p <- p+yell+rdd+
      # GREEN
      geom_polygon(data=data.frame(
        x=c(bthresh, xlim, xlim,btrigger, bthresh, bthresh),
        y=c(0, 0, ftrg,ytarget, ytarget, ythresh)),
        aes(x=x, y=y), fill="green", alpha=alpha) +
      
      # Orange
      geom_polygon(data=data.frame(
        x=c(bthresh, btrigger, xlim, xlim, bthresh,bthresh),
        y=c(ythresh , ftrg, ftrg, ylim, ylim, ythresh)),
        aes(x=x, y=y), fill="orange", alpha=alpha)+geom_line()

  
  if(critical){
    p<-p+geom_polygon(data=data.frame(
      x=c(0, blim, blim,  0,0),
      y=c(0, 0, ylim, ylim,0)),
      aes(x=x, y=y), fill="red4", alpha=0.7)  
  }   
  }
  
  # Btrg
  p <- p+geom_segment(aes(x=0, xend=btrigger * 1.25, y=ftrg, yend=ftrg), linetype=2) +
    # Btrigger
    geom_segment(aes(x=btrigger, xend=btrigger, y=0, yend=ftrg), linetype=2)+
    annotate("text", x=btrigger*1.02, y=ftrg*1.03, label=paste0("B[trigger]"), 
             vjust="bottom",parse=TRUE) +
    # Btrg
    geom_segment(aes(x=btrg, xend=btrg, y=0, yend=ytarget), linetype=1,color="blue",cex=0.9)+
    annotate("text", x=btrg*1.02, y=ftrg*0.5, label=paste0("B[trg]"), 
             hjust="left",parse=TRUE)+  
    # Btresh
    geom_segment(aes(x=bthresh, xend=bthresh, y=0, yend=ythresh), linetype=1)+
    annotate("text", x=bthresh*0.98, y=ftrg*0.4, label=paste0("B[thresh]"), 
             hjust="right",parse=TRUE) +
    
    # Blim
    geom_segment(aes(x=blim, xend=blim, y=0, yend=ftrg), linetype=1,cex=0.9)+
    annotate("text", x=blim*0.98, y=ftrg*1.03, label=paste0("B[lim]"), 
             vjust="bottom",hjust="left",parse=TRUE)+  
    
    geom_line()
    if(labels){
    p <- p+annotate("text", x=xlim*0.8, y=ftrg + ylim / 30*1.08, label="F[trg]",parse=TRUE, hjust="left") +
    annotate("text", x=btrigger*1.02, y=ftrg*1.03, label=paste0("B[trigger]"), 
               vjust="bottom",parse=TRUE) +
    # Btrg
    annotate("text", x=btrg*1.02, y=ftrg*0.5, label=paste0("B[trg]"), 
               hjust="left",parse=TRUE)+  
    # Btresh
    annotate("text", x=bthresh*0.98, y=ftrg*0.4, label=paste0("B[thresh]"), 
               hjust="right",parse=TRUE) +
    # Blim
    annotate("text", x=blim*0.98, y=ftrg*1.03, label=paste0("B[lim]"), 
               vjust="bottom",hjust="left",parse=TRUE)  
    }  

  if(!missing(obs)) {
    p=p+geom_path(data=obs,cex=0.2,col="blue",linetype=1,alpha=0.5)+
      geom_point(data=obs,cex=c(rep(1.5,nrow(obs)-1),2.5),pch=21,fill=c(rep("grey",nrow(obs)-1),"blue"))
    
  }
  return(p)
}  

#{{{
#' plotICES
#
#' Plots the current ICES advice rule
#'
#' @param fmsy Target F = min(Fbrp,Fp0.5) 
#' @param btrigger is also btrg=bthresh 
#' @param blim biomass limit
#' @param bclose biomass that invokes fishing closure 
#' @param fmin minimum allowable (bycatch) fishing mortality under closure 
#' @param obs obtion to show observation with input class `FLStock`
#' @param kobe add kobe colour-coding
#' @param alpha transparency of shading
#' @param xmax multiplier for upper default xlim
#' @param ymax multiplier for upper default ylim
#' @param xlab option customize xlab
#' @param ylab option customize ylab
#' @param rel option to denote x,y labs as relative B/Btrg and F/Ftrg
#' @param expand option to expand the plot area to border - default TRUE
#' @param labels annotate reference point labels
#' @param critical option to highlight critical zone below blim
#' @return ggplot  
#' @export
#' @examples
#' plotICES()
#' # Close fishery at Blim and adjust axis labels to relative
#' plotICES(btrigger=1.4,bclose=1/1.4,rel=TRUE)
#' # Close fishery at Blim, but allow fmin (e.g. bycatch)
#' plotICES(btrigger=1.4,bclose=1/1.4,fmin=0.05,rel=TRUE)
#' # Plot stock data
#' data(ple4)
#' plotICES(fmsy=0.25,blim=2e5,,btrigger = 1.4*2e5,obs=ple4)

plotICES <- function(fmsy = 1,blim=1/1.4,btrigger=1,bclose=0,fmin=0, obs="missing", kobe=TRUE,
                      alpha=1,xmax=3,ymax=1.5,ylab="missing",xlab="missing",rel=FALSE,expand=TRUE,labels=TRUE,critical=FALSE) {
  # Define 
  btrg=bthresh=btrigger
  ftrg=fmsy
  #Define axis
  metric="ssb"
  output="fbar"
  
  # SET args
  xlim <- btrg * xmax
  ylim <- ftrg * ymax
  
  # GET observations
  if(!missing(obs)) {
    if(class(obs)=="FLStock")
         obs <- model.frame(metrics(obs, list(met=get(metric), out=get(output))))
    
    if(class(obs)=="FLQuants"){
      obs=obs[1:2]
      names(obs) = c("met","out")
    }
    
    xlim <- max(c(obs$met,btrg*1.5)) * 1.05
    ylim <- max(c(obs$out,ftrg*1.3)) * 1.05
  }
  
  # SET met values
  met <- seq(0, xlim, length=200)
  
  # BELOW lim
  out <- ifelse(met <= bclose, fmin,
                # Revised
                ifelse(met < btrigger,
                       pmax(c(((ftrg-fmin)/(btrigger-bclose))*(met - bclose) +fmin,fmin)),
                       # ABOVE btrigger
                       ftrg))
  # DATA
  dat <- data.frame(met=met, out=out)
  if(missing(xlab)){ 
    xlab="SSB"
    if(rel) xlab = expression(B/B[trigger])
  }
  if(missing(ylab)){ 
    ylab="Fishing Mortality"
    if(rel) ylab = expression(F/F[MSY])
  }
  p <- ggplot(dat, aes(x=met, y=out))+theme_bw()  +
    xlab(xlab) + ylab(ylab)  
  if(expand){
    p <- p + scale_x_continuous(expand = c(0, 0), limits = c(0, xlim)) + 
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim))+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.line = element_blank())
    
  }
  
  # PREDICT for threshold
  ythresh <- ifelse(bthresh < btrigger,
                    pmax(c(((ftrg-fmin)/(btrigger-bclose))*(bthresh - bclose) +fmin,fmin)),
                    ftrg)
  ytarget <- ifelse(btrg < btrigger,
                    pmax(c(((ftrg-fmin)/(btrigger-bclose))*(btrg - bclose) +fmin,fmin)),
                    ftrg)
  if(kobe) {
    ylwmin = ifelse(fmin>0,0,bclose)
  
      p = p + 
      # Red
      geom_polygon(data=data.frame(
      x=c(0,blim,blim, 0),
      y=c(0,0,ylim, ylim)),
      aes(x=x, y=y), fill="red", alpha=alpha)+
      # Yellow
      geom_polygon(data=data.frame(x=c(blim,btrigger,btrigger,blim),
                                   y=c(0, 0, ylim,ylim)),
                   aes(x=x, y=y), fill="yellow", alpha=alpha)+
      # GREEN
      geom_polygon(data=data.frame(
        x=c(btrigger, xlim, xlim,btrigger),
        y=c(0, 0, ylim,ylim)),
        aes(x=x, y=y), fill="green", alpha=alpha)
    
    
    if(critical){
      p<-p+geom_polygon(data=data.frame(
        x=c(0, blim, blim,  0,0),
        y=c(0, 0, ylim, ylim,0)),
        aes(x=x, y=y), fill="red4", alpha=0.7)  
    }   
  }
  
  # Btrg
  p <- p+geom_segment(aes(x=0, xend=btrigger * 1.25, y=ftrg, yend=ftrg), linetype=2) +
    # Btrigger
    geom_segment(aes(x=btrigger, xend=btrigger, y=0, yend=ftrg), linetype=2)+
    # Blim
    geom_segment(aes(x=blim, xend=blim, y=0, yend=ftrg), linetype=1)+
    
    geom_line()
  if(labels){
    p <- p+annotate("text", x=btrigger*1.6, y=ftrg + ylim / 30*1.09, label="F[MSY]",parse=TRUE, hjust="left") +
      # Btrg
      annotate("text", x=btrigger*1.02, y=ftrg*0.5, label=paste0("B[trigger]"), 
               hjust="left",parse=TRUE)+  
      # Blim
      annotate("text", x=blim*0.98, y=ftrg*1.03, label=paste0("B[lim]"), 
               vjust="bottom",hjust="right",parse=TRUE)  
  }  
  
  if(!missing(obs)) {
    p=p+geom_path(data=obs,cex=0.2,col="blue",linetype=1,alpha=0.5)+
      geom_point(data=obs,cex=c(rep(1.5,nrow(obs)-1),2.5),pch=21,fill=c(rep("grey",nrow(obs)-1),"blue"))
  }
  return(p)
}  



#{{{
#' plotMAjuro
#
#' Plots the new proposed ICES advice rule
#'
#' @param ftrg Target F = min(Fbrp,Fp0.5) 
#' @param btrg Biomass target corresponding to Fbrp
#' @param blim biomass limit
#' @param btrigger biomass trigger below which F is linearly reduced   
#' @param bthresh biomass threshold beyond which biomass is classified sustainable
#' @param bclose biomass that invokes fishing closure 
#' @param fmin minimum allowable (bycatch) fishing mortality under closure 
#' @param obs obtion to show observation with input class `FLStock`
#' @param kobe add kobe colour-coding
#' @param alpha transparency of shading
#' @param xmax multiplier for upper default xlim
#' @param ymax multiplier for upper default ylim
#' @param xlab option customize xlab
#' @param ylab option customize ylab
#' @param rel option to denote x,y labs as relative B/Btrg and F/Ftrg
#' @param expand option to expand the plot area to border - default TRUE
#' @param labels annotate reference point labels
#' @param critical option to highlight critical zone below blim
#' @return ggplot  
#' @export
#' @examples
#' plotWKREF()
#' # Close fishery at Blim and adjust axis labels to relative
#' plotWKREF(blim=0.2,bclose=0.2,rel=TRUE)
#' # Close fishery at Blim, but allow fmin (e.g. bycatch)
#' plotWKREF(blim=0.2,bclose=0.2,fmin=0.1,rel=TRUE)
#' # Change Btrigger above Btrg
#' plotWKREF(blim=0.2,bclose=0.2,fmin=0.1,btrigger=1.1,rel=TRUE)
#' # Plot stock data
#' data(ple4)
#' plotWKREF(ftrg=0.25,btrg=8e+05,btrigger = 0.9*8e+05, blim=2e5,bclose=3e5,fmin=0.03,obs=ple4)



plotMajuro<- function(ftrg = 1,fthresh=1.1,btrg=1,blim=0.1,btrigger=0.8*btrg,bthresh=0.5*btrg,bclose=0,fmin=0, obs="missing", kobe=TRUE,
                      alpha=1,xmax=1.5,ymax=1.5,ylab="missing",xlab="missing",rel=FALSE,expand=TRUE,labels=TRUE,critical=kobe) {
  #Define axis
  metric="ssb"
  output="fbar"
  
  # SET args
  xlim <- btrg * xmax
  ylim <- ftrg * ymax
  
  # GET observations
  if(!missing(obs)) {
    if(class(obs)=="FLStock")
      obs <- model.frame(metrics(obs, list(met=get(metric), out=get(output))))
    
    if(class(obs)=="FLQuants"){
      obs=obs[1:2]
      names(obs) = c("met","out")
    }
    
    xlim <- max(c(obs$met,btrg*1.5)) * 1.05
    ylim <- max(c(obs$out,ftrg*1.3)) * 1.05
  }
  
  # SET met values
  met <- seq(0, xlim, length=200)
  
  # BELOW lim
  out <- ifelse(met <= bclose, fmin,
                # Revised
                ifelse(met < btrigger,
                       pmax(c(((ftrg-fmin)/(btrigger-bclose))*(met - bclose) +fmin,fmin)),
                       # ABOVE btrigger
                       ftrg))
  # DATA
  dat <- data.frame(met=met, out=out)
  
  
  
  if(missing(xlab)){ 
    xlab="SSB"
    if(rel) xlab = expression(B/B[trg])
  }
  if(missing(ylab)){ 
    ylab="Fishing Mortality"
    if(rel) ylab = expression(F/F[trg])
  }
  p <- ggplot(dat, aes(x=met, y=out))+theme_bw()  +
    xlab(xlab) + ylab(ylab)  
  if(expand){
    p <- p + scale_x_continuous(expand = c(0, 0), limits = c(0, xlim)) + 
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim))+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), panel.background = element_blank(), 
            axis.line = element_blank())
    
  }
  
  # PREDICT for threshold
  yfthresh <- ifelse(bthresh < btrigger,
                    pmax(c(((fthresh-fmin)/(btrigger-bclose))*(bthresh - bclose) +fmin,fmin)),
                    fthresh)
  
  
  ythresh <- ifelse(bthresh < btrigger,
                    pmax(c(((ftrg-fmin)/(btrigger-bclose))*(bthresh - bclose) +fmin,fmin)),
                    ftrg)
  ytarget <- ifelse(btrg < btrigger,
                    pmax(c(((ftrg-fmin)/(btrigger-bclose))*(btrg - bclose) +fmin,fmin)),
                    ftrg)
  if(kobe) {
    ylwmin = ifelse(fmin>0,0,bclose)
    
    
    p <- p+
      # GREEN
      geom_polygon(data=data.frame(
        x=c(bthresh, xlim, xlim,bthresh),
        y=c(0, 0, fthresh,fthresh)),
        aes(x=x, y=y), fill="green", alpha=alpha) +
      # Orange
      geom_polygon(data=data.frame(
        x=c(bthresh, xlim, xlim, bthresh),
        y=c(fthresh, fthresh, ylim, ylim)),
        aes(x=x, y=y), fill="orange", alpha=alpha)+
      
      # Yellow
      geom_polygon(data=data.frame(
        x=c(blim, bthresh, bthresh, blim),
        y=c(0, 0, fthresh,fthresh)),
        aes(x=x, y=y), fill="yellow", alpha=alpha)+
      # red
      geom_polygon(data=data.frame(
        x=c(0,blim, blim, bthresh, bthresh,0,0),
        y=c(0, 0, fthresh,fthresh,ylim,ylim,0)),
        aes(x=x, y=y), fill="red", alpha=alpha)+
    
    
    geom_line()
    
    
    if(critical){
      p<-p+geom_polygon(data=data.frame(
        x=c(0, blim, blim,  0,0),
        y=c(0, 0, ylim, ylim,0)),
        aes(x=x, y=y), fill="red4", alpha=0.7)  
    }   
  }
  
  # Btrg
  p <- p+geom_segment(aes(x=0, xend=btrigger * 1.25, y=ftrg, yend=ftrg), linetype=2) +
    # Btrigger
    geom_segment(aes(x=btrigger, xend=btrigger, y=0, yend=ftrg), linetype=2)+
   
    # Btrg
    #geom_segment(aes(x=btrg, xend=btrg, y=0, yend=ytarget), linetype=1,color="blue",cex=0.9)+
    
    # Btresh
    geom_segment(aes(x=bthresh, xend=bthresh, y=0, yend=fthresh), linetype=1)+
  
    
    # Blim
    geom_segment(aes(x=blim, xend=blim, y=0, yend=ylim), linetype=1,cex=0.9)+
    
    geom_segment(aes(x=blim, xend=xlim, y=fthresh, yend=fthresh),cex=0.2)+
    
    geom_line()
  
  
  if(labels){
    p <- p+annotate("text", x=xlim*0.7, y=ftrg + ylim / 30*1.08, label="F[advive]/F[MSYpa]",parse=TRUE, hjust="left") +
      annotate("text", x=btrigger*1.02, y=ftrg*1.03, label=paste0("B[trigger]"), 
               vjust="bottom",parse=TRUE) +
      # Btrg
      #annotate("text", x=btrg*1.02, y=ftrg*0.5, label=paste0("B[MSY]"), 
      #         hjust="left",parse=TRUE)+  
      # Btresh
      annotate("text", x=bthresh, y=fthresh*1.01, label=paste0("B[safe]"), 
               vjust="bottom",parse=TRUE) +
      # Fthresh
      annotate("text", x=btrg, y=fthresh*1.01, label=paste0("F[MSY]"), 
               vjust="bottom",parse=TRUE) +
      # Blim
      annotate("text", x=blim*1.05, y=ftrg*0.7, label=paste0("B[lim]"), 
               vjust="bottom",hjust="left",parse=TRUE)  
  }  
  
  if(!missing(obs)) {
    p=p+geom_path(data=obs,cex=0.2,col="blue",linetype=1,alpha=0.5)+
      geom_point(data=obs,cex=c(rep(1.5,nrow(obs)-1),2.5),pch=21,fill=c(rep("grey",nrow(obs)-1),"blue"))
    
  }
  return(p)
}  

