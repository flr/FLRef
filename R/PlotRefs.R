
#{{{
#' ploteq()
#
#' Modification of method plot(`FLBRP`) to plot equilibrium output of computeFbrp()  
#'
#' @param brp output object from computeFbrp of class FLBRP  
#' @param refpts Reference points, defaults are computed refpts from computeFbrp()  
#' \itemize{
#'   \item Fbrp  
#'   \item Blim 
#'   \item B0  
#'   \item Btri  
#' }    
#' @param obs Should observations be plotted? Defaults to `FALSE`.
#' @param labels plot refpts label, default to `FALSE`
#' @param colours refpts colours, default is designed for computeFbrp() output
#' @param shapes refpts symbols, default is designed for computeFbrp() output
#' @param panels plot panel option 1:4 
#' @param ncol number of plot panel columns
#' @return ggplot  
#' @export

ploteq <- function(brp, refpts=c("Fbrp","Btri","Blim","B0"), obs=FALSE, labels=FALSE,
                   shapes="missing", colours="missing", panels=NULL, ncol=2) {
            
            x = brp
            # EXTRACT metrics
            df <- model.frame(metrics(x,
                                      list(ssb=ssb, harvest=fbar, rec=rec, yield=landings, profit=profit)),
                              drop=FALSE)
            # refpts
            drps <- dimnames(refpts(x))$refpt
            rps <- refpts(x)[drps %in% refpts,]
            
            # estimated?
            rpf <- !all(is.na(rps))
            
            # SUBSET df IF rpf
            if(rpf && "crash" %in% dimnames(rps)$refpt)
              df <- df[df$harvest <= c(rps['crash', 'harvest']),]
            
            # NO economics
            plots <- list(
              P1=c(x="harvest", y="ssb", panel="SSB ~ F", pos=1),
              P2=c(x="ssb", y="rec", panel="Recruitment ~ SSB", pos=2),
              P3=c(x="harvest", y="yield", panel="Yield ~ F", pos=3),
              P4=c(x="ssb", y="yield", panel="Yield ~ SSB", pos=4))
            
            # WITH economics
            if(!all(is.na(rps[, 'profit']))) {
              plots <- c(plots, list(
                P5=c(x="harvest", y="profit", panel="Equilibrium Profit v. F", pos=5),
                P6=c(x="ssb", y="profit", panel="Equilibrium Profit v. SSB", pos=6)))
            } else {
              dms <- dimnames(rps)
              rps <- rps[!dms$refpt %in% "mey",
                         !dms$quant %in% c("revenue", "cost", "profit")]
            }
            
            # SUBSET panels if not NULL
            if(!is.null(panels))
              plots <- plots[panels]
            
            # APPLY over plots to extract x, y and panel for each element
            dat <- lapply(plots, function(p) {
              data.frame(x=df[,p['x']], y=df[,p['y']], iter=df[,'iter'],
                         panel=p['panel'], pos=p['pos'], row.names=NULL)
            })
            
            # RBIND into single df
            dat <- do.call(rbind, c(dat, list(make.row.names = FALSE)))
            
            # Limit y to 0 in panels 1:4
            dat[dat$pos %in% 1:4 & dat$y<0, "y"] <- 0
            
            # CREATE facet labels vector
            facl <- setNames(unique(dat$panel), nm=unique(dat$pos))
            
            # PLOT
            p <- ggplot(dat, aes_(x=~x, y=~y, group=~iter)) + geom_line() + theme_bw()+
              facet_wrap(~pos, scales="free", ncol=ncol, labeller=labeller(pos=facl)) +
              xlab("") + ylab("") +
              scale_x_continuous(expand = expansion(mult = c(0, .05)),labels=human_numbers, limits=c(0, NA))+
              scale_y_continuous(expand = expansion(mult = c(0, .1)))+
              theme(legend.title = element_blank())
            
            # PLOT observations
            if(obs) {
              
              dfo <- model.frame(metrics(x,
                                         list(ssb=ssb.obs, harvest=fbar.obs, rec=rec.obs, yield=landings.obs,
                                              profit=profit.obs)), drop=FALSE)
              
              # APPLY over plots to extract x, y and panel for each element
              dato <- lapply(plots, function(p)
                data.frame(x=dfo[,p['x']], y=dfo[,p['y']], iter=dfo[,'iter'],
                           pos=p['pos'], row.names=NULL))
              
              # REMOVE if NA
              idx <- unlist(lapply(dato, function(x) all(is.na(x$y))))
              
              dato <- do.call(rbind, c(dato[!idx], list(make.row.names = FALSE)))
              
              p <- p + geom_point(data=dato,pch=21,bg="grey")
            }
            
            # PLOT refpts
            if(rpf) {
              rpdat <- lapply(plots, function(p) {
                # CBIND x, + refpt, iter ...
                cbind(as(rps[, p['x']], 'data.frame')[, -2],
                      # ... y, panel
                      y=c(rps[,p['y']]), pos=unname(p['pos']))
              })
              rpdat <- do.call(rbind, c(rpdat, list(make.row.names = FALSE)))
              
              # CALCULATE ymin per panel
              rpdat$ymin <- ave(rpdat$y, rpdat$pos, FUN=function(x) pmin(min(x), 0))
              
              # SET shapes and colors
              if(missing(shapes))
                shapes <- rep(21,length(refpts))
              if(missing(colours))
                colours <- c(rainbow(3)[c(2,1,3)],"orange", rep(c("#e69f00",
                                                                     "#56b4e9", "#f0e442", "#0072b2", "#d55e00", "#cc79a7"), 4))
              
              # ADD rps points
              p <- p + geom_point(data=rpdat, size=2.5,
                                  aes_(x=~data, y=~y, group=~refpt, fill=~refpt, shape=~refpt),alpha=0.7) +
                scale_shape_manual(values=shapes) +
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
                
                # LABEL
                if(labels) {
                  p <- p + geom_text(data=rpdat,
                                     aes_(x=~data, y=~ymin, label=~refpt), angle = 90, size=3, vjust="left") +
                    # LINES
                    geom_segment(data=rpdat, aes_(x=~data, y=~ystart, xend=~data, yend=~yend),
                                 colour="grey")
                } else {
                  p<-p+geom_segment(data=rpdat, aes_(x=~data, y=~ystart, xend=~data, yend=~yend),linetype="dashed",color="darkgrey") 
                   
                }
              }
            }
            return(p)
          }
 # }}}

#{{{
#' plotFbrps()
#
#' Modification of method plot(`FLBRP`) to plot equilibrium output of computeFbrp()  
#'
#' @param brp output object from computeFbrp of class FLBRP  
#' @param refpts Reference points, defaults are computed refpts from computeFbrp()  
#' \itemize{
#'   \item all  
#'   \item sprx 
#'   \item bx  
#' }    
#' @param obs Should observations be plotted? Defaults to `FALSE`.
#' @param labels plot refpts label, default to `FALSE`
#' @param colours refpts colours, default is designed for computeFbrp() output
#' @param shapes refpts symbols, default is designed for computeFbrp() output
#' @param panels plot panel option 1:4 
#' @param ncol number of plot panel columns
#' @return ggplot  
#' @export

plotFbrps <- function(brp, proxies=c("all","sprx","bx"), obs=FALSE, labels=FALSE,
                   shapes="missing", colours="missing", panels=NULL, ncol=2) {
  
  x = brp
  # EXTRACT metrics
  df <- model.frame(metrics(x,
                            list(ssb=ssb, harvest=fbar, rec=rec, yield=landings, profit=profit)),
                    drop=FALSE)
  
  if(SRModelName(model(sr))%in%c("segregA1","segreg","geomean"))
  
    
  # refpts
  drps <- dimnames(refpts(x))$refpt
  if(proxies=="sprx") refpts = c("msy","virgin",drps[8:11])
  
  rps <- refpts(x)[drps %in% refpts,]
  
  # estimated?
  rpf <- !all(is.na(rps))
  
  # SUBSET df IF rpf
  if(rpf && "crash" %in% dimnames(rps)$refpt)
    df <- df[df$harvest <= c(rps['crash', 'harvest']),]
  
  # NO economics
  plots <- list(
    P1=c(x="harvest", y="ssb", panel="SSB ~ F", pos=1),
    P2=c(x="ssb", y="rec", panel="Recruitment ~ SSB", pos=2),
    P3=c(x="harvest", y="yield", panel="Yield ~ F", pos=3),
    P4=c(x="ssb", y="yield", panel="Yield ~ SSB", pos=4))
  
  # WITH economics
  if(!all(is.na(rps[, 'profit']))) {
    plots <- c(plots, list(
      P5=c(x="harvest", y="profit", panel="Equilibrium Profit v. F", pos=5),
      P6=c(x="ssb", y="profit", panel="Equilibrium Profit v. SSB", pos=6)))
  } else {
    dms <- dimnames(rps)
    rps <- rps[!dms$refpt %in% "mey",
               !dms$quant %in% c("revenue", "cost", "profit")]
  }
  
  # SUBSET panels if not NULL
  if(!is.null(panels))
    plots <- plots[panels]
  
  # APPLY over plots to extract x, y and panel for each element
  dat <- lapply(plots, function(p) {
    data.frame(x=df[,p['x']], y=df[,p['y']], iter=df[,'iter'],
               panel=p['panel'], pos=p['pos'], row.names=NULL)
  })
  
  # RBIND into single df
  dat <- do.call(rbind, c(dat, list(make.row.names = FALSE)))
  
  # Limit y to 0 in panels 1:4
  dat[dat$pos %in% 1:4 & dat$y<0, "y"] <- 0
  
  # CREATE facet labels vector
  facl <- setNames(unique(dat$panel), nm=unique(dat$pos))
  
  # PLOT
  p <- ggplot(dat, aes_(x=~x, y=~y, group=~iter)) + geom_line() + theme_bw()+
    facet_wrap(~pos, scales="free", ncol=ncol, labeller=labeller(pos=facl)) +
    xlab("") + ylab("") +
    scale_x_continuous(expand = expansion(mult = c(0, .05)),labels=human_numbers, limits=c(0, NA))+
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    theme(legend.title = element_blank())
  
  # PLOT observations
  if(obs) {
    
    dfo <- model.frame(metrics(x,
                               list(ssb=ssb.obs, harvest=fbar.obs, rec=rec.obs, yield=landings.obs,
                                    profit=profit.obs)), drop=FALSE)
    
    # APPLY over plots to extract x, y and panel for each element
    dato <- lapply(plots, function(p)
      data.frame(x=dfo[,p['x']], y=dfo[,p['y']], iter=dfo[,'iter'],
                 pos=p['pos'], row.names=NULL))
    
    # REMOVE if NA
    idx <- unlist(lapply(dato, function(x) all(is.na(x$y))))
    
    dato <- do.call(rbind, c(dato[!idx], list(make.row.names = FALSE)))
    
    p <- p + geom_point(data=dato,pch=21,bg="grey")
  }
  
  # PLOT refpts
  if(rpf) {
    rpdat <- lapply(plots, function(p) {
      # CBIND x, + refpt, iter ...
      cbind(as(rps[, p['x']], 'data.frame')[, -2],
            # ... y, panel
            y=c(rps[,p['y']]), pos=unname(p['pos']))
    })
    rpdat <- do.call(rbind, c(rpdat, list(make.row.names = FALSE)))
    
    # CALCULATE ymin per panel
    rpdat$ymin <- ave(rpdat$y, rpdat$pos, FUN=function(x) pmin(min(x), 0))
    
    # SET shapes and colors
    if(missing(shapes))
      shapes <- rep(21,length(refpts))
    if(missing(colours))
      colours <- c(rainbow(3)[c(2,1,3)],"orange", rep(c("#e69f00",
                                                        "#56b4e9", "#f0e442", "#0072b2", "#d55e00", "#cc79a7"), 4))
    
    # ADD rps points
    p <- p + geom_point(data=rpdat, size=2.5,
                        aes_(x=~data, y=~y, group=~refpt, fill=~refpt, shape=~refpt),alpha=0.7) +
      scale_shape_manual(values=shapes) +
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
      
      # LABEL
      if(labels) {
        p <- p + geom_text(data=rpdat,
                           aes_(x=~data, y=~ymin, label=~refpt), angle = 90, size=3, vjust="left") +
          # LINES
          geom_segment(data=rpdat, aes_(x=~data, y=~ystart, xend=~data, yend=~yend),
                       colour="grey")
      } else {
        p<-p+geom_segment(data=rpdat, aes_(x=~data, y=~ystart, xend=~data, yend=~yend),linetype="dashed",color="darkgrey") 
        
      }
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
#' @return ggplot  
#' @export
 
plotFsim <- function(object,worms=TRUE,thinning = 10,probs=c(0.05,0.2,0.50,0.8,0.95),plotrefs=TRUE,
                     colour="dodgerblue",ncol=2,yrs.eval=NULL){
stock = object$stock
brp = object$brp

nyears=dims(stock)$maxyear-dims(stock)$minyear+1
if(is.null(yrs.eval)) yrs.eval=ceiling(nyears/2)
iters = dims(stock)$iter
if(!worms){  
p <- ggplotFL::plot(stock)+theme_bw()+xlab("Year")+theme(legend.position = "none")
} else {
p <- ggplotFL::plot(stock,iter=seq(1,iters,thinning))+scale_color_manual(values=c(grey.colors(length(seq(1,iters,thinning)))))+
  theme_bw()+xlab("Year")+theme(legend.position = "none")
}  

p = p +ggplotFL::geom_flquantiles(fill=colour, probs=probs[c(1,3,5)], alpha=0.4) +
  ggplotFL::geom_flquantiles(fill=colour, probs=probs[c(2,3,4)], alpha=0.5)+
  facet_wrap(~qname, scales="free",ncol=ncol)


if(plotrefs){
p = p + ggplotFL::geom_flpar(data=FLPars(SSB  =FLPar(Btrg=refpts(brp)["Fbrp","ssb"],Blim=refpts(brp)["Blim","ssb"],B0=refpts(brp)["virgin","ssb"]),
                         F    =FLPar(Fbrp = refpts(brp)["Fbrp","harvest"],Flim = refpts(brp)["Blim","harvest"]),
                         Catch    =FLPar(Yeq = refpts(brp)["Fbrp","yield"]),
                         Rec=FLPar(R0=refpts(brp)["virgin","rec"])),
             x=c(0.2*nyears),colour = c(c("darkgreen","red","blue"),c("darkgreen","red"),c("darkgreen","blue")))+
  geom_vline(xintercept = nyears-yrs.eval+0.5,linetype="dashed")
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
   if(rel) xlab = expression(F/F[trg])
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
    p <- p+annotate("text", x=xlim*0.8, y=ftrg + ylim / 30*1.05, label="F[trg]",parse=TRUE, hjust="left") +
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
    p <- p + geom_point(data=obs)
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

plotICES <- function(fmsy = 1,blim=1/1.4,btrigger=1.4,bclose=0,fmin=0, obs="missing", kobe=TRUE,
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
    if(rel) xlab = expression(B/B[trg])
  }
  if(missing(ylab)){ 
    ylab="Fishing Mortality"
    if(rel) xlab = expression(F/F[trg])
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
    p <- p+annotate("text", x=btrigger*1.6, y=ftrg + ylim / 30*0.9, label="F[MSY]",parse=TRUE, hjust="left") +
      # Btrg
      annotate("text", x=btrigger*1.02, y=ftrg*0.5, label=paste0("B[trigger]"), 
               hjust="left",parse=TRUE)+  
      # Blim
      annotate("text", x=blim*0.98, y=ftrg*1.03, label=paste0("B[lim]"), 
               vjust="bottom",hjust="right",parse=TRUE)  
  }  
  
  if(!missing(obs)) {
    p <- p + geom_point(data=obs)
  }
  return(p)
}  
