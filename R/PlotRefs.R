
#{{{
#' ploteq()
#
#' Modification of method plot(FLBRP) to plot equilibrium output of computeFbrp()  
#'
#' @param brp output object from computeFbrp of class FLBRP  
#' @param refpts Reference points, defaults are computed refpts from computeFbrp()  
#' \itemize{
#'   \item Fbrp Fmsy proxy 
#'   \item Blim 
#'   \item B0  unfished biomass for reference period
#'   \item Btri Btrigger as fraction of Btrg corresponding Fbrp 
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
                   shapes="missing", colours="missing", panels=NULL, ncol=2, ...) {
            
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
                                  aes_(x=~data, y=~y, group=~refpt, fill=~refpt, shape=~refpt)) +
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
