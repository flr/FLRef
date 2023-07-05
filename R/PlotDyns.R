#{{{
#' plotdyn()
#
#' Plots stock trajectories at age 
#'
#' @param stk stock object class FLStock  
#' @param ncol number of columns in multiplot 
#' @return ggplot  
#' @export
#' @examples 
#' data(ple4)
#' plotdyn(ple4)

plotdyn <- function(stk,ncol=2){
  if(dim(stk)[6]>1) stk = stockMedians(stk)
  if(any(dim(stk)[3:5]>1)){
  stk = simplify(stk)
  }
  dat = as.data.frame(
    FLQuants(Biomass=stock.n(stk)*stock.wt(stk), Vuln.Bio=stock.n(stk)*stock.wt(stk)*catch.sel(stk),SSB=stock.n(stk)*stock.wt(stk)*mat(stk),
             Catch = catch.n(stk),F=harvest(stk),
             "Selectivity"=catch.sel(stk)))
  dat$Age = factor(dat$age)
  # Plotting dynamics at age
  p =ggplot2::ggplot(dat)+
    geom_line(aes(year,data,group=age,col=Age))+
    facet_wrap(~qname,scale="free",ncol=ncol)+
    theme_bw()+theme(legend.position="right")+
    xlab("Year")+ylab("Quantities")+
    theme(legend.position="right",legend.text = element_text(size=8),
          legend.title = element_text(size=9),legend.key.height = unit(.5, 'cm'))
  
  return(p)
}

#{{{
#' plotbioyr()
#
#' Plots stock N_a, W_a, M_a and Mat_a across years 
#'
#' @param stk stock object class FLStock  
#' @param ncol number of columns in multiplot  
#' @return ggplot  
#' @export
#' @examples 
#' data(ple4)
#' plotbioyr(ple4)

plotbioyr <- function(stk,ncol=2){
  if(dim(stk)[6]>1) stk = stockMedians(stk)
  if(any(dim(stk)[3:5]>1)){
    stk = simplify(stk)
  }
  flqs = metrics(stk,metrics=list(Numbers=stock.n, Weight=catch.wt,M=m,Maturity=mat))
  dat=as.data.frame(flqs)

  dat$Age = factor(dat$age)
  dat
  # Plotting dynamics at age
  
  p <- ggplot(dat,aes(year,data,group=age,col=Age))+
    geom_line()+
    facet_wrap(~qname,scale="free",ncol=ncol)+
    theme_bw()+
    xlab("Year")+ylab("Quantities")#+
    theme(legend.position="right",legend.text = element_text(size=5),
          legend.title = element_text(size=6),legend.key.height = unit(.3, 'cm'))
  return(p)
}

#{{{
#' plotbioage()
#
#' Plots stock N_a, W_a, M_a and Mat_a by year 
#'
#' @param stk stock object class FLStock  
#' @param ncol number of columns in multiplot  
#' @return ggplot  
#' @export
#' @examples 
#' data(ple4)
#' plotbioage(ple4)

plotbioage = function(stk,ncol=2){
  if(dim(stk)[6]>1) stk = stockMedians(stk)
  if(any(dim(stk)[3:5]>1)){
    stk = simplify(stk)
  }
 
  dat = as.data.frame(
    FLQuants(Weight=stock.wt(stk),Maturity=mat(stk),M=m(stk),Selectivity=catch.sel(stk))
  )
  dat$Year = factor(dat$year)
  p =ggplot2::ggplot(dat)+
    geom_line(aes(age,data,group=Year,col=Year))+
    facet_wrap(~qname,scale="free",ncol=ncol)+
    theme_bw()+
    theme(legend.position="right",legend.text = element_text(size=5),
          legend.title = element_text(size=6),legend.key.height = unit(.3, 'cm'))+
    xlab("Age")+ylab("Value-at-age")
  return(p)
}

#{{{
#' plotspr()
#
#' Plots current vs unfished spawning biomass per recruit at age
#'
#' @param stk stock object class FLStock  
#' @param ncol number of columns in multiplot
#' @param nyears number of current last years, default is 3  
#' @return ggplot  
#' @export
#' @examples 
#' data(ple4)
#' plotbioage(ple4)

plotspr = function(stk,nyears=3){ 
  if(dim(stk)[6]>1) stk = stockMedians(stk)
  if(any(dim(stk)[3:5]>1)){
    stk = simplify(stk)
  }
  flqs = sprFy(stk[,,1])/spr0y(stk[,,1])
  df = as.data.frame(flqs)
  Fcur = round(mean(tail(fbar(stk),nyears)),2)
  
  df= FLQuants(
    "Current"=yearMeans(tail(sprFy(stk,byage=T),nyears)),
    "Unfished"=yearMeans(tail(spr0y(stk,byage=T),nyears)))
  names(df) = c(paste0("Current(F = ",Fcur,")"),paste("Unfished (F = 0)"))
  
  p=ggplot2::ggplot(as.data.frame(df),aes(x=age,y=data,fill=qname))+
    geom_bar(position="dodge", stat="identity")+theme_bw()+
    theme(legend.position="bottom",legend.title = element_blank(),legend.text = element_text(size=6),legend.key.height = unit(.4, 'cm'))+
    xlab("Age")+ylab("SPR-at-age")
 return(p)
}
