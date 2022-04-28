

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
