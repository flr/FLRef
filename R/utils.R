

# rGclass {{{ 

#' Function to characterize Productivity and refpts based on r and Generation
#'
#' @param r value of the intrinsic rate of population increase
#' @param gt generation time G
#'
#' @return list with Productivity category and suggest Fbrps

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
