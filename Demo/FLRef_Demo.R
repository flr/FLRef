#----------------------------------------------------------------
# Demo: FLRef a proto type to support reference point estimation
#----------------------------------------------------------------

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Installation: Latest Versions of FLR
# library(devtools)
# install_github("flr/FLCore")
# install_github("flr/FLBRP")
# install_github("flr/FLFlasher")
# install_github("flr/mse")
# install_github("flr/FLSRTMB")
# install_github("henning-winker/FLRef")

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


library(FLCore)
library(mse)
library(FLSRTMB)
library(FLRef)
data(ple4)
stk = ple4

# Fit SSR BevHolt with FishLife priors
bh= srrTMB(as.FLSR(stk,model=bevholtSV),spr0=spr0y(stk),nyears=3)
# Fit Hockey-Stick
hs= srrTMB(as.FLSR(stk,model=segreg),spr0=spr0y(stk),plim=0.1,pmax=.3,nyears = 3)
# Estimate Blim Type 1 from Hockey-Stick
blim1 = params(hs)[[2]]
plot(FLSRs(bh=bh,hs=hs))+theme_bw()+theme(legend.position = "right")

plot(FLSRs(bh=bh))+theme_bw()+theme(legend.position = "none")


# Check Productivity 
rG = productivity(stk)
plot(rG)

# Get suggested Fbrp
Fref = rGclass(r=mean(rG$r),gt=mean(rG$gt))
Fref

# Compute Fspr40 with Blim type 1 with a BevHolt SRR
brp.hs = computeFbrp(stk,sr=hs,proxy="sprx",x=Fref$Fspr,blim=blim1)
# Compute Fb40 with Blim type 2 set to Blim = 0.1B0 with a Hockey-Stick
brp.bh = computeFbrp(stk,sr=bh,proxy="bx",x=Fref$Fsb,blim=blim1)

# Check Refpoints 
Fbrp(brp.hs)
Fbrp(brp.bh)
# Plot
ploteq(brp.hs)
ploteq(brp.bh)


p1= plotWKREF(btrigger=0.6,kobe=F,rel=T)
p2= plotWKREF(btrigger=0.6,kobe=F)


# Run Fsim with 500 iterations, using FishLife Priors
sigmaR = 0.6
rho = 0.3
sim.bh = Fsim(brp.bh,sigmaR=sigmaR,rho=rho,iters=500)
sim.hs = Fsim(brp.hs,sigmaR=sigmaR,rho=rho,iters=500)

# Plot 
plotFsim(sim.bh)
plotFsim(sim.hs)

# Compute Fp.05 for Hockey Stpck example
Fp.05 = Fp05(sim.hs)

# ICES Advice Rule
plotICES(fmsy=stk@benchmark[["Fmsy"]],
         blim=stk@benchmark[["Blim"]],
         btrigger=stk@benchmark[["Btrigger"]],
         obs=stk)

# WKREF Advice Rule
plotWKREF(ftrg=an(Fbrp(brp.bh)["Fbrp"]),
         blim=an(Fbrp(brp.bh)["Blim"]),
         btrg=an(Fbrp(brp.bh)["Btrg"]),
         obs=stk)

plotWKREF(ftrg=an(Fbrp(brp.bh)["Fbrp"]),
          blim=an(Fbrp(brp.bh)["Blim"]),
          btrg=an(Fbrp(brp.bh)["Btrg"]),
          obs=stk)



# Illustrate new advice rule
plotWKREF(xmax=1.8)+
  annotate("text", x=0.4, y=1.35, label="Overfished", 
           hjust="center",cex=4)+
  annotate("text", x=1.25, y=1.35, label="Overfishing", 
           hjust="center",cex=4)+
  annotate("text", x=0.4, y=0.2, label="Rebuilding", 
           hjust="center",cex=4)+
  annotate("text", x=1.25, y=0.2, label="Sustainable", 
           hjust="center",cex=4)+
  annotate("text", x=0.1, y=0.7, label="Critical", 
           hjust="center",cex=4,angle = 90)





# Close fishery at Blim and adjust axis labels to relative
plotWKREF(blim=0.2,bclose=0.2,rel=TRUE)
# Close fishery at Blim, but allow fmin (e.g. bycatch)
plotWKREF(blim=0.2,bclose=0.2,fmin=0.1,rel=TRUE)
# Change Btrigger above Btrg
plotWKREF(blim=0.2,bclose=0.2,fmin=0.1,btrigger=1.1,rel=TRUE)
# Plot stock data
data(ple4)
plotWKREF(ftrg=0.25,btrg=8e+05,btrigger = 0.9*8e+05, blim=2e5,bclose=3e5,fmin=0.03,obs=ple4)

# Figure 1
plotWKREF(xmax=1.8)+
  annotate("text", x=0.4, y=1.35, label="Overfished", 
           hjust="center",cex=4)+
  annotate("text", x=1.25, y=1.35, label="Overfishing", 
           hjust="center",cex=4)+
  annotate("text", x=0.4, y=0.2, label="Rebuilding", 
           hjust="center",cex=4)+
  annotate("text", x=1.25, y=0.2, label="Sustainable", 
           hjust="center",cex=4)+
  annotate("text", x=0.1, y=0.7, label="Critical", 
           hjust="center",cex=4,angle = 90)
