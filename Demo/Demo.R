
library(FLCore)
library(mse)
library(FLSRTMB)
library(FLRef)

data(ple4)
# Fit SSR BevHolt
bh= srrTMB(as.FLSR(ple4,model=bevholtSV),spr0=spr0y(ple4))
# Fit Hockey-Stick
hs= srrTMB(as.FLSR(ple4,model=segreg),spr0=spr0y(ple4),plim=0.05,pmax=.3)
blim2 = params(hs)[[2]]
plot(FLSRs(bh=bh,hs=hs))



# Check Productivity 
rG = productivity(ple4)
plot(rG)

# Get suggest Fbrp
Fref = rGclass(r=mean(rG$r),gt=mean(rG$gt))
Fref
# Compute Fspr40 with Blim type 2
brp.bh = computeFbrp(ple4,sr=bh,proxy="bx",x=Fref$Fsb,blim=0.1)
# Compute Fb35 with Blim type 1
brp.hs = computeFbrp(ple4,sr=bh,proxy="sprx",x=Fref$Fspr,blim=blim2)

# Check 
Fbrp(brp.bh)
Fbrp(brp.hs)

# Plot
ploteq(brp.bh)
ploteq(brp.hs)

# Run Fsim
sim.bh = Fsim(brp.bh,sigmaR=0.6,rho=0.4,iters=500)
sim.hs = Fsim(brp.hs,sigmaR=0.6,rho=0.4,iters=500)

# Plot 
plotFsim(sim.bh)
plotFsim(sim.hs)

# Compute Fp.05
Fp.05 = Fp05(sim.hs)

# TODO Advice plot with input stock, brp, Fp.05

# Proposed advice
plotWKREF()
# Add zones (Figure 1)

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
