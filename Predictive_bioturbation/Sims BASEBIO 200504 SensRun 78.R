# Author: Sebastiaan van de Velde
###############################################################################
is.local <- FALSE
FOLDER   <- "Sensitivityrun"
NUMBER <- 78
BASENAME <- paste(NUMBER,"Sens BASEBIO Transient")
#=============================================================================
# Compile packages for current node
if (!is.local){
  
  dir.create(BASENAME)
  
  install.packages("../packages/gsw_1.0-5.tar.gz", repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/oce_1.2-0.tar.gz", repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/rootSolve_1.8.2.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/deSolve_1.27.1.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/shape_1.4.4.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/ReacTran_1.4.3.1.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/AquaEnv_1.0-4.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/seacarb_3.2.12.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  install.packages("../packages/marelac_2.1.10.tar.gz",repos=NULL, type="source", lib=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  
  # public packages
  require(AquaEnv,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(rootSolve,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(deSolve,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(shape,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(ReacTran,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(gsw,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(oce,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(seacarb,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  require(marelac,lib.loc=paste("/theia/home/brussel/102/vsc10244/",FOLDER,"/",BASENAME,"/",sep=""))
  
  # non-public packages
  require(diagenesis,lib.loc=getwd())
}
#=============================================================================
# Source file containing model function
if (!is.local){source("../Model_CSFe_FeisoBASE_HPCversion_final.R")}
if (is.local){source("../Model_CSFe_FeisoBASE_localversion_v07.R")}
#=============================================================================
# Aux function (blame Dale)
loc.erf <- function (x) 
{
  pchisq(2 * x^2, 1) * sign(x)
}
p.Dale <- function(x, Db.0 = 0., x.att = 1){
  return(Db.0*exp(-(x^2/(2*x.att^2))))
}
#=============================================================================
# parameters to be tested
F.POC.sens   <- 12
O2.ow.sens   <- 100
F.FeOOH.sens <- 1110
# Number of simulations
model <- CSFe.model
#=============================================================================
# Dynamic simulation 1 (start from roughs initial conditions)
# set-up OMM
#=============================================================================
# Initialisation simulation type 
sim.info$index <- 1
sim.info$code <- BASENAME
sim.info$name <- paste(sim.info$code,sim.info$index)
PL <- initialise.parameters(PL)
# Adapt parameter list
# PL$simulation.type <- "direct.steady.state"
PL$simulation.type <- "time.dependent"
# OMM dynamics 
# OMM dynamics (deconvolution following Dale et al., 2015)
PL$F.POC    <- F.POC.sens*365.25/10
PL$F.POC.1  <- PL$F.POC*0.021746
PL$F.POC.2  <- PL$F.POC*0.00725275
PL$F.POC.3  <- PL$F.POC*0.0096717
PL$F.POC.4  <- PL$F.POC*0.0128974
PL$F.POC.5  <- PL$F.POC*0.017199
PL$F.POC.6  <- PL$F.POC*0.0229352
PL$F.POC.7  <- PL$F.POC*0.0305846
PL$F.POC.8  <- PL$F.POC*0.0407852
PL$F.POC.9  <- PL$F.POC*0.0543879
PL$F.POC.10 <- PL$F.POC*0.0725265
PL$F.POC.11 <- PL$F.POC*0.0967046
PL$F.POC.12 <- PL$F.POC*0.12881
PL$F.POC.13 <- PL$F.POC*0.169822
PL$F.POC.14 <- PL$F.POC*0.314677
# Bottom water conditions
PL$TC    <- 10.            # temperature
PL$P     <- 36.           # pressure relation to depth? (1 bar per 10 m + 1 bar atmospheric pressure) Dale et al. 2015
c.fac <- 10000*(3600*24*365.25)
PL$por.0     <- 0.9  # porosity at the sediment-water interface
PL$por.inf   <- 0.7  # asymptotic porosity at depth
PL$por.x.att <- 10.  # attenuation depth [cm]
PL$por.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$por.0,y.inf=PL$por.inf,x.L=0,x.att=PL$por.x.att)
PL$D.HCO3.grid <- init.diffusion.coefficient(species="HCO3",   S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.NH4.grid  <- init.diffusion.coefficient(species="NH4",    S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.O2.grid   <- init.diffusion.coefficient(species="O2",     S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.NO3.grid  <- init.diffusion.coefficient(species="NO3",    S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.NO2.grid  <- init.diffusion.coefficient(species="NO2",    S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.SO4.grid  <- init.diffusion.coefficient(species="SO4",    S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.N2.grid   <- init.diffusion.coefficient(species="N2",     S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.Mn.grid   <- init.diffusion.coefficient(species="Mn",     S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.Fe.grid   <- init.diffusion.coefficient(species="Fe",     S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.HS.grid   <- init.diffusion.coefficient(species="HS",     S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.CH4.grid  <- init.diffusion.coefficient(species="CH4",    S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.H2O.grid  <- init.diffusion.coefficient(species="H2O",    S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.H.grid    <- init.diffusion.coefficient(species="H",      S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$D.H2.grid   <- init.diffusion.coefficient(species="H2",     S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
PL$O2.ow <- O2.ow.sens/1000 #  
# sedimentation determined 
PL$svf.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=(1-PL$por.0),y.inf=(1-PL$por.inf),x.L=0,x.att=PL$por.x.att)
PL$u.inf   <- 0.06  # sedimentation velocity pore water [cm yr-1] Dale et al. 2015
PL$v.inf   <- 0.06  # sedimentation velocity solids [cm yr-1] Dale et al. 2015
PL$v.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$v
PL$u.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$u
# Meysman et al., 2003
PL$F.FeOOH.HR      <- 0.0*F.FeOOH.sens/6*365.25*1e-4
PL$F.Fe56_FeOOH.HR <- 0.0*((PL$d56Fe_FeOOH.HR*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.HR*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.HR
PL$F.FeOOH.MR      <- 0.0*F.FeOOH.sens/6*365.25*1e-4
PL$F.Fe56_FeOOH.MR <- 0.0*((PL$d56Fe_FeOOH.MR*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.MR*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.MR
PL$F.FeOOH.PR      <- 0.0*F.FeOOH.sens/6*365.25*1e-4
PL$F.Fe56_FeOOH.PR <- 0.0*((PL$d56Fe_FeOOH.PR*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.PR*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.PR
PL$F.FeOOH.U       <- 0.0*F.FeOOH.sens/2*365.25*1e-4
PL$F.Fe56_FeOOH.U  <- 0.0*((PL$d56Fe_FeOOH.U*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.U*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.U
#PL$Db.0   <- 23.0*(0.5+0.5*loc.erf((O2.ow.sens-20.0)/12.0))
#PL$x.L    <- 4.0 # to be near-identical as Dale et al. 2015
#PL$x.att  <- 3.0 # to be near-identical as Dale et al. 2015
#PL$y.inf  <- 0.0 # to be near-identical as Dale et al. 2015
#PL$Db.grid <- setup.prop.1D(func=p.Dale, Db.0=PL$Db.0, x.att=PL$x.att, grid = PL$grid)
PL$Db.0   <- 10.0*(0.5+0.5*loc.erf((O2.ow.sens-20.0)/12.0))
PL$x.L    <- 1.0 + 9.0*(1-exp(-PL$Db.0/3.0))
PL$x.att  <- 2.0
PL$y.inf  <- 0.0
PL$Db.grid <- setup.prop.1D(func=p.sig,y.0=PL$Db.0,y.inf=PL$y.inf,x.att=PL$x.att,x.L=PL$x.L, grid = PL$grid)
PL$irr.0   <- 290.0*(0.5+0.5*loc.erf((O2.ow.sens-20.0)/12.0)) # irrigation rate at SWI [yr-1]
PL$irr.att <- 1.4 
PL$irr.inf <- 0.0
PL$irr.grid <- setup.prop.1D(func=p.exp,y.0=PL$irr.0,y.inf=PL$irr.inf,x.L=PL$x.irr,x.att=PL$irr.att ,grid=PL$grid)
PL$irr.Fe   <- 0.2
PL$irr.HS   <- 1.0
PL$k.FIO     <- 0.0*1E+07
PL$k.NFO     <- 0.0*1E+02
PL$k.ISP     <- 0.0*1E+03
PL$k.PyP     <- 0.0*3.25
PL$K.FIS     <- 0.0*696.8
PL$K_MnO2  <- 8.0*PL$rho.sed 
PL$K_FeOOH <- 100.0*PL$rho.sed
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------
sim.info$time.seq <- c(0,1,6,12,24,1*24*365.25,1000*24*365.25,10000*24*365.25)/(24*365.25) 
#sim.info$time.seq <- c(0,1,2,3,4,5,6,10)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","1 yr","1000 yr","10000 yr")
sim.info$N.out <- length(sim.info$time.seq)
#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------
# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix
# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV.init, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)
# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# add FeOOH
#=============================================================================
load(paste(BASENAME,"1.Rdata"))
unlink(paste(BASENAME,"1.Rdata"))
# Initialisation simulation type 
sim.info$index <- 2
sim.info$code <- BASENAME
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL
# Adapt parameter list
PL$simulation.type <- "time.dependent"
#PL$simulation.type <- "direct.steady.state"
# Iron input
PL$F.FeOOH.HR      <- F.FeOOH.sens/6*365.25*1e-4
PL$F.Fe56_FeOOH.HR <- ((PL$d56Fe_FeOOH.HR*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.HR*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.HR
PL$F.FeOOH.MR      <- F.FeOOH.sens/6*365.25*1e-4
PL$F.Fe56_FeOOH.MR <- ((PL$d56Fe_FeOOH.MR*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.MR*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.MR
PL$F.FeOOH.PR      <- F.FeOOH.sens/6*365.25*1e-4
PL$F.Fe56_FeOOH.PR <- ((PL$d56Fe_FeOOH.PR*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.PR*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.PR
PL$F.FeOOH.U       <- F.FeOOH.sens/2*365.25*1e-4
PL$F.Fe56_FeOOH.U  <- ((PL$d56Fe_FeOOH.U*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.U*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.U
# Isotope dynamics
PL$a56Fe_DIR     <- (1.0 - 1.3/1000.0) # alpha for Dissimilatory Iron Reduction
PL$a56Fe_ISO     <- (1.0 - 1.3/1000.0) # alpha for Dissimilatory Iron Reduction
PL$a56Fe_FIO     <- (1.0 + 0.4/1000.0) # alpha for Ferrous Iron oxidation
PL$a56Fe_ISP     <- (1.0 + 0.5/1000.0) # alpha for Iron sulphide precipitation
PL$a56Fe_AIO     <- (1.0 + 0.4/1000.0) # alpha for Ferrous Iron oxidation
PL$a56Fe_ISD     <- (1.0 - 0.5/1000.0) # alpha for Ferrous Iron oxidation
PL$a56Fe_PyP     <- (1.0 - 0.7/1000.0) # alpha for Pyrite precipitation
PL$a56Fe_PyO     <- (1.0 + 0.5/1000.0) # alpha for Pyrite precipitation
PL$a56Fe_IOA     <- (1.0 + 0.0/1000.0) # alpha for Iron oxide aging
PL$a56Fe_FIS     <- (1.0 + 0.4/1000.0) # alpha for Ferrous Iron Sorption
PL$k.FIO     <- 5E+05
PL$k.NFO     <- 1E+02
PL$k.ISP     <- 0.0*1E+03
PL$k.PyP     <- 0.0*3.25
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------
sim.info$time.seq <- c(0,1,6,12,24,1*24*365.25,1000*24*365.25,10000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","1 yr","1000 yr","10000 yr")
sim.info$N.out <- length(sim.info$time.seq)
#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------
# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix
# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV.init, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)
# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# Precipitate FeS
#=============================================================================
load(paste(BASENAME,"2.Rdata"))
# Initialisation simulation type 
sim.info$index <- 3
sim.info$code <- BASENAME
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL
# Adapt parameter list
PL$simulation.type <- "time.dependent"
#PL$simulation.type <- "direct.steady.state"
PL$k.ISP     <- 1E+03
PL$k.PyP.1   <- 3.25 #1E+02 # -> Dale et al., 2015 - 3.25 (is from myself)
PL$k.PyP.2   <- 3.25 #1E+02
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------
sim.info$time.seq <- c(0,1,6,12,24,1*24*365.25,100*24*365.25,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","1 yr","100 yr","500 yr")
sim.info$N.out <- length(sim.info$time.seq)
#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------
# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix
# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV.init, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)
# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
unlink(paste(BASENAME,"2.Rdata"))
#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# Precipitate FeS2
#=============================================================================
load(paste(BASENAME,"3.Rdata"))
# Initialisation simulation type 
sim.info$index <- 4
sim.info$code <- BASENAME
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[length(simulation.list)]]
PL <- sim.start$PL
# Adapt parameter list
PL$simulation.type <- "time.dependent"
#PL$simulation.type <- "direct.steady.state"
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------
sim.info$time.seq <- c(0,1,6,24,1*24*365.25,10*24*365.25,100*24*365.25,500*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","1 d","1 year","10 yr","100 yr","500 yr")
sim.info$N.out <- length(sim.info$time.seq)
#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------
# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix
# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV.init, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 300000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)
# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
unlink(paste(BASENAME,"3.Rdata"))
#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# Precipitate FeS2
#=============================================================================
load(paste(BASENAME,"4.Rdata"))
# Initialisation simulation type 
sim.info$index <- 5
sim.info$code <- BASENAME
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[length(simulation.list)]]
PL <- sim.start$PL
# Adapt parameter list
PL$simulation.type <- "time.dependent"
#PL$simulation.type <- "direct.steady.state"
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------
sim.info$time.seq <- c(0,6,24,100*24*365.25,500*24*365.25,1000*24*365.25,2000*24*365.25,5000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","6 hr","1 d","100 yr","500 yr","1000 yr","2000 yr","5000 yr")
sim.info$N.out <- length(sim.info$time.seq)
#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------
# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix
# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV.init, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)
# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
unlink(paste(BASENAME,"4.Rdata"))
#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# Precipitate FeS2
#=============================================================================
load(paste(BASENAME,"5.Rdata"))
# Initialisation simulation type 
sim.info$index <- 6
sim.info$code <- BASENAME
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[length(simulation.list)]]
PL <- sim.start$PL
# Adapt parameter list
PL$simulation.type <- "time.dependent"
#PL$simulation.type <- "direct.steady.state"
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------
sim.info$time.seq <- c(0,6,24,100*24*365.25,500*24*365.25,1000*24*365.25,2000*24*365.25,5000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","6 hr","1 d","100 yr","500 yr","1000 yr","2000 yr","5000 yr")
sim.info$N.out <- length(sim.info$time.seq)
#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------
# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix
# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV.init, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)
# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
unlink(paste(BASENAME,"5.Rdata"))
#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# Precipitate FeS2
#=============================================================================
load(paste(BASENAME,"6.Rdata"))
# Initialisation simulation type 
sim.info$index <- 7
sim.info$code <- BASENAME
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[length(simulation.list)]]
PL <- sim.start$PL
# Adapt parameter list
PL$simulation.type <- "time.dependent"
#PL$simulation.type <- "direct.steady.state"
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------
sim.info$time.seq <- c(0,6,24,100*24*365.25,500*24*365.25,1000*24*365.25,2000*24*365.25,5000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","6 hr","1 d","100 yr","500 yr","1000 yr","2000 yr","5000 yr")
sim.info$N.out <- length(sim.info$time.seq)
#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------
# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix
# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV.init, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)
# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
unlink(paste(BASENAME,"6.Rdata"))
#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# Precipitate FeS2
#=============================================================================
load(paste(BASENAME,"7.Rdata"))
# Initialisation simulation type 
sim.info$index <- 8
sim.info$code <- BASENAME
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[length(simulation.list)]]
PL <- sim.start$PL
# Adapt parameter list
PL$simulation.type <- "time.dependent"
#PL$simulation.type <- "direct.steady.state"
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------
sim.info$time.seq <- c(0,6,24,100*24*365.25,500*24*365.25,1000*24*365.25,2000*24*365.25,5000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","6 hr","1 d","100 yr","500 yr","1000 yr","2000 yr","5000 yr")
sim.info$N.out <- length(sim.info$time.seq)
#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------
# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix
# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV.init, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)
# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
unlink(paste(BASENAME,"7.Rdata"))
#=============================================================================
# Dynamic simulation 3 (start from previous initial conditions)
# check steady state
#=============================================================================
load(paste(BASENAME,"8.Rdata"))
# Initialisation simulation type 
sim.info$index <- 1
sim.info$code <- paste(NUMBER,"Sens BASEBIO 200504 SteadySta")
sim.info$name <- paste(sim.info$code,sim.info$index)
# Initialisation parameter list
sim.start <- simulation.list[[length(simulation.list)]]
PL <- sim.start$PL
# Adapt parameter list
#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"
PL <- initialise.parameters(PL)
# Initialise set of state variables
PL <- initialise.state.variables(PL)
# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)
# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)
# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)
# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)
#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------
#Actual steady state simulation
sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)
#prepare matrix out
out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2
#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------
# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation
save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
unlink(paste(BASENAME,"8.Rdata"))
#=============================================================================
# Send email to say you are finished (and clean up after yourself)
#=============================================================================
if (!is.local){
  unlink(BASENAME,recursive=T)
}
if (is.local){
  
  require(mailR)
  
  send.mail(from = "rsebsebr@gmail.com",
            to = c("sebastiv@ucr.edu"),
            subject = paste("I think ", BASENAME, "is done!"),
            body = "Hello, I think the code you were running has finished. Best come check what mistakes you made ...",
            smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "rsebsebr@gmail.com", passwd = "Rmail888", ssl = TRUE),
            authenticate = TRUE,  send = TRUE)
}
