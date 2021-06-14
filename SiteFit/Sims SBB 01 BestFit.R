###############################################################################
# Simulation file - try to fit Monterey Bay site from Severmann et al., 2006, GCA
# CSFe - Fe isotope model 
# Author: Sebastiaan van de Velde
###############################################################################

is.local <- FALSE
FOLDER   <- "SBBrun"
BASENAME <- "01 SBB Transient FPOC4_6 FFeOOH20_25"

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

if (!is.local){source("../Model_CSFe_Feiso_HPCversion_v02.R")}
if (is.local){source("../../Model_CSFe_Feiso_localversion_v03.R")}

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

# OMM dynamics (Severmann et al., 2006 GCA) -> 6 to 12 mmol m-2 d-1
# OMM dynamics (deconvolution following Dale et al., 2015)
PL$F.POC    <- 4.6*365.25/10
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

#0.66
PL$K_O2    <- 0.002
PL$K_FeOOH <- 12*2.6

# anoxic, full marine site

PL$pH      <- 8.0  # After 2 cm in sediment (Reimers et al., 1996, GCA)
PL$pH.grid <- setup.prop.1D(value=PL$pH,grid=PL$grid)
PL$O2.ow   <- 0.01 # Severmann et al., 2006, GCA
PL$SO4.ow  <- 28.0 # Meysman et al., 2003
PL$HCO3.ow <- 2.45 # Meysman et al., 2003

# sedimentation determined by Pb dating (Severmann et al., 2006, GCA)

PL$u.inf   <- 0.25  # sedimentation velocity pore water [cm yr-1]
PL$v.inf   <- 0.25  # sedimentation velocity solids [cm yr-1]
PL$v.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$v
PL$u.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$u

# Meysman et al., 2003
PL$F.FeOOH.f      <- 0.0*20.25
PL$F.Fe56_FeOOH.f <- 0.0*((PL$d56Fe_FeOOH.f*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.f*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.f
PL$F.FeOOH.a      <- 0.0*0.0
PL$F.Fe56_FeOOH.a <- ((PL$d56Fe_FeOOH.a*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.a*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.a

PL$k.FIO     <- 0.0*1E+07
PL$k.ISP     <- 0.0*1E+03
PL$k.PyP     <- 0.0*3.25
PL$frac.limit <- 1e-6
PL$K.FIS      <- 0.1*PL$rho.sed*268

# Tuning
PL$k.IOA   <- 0.6#1.0

PL$x.L    <- 1.0
PL$x.att  <- 0.5
PL$y.inf  <- 0.0
PL$Db.0   <- 0.1
PL$Db.grid <- setup.prop.1D(func=p.sig,y.0=PL$Db.0,y.inf=PL$y.inf,x.att=PL$x.att,x.L=PL$x.L, grid = PL$grid)

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

# Meysman et al., 2003
PL$F.FeOOH.f      <- 20.25 
PL$F.Fe56_FeOOH.f <- ((PL$d56Fe_FeOOH.f*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.f*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.f
PL$F.FeOOH.a      <- 0.0
PL$F.Fe56_FeOOH.a <- ((PL$d56Fe_FeOOH.a*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.a*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.a

PL$k.FIO     <- 1E+07
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
PL$k.PyP     <- 3.25

PL$F.FeS2      <- 1.0

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
sim.start <- simulation.list[[sim.info$N.out]]
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


#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# Precipitate FeS2
#=============================================================================

#load(paste(BASENAME,"8.Rdata"))
load("01 SBB Transient BestFit 1.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "01 SBB Transient BestFit"
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

sim.info$time.seq <- c(0,6,24,100*24*365.25,500*24*365.25,1000*24*365.25,2000*24*365.25,100000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","6 hr","1 d","100 yr","500 yr","1000 yr","2000 yr","100000 yr")
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
# Dynamic simulation 3 (start from previous initial conditions)
# check steady state
#=============================================================================

load(paste(BASENAME,"8.Rdata"))
substr(BASENAME,8,16) <- "SteadySta"

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- BASENAME
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

#=============================================================================
# Send email to say you are finished (and clean up after yourself)
#=============================================================================

if (!is.local){
  substr(BASENAME,8,16) <- "Transient"
  unlink(BASENAME,recursive=T)
  
  for (i in 1:8){
    unlink(paste(BASENAME,i))
  }
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
