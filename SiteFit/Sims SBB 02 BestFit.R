###############################################################################
# Simulation file - try to fit Monterey Bay site from Severmann et al., 2006, GCA
# CSFe - Fe isotope model 
# Author: Sebastiaan van de Velde
###############################################################################

is.local <- TRUE
FOLDER   <- "SBBrun"
BASENAME <- "01 SBB Transient BestFit"

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
if (is.local){source("../Model_CSFe_Feiso_localversion_v03.R")}

# Number of simulations

model <- CSFe.model

#=============================================================================
# Dynamic simulation 1 (start from previous initial conditions)
# add isotope expressions 
#=============================================================================

load("02 SBB Transient BestFit 1.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "02 SBB Transient BestFit"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"
#PL$simulation.type <- "direct.steady.state"

# Isotope expression
PL$a56Fe_DIR     <- (1.0 - 1.3/1000.0) # alpha for Dissimilatory Iron Reduction
PL$a56Fe_CSFO    <- (1.0 - 1.3/1000.0) # alpha for Dissimilatory Iron Reduction
PL$a56Fe_FIO     <- (1.0 + 0.4/1000.0) # alpha for Ferrous Iron oxidation
PL$a56Fe_ISP     <- (1.0 + 0.5/1000.0) # alpha for Iron sulphide precipitation
PL$a56Fe_ISO     <- (1.0 + 0.4/1000.0) # alpha for Ferrous Iron oxidation
PL$a56Fe_ISD     <- (1.0 - 0.5/1000.0) # alpha for Ferrous Iron oxidation
PL$a56Fe_PyP     <- (1.0 - 0.7/1000.0) # alpha for Pyrite precipitation
PL$a56Fe_PyO     <- (1.0 + 0.5/1000.0) # alpha for Pyrite precipitation
PL$a56Fe_IOA     <- (1.0 + 0.0/1000.0) # alpha for Iron oxide aging
PL$a56Fe_FIS     <- (1.0 + 0.4/1000.0) # alpha for Ferrous Iron Sorption

PL$F.FeOOH.f      <- 20.525
PL$d56Fe_FeOOH.f  <- -1.5
PL$d56Fe_FeOOH.a  <- -1.5
PL$F.Fe56_FeOOH.f <- ((PL$d56Fe_FeOOH.f*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.f*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.f
PL$F.Fe56_FeOOH.a <- ((PL$d56Fe_FeOOH.a*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.a*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.a

PL$F.FeS2      <- 1.0
PL$d56Fe_FeS2  <- -0.4
PL$F.Fe56_FeS2 <- ((PL$d56Fe_FeS2*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeS2*1e-3+1.0)*PL$IRMM014))*PL$F.FeS2

PL$P <- 51

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

# print mass balance
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

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

load("02 SBB Transient BestFit 1.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "02 SBB SteadySta BestFit"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
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

# print mass balance
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Send email to say you are finished
#=============================================================================

if (is.local){
  
  require(mailR)
  
  send.mail(from = "rsebsebr@gmail.com",
            to = c("sebastiv@ucr.edu"),
            subject = paste("I think ", sim.info$code,"is done!"),
            body = "Hello, I think the code you were running has finished. Best come check what mistakes you made ...",
            smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "rsebsebr@gmail.com", passwd = "Rmail888", ssl = TRUE),
            authenticate = TRUE,  send = TRUE)
}
