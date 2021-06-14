###############################################################################
# Model of early diagenesis
# Type: Carbon,sulfur and iron cycling  
# Author: Filip Meysman (filip.meysman@nioz.nl)
#         Sebastiaan van de Velde (sevdevel@vub.ac.be) 
# Affiliation: NIOZ, Korringaweg 7, 4401 NT, Yerseke 
#              VUB, Pleinlaan 2, 1050 BE, Brussel
###############################################################################

# public packages
require(AquaEnv)
require(rootSolve)
require(deSolve)
require(shape)
require(ReacTran)
require(gsw)
require(oce)
require(seacarb)
require(marelac)

# non-public packages
require(diagenesis)

#=============================================================================
# Definition of specific functions
#=============================================================================

#-----------------------------------------------------------------------------
# Saturation function
#-----------------------------------------------------------------------------

FSAT <- function(C,K,n) (C/K)^n/((C/K)^n + 1)

#-----------------------------------------------------------------------------
# Function: tortuosity 
#-----------------------------------------------------------------------------

tortuosity <- function(por) 1-2*log(por)

#-----------------------------------------------------------------------------
# Function: init.diffusion.coefficient 
#-----------------------------------------------------------------------------

init.diffusion.coefficient <- function (species,S,TC,P,conv.fac,grid,tort.grid)
{
  Dmol <- conv.fac*diffcoeff(S=S,t=TC,P=P,species=species)[[species]]
  D.grid <- setup.prop.1D(value=Dmol,grid=grid)
  D.grid$mid <- D.grid$mid/tort.grid$mid
  D.grid$int <- D.grid$int/tort.grid$int
  return(D.grid)
}

#=============================================================================
# Definition of objects: species, reactions, elements, simulations
#=============================================================================

# SL = list of species incorporated in the model 
SL <- list()

#--------------------------------------------------
# chemical species incorporated in the model 
#--------------------------------------------------
#---------- C -------------------------------------

# Organic matter -> multi-G approx. of reactive continuum (following Dale et al. 2015) 

for (i in 1:14){
  SL <- initialise.list.species(SL,"1D",name=paste("POC.",i,sep=""),phase="solid",type="chemical")
}

SL <- initialise.list.species(SL,"1D",name="HCO3",phase="solute",type="chemical")

SL <- initialise.list.species(SL,"1D",name="CH4", phase="solute",type="chemical")

#---------- O -------------------------------------

SL <- initialise.list.species(SL,"1D",name="O2",phase="solute",type="chemical")

#---------- N -------------------------------------

SL <- initialise.list.species(SL,"1D",name="NO3",  phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="NO2",  phase="solute",type="chemical")

SL <- initialise.list.species(SL,"1D",name="NH4",  phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="X.NH4",phase="solid",type="chemical")

SL <- initialise.list.species(SL,"1D",name="N2",   phase="solute",type="chemical")

#---------- Mn ------------------------------------

SL <- initialise.list.species(SL,"1D",name="MnO2.HR",phase="solid", type="chemical")
SL <- initialise.list.species(SL,"1D",name="MnO2.MR",phase="solid", type="chemical")

SL <- initialise.list.species(SL,"1D",name="Mn",     phase="solute",type="chemical")

#---------- Fe ------------------------------------

SL <- initialise.list.species(SL,"1D",name="FeOOH.HR",     phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeOOH.HR",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeOOH.MR",     phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeOOH.MR",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeOOH.PR",     phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeOOH.PR",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeOOH.U",      phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeOOH.U", phase="solid",type="chemical")

SL <- initialise.list.species(SL,"1D",name="Fe",           phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_Fe",      phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="X.Fe",         phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="X.Fe56_Fe",    phase="solid",type="chemical")

SL <- initialise.list.species(SL,"1D",name="FeS",          phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeS",     phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeS2",         phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeS2",    phase="solid",type="chemical")

#---------- S -------------------------------------

SL <- initialise.list.species(SL,"1D",name="SO4",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="S0", phase="solid", type="chemical")
SL <- initialise.list.species(SL,"1D",name="HS", phase="solute",type="chemical")

#---------- other----------------------------------

SL <- initialise.list.species(SL,"1D",name="H",  phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="H2", phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="H2O",phase="solute",type="chemical")

#--------------------------------------------------
# Composite species incorporated in the model 
#--------------------------------------------------

SL <- initialise.list.species(SL,"1D",name="POC",       phase="solid",type="composite")
SL <- initialise.list.species(SL,"1D",name="MnO2",      phase="solid",type="composite")
SL <- initialise.list.species(SL,"1D",name="FeOOH",     phase="solid",type="composite")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeOOH",phase="solid",type="composite")

#--------------------------------------------------
# Artificial species incorporated in the model 
#--------------------------------------------------

SL <- initialise.list.species(SL,"1D",name="Omega.FeS",phase="solute",type="artificial")
SL <- initialise.list.species(SL,"1D",name="Omega.NH4",phase="solute",type="artificial")

#--------------------------------------------------
# Initialise species matrix 
#--------------------------------------------------

spec.mat <- initialise.species.matrix(SL)
spec.mat["POC",paste("POC.",1:14,sep="")] <- rep(1,14)
spec.mat["MnO2",c("MnO2.HR","MnO2.MR")] <- c(1,1)
spec.mat["FeOOH",c("FeOOH.HR","FeOOH.MR","FeOOH.PR","FeOOH.U")] <- c(1,1,1,1)
spec.mat["Fe56_FeOOH",c("Fe56_FeOOH.HR","Fe56_FeOOH.MR","Fe56_FeOOH.PR","Fe56_FeOOH.U")] <- c(1,1,1,1)

spec.table <- as.data.frame(spec.mat)
remove(spec.mat)

#--------------------------------------------------
# List of reactions included 
#--------------------------------------------------

RL <- list()

#---------- mineralization terms ------------------

for (i in 1:14){
  RL <- initialise.list.reaction(RL,"1D",name=paste("Cmin.",i,sep=""),type="kinetic")
}

RL <- initialise.list.reaction(RL,"1D",name="Cmin",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="H2min",type="kinetic")

#---------- individual mineralization terms -------

RL <- initialise.list.reaction(RL,"1D",name="AR",      type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DN.1",    type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DN.2",    type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MR",      type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DIR",     type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DIR_56Fe",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SR",      type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MG",      type="kinetic")

#---------- individual H2 oxidation terms ---------

RL <- initialise.list.reaction(RL,"1D",name="AR.H2",      type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DN.1.H2",    type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DN.2.H2",    type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MR.H2",      type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DIR.H2",     type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DIR_56Fe.H2",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SR.H2",      type="kinetic")

#---------- aerobic oxidation terms ---------------

RL <- initialise.list.reaction(RL,"1D",name="NIT.1",   type="kinetic") # NH4 + O2
RL <- initialise.list.reaction(RL,"1D",name="NIT.2",   type="kinetic") # NO2 + O2
RL <- initialise.list.reaction(RL,"1D",name="MnO",     type="kinetic") # Mn + O2
RL <- initialise.list.reaction(RL,"1D",name="FIO",     type="kinetic") # Fe + O2
RL <- initialise.list.reaction(RL,"1D",name="FIO_56Fe",type="kinetic") # Fe + O2
RL <- initialise.list.reaction(RL,"1D",name="SIO",     type="kinetic") # XFe + O2
RL <- initialise.list.reaction(RL,"1D",name="SIO_56Fe",type="kinetic") # XFe + O2
RL <- initialise.list.reaction(RL,"1D",name="CSO",     type="kinetic") # H2S + O2
RL <- initialise.list.reaction(RL,"1D",name="AMO",     type="kinetic") # CH4 + O2

RL <- initialise.list.reaction(RL,"1D",name="AIO",     type="kinetic") # FeS + O2
RL <- initialise.list.reaction(RL,"1D",name="AIO_56Fe",type="kinetic") # FeS + O2
RL <- initialise.list.reaction(RL,"1D",name="PyO",     type="kinetic") # FeS2 + O2
RL <- initialise.list.reaction(RL,"1D",name="PyO_56Fe",type="kinetic") # FeS2 + O2

#---------- anaerobic oxidation terms -------------

RL <- initialise.list.reaction(RL,"1D",name="ANA",type="kinetic") # ANAMOX

RL <- initialise.list.reaction(RL,"1D",name="NFO",     type="kinetic") # Fe2+ + NO3
RL <- initialise.list.reaction(RL,"1D",name="NFO_56Fe",type="kinetic") # Fe2+ + NO3
RL <- initialise.list.reaction(RL,"1D",name="NSO",     type="kinetic") # H2S + NO3

RL <- initialise.list.reaction(RL,"1D",name="MFO.HR",     type="kinetic") # Fe2+ + MnHR
RL <- initialise.list.reaction(RL,"1D",name="MFO.HR_56Fe",type="kinetic") # Fe2+ + MnHR
RL <- initialise.list.reaction(RL,"1D",name="MFO.MR",     type="kinetic") # Fe2+ + MnMR
RL <- initialise.list.reaction(RL,"1D",name="MFO.MR_56Fe",type="kinetic") # Fe2+ + MnMR
RL <- initialise.list.reaction(RL,"1D",name="MSO.HR",     type="kinetic") # H2S + MnHR
RL <- initialise.list.reaction(RL,"1D",name="MSO.MR",     type="kinetic") # H2S + MnMR

RL <- initialise.list.reaction(RL,"1D",name="ISO.HR",     type="kinetic") # H2S + FeHR
RL <- initialise.list.reaction(RL,"1D",name="ISO.HR_56Fe",type="kinetic") # H2S + FeHR
RL <- initialise.list.reaction(RL,"1D",name="ISO.MR",     type="kinetic") # H2S + FeMR
RL <- initialise.list.reaction(RL,"1D",name="ISO.MR_56Fe",type="kinetic") # H2S + FemR
RL <- initialise.list.reaction(RL,"1D",name="ISO.PR",     type="kinetic") # H2S + FePR
RL <- initialise.list.reaction(RL,"1D",name="ISO.PR_56Fe",type="kinetic") # H2S + FePR
RL <- initialise.list.reaction(RL,"1D",name="ISO.U",      type="kinetic") # H2S + FeU
RL <- initialise.list.reaction(RL,"1D",name="ISO.U_56Fe", type="kinetic") # H2S + FeU

RL <- initialise.list.reaction(RL,"1D",name="SMO",type="kinetic") # CH4 + SO4

#---------- mineral reactions -------------------

RL <- initialise.list.reaction(RL,"1D",name="AmS",type="kinetic") # NH4 -> ads
RL <- initialise.list.reaction(RL,"1D",name="SDP",type="kinetic") # S0

RL <- initialise.list.reaction(RL,"1D",name="ISP",     type="kinetic") # Fe + H2S 
RL <- initialise.list.reaction(RL,"1D",name="ISP_56Fe",type="kinetic") # Fe + H2S 
RL <- initialise.list.reaction(RL,"1D",name="ISD",     type="kinetic") # FeS 
RL <- initialise.list.reaction(RL,"1D",name="ISD_56Fe",type="kinetic") # FeS 

RL <- initialise.list.reaction(RL,"1D",name="PyP.1",     type="kinetic") # FeS + H2S   
RL <- initialise.list.reaction(RL,"1D",name="PyP.1_56Fe",type="kinetic") # FeS + H2S  
RL <- initialise.list.reaction(RL,"1D",name="PyP.2",     type="kinetic") # FeS + S0   
RL <- initialise.list.reaction(RL,"1D",name="PyP.2_56Fe",type="kinetic") # FeS + S0  

#---------- metal oxide ageing -----------------

RL <- initialise.list.reaction(RL,"1D",name="MOA",     type="kinetic") # MnHR
RL <- initialise.list.reaction(RL,"1D",name="IOA",     type="kinetic") # FeHR
RL <- initialise.list.reaction(RL,"1D",name="IOA_56Fe",type="kinetic") # FeHR

#-----------------------------------------------
# List of elements included 
#-----------------------------------------------

EL <- list() # create empty element object list

EL <- initialise.list.element(EL,"1D",name="C")
EL <- initialise.list.element(EL,"1D",name="N")
EL <- initialise.list.element(EL,"1D",name="O")
EL <- initialise.list.element(EL,"1D",name="H")
EL <- initialise.list.element(EL,"1D",name="S")
EL <- initialise.list.element(EL,"1D",name="Mn")
EL <- initialise.list.element(EL,"1D",name="Fe")
EL <- initialise.list.element(EL,"1D",name="56Fe")
EL <- initialise.list.element(EL,"1D",name="e")

# Initialise element matrix 

elt.mat <- initialise.element.matrix(SL, EL)

#---------- C -------------------------------------

for (i in 1:14){
  elt.mat[paste("POC.",i,sep=""),c("C","H","O")] <- c(1,2,1)
}
elt.mat["HCO3",c("H","C","O","e")] <- c(1,1,3,1)

elt.mat["CH4",c("C","H")] <- c(1,4)

#---------- O -------------------------------------

elt.mat["O2",c("O")] <- c(2)

#---------- N -------------------------------------

elt.mat["NO3",c("N","O","e")] <- c(1,3,1)
elt.mat["NO2",c("N","O","e")] <- c(1,3,1)

elt.mat["NH4",c("N","H","e")]   <- c(1,4,-1)
elt.mat["X.NH4",c("N","H","e")] <- c(1,4,-1)

elt.mat["N2",c("N")] <- c(2)

#---------- Mn ------------------------------------

elt.mat["MnO2.HR",c("Mn","O")] <- c(1,2)
elt.mat["MnO2.MR",c("Mn","O")] <- c(1,2)

elt.mat["Mn",c("Mn","e")] <- c(1,-2)

#---------- Fe ------------------------------------

elt.mat["FeOOH.HR",     c("Fe","O","H")]   <- c(1,2,1)
elt.mat["Fe56_FeOOH.HR",c("56Fe","O","H")] <- c(1,2,1)
elt.mat["FeOOH.MR",     c("Fe","O","H")]   <- c(1,2,1)
elt.mat["Fe56_FeOOH.MR",c("56Fe","O","H")] <- c(1,2,1)
elt.mat["FeOOH.PR",     c("Fe","O","H")]   <- c(1,2,1)
elt.mat["Fe56_FeOOH.PR",c("56Fe","O","H")] <- c(1,2,1)
elt.mat["FeOOH.U",      c("Fe","O","H")]   <- c(1,2,1)
elt.mat["Fe56_FeOOH.U", c("56Fe","O","H")] <- c(1,2,1)

elt.mat["Fe",       c("Fe","e")]   <- c(1,-2)
elt.mat["Fe56_Fe",  c("56Fe","e")] <- c(1,-2)
elt.mat["X.Fe",     c("Fe","e")]   <- c(1,-2)
elt.mat["X.Fe56_Fe",c("56Fe","e")] <- c(1,-2)

elt.mat["FeS",      c("S","Fe")]   <- c(1,1)
elt.mat["Fe56_FeS", c("S","56Fe")] <- c(1,1)
elt.mat["FeS2",     c("S","Fe")]   <- c(2,1)
elt.mat["Fe56_FeS2",c("S","56Fe")] <- c(2,1)

#---------- S -------------------------------------

elt.mat["SO4",c("S","O","e")] <- c(1,4,2)
elt.mat["HS", c("H","S","e")] <- c(1,1,1)
elt.mat["S0", c("S")]         <- c(1)

#---------- other----------------------------------

elt.mat["H2", c("H")]     <- c(2)
elt.mat["H",  c("H","e")] <- c(1,-1)
elt.mat["H2O",c("H","O")] <- c(2,1)

elt.table <- as.data.frame(elt.mat)
remove(elt.mat)

################################################################################
# initialise.parameters function
# This function initialises a default set of parameters (baseline simulation),
# but this function is also used when updating the parameter set.
################################################################################

initialise.parameters <- function(PL=NULL)
  
  # Model parameters are intialised on the first call (PL = NULL), and updated 
  # when a parmeter list is supplied (PL not NULL)
  
  # begin initialise.parameters
{
  
  #=============================================================================
  # Physical units 
  #=============================================================================
  
  if (is.null(PL$units$M)) PL$units$M <- "umol" # Mass unit
  if (is.null(PL$units$T)) PL$units$T <- "yr"   # Time unit
  if (is.null(PL$units$L)) PL$units$L <- "cm"   # Length unit
  
  # Limiting concentration: numerical limit on consumption in kinetic reactions
  PL$C.lim <- 1E-3
  
  #=============================================================================
  # Model domain and grid definition
  #=============================================================================
  
  if (is.null (PL$L)) PL$L <- 30#150   # depth of sediment domain [cm]
  if (is.null (PL$N)) PL$N <- 200  # number of grid layers
  
  # Call "setup.grid.1D" from ReacTran
  if (is.null (PL$grid.type)) PL$grid.type <- "uniform"  # number of grid layers
  if (is.null (PL$dx.1))      PL$dx.1      <- PL$L/(10*PL$N)
  #if (PL$grid.type == "uniform")
  #{
    # even grid - all depth layers have same thickness 
    #PL$grid <- setup.grid.1D(x.up=0,L=PL$L,N=PL$N)
  #} else {
    # uneven grid - higher resolution near the sediment water interface
    PL$grid <- setup.grid.1D(x.up=0,L=PL$L,N=PL$N,dx.1=PL$dx.1)
  #}
  
  #=============================================================================
  # Sediment parameters: 
  #=============================================================================
  
  # Environmental parameters
  
  if (is.null(PL$S))  PL$S  <- 35  # salinity
  if (is.null(PL$TC)) PL$TC <- 6   # temperature [deg C]
  if (is.null(PL$P))  PL$P  <- 36  # pressure [bar] (350m = 35 bar + 1 bar atmosphere)
    
  # Porosity profile 
  
  if (is.null(PL$por.0))     PL$por.0     <- 0.9 # porosity at the sediment-water interface
  if (is.null(PL$por.inf))   PL$por.inf   <- 0.7 # asymptotic porosity at depth
  if (is.null(PL$por.x.att)) PL$por.x.att <- 10  # attenuation depth [cm]
    
  PL$por.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$por.0,y.inf=PL$por.inf,x.L=0,x.att=PL$por.x.att)
  PL$svf.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=(1-PL$por.0),y.inf=(1-PL$por.inf),x.L=0,x.att=PL$por.x.att)
  
  # Initialisation tortuosity 

  PL$tort.grid <- setup.prop.1D(value=0,grid=PL$grid)
  PL$tort.grid$mid <- tortuosity(PL$por.grid$mid)
  PL$tort.grid$int <- tortuosity(PL$por.grid$int)

  # Fixed profiles: pH
  
  if (is.null(PL$pH)) PL$pH      <- 8.0
  PL$pH.grid <-  setup.prop.1D(value=PL$pH,grid=PL$grid)
  
  #=============================================================================
  # Diffusion coefficients 
  # Uses routine 'diffcoeff' from 'marelac' package
  #=============================================================================
  
  c.fac <- 10000*(3600*24*365.25) # conversion from [m2 s-1] to [cm2 yr-1]
  
  # Diffusion coefficients from marelac package
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
  
  #=============================================================================
  # Transport parameters
  #=============================================================================
  
  # Diffusive boundary layer
  
  if (is.null (PL$x.DBL)) PL$x.DBL <- 0 # thickness of DBL [cm]
  
  # Advective velocities 
  
  if (is.null(PL$rho.sed)) PL$rho.sed <- 2.6  # density solid sediment [g cm-3]
  if (is.null(PL$u.inf))   PL$u.inf   <- 0.06  # sedimentation velocity pore water [cm yr-1]
  if (is.null(PL$v.inf))   PL$v.inf   <- 0.06  # sedimentation velocity solids [cm yr-1]
  
  #PL$SedFlux <- 0.1  # sedimentation flux [g cm-2 yr-1]
  
  PL$v.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$v # advection velocity solid
  PL$u.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$u # advection velocity pore water
  
  # Bioturbation profile
  
  p.Dale <- function(x, Db.0 = 0., x.att = 1){
    return(Db.0*exp(-(x^2/(2*x.att^2))))
  }
  
  if (is.null(PL$x.L))  PL$x.L    <- 0.0
  if (is.null(PL$x.att))PL$x.att  <- 3.0
  if (is.null(PL$y.inf))PL$y.inf  <- 0.0
  if (is.null(PL$Db.0)) PL$Db.0   <- 23.0     # biodiffusion coefficient Db at SWI [cm2 yr-1]
  #PL$Db.grid <- setup.prop.1D(func=p.sig,y.0=PL$Db.0,y.inf=PL$y.inf,x.att=PL$x.att,x.L=PL$x.L, grid = PL$grid)
  PL$Db.grid <- setup.prop.1D(func=p.Dale, Db.0=PL$Db.0, x.att=PL$x.att, grid = PL$grid)
  
  # Irrigation profile
  
  if (is.null(PL$irr.0))   PL$irr.0   <- 290.0   # irrigation rate at SWI [yr-1]
  if (is.null(PL$x.irr))   PL$x.irr   <- 0.0     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.at))  PL$irr.att <- 1.4 
  if (is.null(PL$irr.inf)) PL$irr.inf <- 0.0
  
  if (is.null(PL$irr.HCO3)) PL$irr.HCO3 <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.NH4))  PL$irr.NH4  <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.O2))   PL$irr.O2   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.NO3))  PL$irr.NO3  <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.NO2))  PL$irr.NO2  <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.SO4))  PL$irr.SO4  <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.N2))   PL$irr.N2   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.Mn))   PL$irr.Mn   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.Fe))   PL$irr.Fe   <- 0.2   # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.HS))   PL$irr.HS   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.CH4))  PL$irr.CH4  <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.H2))   PL$irr.H2   <- 1     # irrigation rate at SWI [yr-1]
 
  if (is.null(PL$irr.H))    PL$irr.H    <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.H2O))  PL$irr.H2O  <- 1     # irrigation rate at SWI [yr-1]
  
  PL$irr.grid <- setup.prop.1D(func=p.exp,y.0=PL$irr.0,y.inf=PL$irr.inf,x.L=PL$x.irr,x.att=PL$irr.att ,grid=PL$grid)
  #PL$irr.grid <- setup.prop.1D(func=p.sig,y.0=PL$irr.0,y.inf=PL$irr.inf,x.att=PL$irr.att,x.L=PL$x.irr, grid = PL$grid)
  
  #=============================================================================
  # Reaction parameters
  #=============================================================================
  #---------- Decay constants organic matter (from Dale et al., 2015) ----------
  
  if (is.null(PL$k.1))    PL$k.1   <- 1e-10     # decay constant of organic matter [yr-1]
  if (is.null(PL$k.2))    PL$k.2   <- 3.16e-10  # decay constant of organic matter [yr-1]
  if (is.null(PL$k.3))    PL$k.3   <- 3.16e-9   # decay constant of organic matter [yr-1]
  if (is.null(PL$k.4))    PL$k.4   <- 3.16e-8   # decay constant of organic matter [yr-1] 
  if (is.null(PL$k.5))    PL$k.5   <- 3.16e-7   # decay constant of organic matter [yr-1] 
  if (is.null(PL$k.6))    PL$k.6   <- 3.16e-6   # decay constant of organic matter [yr-1] 
  if (is.null(PL$k.7))    PL$k.7   <- 3.16e-5   # decay constant of organic matter [yr-1] 
  if (is.null(PL$k.8))    PL$k.8   <- 3.16e-4   # decay constant of organic matter [yr-1] 
  if (is.null(PL$k.9))    PL$k.9   <- 3.16e-3   # decay constant of organic matter [yr-1] 
  if (is.null(PL$k.10))   PL$k.10  <- 3.16e-2   # decay constant of organic matter [yr-1] 
  if (is.null(PL$k.11))   PL$k.11  <- 3.16e-1   # decay constant of organic matter [yr-1] 
  if (is.null(PL$k.12))   PL$k.12  <- 3.16e-0   # decay constant of organic matter [yr-1] 
  if (is.null(PL$k.13))   PL$k.13  <- 3.16e+1   # decay constant of organic matter [yr-1] 
  if (is.null(PL$k.14))   PL$k.14  <- 100       # decay constant of organic matter [yr-1] 
  
  if (is.null(PL$k.H2))   PL$k.H2  <- 1E+03     # decay constant of hydrogen [umol-1 dm3 yr-1] 
  
  if (is.null(PL$K_O2))   PL$K_O2    <- 0.001            # Monod constant O2 consumption [umol cm-3 or mM] (Meysman et al., 2003)
  if (is.null(PL$K_NO3))  PL$K_NO3   <- 0.010            # Monod constant NO3 reduction [umol cm-3 or mM] (Van Cappellen and Wang, 1996)
  if (is.null(PL$K_NO2))  PL$K_NO2   <- 0.010            # Monod constant NO2 reduction [umol cm-3 or mM] (Dale et al., 2015)
  if (is.null(PL$K_MnO2)) PL$K_MnO2  <- 8.0*PL$rho.sed   # Monod constant MnIV reduction [umol cm-3 or mM] (???Van Cappellen and Wang, 1996)
  if (is.null(PL$K_FeOOH))PL$K_FeOOH <- 100.0*PL$rho.sed # Monod constant FeIII reduction [umol cm-3 or mM] (???Dale et al., 2015)
  if (is.null(PL$K_SO4))  PL$K_SO4   <- 0.5              # Monod constant SO4 reduction [umol cm-3 or mM] (Dale et al., 2015)
  
  if (is.null(PL$CN.ratio))  PL$CN.ratio <- 106/16 # CN ratio organic matter (Redfield, 1984)
  
  #---------- aerobic oxidation reactions ----------------------------------------
  
  if (is.null(PL$k.NIT.1)) PL$k.NIT.1 <- 1E+04  # kinetic constant ammonium oxidation [umol-1 cm3 yr-1] (Dale et al., 2015)
  if (is.null(PL$k.NIT.2)) PL$k.NIT.2 <- 1E+04  # kinetic constant nitrite oxidation [umol-1 cm3 yr-1] (Dale et al., 2015)
  if (is.null(PL$k.MnO))   PL$k.MnO   <- 5E+03  # kinetic constant MnII oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.FIO))   PL$k.FIO   <- 5E+05  # kinetic constant Ferrous iron oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.SIO))   PL$k.SIO   <- 5E+05  # kinetic constant Ferrous iron oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.CSO))   PL$k.CSO   <- 1E+02  # kinetic constant sulfide oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.AMO))   PL$k.AMO   <- 1E+02  # kinetic constant aerobic methane oxidation [µmol-1 cm3 yr-1] (Meysman et al., 2003)

  if (is.null(PL$k.AIO))   PL$k.AIO   <- 1E+02  # kinetic constant iron sulfide oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996)
  if (is.null(PL$k.PyO))   PL$k.PyO   <- 1E+00  # kinetic constant Pyrite oxidation [umol-1 cm3 yr-1] (Berg et al., 2003)
  
  #---------- anaerobic oxidation reactions --------------------------------------
  
  if (is.null(PL$k.ANA)) PL$k.ANA <- 1E+05  # kinetic constant ammonium oxidation [umol-1 cm3 yr-1] (Dale et al., 2015)
  
  if (is.null(PL$k.NFO))   PL$k.NFO   <- 1E+02   # kinetic constant ferrous iron oxidation by NO3- [umol-1 cm3 yr-1] (Dale et al., 2015)
  if (is.null(PL$k.NSO))   PL$k.NSO   <- 0.0     # kinetic constant sulphide oxidation by NO3- [umol-1 cm3 yr-1] (Dale et al., 2015)
  if (is.null(PL$k.MFOHR)) PL$k.MFOHR <- 1E+04   # kinetic constant ferrous iron oxidation by fresh MnIV [umol-1 cm3 yr-1] (Dale et al., 2015)
  if (is.null(PL$k.MFOMR)) PL$k.MFOMR <- 1E+02   # kinetic constant ferrous iron oxidation by aged MnIV [umol-1 cm3 yr-1] (Dale et al., 2015)
  if (is.null(PL$k.MSOHR)) PL$k.MSOHR <- 1E+02   # kinetic constant sulphide oxidation by fresh MnIV [umol-1 cm3 yr-1] (Dale et al., 2015)
  if (is.null(PL$k.MSOMR)) PL$k.MSOMR <- 1E+00   # kinetic constant sulphide oxidation by aged MnIV [umol-1 cm3 yr-1] (Dale et al., 2015)
 
  if (is.null(PL$k.ISOHR)) PL$k.ISOHR <- 3.   # kinetic constant sulfide oxidation by fresh FeIII [umol-1 cm3 yr-1] (Poulton et al., 2004)
  if (is.null(PL$k.ISOMR)) PL$k.ISOMR <- 3E-3 # kinetic constant sulfide oxidation by aged FeIII [umol-1 cm3 yr-1] (Poulton et al., 2004)
  if (is.null(PL$k.ISOPR)) PL$k.ISOPR <- 1E-5 # kinetic constant sulfide oxidation by aged FeIII [umol-1 cm3 yr-1] (Poulton et al., 2004)
  if (is.null(PL$k.ISOU))  PL$k.ISOU  <- 0.   # kinetic constant sulfide oxidation by aged FeIII [umol-1 cm3 yr-1] (Poulton et al., 2004)
  
  if (is.null(PL$k.SMO))   PL$k.SMO   <- 1E+02   # kinetic constant Sulfate methane oxidation [µmol-1 cm3 yr-1] (-)
  
  #---------- mineral reactions ---------------------------------------------------
  
  if (is.null(PL$K.FIS)) PL$K.FIS <- PL$rho.sed*268   # Equilibrium constant Ferrous iron sorption [] (Berg et al., 2003)
  
  if (is.null(PL$k.AmS)) PL$k.AmS  <- 1E-04    # kinetic constant ammonium adsorption [umol-1 cm3 yr-1] (this paper)
  if (is.null(PL$K.NH4)) PL$K.NH4  <- 1.6*PL$rho.sed  # K ammonium adsorption [cm3 g-1] (Bohlen, 2011)
  if (is.null(PL$k.SDP)) PL$k.SDP  <- 1        # kinetic constant sulfur disproportionation [yr-1] (Dale et al., 2015)
  
  if (is.null(PL$k.ISP)) PL$k.ISP <- 1E+03   # kinetic constant FeS precipitation [umol-1 cm3 yr-1] (Meysman et al., 2015)
  if (is.null(PL$k.ISD)) PL$k.ISD <- 3       # kinetic constant FeS dissolution [yr-1] (Meysman et al., 2015)
  if (is.null(PL$n.ISP)) PL$n.ISP <- 1       # kinetic exponent FeS precipitation [] (Meysman et al., 2003)
  if (is.null(PL$n.ISD)) PL$n.ISD <- 1       # kinetic exponent FeS dissolution [] (Meysman et al., 2003)
  if (is.null(PL$K_IS))  PL$K_IS  <- (10^-PL$pH)*1000*3.16    # Saturation constant [umol-1 cm3] (Rickard, 2006)

  if (is.null(PL$k.PyP.1))  PL$k.PyP.1  <- 3.25    # kinetic constant Pyrite precipitation [umol-1 cm3 yr-1] (vdV et al., 2020 + Rickard)
  if (is.null(PL$k.PyP.2))  PL$k.PyP.2  <- 3.25    # kinetic constant Pyrite precipitation [umol-1 cm3 yr-1] (vdV et al., 2020 + Rickard)
  
  #---------- metal oxide ageing -------------------------------------------------- 
  
  if (is.null(PL$k.IOA)) PL$k.IOA <- 0.7  # kinetic constant iron oxide ageing [yr-1] (Berg et al., 2003)
  if (is.null(PL$k.MOA)) PL$k.MOA <- 1.7  # kinetic constant iron oxide ageing [yr-1] (Berg et al., 2003)
  
  #=============================================================================
  # Iron isotope system
  #=============================================================================
  
  if (is.null(PL$IRMM014))        PL$IRMM014        <- 15.69786  # Fe isotope reference (Dauphas et al., Reviews in Mineralogy & Geochemistry, 2017)
  if (is.null(PL$frac.limit))     PL$frac.limit     <- 1e-6      # concentration below which no fractionation happens
  if (is.null(PL$d56Fe_Fe))       PL$d56Fe_Fe       <- 0.0       # overlying water Fe isotope signature
  if (is.null(PL$d56Fe_FeOOH.HR)) PL$d56Fe_FeOOH.HR <- 0.0       # influx FeOOH.HR isotope signature
  if (is.null(PL$d56Fe_FeOOH.MR)) PL$d56Fe_FeOOH.MR <- 0.0       # influx FeOOH.MR isotope signature
  if (is.null(PL$d56Fe_FeOOH.PR)) PL$d56Fe_FeOOH.PR <- 0.0       # influx FeOOH.PR isotope signature
  if (is.null(PL$d56Fe_FeOOH.U))  PL$d56Fe_FeOOH.U  <- 0.0       # influx FeOOH.U isotope signature
  if (is.null(PL$d56Fe_FeS))      PL$d56Fe_FeS      <- 0.0       # influx FeS isotope signature
  if (is.null(PL$d56Fe_FeS2))     PL$d56Fe_FeS2     <- 0.0       # influx FeS2 isotope signature
  if (is.null(PL$d56Fe_X.Fe))     PL$d56Fe_X.Fe     <- 0.0       # influx sorbed Fe isotope signature
  
  if (is.null(PL$a56Fe_DIR))     PL$a56Fe_DIR     <- (1.0 - 1.3/1000.0) # alpha for Dissimilatory Iron Reduction
  if (is.null(PL$a56Fe_ISO))     PL$a56Fe_ISO     <- (1.0 - 1.3/1000.0) # alpha for Dissimilatory Iron Reduction
  if (is.null(PL$a56Fe_FIO))     PL$a56Fe_FIO     <- (1.0 + 0.4/1000.0) # alpha for Ferrous Iron oxidation
  if (is.null(PL$a56Fe_ISP))     PL$a56Fe_ISP     <- (1.0 + 0.5/1000.0) # alpha for Iron sulphide precipitation
  if (is.null(PL$a56Fe_AIO))     PL$a56Fe_AIO     <- (1.0 + 0.4/1000.0) # alpha for Ferrous Iron oxidation
  if (is.null(PL$a56Fe_ISD))     PL$a56Fe_ISD     <- (1.0 - 0.5/1000.0) # alpha for Ferrous Iron oxidation
  if (is.null(PL$a56Fe_PyP))     PL$a56Fe_PyP     <- (1.0 - 0.7/1000.0) # alpha for Pyrite precipitation
  if (is.null(PL$a56Fe_PyO))     PL$a56Fe_PyO     <- (1.0 + 0.5/1000.0) # alpha for Pyrite precipitation
  if (is.null(PL$a56Fe_IOA))     PL$a56Fe_IOA     <- (1.0 + 0.0/1000.0) # alpha for Dissimilatory Iron Reduction
  if (is.null(PL$a56Fe_FIS))     PL$a56Fe_FIS     <- (1.0 + 0.4/1000.0) # alpha for Ferrous Iron Sorption
  
  #=============================================================================
  # Boundary conditions
  #=============================================================================
  #---------- Upper boundary fluxes --------------------------------------------
  
  F.POC <- 15 # mmol m-2 d-1
  if (is.null(PL$F.POC)) PL$F.POC <- F.POC*365.25/10      # conversion to umol cm-2 yr-1 
  
  # deconvolution of POC flux in fraction following Dale et al., 2015
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
  
  F.X.NH4 <- 0 # mmol m-2 d-1
  if (is.null(PL$F.X.NH4)) PL$F.X.NH4 <- F.X.NH4*365.25/10         # conversion to umol cm-2 yr-1 
  
  F.MnO2 <- 108e-3 # mmol m-2 d-1
  if (is.null(PL$F.MnO2))    PL$F.MnO2    <- F.MnO2*365.25/10      # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.MnO2.HR)) PL$F.MnO2.HR <- PL$F.MnO2/2
  if (is.null(PL$F.MnO2.MR)) PL$F.MnO2.MR <- PL$F.MnO2/2
  
  F.FeOOH <- 1110e-3 # mmol m-2 d-1
  if (is.null(PL$F.FeOOH))         PL$F.FeOOH         <- F.FeOOH*365.25/10 # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.FeOOH.HR))      PL$F.FeOOH.HR      <- PL$F.FeOOH/6      # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.Fe56_FeOOH.HR)) PL$F.Fe56_FeOOH.HR <- ((PL$d56Fe_FeOOH.HR*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.HR*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.HR # conversion to amount of 56Fe 
  if (is.null(PL$F.FeOOH.MR))      PL$F.FeOOH.MR      <- PL$F.FeOOH/6      # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.Fe56_FeOOH.MR)) PL$F.Fe56_FeOOH.MR <- ((PL$d56Fe_FeOOH.MR*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.MR*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.MR # conversion to amount of 56Fe 
  if (is.null(PL$F.FeOOH.PR))      PL$F.FeOOH.PR      <- PL$F.FeOOH/6      # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.Fe56_FeOOH.PR)) PL$F.Fe56_FeOOH.PR <- ((PL$d56Fe_FeOOH.PR*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.PR*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.PR # conversion to amount of 56Fe 
  if (is.null(PL$F.FeOOH.U))       PL$F.FeOOH.U       <- PL$F.FeOOH/2      # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.Fe56_FeOOH.U))  PL$F.Fe56_FeOOH.U  <- ((PL$d56Fe_FeOOH.U*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.U*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.U # conversion to amount of 56Fe 
  
  F.X.Fe <- 0 # mmol m-2 d-1
  if (is.null(PL$F.X.Fe))      PL$F.X.Fe      <- F.X.Fe*365.25/10         # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.X.Fe56_Fe)) PL$F.X.Fe56_Fe <- ((PL$d56Fe_X.Fe*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_X.Fe*1e-3+1.0)*PL$IRMM014))*PL$F.X.Fe  # conversion to amount of 56Fe 
  F.FeS <- 0 # mmol m-2 d-1
  if (is.null(PL$F.FeS))       PL$F.FeS       <- F.FeS*365.25/10       # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.Fe56_FeS))  PL$F.Fe56_FeS  <- ((PL$d56Fe_FeS*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeS*1e-3+1.0)*PL$IRMM014))*PL$F.FeS # conversion to amount of 56Fe 
  F.FeS2 <- 0 # mmol m-2 d-1
  if (is.null(PL$F.FeS2))      PL$F.FeS2      <- F.FeS2*365.25/10      # conversion to umol cm-2 yr-1
  if (is.null(PL$F.Fe56_FeS2)) PL$F.Fe56_FeS2 <- ((PL$d56Fe_FeS2*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeS2*1e-3+1.0)*PL$IRMM014))*PL$F.FeS2 # conversion to amount of 56Fe 
  
  F.S0 <- 0 # mmol m-2 d-1
  if (is.null(PL$F.FeS0)) PL$F.S0 <- F.S0*365.25/10      # conversion to umol cm-2 yr-1
  
  #---------- Upper boundary concentrations ----------------------------------------
  
  if (is.null(PL$HCO3.ow))     PL$HCO3.ow    <- 2.2  # SumCO2 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$CH4.ow))      PL$CH4.ow     <- 0.0
  
  if (is.null(PL$O2.ow))       PL$O2.ow      <- 0.28  # O2 concentration bottom water [umol cm-3 or mM]
  
  if (is.null(PL$NH4.ow))      PL$NH4.ow     <- 0.001 # NH4 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$NO3.ow))      PL$NO3.ow     <- 0.035 # NO3 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$NO2.ow))      PL$NO2.ow     <- 0.0   # NO3 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$N2.ow))       PL$N2.ow      <- 0.0
  
  if (is.null(PL$Mn.ow))       PL$Mn.ow      <- 0.0
  
  if (is.null(PL$Fe.ow))       PL$Fe.ow      <- 0.0
  if (is.null(PL$Fe56_Fe.ow))  PL$Fe56_Fe.ow <- ((PL$d56Fe_Fe*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_Fe*1e-3+1.0)*PL$IRMM014))*PL$Fe.ow
  
  if (is.null(PL$SO4.ow))      PL$SO4.ow     <- 28.2  # SO4 concentration bottom water
  if (is.null(PL$HS.ow))       PL$HS.ow      <- 0.0   # SumH2S concentration bottom water [umol cm-3 or mM]
  
  if (is.null(PL$H2O.ow))      PL$H2O.ow     <- 0.0
  if (is.null(PL$H.ow))        PL$H.ow       <- 0.0 # Proton concentration bottom water, without equilibrium with H2O
  if (is.null(PL$H2.ow))       PL$H2.ow      <- 0.0 # Proton concentration bottom water, without equilibrium with H2O
  
  #---------- Lower boundary conditions ----------------------------------------
  
  # All species: no gradient
  
  if (is.null(PL$HCO3.ds))     PL$HCO3.ds     <- NA
  if (is.null(PL$CH4.ds))      PL$CH4.ds      <- NA
  
  if (is.null(PL$O2.ds))       PL$O2.ds       <- NA
  
  if (is.null(PL$NH4.ds))      PL$NH4.ds      <- NA
  if (is.null(PL$NO3.ds))      PL$NO3.ds      <- NA
  if (is.null(PL$NO2.ds))      PL$NO2.ds      <- NA
  if (is.null(PL$N2.ds))       PL$N2.ds       <- NA
  
  if (is.null(PL$Mn.ds))       PL$Mn.ds       <- NA
  
  if (is.null(PL$Fe.ds))       PL$Fe.ds       <- NA
  if (is.null(PL$Fe56_Fe.ds))  PL$Fe56_Fe.ds  <- NA
  
  if (is.null(PL$SO4.ds))      PL$SO4.ds      <- NA
  if (is.null(PL$HS.ds))       PL$HS.ds       <- NA
  
  if (is.null(PL$H2O.ds))      PL$H2O.ds      <- NA
  if (is.null(PL$H.ds))        PL$H.ds        <- NA
  if (is.null(PL$H2.ds))       PL$H2.ds       <- NA
  
  #=============================================================================
  # Flags 
  #=============================================================================
  
  if (is.null(PL$simulation.type)) PL$simulation.type <- "direct.steady.state"
  #simulation.type <- "dynamic.steady.state"
  #simulation.type <- "time.dependent"
  
  ##############################################################################
  return(PL)
  # end initialise.parameters
}
##############################################################################

#=============================================================================
# Model formulation
#=============================================================================

CSFe.model <- function (t,state,parameters,summary.call=FALSE) 
{
  
  with(as.list(c(parameters)),{
    
    #---------------------------------------------------------------------------
    # Initialisation of 
    # SV : matrix of state variables
    # OL : Object List
    #---------------------------------------------------------------------------
    
    SV <- matrix(nrow=N,ncol=N.var,data=state)
    for (i in 1:N.var) OL[[var.names[i]]]$C <- SV[,i]
    
    #---------------------------------------------------------------------------
    OL <- within(OL,{
    #---------------------------------------------------------------------------
      
    #=========================================================================
    # Initialisation of concentration depth profiles
    #=========================================================================
    
    #-------------------------------------------------------------------------
    # Initialisation of species with fixed concentration depth profiles 
    #-------------------------------------------------------------------------
    
    X.Fe$C <- K.FIS*Fe$C
    
    #=========================================================================
    # Isotope ratios 
    #=========================================================================
    
    R_56FeOOH.HR                                            <- Fe56_FeOOH.HR$C/(FeOOH.HR$C-Fe56_FeOOH.HR$C)
    R_56FeOOH.HR[(FeOOH.HR$C-Fe56_FeOOH.HR$C)<(C.lim*1e-6)] <- 0.0
    R_56FeOOH.HR[FeOOH.HR$C<(C.lim*1e-6)]                   <- 0.0
    R_56FeOOH.HR[Fe56_FeOOH.HR$C<(C.lim*1e-6)]              <- 0.0
    R_56FeOOH.HR[R_56FeOOH.HR<0.0]                          <- 0.0
    r_56FeOOH.HR                                            <- Fe56_FeOOH.HR$C/FeOOH.HR$C
    r_56FeOOH.HR[FeOOH.HR$C<(C.lim*1e-6)]                   <- 0.0
    r_56FeOOH.HR[Fe56_FeOOH.HR$C<(C.lim*1e-6)]              <- 0.0
    r_56FeOOH.HR[r_56FeOOH.HR<0.0]                          <- 0.0
    
    R_56FeOOH.MR                                            <- Fe56_FeOOH.MR$C/(FeOOH.MR$C-Fe56_FeOOH.MR$C)
    R_56FeOOH.MR[(FeOOH.MR$C-Fe56_FeOOH.MR$C)<(C.lim*1e-6)] <- 0.0
    R_56FeOOH.MR[(FeOOH.MR$C)<(C.lim*1e-6)]                 <- 0.0
    R_56FeOOH.MR[(Fe56_FeOOH.MR$C)<(C.lim*1e-6)]            <- 0.0
    R_56FeOOH.MR[R_56FeOOH.MR<0.0]                          <- 0.0
    r_56FeOOH.MR                                            <- Fe56_FeOOH.MR$C/FeOOH.MR$C
    r_56FeOOH.MR[FeOOH.MR$C<(C.lim*1e-6)]                   <- 0.0
    r_56FeOOH.MR[Fe56_FeOOH.MR$C<(C.lim*1e-6)]              <- 0.0
    r_56FeOOH.MR[r_56FeOOH.MR<0.0]                          <- 0.0
    
    R_56FeOOH.PR                                            <- Fe56_FeOOH.PR$C/(FeOOH.PR$C-Fe56_FeOOH.PR$C)
    R_56FeOOH.PR[(FeOOH.PR$C-Fe56_FeOOH.PR$C)<(C.lim*1e-6)] <- 0.0
    R_56FeOOH.PR[(FeOOH.PR$C)<(C.lim*1e-6)]                 <- 0.0
    R_56FeOOH.PR[(Fe56_FeOOH.PR$C)<(C.lim*1e-6)]            <- 0.0
    R_56FeOOH.PR[R_56FeOOH.PR<0.0]                          <- 0.0
    r_56FeOOH.PR                                            <- Fe56_FeOOH.PR$C/FeOOH.PR$C
    r_56FeOOH.PR[FeOOH.PR$C<(C.lim*1e-6)]                   <- 0.0
    r_56FeOOH.PR[Fe56_FeOOH.PR$C<(C.lim*1e-6)]              <- 0.0
    r_56FeOOH.PR[r_56FeOOH.PR<0.0]                          <- 0.0
    
    R_56FeOOH.U                                             <- Fe56_FeOOH.U$C/(FeOOH.U$C-Fe56_FeOOH.U$C)
    R_56FeOOH.U[(FeOOH.U$C-Fe56_FeOOH.U$C)<(C.lim*1e-6)]    <- 0.0
    R_56FeOOH.U[(FeOOH.U$C)<(C.lim*1e-6)]                   <- 0.0
    R_56FeOOH.U[(Fe56_FeOOH.U$C)<(C.lim*1e-6)]              <- 0.0
    R_56FeOOH.U[R_56FeOOH.U<0.0]                            <- 0.0
    r_56FeOOH.U                                             <- Fe56_FeOOH.U$C/FeOOH.U$C
    r_56FeOOH.U[FeOOH.U$C<(C.lim*1e-6)]                     <- 0.0
    r_56FeOOH.U[Fe56_FeOOH.U$C<(C.lim*1e-6)]                <- 0.0
    r_56FeOOH.U[r_56FeOOH.U<0.0]                            <- 0.0
    
    R_56FeS                                  <- Fe56_FeS$C/(FeS$C-Fe56_FeS$C)
    R_56FeS[(FeS$C-Fe56_FeS$C)<(C.lim*1e-6)] <- 0.0
    R_56FeS[(FeS$C)<(C.lim*1e-6)]            <- 0.0
    R_56FeS[(Fe56_FeS$C)<(C.lim*1e-6)]       <- 0.0
    R_56FeS[R_56FeS<0.0]                     <- 0.0
    r_56FeS                                  <- Fe56_FeS$C/FeS$C
    r_56FeS[FeS$C<(C.lim*1e-6)]              <- 0.0
    r_56FeS[Fe56_FeS$C<(C.lim*1e-6)]         <- 0.0
    r_56FeS[r_56FeS<0.0]                     <- 0.0
    
    R_56FeS2                                    <- Fe56_FeS2$C/(FeS2$C-Fe56_FeS2$C)
    R_56FeS2[(FeS2$C-Fe56_FeS2$C)<(C.lim*1e-6)] <- 0.0
    R_56FeS2[(FeS2$C)<(C.lim*1e-6)]             <- 0.0
    R_56FeS2[(Fe56_FeS2$C)<(C.lim*1e-6)]        <- 0.0
    R_56FeS2[R_56FeS2<0.0]                      <- 0.0
    r_56FeS2                                    <- Fe56_FeS2$C/FeS2$C
    r_56FeS2[FeS2$C<(C.lim*1e-6)]               <- 0.0
    r_56FeS2[Fe56_FeS2$C<(C.lim*1e-6)]          <- 0.0
    r_56FeS2[r_56FeS2<0.0]                      <- 0.0
    
    R_56Fe                                <- Fe56_Fe$C/(Fe$C-Fe56_Fe$C)
    R_56Fe[(Fe$C-Fe56_Fe$C)<(C.lim*1e-6)] <- 0.0
    R_56Fe[(Fe$C)<(C.lim*1e-6)]           <- 0.0
    R_56Fe[(Fe56_Fe$C)<(C.lim*1e-6)]      <- 0.0
    R_56Fe[R_56Fe<0.0]                    <- 0.0
    r_56Fe                                <- Fe56_Fe$C/Fe$C
    r_56Fe[Fe$C<(C.lim*1e-6)]             <- 0.0
    r_56Fe[Fe56_Fe$C<(C.lim*1e-6)]        <- 0.0
    r_56Fe[r_56Fe<0.0]                    <- 0.0
    
    X.Fe56_Fe$C[Fe$C>=frac.limit] <- K.FIS*((a56Fe_FIS+a56Fe_FIS*R_56Fe[Fe$C>=frac.limit])/(1.0+a56Fe_FIS*R_56Fe[Fe$C>=frac.limit]))*Fe56_Fe$C[Fe$C>=frac.limit]
    X.Fe56_Fe$C[Fe$C<frac.limit]  <- K.FIS*Fe56_Fe$C[Fe$C<frac.limit]
    
    R_56XFe                                    <- X.Fe56_Fe$C/(X.Fe$C-X.Fe56_Fe$C)
    R_56XFe[(X.Fe$C-X.Fe56_Fe$C)<(C.lim*1e-6)] <- 0.0
    R_56XFe[(X.Fe$C)<(C.lim*1e-6)]             <- 0.0
    R_56XFe[(X.Fe56_Fe$C)<(C.lim*1e-6)]        <- 0.0
    R_56XFe[R_56XFe<0.0]                       <- 0.0
    r_56XFe                                    <- X.Fe56_Fe$C/X.Fe$C
    r_56XFe[X.Fe$C<(C.lim*1e-6)]               <- 0.0
    r_56XFe[X.Fe56_Fe$C<(C.lim*1e-6)]          <- 0.0
    r_56XFe[r_56XFe<0.0]                       <- 0.0
    
    #=========================================================================
    # Upper and lower boundary conditions
    #=========================================================================

    HCO3$C.up    <- HCO3.ow
    CH4$C.up     <- CH4.ow
    
    O2$C.up      <- O2.ow
    
    NH4$C.up     <- NH4.ow
    NO3$C.up     <- NO3.ow
    NO2$C.up     <- NO2.ow
    N2$C.up      <- N2.ow
    
    Mn$C.up      <- Mn.ow
    
    Fe$C.up      <- Fe.ow
    Fe56_Fe$C.up <- Fe56_Fe.ow
    
    SO4$C.up     <- SO4.ow
    HS$C.up      <- HS.ow
    
    H2O$C.up     <- H2O.ow
    H$C.up       <- H.ow
    H2$C.up      <- H2.ow
    
    if (is.na(HCO3.ds))     HCO3$C.down     <- HCO3$C[N]    else HCO3$C.down    <- HCO3.ds 
    if (is.na(CH4.ds))      CH4$C.down      <- CH4$C[N]     else CH4$C.down     <- CH4.ds 
    
    if (is.na(O2.ds))       O2$C.down       <- O2$C[N]      else O2$C.down      <- O2.ds 
    
    if (is.na(NH4.ds))      NH4$C.down      <- NH4$C[N]     else NH4$C.down     <- NH4.ds 
    if (is.na(NO3.ds))      NO3$C.down      <- NO3$C[N]     else NO3$C.down     <- NO3.ds 
    if (is.na(NO2.ds))      NO2$C.down      <- NO2$C[N]     else NO2$C.down     <- NO2.ds 
    if (is.na(N2.ds))       N2$C.down       <- N2$C[N]      else N2$C.down      <- N2.ds 
    
    if (is.na(Mn.ds))       Mn$C.down       <- Mn$C[N]      else Mn$C.down      <- Mn.ds 
    
    if (is.na(Fe.ds))       Fe$C.down       <- Fe$C[N]      else Fe$C.down      <- Fe.ds 
    if (is.na(Fe56_Fe.ds))  Fe56_Fe$C.down  <- Fe56_Fe$C[N] else Fe56_Fe$C.down <- Fe56_Fe.ds 
    
    if (is.na(SO4.ds))      SO4$C.down      <- SO4$C[N]     else SO4$C.down     <- SO4.ds 
    if (is.na(HS.ds))       HS$C.down       <- HS$C[N]      else HS$C.down      <- HS.ds 
    
    if (is.na(H2O.ds))      H2O$C.down      <- H2O$C[N]     else H2O$C.down     <- H2O.ds 
    if (is.na(H.ds))        H$C.down        <- H$C[N]       else H$C.down       <- H.ds 
    if (is.na(H2.ds))       H2$C.down       <- H2$C[N]      else H2$C.down      <- H2.ds 
    
    #=========================================================================
    # TRANSPORT terms using tran.1D from ReacTran
    #=========================================================================
    
    if (summary.call) {ind <- c("tran","C.up","C.down","dif.flux","adv.flux","flux","flux.up","flux.down")} else {ind <- c("tran","flux.up","flux.down")}
    
    #------------ Transport terms: Solids [umol cm-3 solid yr-1] -------------
    
    POC.1[ind]  <- tran.1D(C=POC.1$C, flux.up=F.POC.1,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.2[ind]  <- tran.1D(C=POC.2$C, flux.up=F.POC.2,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.3[ind]  <- tran.1D(C=POC.3$C, flux.up=F.POC.3,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.4[ind]  <- tran.1D(C=POC.4$C, flux.up=F.POC.4,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.5[ind]  <- tran.1D(C=POC.5$C, flux.up=F.POC.5,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.6[ind]  <- tran.1D(C=POC.6$C, flux.up=F.POC.6,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.7[ind]  <- tran.1D(C=POC.7$C, flux.up=F.POC.7,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.8[ind]  <- tran.1D(C=POC.8$C, flux.up=F.POC.8,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.9[ind]  <- tran.1D(C=POC.9$C, flux.up=F.POC.9,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.10[ind] <- tran.1D(C=POC.10$C,flux.up=F.POC.10,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.11[ind] <- tran.1D(C=POC.11$C,flux.up=F.POC.11,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.12[ind] <- tran.1D(C=POC.12$C,flux.up=F.POC.12,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.13[ind] <- tran.1D(C=POC.13$C,flux.up=F.POC.13,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    POC.14[ind] <- tran.1D(C=POC.14$C,flux.up=F.POC.14,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    X.NH4[ind]         <- tran.1D(C=X.NH4$C,        flux.up=F.X.NH4,        v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    MnO2.HR[ind]       <- tran.1D(C=MnO2.HR$C,      flux.up=F.MnO2.HR,      v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    MnO2.MR[ind]       <- tran.1D(C=MnO2.MR$C,      flux.up=F.MnO2.MR,      v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    FeOOH.HR[ind]      <- tran.1D(C=FeOOH.HR$C,     flux.up=F.FeOOH.HR,     v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    Fe56_FeOOH.HR[ind] <- tran.1D(C=Fe56_FeOOH.HR$C,flux.up=F.Fe56_FeOOH.HR,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeOOH.MR[ind]      <- tran.1D(C=FeOOH.MR$C,     flux.up=F.FeOOH.MR,     v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    Fe56_FeOOH.MR[ind] <- tran.1D(C=Fe56_FeOOH.MR$C,flux.up=F.Fe56_FeOOH.MR,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeOOH.PR[ind]      <- tran.1D(C=FeOOH.PR$C,     flux.up=F.FeOOH.PR,     v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    Fe56_FeOOH.PR[ind] <- tran.1D(C=Fe56_FeOOH.PR$C,flux.up=F.Fe56_FeOOH.PR,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeOOH.U[ind]       <- tran.1D(C=FeOOH.U$C,      flux.up=F.FeOOH.U,      v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    Fe56_FeOOH.U[ind]  <- tran.1D(C=Fe56_FeOOH.U$C, flux.up=F.Fe56_FeOOH.U, v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    X.Fe[ind]          <- tran.1D(C=X.Fe$C,         flux.up=F.X.Fe,         v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    X.Fe56_Fe[ind]     <- tran.1D(C=X.Fe56_Fe$C,    flux.up=F.X.Fe56_Fe,    v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    FeS[ind]           <- tran.1D(C=FeS$C,          flux.up=F.FeS,          v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    Fe56_FeS[ind]      <- tran.1D(C=Fe56_FeS$C,     flux.up=F.Fe56_FeS,     v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeS2[ind]          <- tran.1D(C=FeS2$C,         flux.up=F.FeS2,         v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    Fe56_FeS2[ind]     <- tran.1D(C=Fe56_FeS2$C,    flux.up=F.Fe56_FeS2,    v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    S0[ind]            <- tran.1D(C=S0$C,           flux.up=F.S0,           v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    #------------ Transport terms: Solutes [umol cm-3 pore water yr-1] -------------------------------------------------

    HCO3[ind]     <- tran.1D(C=HCO3$C,   C.up=HCO3$C.up,   C.down=HCO3$C.down,   v=u.grid,D=D.HCO3.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    CH4[ind]      <- tran.1D(C=CH4$C,    C.up=CH4$C.up,    C.down=CH4$C.down,    v=u.grid,D=D.CH4.grid, AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    O2[ind]       <- tran.1D(C=O2$C,     C.up=O2$C.up,     C.down=O2$C.down,     v=u.grid,D=D.O2.grid,  AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    NH4[ind]      <- tran.1D(C=NH4$C,    C.up=NH4$C.up,    C.down=NH4$C.down,    v=u.grid,D=D.NH4.grid, AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    NO3[ind]      <- tran.1D(C=NO3$C,    C.up=NO3$C.up,    C.down=NO3$C.down,    v=u.grid,D=D.NO3.grid, AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    NO2[ind]      <- tran.1D(C=NO2$C,    C.up=NO2$C.up,    C.down=NO2$C.down,    v=u.grid,D=D.NO2.grid, AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    N2[ind]       <- tran.1D(C=N2$C,     C.up=N2$C.up,     C.down=N2$C.down,     v=u.grid,D=D.N2.grid,  AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    Mn[ind]       <- tran.1D(C=Mn$C,     C.up=Mn$C.up,     C.down=Mn$C.down,     v=u.grid,D=D.Mn.grid,  AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    Fe[ind]       <- tran.1D(C=Fe$C,     C.up=Fe$C.up,     C.down=Fe$C.down,     v=u.grid,D=D.Fe.grid,  AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    Fe56_Fe[ind]  <- tran.1D(C=Fe56_Fe$C,C.up=Fe56_Fe$C.up,C.down=Fe56_Fe$C.down,v=u.grid,D=D.Fe.grid,  AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    SO4[ind]      <- tran.1D(C=SO4$C,    C.up=SO4$C.up,    C.down=SO4$C.down,    v=u.grid,D=D.SO4.grid, AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    HS[ind]       <- tran.1D(C=HS$C,     C.up=HS$C.up,     C.down=HS$C.down,     v=u.grid,D=D.HS.grid,  AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    H2O[ind]      <- tran.1D(C=H2O$C,    C.up=H2O$C.up,    C.down=H2O$C.down,    v=u.grid,D=D.H2O.grid, AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    H[ind]        <- tran.1D(C=H$C,      C.up=H$C.up,      C.down=H$C.down,      v=u.grid,D=D.H.grid,   AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    H2[ind]       <- tran.1D(C=H2$C,     C.up=H2$C.up,     C.down=H2$C.down,     v=u.grid,D=D.H.grid,   AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    #=========================================================================
    # Saturation state calculation
    #=========================================================================
    
    # Saturation state FeS 
    
    Omega.FeS$C      <- Fe$C*HS$C/(K_IS)
    Omega.FeS$C.up   <- Fe.ow*HS$C.up/K_IS
    Omega.FeS$C.down <- Fe$C.down*HS$C.down/K_IS
    
    Omega.NH4$C      <- X.NH4$C/(NH4$C*K.NH4)
    Omega.NH4$C.up   <- X.NH4$C.up/(NH4$C.up*K.NH4)
    Omega.NH4$C.down <- X.NH4$C.down/(NH4$C.down*K.NH4)
    
    #=========================================================================
    # Rates of Kinetic reactions 
    # Expressed per volume of bulk sediment [umol cm-3 yr-1]
    #=========================================================================
    #---------------- Mineralisation reactions -------------------------------
    
    Cmin.1$R  <- k.1*svf.grid$mid*POC.1$C*FSAT(POC.1$C,C.lim,5)
    Cmin.2$R  <- k.2*svf.grid$mid*POC.2$C*FSAT(POC.2$C,C.lim,5) 
    Cmin.3$R  <- k.3*svf.grid$mid*POC.3$C*FSAT(POC.3$C,C.lim,5)
    Cmin.4$R  <- k.4*svf.grid$mid*POC.4$C*FSAT(POC.4$C,C.lim,5)
    Cmin.5$R  <- k.5*svf.grid$mid*POC.5$C*FSAT(POC.5$C,C.lim,5)
    Cmin.6$R  <- k.6*svf.grid$mid*POC.6$C*FSAT(POC.6$C,C.lim,5)
    Cmin.7$R  <- k.7*svf.grid$mid*POC.7$C*FSAT(POC.7$C,C.lim,5)
    Cmin.8$R  <- k.8*svf.grid$mid*POC.8$C*FSAT(POC.8$C,C.lim,5)
    Cmin.9$R  <- k.9*svf.grid$mid*POC.9$C*FSAT(POC.9$C,C.lim,5)
    Cmin.10$R <- k.10*svf.grid$mid*POC.10$C*FSAT(POC.10$C,C.lim,5)
    Cmin.11$R <- k.11*svf.grid$mid*POC.11$C*FSAT(POC.11$C,C.lim,5)
    Cmin.12$R <- k.12*svf.grid$mid*POC.12$C*FSAT(POC.12$C,C.lim,5)
    Cmin.13$R <- k.13*svf.grid$mid*POC.13$C*FSAT(POC.13$C,C.lim,5)
    Cmin.14$R <- k.14*svf.grid$mid*POC.14$C*FSAT(POC.14$C,C.lim,5)
    
    Cmin$R  <- Cmin.1$R + Cmin.2$R + Cmin.3$R + Cmin.4$R + Cmin.5$R + Cmin.6$R + Cmin.7$R + Cmin.8$R + Cmin.9$R + Cmin.10$R + Cmin.11$R + Cmin.12$R + Cmin.13$R + Cmin.14$R  
    H2min$R <- k.H2*por.grid$mid*H2$C*FSAT(H2$C,C.lim/100,5)
    
    O2.lim    <- O2$C/(O2$C+K_O2)*FSAT(O2$C,C.lim,5) #*(O2$C>sqrt(.Machine$double.eps))
    O2.inh    <- K_O2/(O2$C+K_O2)
    NO3.lim   <- NO3$C/(NO3$C+K_NO3)*FSAT(NO3$C,C.lim,5)
    NO3.inh   <- K_NO3/(NO3$C+K_NO3)
    NO2.lim   <- NO2$C/(NO2$C+K_NO2)*FSAT(NO2$C,C.lim,5)
    NO2.inh   <- K_NO2/(NO2$C+K_NO2)
    MnO2.lim  <- MnO2.HR$C/(MnO2.HR$C+K_MnO2)*FSAT(MnO2.HR$C,C.lim,5)
    MnO2.inh  <- K_MnO2/(MnO2.HR$C+K_MnO2)
    FeOOH.lim <- FeOOH.HR$C/(FeOOH.HR$C+K_FeOOH)*FSAT(FeOOH.HR$C,C.lim,5)
    FeOOH.inh <- K_FeOOH/(FeOOH.HR$C+K_FeOOH)
    SO4.lim   <- SO4$C/(SO4$C+K_SO4)*FSAT(SO4$C,C.lim,5)
    SO4.inh   <- K_SO4/(SO4$C+K_SO4)
    
    a.f   <- O2.lim/(O2.lim + NO3.lim*O2.inh + NO2.lim*NO3.inh*O2.inh + MnO2.lim*O2.inh*NO3.inh*NO2.inh + FeOOH.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh + SO4.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh + O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh*SO4.inh)
    n1.f  <- (NO3.lim*O2.inh)/(O2.lim + NO3.lim*O2.inh + NO2.lim*NO3.inh*O2.inh + MnO2.lim*O2.inh*NO3.inh*NO2.inh + FeOOH.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh + SO4.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh + O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh*SO4.inh)
    n2.f  <- (NO2.lim*O2.inh*NO3.inh)/(O2.lim + NO3.lim*O2.inh + NO2.lim*NO3.inh*O2.inh + MnO2.lim*O2.inh*NO3.inh*NO2.inh + FeOOH.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh + SO4.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh + O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh*SO4.inh)
    mn.f  <- (MnO2.lim*O2.inh*NO3.inh*NO2.inh)/(O2.lim + NO3.lim*O2.inh + NO2.lim*NO3.inh*O2.inh + MnO2.lim*O2.inh*NO3.inh*NO2.inh + FeOOH.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh + SO4.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh + O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh*SO4.inh)
    f.f   <- (FeOOH.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh)/(O2.lim + NO3.lim*O2.inh + NO2.lim*NO3.inh*O2.inh + MnO2.lim*O2.inh*NO3.inh*NO2.inh + FeOOH.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh + SO4.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh + O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh*SO4.inh)
    s.f   <- (SO4.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh)/(O2.lim + NO3.lim*O2.inh + NO2.lim*NO3.inh*O2.inh + MnO2.lim*O2.inh*NO3.inh*NO2.inh + FeOOH.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh + SO4.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh + O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh*SO4.inh)
    m.f   <- (O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh*SO4.inh)/(O2.lim + NO3.lim*O2.inh + NO2.lim*NO3.inh*O2.inh + MnO2.lim*O2.inh*NO3.inh*NO2.inh + FeOOH.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh + SO4.lim*O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh + O2.inh*NO3.inh*NO2.inh*MnO2.inh*FeOOH.inh*SO4.inh)
    
    AR$R   <- a.f*Cmin$R
    DN.1$R <- n1.f*Cmin$R 
    DN.2$R <- n2.f*Cmin$R 
    MR$R   <- mn.f*Cmin$R 
    DIR$R  <- f.f*Cmin$R 
    DIR_56Fe$R[FeOOH.HR$C>=frac.limit] <- (a56Fe_DIR*R_56FeOOH.HR[FeOOH.HR$C>=frac.limit])/(1.0+a56Fe_DIR*R_56FeOOH.HR[FeOOH.HR$C>=frac.limit])*DIR$R[FeOOH.HR$C>=frac.limit]# FeIII reduction
    DIR_56Fe$R[FeOOH.HR$C<frac.limit]  <- r_56FeOOH.HR[FeOOH.HR$C<frac.limit]*DIR$R[FeOOH.HR$C<frac.limit]# FeIII reduction
    SR$R   <- s.f*Cmin$R
    MG$R   <- m.f*Cmin$R 
    
    AR.H2$R   <- a.f*H2min$R
    DN.1.H2$R <- n1.f*H2min$R 
    DN.2.H2$R <- n2.f*H2min$R 
    MR.H2$R   <- mn.f*H2min$R 
    DIR.H2$R  <- f.f*H2min$R 
    DIR_56Fe.H2$R[FeOOH.HR$C>=frac.limit] <- (a56Fe_DIR*R_56FeOOH.HR[FeOOH.HR$C>=frac.limit])/(1.0+a56Fe_DIR*R_56FeOOH.HR[FeOOH.HR$C>=frac.limit])*DIR.H2$R[FeOOH.HR$C>=frac.limit]# FeIII reduction
    DIR_56Fe.H2$R[FeOOH.HR$C<frac.limit]  <- r_56FeOOH.HR[FeOOH.HR$C<frac.limit]*DIR.H2$R[FeOOH.HR$C<frac.limit]# FeIII reduction
    SR.H2$R   <- s.f*H2min$R
    
    #---------------- Aerobic oxidation reactions ----------------------------
    
    NIT.1$R  <- por.grid$mid*k.NIT.1*O2$C*FSAT(O2$C,C.lim/100,5)*NH4$C*FSAT(NH4$C,C.lim/100,5) # aerobic ammonium oxidation
    NIT.2$R  <- por.grid$mid*k.NIT.2*O2$C*FSAT(O2$C,C.lim/100,5)*NO2$C*FSAT(NO2$C,C.lim/100,5) # aerobic nitrite oxidation
    
    MnO$R  <- por.grid$mid*k.MnO*O2$C*FSAT(O2$C,C.lim/100,5)*Mn$C*FSAT(Mn$C,C.lim/100,5)   # aerobic manganese re-oxidation
    
    FIO$R                         <- por.grid$mid*k.FIO*O2$C*FSAT(O2$C,C.lim/100,5)*Fe$C*FSAT(Fe$C,C.lim/100,5) # iron re-oxidation
    FIO_56Fe$R[Fe$C>=frac.limit]  <- (a56Fe_FIO*R_56Fe[Fe$C>=frac.limit])/(1+a56Fe_FIO*R_56Fe[Fe$C>=frac.limit])*FIO$R[Fe$C>=frac.limit]# FeIII reduction
    FIO_56Fe$R[Fe$C<frac.limit]   <- r_56Fe[Fe$C<frac.limit]*FIO$R[Fe$C<frac.limit]
    SIO$R                         <- svf.grid$mid*k.SIO*O2$C*FSAT(O2$C,C.lim,5)*X.Fe$C*FSAT(X.Fe$C,C.lim,5) # sorbed iron re-oxidation
    SIO_56Fe$R[X.Fe$C>=frac.limit] <- (a56Fe_FIO*R_56XFe[X.Fe$C>=frac.limit])/(1+a56Fe_FIO*R_56XFe[X.Fe$C>=frac.limit])*SIO$R[X.Fe$C>=frac.limit]# FeIII reduction
    SIO_56Fe$R[X.Fe$C<frac.limit]  <- r_56XFe[X.Fe$C<frac.limit]*SIO$R[X.Fe$C<frac.limit]
    
    CSO$R  <- por.grid$mid*k.CSO*O2$C*FSAT(O2$C,C.lim/100,5)*HS$C*FSAT(HS$C,C.lim/100,5)   # aerobic Sulfide re-oxidation
    
    AMO$R  <- por.grid$mid*k.AMO*O2$C*FSAT(O2$C,C.lim/100,5)*CH4$C*FSAT(CH4$C,C.lim/100,5) # aerobic methane oxidation
 
    AIO$R                         <- svf.grid$mid*k.AIO*O2$C*FSAT(O2$C,C.lim,5)*FeS$C*FSAT(FeS$C,C.lim,5)*(FeS$C>0)*(O2$C>0) # iron sulfide re-oxidation
    AIO_56Fe$R[FeS$C>=frac.limit] <- (a56Fe_AIO*R_56FeS[FeS$C>=frac.limit])/(1+a56Fe_AIO*R_56FeS[FeS$C>=frac.limit])*AIO$R[FeS$C>=frac.limit]# FeIII reduction
    AIO_56Fe$R[FeS$C<frac.limit]  <- r_56FeS[FeS$C<frac.limit]*AIO$R[FeS$C<frac.limit]
    
    PyO$R                          <- svf.grid$mid*k.PyO*O2$C*FSAT(O2$C,C.lim,5)*FeS2$C*FSAT(FeS2$C,C.lim,5)*(FeS2$C>0)*(O2$C>0) # Pyrite re-oxidation
    PyO_56Fe$R[FeS2$C>=frac.limit] <- (a56Fe_PyO*R_56FeS2[FeS2$C>=frac.limit])/(1+a56Fe_PyO*R_56FeS2[FeS2$C>=frac.limit])*PyO$R[FeS2$C>=frac.limit]# FeIII reduction
    PyO_56Fe$R[FeS2$C<frac.limit]  <- r_56FeS2[FeS2$C<frac.limit]*PyO$R[FeS2$C<frac.limit]
    
    #---------------- Anaerobic oxidation reactions --------------------------
    
    ANA$R <- por.grid$mid*k.ANA*NH4$C*FSAT(NH4$C,C.lim/100,5)*NO2$C*FSAT(NO2$C,C.lim/100,5) # anamox                            # ferrous iron oxidation with nitrate
    
    NFO$R                        <- por.grid$mid*k.NFO*NO3$C*FSAT(NO3$C,C.lim/100,5)*Fe$C*FSAT(Fe$C,C.lim/100,5)                            # ferrous iron oxidation with nitrate
    NFO_56Fe$R[Fe$C>=frac.limit] <- (a56Fe_FIO*R_56Fe[Fe$C>=frac.limit])/(1+a56Fe_FIO*R_56Fe[Fe$C>=frac.limit])*NFO$R[Fe$C>=frac.limit] # ferrous iron oxidation with nitrate
    NFO_56Fe$R[Fe$C<frac.limit]  <- r_56Fe[Fe$C<frac.limit]*NFO$R[Fe$C<frac.limit]                                                      # ferrous iron oxidation with nitrate
    NSO$R                        <- por.grid$mid*k.NSO*NO3$C*FSAT(NO3$C,C.lim/100,5)*HS$C*FSAT(HS$C,C.lim/100,5)    # sulphide oxidation with nitrate
   
    MFO.HR$R                              <- svf.grid$mid*k.MFOHR*MnO2.HR$C*FSAT(MnO2.HR$C,C.lim,5)*Fe$C*FSAT(Fe$C,C.lim/100,5)                       # ferrous iron oxidation with MnIV
    MFO.HR_56Fe$R[Fe$C>=frac.limit]       <- (a56Fe_FIO*R_56Fe[Fe$C>=frac.limit])/(1+a56Fe_FIO*R_56Fe[Fe$C>=frac.limit])*MFO.HR$R[Fe$C>=frac.limit] # ferrous iron oxidation with MnIV
    MFO.HR_56Fe$R[Fe$C<frac.limit]        <- r_56Fe[Fe$C<frac.limit]*MFO.HR$R[Fe$C<frac.limit]                                                        # ferrous iron oxidation with MnIV
    MFO.MR$R                              <- svf.grid$mid*k.MFOMR*MnO2.MR$C*FSAT(MnO2.MR$C,C.lim,5)*Fe$C*FSAT(Fe$C,C.lim/100,5)                       # ferrous iron oxidation with MnIV
    MFO.MR_56Fe$R[Fe$C>=frac.limit]       <- (a56Fe_FIO*R_56Fe[Fe$C>=frac.limit])/(1+a56Fe_FIO*R_56Fe[Fe$C>=frac.limit])*MFO.MR$R[Fe$C>=frac.limit] # ferrous iron oxidation with MnIV
    MFO.MR_56Fe$R[Fe$C<frac.limit]        <- r_56Fe[Fe$C<frac.limit]*MFO.MR$R[Fe$C<frac.limit]                                                        # ferrous iron oxidation with MnIV
    MSO.HR$R                              <- svf.grid$mid*k.MSOHR*MnO2.HR$C*FSAT(MnO2.HR$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim/100,5) # sulphide oxidation with MnIV
    MSO.MR$R                              <- svf.grid$mid*k.MSOMR*MnO2.MR$C*FSAT(MnO2.MR$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim/100,5) # sulphide oxidation with MnIV
    
    ISO.HR$R                              <- svf.grid$mid*k.ISOHR*FeOOH.HR$C*FSAT(FeOOH.HR$C,C.lim,5)*(HS$C^(1./2.))*FSAT(HS$C,C.lim/100,5)# Sulfide re-oxidation by FeIII
    ISO.HR_56Fe$R[FeOOH.HR$C>=frac.limit] <- (a56Fe_ISO*R_56FeOOH.HR[FeOOH.HR$C>=frac.limit])/(1+a56Fe_ISO*R_56FeOOH.HR[FeOOH.HR$C>=frac.limit])*ISO.HR$R[FeOOH.HR$C>=frac.limit]# FeIII reduction
    ISO.HR_56Fe$R[FeOOH.HR$C<frac.limit]  <- r_56FeOOH.HR[FeOOH.HR$C<frac.limit]*ISO.HR$R[FeOOH.HR$C<frac.limit]
    ISO.MR$R                              <- svf.grid$mid*k.ISOMR*FeOOH.MR$C*FSAT(FeOOH.MR$C,C.lim,5)*(HS$C^(1./2.))*FSAT(HS$C,C.lim/100,5)# Sulfide re-oxidation by FeIII
    ISO.MR_56Fe$R[FeOOH.MR$C>=frac.limit] <- (a56Fe_ISO*R_56FeOOH.MR[FeOOH.MR$C>=frac.limit])/(1+a56Fe_ISO*R_56FeOOH.MR[FeOOH.MR$C>=frac.limit])*ISO.MR$R[FeOOH.MR$C>=frac.limit]# FeIII reduction
    ISO.MR_56Fe$R[FeOOH.MR$C<frac.limit]  <- r_56FeOOH.MR[FeOOH.MR$C<frac.limit]*ISO.MR$R[FeOOH.MR$C<frac.limit]
    ISO.PR$R                              <- svf.grid$mid*k.ISOPR*FeOOH.PR$C*FSAT(FeOOH.PR$C,C.lim,5)*(HS$C^(1./2.))*FSAT(HS$C,C.lim/100,5)# Sulfide re-oxidation by FeIII
    ISO.PR_56Fe$R[FeOOH.PR$C>=frac.limit] <- (a56Fe_ISO*R_56FeOOH.PR[FeOOH.PR$C>=frac.limit])/(1+a56Fe_ISO*R_56FeOOH.PR[FeOOH.PR$C>=frac.limit])*ISO.PR$R[FeOOH.PR$C>=frac.limit]# FeIII reduction
    ISO.PR_56Fe$R[FeOOH.PR$C<frac.limit]  <- r_56FeOOH.PR[FeOOH.PR$C<frac.limit]*ISO.PR$R[FeOOH.PR$C<frac.limit]
    ISO.U$R                               <- svf.grid$mid*k.ISOU*FeOOH.U$C*FSAT(FeOOH.U$C,C.lim,5)*(HS$C^(1./2.))*FSAT(HS$C,C.lim/100,5)# Sulfide re-oxidation by FeIII
    ISO.U_56Fe$R[FeOOH.U$C>=frac.limit]   <- (a56Fe_ISO*R_56FeOOH.U[FeOOH.U$C>=frac.limit])/(1+a56Fe_ISO*R_56FeOOH.U[FeOOH.U$C>=frac.limit])*ISO.U$R[FeOOH.U$C>=frac.limit]# FeIII reduction
    ISO.U_56Fe$R[FeOOH.U$C<frac.limit]    <- r_56FeOOH.U[FeOOH.U$C<frac.limit]*ISO.U$R[FeOOH.U$C<frac.limit]
    
    SMO$R <- por.grid$mid*k.SMO*SO4$C*FSAT(SO4$C,C.lim/100,5)*CH4$C*FSAT(CH4$C,C.lim/100,5)  # Sulfate methane oxidation
    
    #---------------- Minerals -----------------------------------------------
   
    AmS$R <- k.AmS*(1.-Omega.NH4$C)
    
    SDP$R <- k.SDP*S0$C*FSAT(S0$C,C.lim,5)
    
    #ISP$R                         <- svf.grid$mid*k.ISP*((Omega.FeS$C-1)^n.ISP)*(Omega.FeS$C>1)*FSAT(HS$C,C.lim/100,5)*FSAT(Fe$C,C.lim/100,5)*(Fe$C>0)*(HS$C>0)              #iron sulfide precipitation
    ISP$R                         <- por.grid$mid*k.ISP*Fe$C*HS$C*FSAT(HS$C,C.lim/100,5)*FSAT(Fe$C,C.lim/100,5)*(Fe$C>0)*(HS$C>0)              #iron sulfide precipitation
    ISP_56Fe$R[Fe$C>=frac.limit]  <- (a56Fe_ISP*R_56Fe[Fe$C>=frac.limit])/(1+a56Fe_ISP*R_56Fe[Fe$C>=frac.limit])*ISP$R[Fe$C>=frac.limit]# FeIII reduction
    ISP_56Fe$R[Fe$C<frac.limit]   <- r_56Fe[Fe$C<frac.limit]*ISP$R[Fe$C<frac.limit]
    
    ISD$R                         <- svf.grid$mid*k.ISD*FeS$C*FSAT(FeS$C,C.lim,5)*(1-Omega.FeS$C)^n.ISD*(Omega.FeS$C<1)*(FeS$C>0) #iron sulfide dissolution
    ISD_56Fe$R[FeS$C>=frac.limit] <- (a56Fe_ISD*R_56FeS[FeS$C>=frac.limit])/(1+a56Fe_ISD*R_56FeS[FeS$C>=frac.limit])*ISD$R[FeS$C>=frac.limit]# FeIII reduction
    ISD_56Fe$R[FeS$C<frac.limit]  <- r_56FeS[FeS$C<frac.limit]*ISD$R[FeS$C<frac.limit]
    
    PyP.1$R                         <- svf.grid$mid*k.PyP.1*FeS$C*FSAT(FeS$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim/100,5)*(FeS$C>0)*(HS$C>0) # pyrite precipitation via H2S
    PyP.1_56Fe$R[FeS$C>=frac.limit] <- (a56Fe_PyP*R_56FeS[FeS$C>=frac.limit])/(1+a56Fe_PyP*R_56FeS[FeS$C>=frac.limit])*PyP.1$R[FeS$C>=frac.limit]
    PyP.1_56Fe$R[FeS$C<frac.limit]  <- r_56FeS[FeS$C<frac.limit]*PyP.1$R[FeS$C<frac.limit]
    
    PyP.2$R                         <- svf.grid$mid*k.PyP.2*FeS$C*FSAT(FeS$C,C.lim,5)*S0$C*FSAT(S0$C,C.lim,5)*(FeS$C>0)*(S0$C>0) # pyrite precipitation via S0
    PyP.2_56Fe$R[FeS$C>=frac.limit] <- (a56Fe_PyP*R_56FeS[FeS$C>=frac.limit])/(1+a56Fe_PyP*R_56FeS[FeS$C>=frac.limit])*PyP.2$R[FeS$C>=frac.limit]
    PyP.2_56Fe$R[FeS$C<frac.limit]  <- r_56FeS[FeS$C<frac.limit]*PyP.2$R[FeS$C<frac.limit]
    
    #---------------- Metal oxide ageing ------------------------------------

    MOA$R <- svf.grid$mid*k.MOA*MnO2.HR$C*FSAT(MnO2.HR$C,C.lim,5) # manganese oxide ageing
    
    IOA$R                              <- svf.grid$mid*k.IOA*FeOOH.HR$C*FSAT(FeOOH.HR$C,C.lim,5) # Iron oxide ageing
    IOA_56Fe$R[FeOOH.HR$C>=frac.limit] <- (a56Fe_IOA*R_56FeOOH.HR[FeOOH.HR$C>=frac.limit])/(1+a56Fe_IOA*R_56FeOOH.HR[FeOOH.HR$C>=frac.limit])*IOA$R[FeOOH.HR$C>=frac.limit]# FeIII reduction
    IOA_56Fe$R[FeOOH.HR$C<frac.limit]  <- r_56FeOOH.HR[FeOOH.HR$C<frac.limit]*IOA$R[FeOOH.HR$C<frac.limit]
    #FIS$R  <- svf.grid$mid*k.FIS*(K.FIS*Fe$C - X.Fe$C)*(Fe$C>0)*(X.Fe$C>0) # ferric iron sorption

    #=========================================================================
    # Consumption terms: kinetic reactions
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================
    
    POC.1$kin.reac  <- - (Cmin.1$R)/svf.grid$mid
    POC.2$kin.reac  <- - (Cmin.2$R)/svf.grid$mid
    POC.3$kin.reac  <- - (Cmin.3$R)/svf.grid$mid
    POC.4$kin.reac  <- - (Cmin.4$R)/svf.grid$mid
    POC.5$kin.reac  <- - (Cmin.5$R)/svf.grid$mid
    POC.6$kin.reac  <- - (Cmin.6$R)/svf.grid$mid
    POC.7$kin.reac  <- - (Cmin.7$R)/svf.grid$mid
    POC.8$kin.reac  <- - (Cmin.8$R)/svf.grid$mid
    POC.9$kin.reac  <- - (Cmin.9$R)/svf.grid$mid
    POC.10$kin.reac <- - (Cmin.10$R)/svf.grid$mid
    POC.11$kin.reac <- - (Cmin.11$R)/svf.grid$mid
    POC.12$kin.reac <- - (Cmin.12$R)/svf.grid$mid
    POC.13$kin.reac <- - (Cmin.13$R)/svf.grid$mid
    POC.14$kin.reac <- - (Cmin.14$R)/svf.grid$mid
    
    HCO3$kin.reac <- (AR$R + DN.1$R + DN.2$R + MR$R + DIR$R + SR$R + 1/2*MG$R + AMO$R + SMO$R)/por.grid$mid
    CH4$kin.reac  <- (                                             + 1/2*MG$R - AMO$R - SMO$R)/por.grid$mid
    
    O2$kin.reac  <- (- AR$R - 1/2*AR.H2$R - 3/2*NIT.1$R - 1/2*NIT.2$R - 1/2*MnO$R - 1/4*FIO$R - 1/4*SIO$R - 2*CSO$R - 2*AMO$R - 9/4*AIO$R - 15/4*PyO$R)/por.grid$mid
    
    NH4$kin.reac   <- (Cmin$R/CN.ratio                                                     - NIT.1$R           - ANA$R             + NSO$R - AmS$R)/por.grid$mid
    X.NH4$kin.reac <- (                                                                                                                    + AmS$R)/svf.grid$mid
    NO3$kin.reac   <- (                - 2*DN.1$R              - DN.1.H2$R                           + NIT.2$R         - 1/5*NFO$R - NSO$R        )/por.grid$mid
    NO2$kin.reac   <- (                + 2*DN.1$R - 4/3*DN.2$R + DN.1.H2$R - 2/3*DN.2.H2$R + NIT.1$R - NIT.2$R - ANA$R                            )/por.grid$mid
    N2$kin.reac    <- (                           + 2/3*DN.2$R             + 1/3*DN.2.H2$R                     + ANA$R + 1/10*NFO$R               )/por.grid$mid
    
    MnO2.HR$kin.reac  <- (- 2*MR$R - MR.H2$R + MnO$R - 1/2*MFO.HR$R                - MSO.HR$R            - MOA$R)/svf.grid$mid  
    MnO2.MR$kin.reac  <- (                                          - 1/2*MFO.MR$R            - MSO.MR$R + MOA$R)/svf.grid$mid  
    Mn$kin.reac       <- (+ 2*MR$R + MR.H2$R - MnO$R + 1/2*MFO.HR$R + 1/2*MFO.MR$R + MSO.HR$R + MSO.MR$R        )/por.grid$mid
   
    FeOOH.HR$kin.reac      <- (- 4*DIR$R      - 2*DIR.H2$R      + FIO$R      + SIO$R      + AIO$R      + PyO$R      + NFO$R      + MFO.HR$R      + MFO.MR$R       - 2*ISO.HR$R                                                                                         - IOA$R     )/svf.grid$mid  
    Fe56_FeOOH.HR$kin.reac <- (- 4*DIR_56Fe$R - 2*DIR_56Fe.H2$R + FIO_56Fe$R + SIO_56Fe$R + AIO_56Fe$R + PyO_56Fe$R + NFO_56Fe$R + MFO.HR_56Fe$R + MFO.MR_56Fe$R  - 2*ISO.HR_56Fe$R                                                                                    - IOA_56Fe$R)/svf.grid$mid  
    FeOOH.MR$kin.reac      <- (                                                                                                                                                     - 2*ISO.MR$R                                                                       + IOA$R     )/svf.grid$mid    
    Fe56_FeOOH.MR$kin.reac <- (                                                                                                                                                     - 2*ISO.MR_56Fe$R                                                                  + IOA_56Fe$R)/svf.grid$mid    
    FeOOH.PR$kin.reac      <- (                                                                                                                                                                       - 2*ISO.PR$R                                                                 )/svf.grid$mid    
    Fe56_FeOOH.PR$kin.reac <- (                                                                                                                                                                       - 2*ISO.PR_56Fe$R                                                            )/svf.grid$mid    
    FeOOH.U$kin.reac       <- (                                                                                                                                                                                         - 2*ISO.U$R                                                )/svf.grid$mid    
    Fe56_FeOOH.U$kin.reac  <- (                                                                                                                                                                                         - 2*ISO.U_56Fe$R                                           )/svf.grid$mid    
    Fe$kin.reac            <- (+ 4*DIR$R      + 2*DIR.H2$R      - FIO$R                                             - NFO$R      - MFO.HR$R      - MFO.MR$R       + 2*ISO.HR$R      + 2*ISO.MR$R      + 2*ISO.PR$R      + 2*ISO.U$R                                                )/por.grid$mid
    Fe56_Fe$kin.reac       <- (+ 4*DIR_56Fe$R + 2*DIR_56Fe.H2$R - FIO_56Fe$R                                        - NFO_56Fe$R - MFO.HR_56Fe$R - MFO.MR_56Fe$R  + 2*ISO.HR_56Fe$R + 2*ISO.MR_56Fe$R + 2*ISO.PR_56Fe$R + 2*ISO.U_56Fe$R                                           )/por.grid$mid
    X.Fe$kin.reac          <- (                                              - SIO$R                                                                                                                                                                                               )/svf.grid$mid
    X.Fe56_Fe$kin.reac     <- (                                              - SIO_56Fe$R                                                                                                                                                                                          )/svf.grid$mid
    FeS$kin.reac           <- (                                                           - AIO$R                                                                                                                                        - PyP.1$R      - PyP.2$R                  )/svf.grid$mid
    Fe56_FeS$kin.reac      <- (                                                           - AIO_56Fe$R                                                                                                                                   - PyP.1_56Fe$R - PyP.2_56Fe$R             )/svf.grid$mid
    FeS2$kin.reac          <- (                                                                        - PyO$R                                                                                                                           + PyP.1$R      + PyP.2$R                  )/svf.grid$mid
    Fe56_FeS2$kin.reac     <- (                                                                        - PyO_56Fe$R                                                                                                                      + PyP.1_56Fe$R + PyP.2_56Fe$R             )/svf.grid$mid
    
    SO4$kin.reac <- (- 1/2*SR$R - 1/4*SR.H2$R + CSO$R + AIO$R + 2*PyO$R + NSO$R                                                                  - SMO$R + 1/4*SDP$R          )/por.grid$mid
    S0$kin.reac  <- (                                                           + MSO.HR$R + MSO.MR$R + ISO.HR$R + ISO.MR$R + ISO.PR$R + ISO.U$R         -     SDP$R - PyP.2$R)/svf.grid$mid
    HS$kin.reac  <- (+ 1/2*SR$R + 1/4*SR.H2$R - CSO$R                   - NSO$R - MSO.HR$R - MSO.MR$R - ISO.HR$R - ISO.MR$R - ISO.PR$R - ISO.U$R + SMO$R + 3/4*SDP$R - PyP.1$R)/por.grid$mid  
     
    H2$kin.reac  <- (- H2min$R + PyP.1$R)/por.grid$mid
    
    H2O$kin.reac <- (                 + 2/3*DN.2$R + 2*MR$R + 6*DIR$R            - 1/2*MG$R + AR.H2$R + DN.1.H2$R + 4/3*DN.2.H2$R + 2*MR.H2$R + 4*DIR.H2$R +     SR.H2$R +   NIT.1$R -   MnO$R - 3/2*FIO$R - 3/2*SIO$R         + AMO$R - 3/2*AIO$R - 5/2*PyO$R + 2*ANA$R - 7/5*NFO$R - NSO$R - MFO.HR$R - MFO.MR$R + 2*MSO.HR$R + 2*MSO.MR$R + 4*ISO.HR$R + 4*ISO.MR$R + 4*ISO.PR$R + 4*ISO.U$R + SMO$R -     SDP$R          )/por.grid$mid
    H$kin.reac   <- (+ AR$R  + DN.1$R - 1/3*DN.2$R - 3*MR$R - 7*DIR$R + 1/2*SR$R + 1/2*MG$R                       - 2/3*DN.2.H2$R - 2*MR.H2$R - 4*DIR.H2$R - 1/4*SR.H2$R + 2*NIT.1$R + 2*MnO$R +   2*FIO$R +   2*SIO$R + CSO$R + AMO$R +   2*AIO$R +   4*PyO$R           + 9/5*NFO$R - NSO$R + MFO.HR$R + MFO.MR$R - 3*MSO.HR$R - 3*MSO.MR$R - 5*ISO.HR$R - 5*ISO.MR$R - 5*ISO.PR$R - 5*ISO.U$R         + 5/4*SDP$R - PyP.1$R)/por.grid$mid

    #=========================================================================
    # Equilibrium reactions
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================
    
    Fe$eq.reac       <- (- ISP$R      + ISD$R     )/por.grid$mid
    Fe56_Fe$eq.reac  <- (- ISP_56Fe$R + ISD_56Fe$R)/por.grid$mid
    FeS$eq.reac      <- (  ISP$R      - ISD$R     )/svf.grid$mid
    Fe56_FeS$eq.reac <- (  ISP_56Fe$R - ISD_56Fe$R)/svf.grid$mid
    
    HS$eq.reac       <- (- ISP$R      + ISD$R)/por.grid$mid
    
    H$eq.reac        <- (  ISP$R      - ISD$R)/por.grid$mid

    #=========================================================================
    # Irrigation 
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================

    HCO3$irr    <- irr.grid$mid*irr.HCO3*(HCO3.ow - HCO3$C)
    CH4$irr     <- irr.grid$mid*irr.CH4*(CH4.ow - CH4$C)
   
    O2$irr      <- irr.grid$mid*irr.O2*(O2.ow - O2$C)
    
    NH4$irr     <- irr.grid$mid*irr.NH4*(NH4.ow - NH4$C)
    NO3$irr     <- irr.grid$mid*irr.NO3*(NO3.ow - NO3$C)
    NO2$irr     <- irr.grid$mid*irr.NO2*(NO2.ow - NO2$C)
    N2$irr      <- irr.grid$mid*irr.N2*(N2.ow - N2$C)
    
    Mn$irr      <- irr.grid$mid*irr.Mn*(Mn.ow - Mn$C)
    
    Fe$irr      <- irr.grid$mid*irr.Fe*(Fe.ow - Fe$C)
    Fe56_Fe$irr <- irr.grid$mid*irr.Fe*(Fe56_Fe.ow - Fe56_Fe$C)
    
    SO4$irr     <- irr.grid$mid*irr.SO4*(SO4.ow - SO4$C)
    HS$irr      <- irr.grid$mid*irr.HS*(HS.ow - HS$C)
    
    H2O$irr     <- irr.grid$mid*(H2O.ow - H2O$C)
    H$irr       <- irr.grid$mid*(H.ow - H$C)
    H2$irr      <- irr.grid$mid*(H2.ow - H2$C)
    
    #=========================================================================
    # Rate of changes composite species 
    #=========================================================================
    
    POC$tran     <- POC.1$tran     + POC.2$tran     + POC.3$tran     + POC.4$tran     + POC.5$tran     + POC.6$tran     + POC.7$tran     + POC.8$tran     + POC.9$tran     + POC.10$tran     + POC.11$tran     + POC.12$tran     + POC.13$tran     + POC.14$tran
    POC$kin.reac <- POC.1$kin.reac + POC.2$kin.reac + POC.3$kin.reac + POC.4$kin.reac + POC.5$kin.reac + POC.6$kin.reac + POC.7$kin.reac + POC.8$kin.reac + POC.9$kin.reac + POC.10$kin.reac + POC.11$kin.reac + POC.12$kin.reac + POC.13$kin.reac + POC.14$kin.reac
    POC$irr      <- POC.1$irr      + POC.2$irr      + POC.3$irr      + POC.4$irr      + POC.5$irr      + POC.6$irr      + POC.7$irr      + POC.8$irr      + POC.9$irr      + POC.10$irr      + POC.11$irr      + POC.12$irr      + POC.13$irr      + POC.14$irr
    POC$eq.reac  <- POC.1$eq.reac  + POC.2$eq.reac  + POC.3$eq.reac  + POC.4$eq.reac  + POC.5$eq.reac  + POC.6$eq.reac  + POC.7$eq.reac  + POC.8$eq.reac  + POC.9$eq.reac  + POC.10$eq.reac  + POC.11$eq.reac  + POC.12$eq.reac  + POC.13$eq.reac  + POC.14$eq.reac

    MnO2$tran       <- MnO2.HR$tran     +  MnO2.MR$tran
    MnO2$kin.reac   <- MnO2.HR$kin.reac +  MnO2.MR$kin.reac
    MnO2$irr        <- MnO2.HR$irr      +  MnO2.MR$irr
    MnO2$eq.reac    <- MnO2.HR$eq.reac  +  MnO2.MR$eq.reac
    
    FeOOH$tran          <- FeOOH.HR$tran          + FeOOH.MR$tran          + FeOOH.PR$tran          + FeOOH.U$tran
    Fe56_FeOOH$tran     <- Fe56_FeOOH.HR$tran     + Fe56_FeOOH.MR$tran     + Fe56_FeOOH.PR$tran     + Fe56_FeOOH.U$tran
    FeOOH$kin.reac      <- FeOOH.HR$kin.reac      + FeOOH.MR$kin.reac      + FeOOH.PR$kin.reac      + FeOOH.U$kin.reac
    Fe56_FeOOH$kin.reac <- Fe56_FeOOH.HR$kin.reac + Fe56_FeOOH.MR$kin.reac + Fe56_FeOOH.PR$kin.reac + Fe56_FeOOH.U$kin.reac
    FeOOH$irr           <- FeOOH.HR$irr           + FeOOH.MR$irr           + FeOOH.PR$irr           + FeOOH.U$irr
    Fe56_FeOOH$irr      <- Fe56_FeOOH.HR$irr      + Fe56_FeOOH.MR$irr      + Fe56_FeOOH.PR$irr      + Fe56_FeOOH.U$irr
    FeOOH$eq.reac       <- FeOOH.HR$eq.reac       + FeOOH.MR$eq.reac       + FeOOH.PR$eq.reac       + FeOOH.U$eq.reac
    Fe56_FeOOH$eq.reac  <- Fe56_FeOOH.HR$eq.reac  + Fe56_FeOOH.MR$eq.reac  + Fe56_FeOOH.PR$eq.reac  + Fe56_FeOOH.U$eq.reac
    
    #=========================================================================
    # Rate of changes of all species
    #=========================================================================
    
    POC.1$ddt  <- POC.1$tran  + POC.1$kin.reac  + POC.1$eq.reac  + POC.1$irr
    POC.2$ddt  <- POC.2$tran  + POC.2$kin.reac  + POC.2$eq.reac  + POC.2$irr
    POC.3$ddt  <- POC.3$tran  + POC.3$kin.reac  + POC.3$eq.reac  + POC.3$irr
    POC.4$ddt  <- POC.4$tran  + POC.4$kin.reac  + POC.4$eq.reac  + POC.4$irr
    POC.5$ddt  <- POC.5$tran  + POC.5$kin.reac  + POC.5$eq.reac  + POC.5$irr
    POC.6$ddt  <- POC.6$tran  + POC.6$kin.reac  + POC.6$eq.reac  + POC.6$irr
    POC.7$ddt  <- POC.7$tran  + POC.7$kin.reac  + POC.7$eq.reac  + POC.7$irr
    POC.8$ddt  <- POC.8$tran  + POC.8$kin.reac  + POC.8$eq.reac  + POC.8$irr
    POC.9$ddt  <- POC.9$tran  + POC.9$kin.reac  + POC.9$eq.reac  + POC.9$irr
    POC.10$ddt <- POC.10$tran + POC.10$kin.reac + POC.10$eq.reac + POC.10$irr
    POC.11$ddt <- POC.11$tran + POC.11$kin.reac + POC.11$eq.reac + POC.11$irr
    POC.12$ddt <- POC.12$tran + POC.12$kin.reac + POC.12$eq.reac + POC.12$irr
    POC.13$ddt <- POC.13$tran + POC.13$kin.reac + POC.13$eq.reac + POC.13$irr
    POC.14$ddt <- POC.14$tran + POC.14$kin.reac + POC.14$eq.reac + POC.14$irr
    
    HCO3$ddt  <- HCO3$tran  + HCO3$kin.reac  + HCO3$eq.reac  + HCO3$irr
    CH4$ddt   <- CH4$tran   + CH4$kin.reac   + CH4$eq.reac   + CH4$irr
    
    O2$ddt    <- O2$tran    + O2$kin.reac    + O2$eq.reac    + O2$irr
    
    NH4$ddt   <- NH4$tran   + NH4$kin.reac   + NH4$eq.reac   + NH4$irr
    X.NH4$ddt <- X.NH4$tran + X.NH4$kin.reac + X.NH4$eq.reac + X.NH4$irr
    NO3$ddt   <- NO3$tran   + NO3$kin.reac   + NO3$eq.reac   + NO3$irr
    NO2$ddt   <- NO2$tran   + NO2$kin.reac   + NO2$eq.reac   + NO2$irr
    N2$ddt    <- N2$tran    + N2$kin.reac    + N2$eq.reac    + N2$irr
    
    MnO2.HR$ddt <- MnO2.HR$tran + MnO2.HR$kin.reac + MnO2.HR$eq.reac + MnO2.HR$irr
    MnO2.MR$ddt <- MnO2.MR$tran + MnO2.MR$kin.reac + MnO2.MR$eq.reac + MnO2.MR$irr
    Mn$ddt      <- Mn$tran      + Mn$kin.reac      + Mn$eq.reac      + Mn$irr 
    
    FeOOH.HR$ddt      <- FeOOH.HR$tran      + FeOOH.HR$kin.reac      + FeOOH.HR$eq.reac      + FeOOH.HR$irr
    Fe56_FeOOH.HR$ddt <- Fe56_FeOOH.HR$tran + Fe56_FeOOH.HR$kin.reac + Fe56_FeOOH.HR$eq.reac + Fe56_FeOOH.HR$irr
    FeOOH.MR$ddt      <- FeOOH.MR$tran      + FeOOH.MR$kin.reac      + FeOOH.MR$eq.reac      + FeOOH.MR$irr
    Fe56_FeOOH.MR$ddt <- Fe56_FeOOH.MR$tran + Fe56_FeOOH.MR$kin.reac + Fe56_FeOOH.MR$eq.reac + Fe56_FeOOH.MR$irr
    FeOOH.PR$ddt      <- FeOOH.PR$tran      + FeOOH.PR$kin.reac      + FeOOH.PR$eq.reac      + FeOOH.PR$irr
    Fe56_FeOOH.PR$ddt <- Fe56_FeOOH.PR$tran + Fe56_FeOOH.PR$kin.reac + Fe56_FeOOH.PR$eq.reac + Fe56_FeOOH.PR$irr
    FeOOH.U$ddt       <- FeOOH.U$tran       + FeOOH.U$kin.reac       + FeOOH.U$eq.reac       + FeOOH.U$irr
    Fe56_FeOOH.U$ddt  <- Fe56_FeOOH.U$tran  + Fe56_FeOOH.U$kin.reac  + Fe56_FeOOH.U$eq.reac  + Fe56_FeOOH.U$irr
    Fe$ddt       <- 1/(por.grid$mid+svf.grid$mid*K.FIS)*(por.grid$mid*(Fe$tran      + Fe$kin.reac      + Fe$eq.reac      + Fe$irr)      + svf.grid$mid*(X.Fe$tran      + X.Fe$kin.reac      + X.Fe$eq.reac      + X.Fe$irr))
    Fe56_Fe$ddt  <- 1/(por.grid$mid+svf.grid$mid*K.FIS)*(por.grid$mid*(Fe56_Fe$tran + Fe56_Fe$kin.reac + Fe56_Fe$eq.reac + Fe56_Fe$irr) + svf.grid$mid*(X.Fe56_Fe$tran + X.Fe56_Fe$kin.reac + X.Fe56_Fe$eq.reac + X.Fe56_Fe$irr))
    FeS$ddt           <- FeS$tran           + FeS$kin.reac           + FeS$eq.reac           + FeS$irr
    Fe56_FeS$ddt      <- Fe56_FeS$tran      + Fe56_FeS$kin.reac      + Fe56_FeS$eq.reac      + Fe56_FeS$irr
    FeS2$ddt          <- FeS2$tran          + FeS2$kin.reac          + FeS2$eq.reac          + FeS2$irr
    Fe56_FeS2$ddt     <- Fe56_FeS2$tran     + Fe56_FeS2$kin.reac     + Fe56_FeS2$eq.reac     + Fe56_FeS2$irr
    
    SO4$ddt   <- SO4$tran + SO4$kin.reac + SO4$eq.reac + SO4$irr
    HS$ddt    <- HS$tran  + HS$kin.reac  + HS$eq.reac  + HS$irr
    S0$ddt    <- S0$tran  + S0$kin.reac  + S0$eq.reac  + S0$irr
    
    H2O$ddt   <- H2O$tran + H2O$kin.reac + H2O$eq.reac  + H2O$irr
    H$ddt     <- H$tran   + H$kin.reac   + H$eq.reac    + H$irr
    H2$ddt    <- H2$tran  + H2$kin.reac  + H2$eq.reac   + H2$irr
    
    POC$ddt        <- POC$tran        + POC$kin.reac        + POC$eq.reac        + POC$irr
    MnO2$ddt       <- MnO2$tran       + MnO2$kin.reac       + MnO2$eq.reac       + MnO2$irr
    FeOOH$ddt      <- FeOOH$tran      + FeOOH$kin.reac      + FeOOH$eq.reac      + FeOOH$irr
    Fe56_FeOOH$ddt <- Fe56_FeOOH$tran + Fe56_FeOOH$kin.reac + Fe56_FeOOH$eq.reac + Fe56_FeOOH$irr
    
    #=========================================================================
    }) # end within call
    #=========================================================================
  
    #-----------------------------------------------------------------------------
    # Assemble matrix with rate of changes (RC) 
    #-----------------------------------------------------------------------------
    
    RC <- 0*SV
    for (i in 1:N.var) RC[,i] <- OL[[var.names[i]]]$ddt

  
  #===========================================================================
  # return statement 
  #===========================================================================
  
  if (!summary.call) {return(list(RC))}
  else {
    
    # Initialisation of lists
    
    SL <- OL[names(SL)]
    RL <- OL[names(RL)]
    EL <- EL
    
    # Updating species information
    
    SL <- updating.chemical.species(SL,PL) 
    SL <- updating.composite.species(spec.table,SL) 
    
    # Updating reaction information
    
    for (i in 1:length(RL)) RL[[i]]$R.int <- sum(RL[[i]]$R*grid$dx) 
    
    # Updating element information
    
    EL <- updating.elements(elt.table,EL,SL) 
    
    return(list(SL=SL,RL=RL,EL=EL,PL=parameters))
    
   }
   
   #-----------------------------------------------------------------------------
 }) # end with call
}
#-----------------------------------------------------------------------------

################################################################################
# Initialisation of state variable matrix
################################################################################

initialise.state.variables <- function(PL)
{
  
  PL$var.names <- c(paste("POC.",1:14,sep=""),"HCO3","CH4","O2","NH4","X.NH4","NO3","NO2","N2","MnO2.HR","MnO2.MR","Mn",
                    "FeOOH.HR","Fe56_FeOOH.HR","FeOOH.MR","Fe56_FeOOH.MR","FeOOH.PR","Fe56_FeOOH.PR","FeOOH.U","Fe56_FeOOH.U",
                    "Fe","Fe56_Fe","X.Fe","X.Fe56_Fe","FeS","Fe56_FeS","FeS2","Fe56_FeS2",
                    "SO4","HS","S0","H2O","H","H2")
  PL$non.var.names <- NULL
  
  PL$N <- PL$grid$N
  
  PL$N.var <- length(PL$var.names)
  PL$N.non.var <- length(PL$non.var.names)
  
  return(PL)
}

initialise.species.arrays <- function(PL,SL)
{
  
  for (i in 1:length(SL))
  {
    SL[[i]]$C <- rep(0,PL$N)
    SL[[i]]$tran <- rep(0,PL$N)
    SL[[i]]$kin.reac <- rep(0,PL$N)
    SL[[i]]$eq.reac <- rep(0,PL$N)
    SL[[i]]$irr <- rep(0,PL$N)
  }
  return(SL)
}

# Function that creates a rough guess for the initial state

crude.initial.state <- function(PL)
{
  with(PL,{
    
    # Initialisation State Variable matrix (SV) 
    
    SV.guess <- matrix(nrow=N,ncol=N.var,dimnames=list(1:N,var.names))
    
    for (i in 1:14){
      if (paste("POC.",i,sep="") %in% PL$var.names) SV.guess[,paste("POC.",i,sep="")] <- rep(0,length.out=N)
    }
    
    if ("HCO3" %in% PL$var.names)    SV.guess[,"HCO3"]    <- rep(HCO3.ow,length.out=N)
    if ("CH4" %in% PL$var.names)     SV.guess[,"CH4"]     <- rep(CH4.ow,length.out=N)
    
    if ("O2" %in% PL$var.names)      SV.guess[,"O2"]      <- rep(O2.ow,length.out=N)
    
    if ("NH4" %in% PL$var.names)     SV.guess[,"NH4"]     <- rep(NH4.ow,length.out=N)
    if ("X.NH4" %in% PL$var.names)   SV.guess[,"X.NH4"]   <- rep(0,length.out=N)
    if ("NO3" %in% PL$var.names)     SV.guess[,"NO3"]     <- rep(NO3.ow,length.out=N)
    if ("NO2" %in% PL$var.names)     SV.guess[,"NO2"]     <- rep(NO2.ow,length.out=N)
    if ("N2" %in% PL$var.names)      SV.guess[,"N2"]      <- rep(N2.ow,length.out=N)
    
    if ("MnO2.HR" %in% PL$var.names)  SV.guess[,"MnO2.HR"]  <- rep(0,length.out=N)
    if ("MnO2.MR" %in% PL$var.names)  SV.guess[,"MnO2.MR"]  <- rep(0,length.out=N)
    if ("Mn" %in% PL$var.names)       SV.guess[,"Mn"]       <- rep(Mn.ow,length.out=N)
  
    if ("FeOOH.HR" %in% PL$var.names)      SV.guess[,"FeOOH.HR"]      <- rep(0,length.out=N)
    if ("Fe56_FeOOH.HR" %in% PL$var.names) SV.guess[,"Fe56_FeOOH.HR"] <- rep(0,length.out=N)
    if ("FeOOH.MR" %in% PL$var.names)      SV.guess[,"FeOOH.MR"]      <- rep(0,length.out=N)
    if ("Fe56_FeOOH.MR" %in% PL$var.names) SV.guess[,"Fe56_FeOOH.MR"] <- rep(0,length.out=N)
    if ("FeOOH.PR" %in% PL$var.names)      SV.guess[,"FeOOH.PR"]      <- rep(0,length.out=N)
    if ("Fe56_FeOOH.PR" %in% PL$var.names) SV.guess[,"Fe56_FeOOH.PR"] <- rep(0,length.out=N)
    if ("FeOOH.U" %in% PL$var.names)       SV.guess[,"FeOOH.U"]       <- rep(0,length.out=N)
    if ("Fe56_FeOOH.U" %in% PL$var.names)  SV.guess[,"Fe56_FeOOH.U"]  <- rep(0,length.out=N)
    if ("Fe" %in% PL$var.names)            SV.guess[,"Fe"]            <- rep(Fe.ow,length.out=N)
    if ("Fe56_Fe" %in% PL$var.names)       SV.guess[,"Fe56_Fe"]       <- rep(Fe56_Fe.ow,length.out=N)
    if ("X.Fe" %in% PL$var.names)          SV.guess[,"X.Fe"]          <- rep(0,length.out=N)
    if ("X.Fe56_Fe" %in% PL$var.names)     SV.guess[,"X.Fe56_Fe"]     <- rep(0,length.out=N)
    if ("FeS" %in% PL$var.names)           SV.guess[,"FeS"]           <- rep(0,length.out=N)
    if ("Fe56_FeS" %in% PL$var.names)      SV.guess[,"Fe56_FeS"]      <- rep(0,length.out=N)
    if ("FeS2" %in% PL$var.names)          SV.guess[,"FeS2"]          <- rep(0,length.out=N)
    if ("Fe56_FeS2" %in% PL$var.names)     SV.guess[,"Fe56_FeS2"]     <- rep(0,length.out=N)
     
    if ("SO4" %in% PL$var.names) SV.guess[,"SO4"] <- rep(SO4.ow,length.out=N)
    if ("HS" %in% PL$var.names)  SV.guess[,"HS"]  <- rep(HS.ow,length.out=N)
    if ("S0" %in% PL$var.names)  SV.guess[,"S0"]  <- rep(0,length.out=N)
    
    if ("H2O" %in% PL$var.names) SV.guess[,"H2O"] <- rep(H2O.ow,length.out=N)
    if ("H" %in% PL$var.names)   SV.guess[,"H"]   <- rep(H.ow,length.out=N)
    if ("H2" %in% PL$var.names)  SV.guess[,"H2"]  <- rep(H2.ow,length.out=N)
    
    return(SV.guess)
    
  })
}

################################################################################
# Test trial model runs
################################################################################

model <- CSFe.model
sim.info <- list(tc= 0,fc=0)

# Initialisation of object list (OL)  

PL <- initialise.parameters(PL=list())
PL$simulation.type <- "direct.steady.state"
SL <- initialise.species.arrays(PL,SL) # default action
OL <- c(SL,RL) # default action

# Test trial 

PL <- initialise.state.variables(PL)
SV <- crude.initial.state(PL)

system.time({
  result <- model(t=0,state=SV,parameters=PL,summary.call=FALSE)
  result <- model(t=0,state=SV,parameters=PL,summary.call=TRUE)
})

remove(result)
remove(SV)

