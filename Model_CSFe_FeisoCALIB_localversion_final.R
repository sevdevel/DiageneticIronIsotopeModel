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

# Organic matter -> multi-G approx. of reactive continuum (following Dale et al. 2015) 

for (i in 1:14){
  SL <- initialise.list.species(SL,"1D",name=paste("POC.",i,sep=""),phase="solid",type="chemical")
}

# chemical species incorporated in the model 

SL <- initialise.list.species(SL,"1D",name="FeOOH.f",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeOOH.f",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeOOH.a",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeOOH.a",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeS",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeS",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="FeS2",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeS2",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="X.Fe",phase="solid",type="chemical")
SL <- initialise.list.species(SL,"1D",name="X.Fe56_Fe",phase="solid",type="chemical")

SL <- initialise.list.species(SL,"1D",name="HCO3",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="O2",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="SO4",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="Fe56_Fe",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="HS",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="CH4",phase="solute",type="chemical")

SL <- initialise.list.species(SL,"1D",name="H",phase="solute",type="chemical")
SL <- initialise.list.species(SL,"1D",name="H2O",phase="solute",type="chemical")

# Composite species incorporated in the model 

SL <- initialise.list.species(SL,"1D",name="POC",phase="solid",type="composite")
SL <- initialise.list.species(SL,"1D",name="FeOOH",phase="solid",type="composite")
SL <- initialise.list.species(SL,"1D",name="Fe56_FeOOH",phase="solid",type="composite")

# Artificial species incorporated in the model 

SL <- initialise.list.species(SL,"1D",name="Omega.FeS",phase="solute",type="artificial")

# Initialise species matrix 

spec.mat <- initialise.species.matrix(SL)
spec.mat["POC",paste("POC.",1:14,sep="")] <- rep(1,14)
spec.mat["FeOOH",c("FeOOH.f","FeOOH.a")] <- c(1,1)
spec.mat["Fe56_FeOOH",c("Fe56_FeOOH.f","Fe56_FeOOH.a")] <- c(1,1)

spec.table <- as.data.frame(spec.mat)
remove(spec.mat)

# List of reactions included 

RL <- list()

#---------------------
# mineralization terms
#---------------------

for (i in 1:14){
  RL <- initialise.list.reaction(RL,"1D",name=paste("Cmin.",i,sep=""),type="kinetic")
}

RL <- initialise.list.reaction(RL,"1D",name="Cmin",type="kinetic")

#for (i in 1:14){
#  RL <- initialise.list.reaction(RL,"1D",name=paste("Cmintot.",i,sep=""),type="kinetic")
#}

#---------------------
# individual mineralization terms
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="AR",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DIR",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="DIR_56Fe",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SR",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="MG",type="kinetic")

#---------------------
# aerobic oxidation terms
# excluding reactions with sulfide minerals
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="AMO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FIO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="FIO_56Fe",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SIO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SIO_56Fe",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="CSO",type="kinetic")

#---------------------
# anaerobic oxdiation terms
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="CSFO.f",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="CSFO.f_56Fe",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="CSFO.a",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="CSFO.a_56Fe",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="SMO",type="kinetic")

#---------------------
# sulfide mineral reactions
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="ISP",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="ISP_56Fe",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="ISD",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="ISD_56Fe",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="ISO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="ISO_56Fe",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="PyP",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="PyP_56Fe",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="PyO",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="PyO_56Fe",type="kinetic")

#---------------------
# iron oxide ageing
#---------------------

RL <- initialise.list.reaction(RL,"1D",name="IOA",type="kinetic")
RL <- initialise.list.reaction(RL,"1D",name="IOA_56Fe",type="kinetic")

#---------------------
# List of elements included 
#---------------------

EL <- list() # create empty element object list

EL <- initialise.list.element(EL,"1D",name="C")
EL <- initialise.list.element(EL,"1D",name="O")
EL <- initialise.list.element(EL,"1D",name="H")
EL <- initialise.list.element(EL,"1D",name="S")
EL <- initialise.list.element(EL,"1D",name="Fe")
EL <- initialise.list.element(EL,"1D",name="56Fe")
EL <- initialise.list.element(EL,"1D",name="e")

# Initialise element matrix 

elt.mat <- initialise.element.matrix(SL, EL)


for (i in 1:14){
  elt.mat[paste("POC.",i,sep=""),c("C","H","O")] <- c(1,2,1)
}

elt.mat["FeOOH.f",c("Fe","O","H")] <- c(1,2,1)
elt.mat["Fe56_FeOOH.f",c("56Fe","O","H")] <- c(1,2,1)
elt.mat["FeOOH.a",c("Fe","O","H")] <- c(1,2,1)
elt.mat["Fe56_FeOOH.a",c("56Fe","O","H")] <- c(1,2,1)
elt.mat["FeS",c("S","Fe")] <- c(1,1)
elt.mat["Fe56_FeS",c("S","56Fe")] <- c(1,1)
elt.mat["FeS2",c("S","Fe")] <- c(2,1)
elt.mat["Fe56_FeS2",c("S","56Fe")] <- c(2,1)

elt.mat["HCO3",c("H","C","O","e")] <- c(1,1,3,1)

elt.mat["O2",c("O")] <- c(2)
elt.mat["SO4",c("S","O","e")] <- c(1,4,2)

elt.mat["Fe",c("Fe","e")] <- c(1,-2)
elt.mat["Fe56_Fe",c("56Fe","e")] <- c(1,-2)
elt.mat["X.Fe",c("Fe","e")] <- c(1,-2)
elt.mat["X.Fe56_Fe",c("56Fe","e")] <- c(1,-2)
elt.mat["HS",c("H","S","e")] <- c(1,1,1)
elt.mat["CH4",c("C","H")] <- c(1,4)

elt.mat["H",c("H","e")] <- c(1,-1)
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
  
  if (is.null (PL$L)) PL$L <- 150   # depth of sediment domain [cm]
  if (is.null (PL$N)) PL$N <- 500  # number of grid layers
  
  # Call "setup.grid.1D" from ReacTran
  if (is.null (PL$grid.type)) PL$grid.type <- "uniform"  # number of grid layers
  if (is.null (PL$dx.1)) PL$dx.1 <- PL$L/(10*PL$N)
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
  
  if (is.null(PL$S))  PL$S  <- 34.2 # salinity
  if (is.null(PL$TC)) PL$TC <- 6    # temperature [deg C]
  if (is.null(PL$P))  PL$P  <- 500  # pressure [bar]
    
  # Porosity profile 
  
  if (is.null(PL$por.0))     PL$por.0     <- 0.948 # porosity at the sediment-water interface
  if (is.null(PL$por.inf))   PL$por.inf   <- 0.824 # asymptotic porosity at depth
  if (is.null(PL$por.x.att)) PL$por.x.att <- 3.6   # attenuation depth [cm]
    
  PL$por.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=PL$por.0,y.inf=PL$por.inf,x.L=0,x.att=PL$por.x.att)
  PL$svf.grid <- setup.prop.1D(func=p.exp,grid=PL$grid,y.0=(1-PL$por.0),y.inf=(1-PL$por.inf),x.L=0,x.att=PL$por.x.att)
  
  # Initialisation tortuosity 

  PL$tort.grid <- setup.prop.1D(value=0,grid=PL$grid)
  PL$tort.grid$mid <- tortuosity(PL$por.grid$mid)
  PL$tort.grid$int <- tortuosity(PL$por.grid$int)

  # Fixed profiles: pH
  
  if (is.null(PL$pH)) PL$pH      <- 7.5
  PL$pH.grid <-  setup.prop.1D(value=PL$pH,grid=PL$grid)
  
  #=============================================================================
  # Diffusion coefficients 
  # Uses routine 'diffcoeff' from 'marelac' package
  #=============================================================================
  
  c.fac <- 10000*(3600*24*365.25) # conversion from [m2 s-1] to [cm2 yr-1]
  
  # Diffusion coefficients from marelac package
  PL$D.HCO3.grid <- init.diffusion.coefficient(species="HCO3",   S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  PL$D.O2.grid   <- init.diffusion.coefficient(species="O2",     S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.SO4.grid  <- init.diffusion.coefficient(species="SO4",    S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  PL$D.Fe.grid   <- init.diffusion.coefficient(species="Fe",     S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.HS.grid   <- init.diffusion.coefficient(species="HS",     S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.CH4.grid  <- init.diffusion.coefficient(species="CH4",    S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  PL$D.H2O.grid  <- init.diffusion.coefficient(species="H2O",    S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  PL$D.H.grid    <- init.diffusion.coefficient(species="H",      S=PL$S,TC=PL$TC,P=PL$P,c.fac,PL$grid,PL$tort.grid)
  
  #=============================================================================
  # Transport parameters
  #=============================================================================
  
  # Diffusive boundary layer
  
  if (is.null (PL$x.DBL)) PL$x.DBL <- 0 # thickness of DBL [cm]
  
  # Advective velocities 
  
  if (is.null(PL$rho.sed)) PL$rho.sed <- 2.6  # density solid sediment [g cm-3]
  if (is.null(PL$u.inf))   PL$u.inf   <- 0.2  # sedimentation velocity pore water [cm yr-1]
  if (is.null(PL$v.inf))   PL$v.inf   <- 0.2  # sedimentation velocity solids [cm yr-1]
  
  #PL$SedFlux <- 0.1  # sedimentation flux [g cm-2 yr-1]
  
  PL$v.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$v # advection velocity solid
  PL$u.grid <- setup.compaction.1D(v.inf=PL$v.inf,por.0=PL$por.0,por.inf=PL$por.inf,por.grid=PL$por.grid)$u # advection velocity pore water
  
  # Bioturbation profile
  
  if (is.null(PL$x.L))  PL$x.L    <- 0.0
  if (is.null(PL$x.att))PL$x.att  <- 1.0
  if (is.null(PL$y.inf))PL$y.inf  <- 0.0
  if (is.null(PL$Db.0)) PL$Db.0   <- 0.0     # biodiffusion coefficient Db at SWI [cm2 yr-1]
  PL$Db.grid <- setup.prop.1D(func=p.sig,y.0=PL$Db.0,y.inf=PL$y.inf,x.att=PL$x.att,x.L=PL$x.L, grid = PL$grid)
  
  # Irrigation profile
  
  if (is.null(PL$irr.0))   PL$irr.0   <- 0.0     # irrigation rate at SWI [yr-1]
  if (is.null(PL$x.irr))   PL$x.irr   <- 0.0     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.at))  PL$irr.att <- 1.0 
  if (is.null(PL$irr.inf)) PL$irr.inf <- 0.0
  
  if (is.null(PL$irr.Fe))   PL$irr.Fe   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.HS))   PL$irr.HS   <- 1     # irrigation rate at SWI [yr-1]
  #if (is.null(PL$irr.NH4))  PL$irr.NH4  <- 1     # irrigation rate at SWI [yr-1]
  #if (is.null(PL$irr.N2))   PL$irr.N2   <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.HCO3)) PL$irr.HCO3 <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.SO4))  PL$irr.SO4  <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.O2))   PL$irr.O2   <- 1     # irrigation rate at SWI [yr-1]
  #if (is.null(PL$irr.H2))   PL$irr.H2   <- 1     # irrigation rate at SWI [yr-1]
  #if (is.null(PL$irr.NO3))  PL$irr.NO3  <- 1     # irrigation rate at SWI [yr-1]
  if (is.null(PL$irr.CH4))  PL$irr.CH4  <- 1     # irrigation rate at SWI [yr-1]
  
  PL$irr.grid <- setup.prop.1D(func=p.exp,y.0=PL$irr.0,y.inf=PL$irr.inf,x.L=PL$x.irr,x.att=PL$irr.att ,grid=PL$grid)
  #PL$irr.grid <- setup.prop.1D(func=p.sig,y.0=PL$irr.0,y.inf=PL$irr.inf,x.att=PL$irr.att,x.L=PL$x.irr, grid = PL$grid)
  
  #=============================================================================
  # Reaction parameters
  #=============================================================================
  
  # Decay constants organic matter (from Dale et al., 2015)
  
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
  
  if (is.null(PL$K_O2))   PL$K_O2    <- 0.001           # Monod constant O2 consumption [umol cm-3 or mM] (Meysman et al., 2003)
  if (is.null(PL$K_FeOOH))PL$K_FeOOH <- 12.0*PL$rho.sed  # Monod constant FeIII reduction [umol cm-3 or mM] (this paper -> 2 fractions, avoid too much overlap)
  if (is.null(PL$K_SO4))  PL$K_SO4   <- 0.9             # Monod constant SO4 reduction [umol cm-3 or mM] (Meysman et al., 2003)
  
  # aerobic oxidation reactions
  # excluding sulfide mineral reactions
  
  if (is.null(PL$k.FIO)) PL$k.FIO <- 1E+07  # kinetic constant Ferrous iron oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.SIO)) PL$k.SIO <- 1E+07  # kinetic constant Ferrous iron oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.CSO)) PL$k.CSO <- 1E+07  # kinetic constant sulfide oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996; Meysman et al., 2003)
  if (is.null(PL$k.AMO)) PL$k.AMO <- 1E+04  # kinetic constant aerobic methane oxidation [µmol-1 cm3 yr-1] (Meysman et al., 2003)

  # anaerobic oxidation reactions
  
  if (is.null(PL$k.CSFOf))PL$k.CSFOf <- 494     # kinetic constant sulfide oxidation by fresh FeIII [umol-1 cm3 yr-1] (Poulton et al., 2004)
  if (is.null(PL$k.CSFOa))PL$k.CSFOa <- 3.6     # kinetic constant sulfide oxidation by aged FeIII [umol-1 cm3 yr-1] (Poulton et al., 2004)
  if (is.null(PL$k.SMO))  PL$k.SMO   <- 1E+04   # kinetic constant Sulfate methane oxidation [µmol-1 cm3 yr-1] (-)
  
  # reactions involving sulfide minerals
  
  if (is.null(PL$k.ISO))  PL$k.ISO  <- 1E+07    # kinetic constant iron sulfide oxidation [umol-1 cm3 yr-1] (Van Cappellen and Wang, 1996)
  if (is.null(PL$k.PyO))  PL$k.PyO  <- 9.47    # kinetic constant Pyrite oxidation [umol-1 cm3 yr-1] (Berg et al., 2003)
  if (is.null(PL$k.PyP))  PL$k.PyP  <- 3.25    # kinetic constant Pyrite precipitation [umol-1 cm3 yr-1] (this paper)
    
  if (is.null(PL$k.ISP)) PL$k.ISP <- 1E+04   # kinetic constant FeS precipitation [umol-1 cm3 yr-1] (Meysman et al., 2015)
  if (is.null(PL$k.ISD)) PL$k.ISD <- 3       # kinetic constant FeS dissolution [yr-1] (Meysman et al., 2015)
  if (is.null(PL$n.ISP)) PL$n.ISP <- 1       # kinetic exponent FeS precipitation [] (Meysman et al., 2003)
  if (is.null(PL$n.ISD)) PL$n.ISD <- 1       # kinetic exponent FeS dissolution [] (Meysman et al., 2003)
  
  if (is.null(PL$K_IS))  PL$K_IS  <- (10^-PL$pH)*1000*3.16    # Saturation constant [umol-1 cm3] (Rickard, 2006)
  
  # iron oxide ageing and adsorption reactions 
  
  if (is.null(PL$k.IOA)) PL$k.IOA <- 0.57  # kinetic constant iron oxide ageing [yr-1] (Berg et al., 2003)
  
  #if (is.null(PL$k.FIS)) PL$k.FIS <- 1    # kinetic constant Ferrous iron sorption  [yr-1] (Berg et al., 2003)
  if (is.null(PL$K.FIS)) PL$K.FIS <- PL$rho.sed*268   # Equilibrium constant Ferrous iron sorption [] (Berg et al., 2003)
  
  #=============================================================================
  # Ison isotope system
  #=============================================================================
  
  if (is.null(PL$IRMM014))       PL$IRMM014       <- 15.81034 # Fe isotope reference (IUPAC report)
  if (is.null(PL$frac.limit))    PL$frac.limit    <- 0.001    # concentration below which no fractionation happens
  if (is.null(PL$d56Fe_Fe))      PL$d56Fe_Fe      <- 0.0      # overlying water Fe isotope signature
  if (is.null(PL$d56Fe_FeOOH.f)) PL$d56Fe_FeOOH.f <- 0.0      # influx FeOOH.f isotope signature
  if (is.null(PL$d56Fe_FeOOH.a)) PL$d56Fe_FeOOH.a <- 0.0      # influx FeOOH.a isotope signature
  if (is.null(PL$d56Fe_FeS))     PL$d56Fe_FeS     <- 0.0      # influx FeS isotope signature
  if (is.null(PL$d56Fe_FeS2))    PL$d56Fe_FeS2    <- 0.0      # influx FeS2 isotope signature
  if (is.null(PL$d56Fe_X.Fe))    PL$d56Fe_X.Fe    <- 0.0      # influx sorbed Fe isotope signature
  
  if (is.null(PL$a56Fe_DIR))     PL$a56Fe_DIR     <- (1.0 + 0.0/1000.0) # alpha for Dissimilatory Iron Reduction
  if (is.null(PL$a56Fe_CSFO))    PL$a56Fe_CSFO    <- (1.0 + 0.0/1000.0) # alpha for Dissimilatory Iron Reduction
  if (is.null(PL$a56Fe_FIO))     PL$a56Fe_FIO     <- (1.0 + 0.0/1000.0) # alpha for Ferrous Iron oxidation
  if (is.null(PL$a56Fe_ISP))     PL$a56Fe_ISP     <- (1.0 + 0.0/1000.0) # alpha for Iron sulphide precipitation
  if (is.null(PL$a56Fe_ISO))     PL$a56Fe_ISO     <- (1.0 + 0.0/1000.0) # alpha for Ferrous Iron oxidation
  if (is.null(PL$a56Fe_ISD))     PL$a56Fe_ISD     <- (1.0 + 0.0/1000.0) # alpha for Ferrous Iron oxidation
  if (is.null(PL$a56Fe_PyP))     PL$a56Fe_PyP     <- (1.0 + 0.0/1000.0) # alpha for Pyrite precipitation
  if (is.null(PL$a56Fe_PyO))     PL$a56Fe_PyO     <- (1.0 + 0.0/1000.0) # alpha for Pyrite precipitation
  if (is.null(PL$a56Fe_IOA))     PL$a56Fe_IOA     <- (1.0 + 0.0/1000.0) # alpha for Dissimilatory Iron Reduction
  
  #=============================================================================
  # Boundary conditions
  #=============================================================================
  
  # Upper boundary fluxes
  F.POC <- 20 # mmol m-2 d-1
  if (is.null(PL$F.POC))      PL$F.POC      <- F.POC*365.25/10      # conversion to umol cm-2 yr-1 
  
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
  
  F.FeS <- 0 # mmol m-2 d-1
  if (is.null(PL$F.FeS))       PL$F.FeS       <- F.FeS*365.25/10       # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.Fe56_FeS))  PL$F.Fe56_FeS  <- ((PL$d56Fe_FeS*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeS*1e-3+1.0)*PL$IRMM014))*PL$F.FeS # conversion to amount of 56Fe 
  F.FeS2 <- 0 # mmol m-2 d-1
  if (is.null(PL$F.FeS2))      PL$F.FeS2      <- F.FeS2*365.25/10      # conversion to umol cm-2 yr-1
  if (is.null(PL$F.Fe56_FeS2)) PL$F.Fe56_FeS2 <- ((PL$d56Fe_FeS2*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeS2*1e-3+1.0)*PL$IRMM014))*PL$F.FeS2 # conversion to amount of 56Fe 
  
  # 1% Fe content of mineral deposition 
  F.FeOOH.f <- 0.5 # mmol m-2 d-1
  if (is.null(PL$F.FeOOH.f))      PL$F.FeOOH.f      <- F.FeOOH.f*365.25/10      # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.Fe56_FeOOH.f)) PL$F.Fe56_FeOOH.f <- ((PL$d56Fe_FeOOH.f*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.f*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.f # conversion to amount of 56Fe 
  F.FeOOH.a <- 0.5 # mmol m-2 d-1
  if (is.null(PL$F.FeOOH.a))      PL$F.FeOOH.a      <- F.FeOOH.a*365.25/10      # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.Fe56_FeOOH.a)) PL$F.Fe56_FeOOH.a <- ((PL$d56Fe_FeOOH.a*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_FeOOH.a*1e-3+1.0)*PL$IRMM014))*PL$F.FeOOH.a  # conversion to amount of 56Fe 
  
  F.X.Fe <- 0 # mmol m-2 d-1
  if (is.null(PL$F.X.Fe))         PL$F.X.Fe         <- F.X.Fe*365.25/10         # conversion to umol cm-2 yr-1 
  if (is.null(PL$F.X.Fe56_Fe))    PL$F.X.Fe56_Fe    <- ((PL$d56Fe_X.Fe*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_X.Fe*1e-3+1.0)*PL$IRMM014))*PL$F.X.Fe  # conversion to amount of 56Fe 

  # Upper boundary concentrations

  if (is.null(PL$HCO3.ow))     PL$HCO3.ow    <- 2.1 # SumCO2 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$O2.ow))       PL$O2.ow      <- 0.28#0.0  # O2 concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$SO4.ow))      PL$SO4.ow     <- 28.2#0.01 #140.0*PL$S/((32.065+4*15.999)*1.80655)
  
  if (is.null(PL$Fe.ow))       PL$Fe.ow      <- 0.0#0.05
  if (is.null(PL$Fe56_Fe.ow))  PL$Fe56_Fe.ow <- ((PL$d56Fe_Fe*1e-3+1.0)*PL$IRMM014)/(1+((PL$d56Fe_Fe*1e-3+1.0)*PL$IRMM014))*PL$Fe.ow
  if (is.null(PL$HS.ow))       PL$HS.ow      <- 0.0   # SumH2S concentration bottom water [umol cm-3 or mM]
  if (is.null(PL$CH4.ow))      PL$CH4.ow     <- 0.0
  if (is.null(PL$H2O.ow))      PL$H2O.ow     <- 0
  if (is.null(PL$H.ow))        PL$H.ow       <- 0 # Proton concentration bottom water, without equilibrium with H2O
    
  # Boundary condition lower boundary 
  # All species: no gradient
  
  if (is.null(PL$HCO3.ds))     PL$HCO3.ds     <- NA
  if (is.null(PL$O2.ds))       PL$O2.ds       <- NA
  if (is.null(PL$SO4.ds))      PL$SO4.ds      <- NA
  if (is.null(PL$Fe.ds))       PL$Fe.ds       <- NA
  if (is.null(PL$Fe56_Fe.ds))  PL$Fe56_Fe.ds  <- NA
  if (is.null(PL$HS.ds))       PL$HS.ds       <- NA
  if (is.null(PL$CH4.ds))      PL$CH4.ds      <- NA
  if (is.null(PL$H2O.ds))      PL$H2O.ds      <- NA
  if (is.null(PL$H.ds))        PL$H.ds        <- NA
  
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
    
    #if (!("POC.1" %in% var.names)) POC.1$C <- POC.1.fix
    #if (!("POC.2" %in% var.names)) POC.2$C <- POC.2.fix
    #if (!("POC.3" %in% var.names)) POC.3$C <- POC.3.fix
    
    X.Fe$C <- K.FIS*Fe$C
    X.Fe56_Fe$C <- K.FIS*Fe56_Fe$C
    
    #-------------------------------------------------------------------------
    # Upper and lower boundary conditions
    #-------------------------------------------------------------------------

    HCO3$C.up    <- HCO3.ow
    O2$C.up      <- O2.ow
    SO4$C.up     <- SO4.ow
    
    Fe$C.up      <- Fe.ow
    Fe56_Fe$C.up <- Fe56_Fe.ow
    HS$C.up      <- HS.ow
    CH4$C.up     <- CH4.ow
    
    H2O$C.up     <- H2O.ow
    H$C.up       <- H.ow
    
    if (is.na(HCO3.ds))     HCO3$C.down     <- HCO3$C[N]    else HCO3$C.down    <- HCO3.ds 
    if (is.na(O2.ds))       O2$C.down       <- O2$C[N]      else O2$C.down      <- O2.ds 
    if (is.na(SO4.ds))      SO4$C.down      <- SO4$C[N]     else SO4$C.down     <- SO4.ds 
    
    if (is.na(Fe.ds))       Fe$C.down       <- Fe$C[N]      else Fe$C.down      <- Fe.ds 
    if (is.na(Fe56_Fe.ds))  Fe56_Fe$C.down  <- Fe56_Fe$C[N] else Fe56_Fe$C.down <- Fe56_Fe.ds 
    if (is.na(HS.ds))       HS$C.down       <- HS$C[N]      else HS$C.down      <- HS.ds 
    if (is.na(CH4.ds))      CH4$C.down      <- CH4$C[N]     else CH4$C.down     <- CH4.ds 
    if (is.na(H2O.ds))      H2O$C.down      <- H2O$C[N]     else H2O$C.down     <- H2O.ds 
    if (is.na(H.ds))        H$C.down        <- H$C[N]       else H$C.down       <- H.ds 
  
    #=========================================================================
    # Isotope ratios 
    #=========================================================================
    
    R_56FeOOH.f                                          <- Fe56_FeOOH.f$C/(FeOOH.f$C-Fe56_FeOOH.f$C)
    R_56FeOOH.f[(FeOOH.f$C-Fe56_FeOOH.f$C)<(C.lim*1e-9)] <- 0.0
    R_56FeOOH.f[FeOOH.f$C<(C.lim*1e-9)]                  <- 0.0
    R_56FeOOH.f[Fe56_FeOOH.f$C<(C.lim*1e-9)]             <- 0.0
    R_56FeOOH.f[R_56FeOOH.f<0.0]                         <- 0.0
    r_56FeOOH.f                                          <- Fe56_FeOOH.f$C/FeOOH.f$C
    r_56FeOOH.f[FeOOH.f$C<(C.lim*1e-9)]                  <- 0.0
    r_56FeOOH.f[Fe56_FeOOH.f$C<(C.lim*1e-9)]             <- 0.0
    r_56FeOOH.f[r_56FeOOH.f<0.0]                         <- 0.0
    
    R_56FeOOH.a                                          <- Fe56_FeOOH.a$C/(FeOOH.a$C-Fe56_FeOOH.a$C)
    R_56FeOOH.a[(FeOOH.a$C-Fe56_FeOOH.a$C)<(C.lim*1e-9)] <- 0.0
    R_56FeOOH.a[(FeOOH.a$C)<(C.lim*1e-9)]                <- 0.0
    R_56FeOOH.a[(Fe56_FeOOH.a$C)<(C.lim*1e-9)]           <- 0.0
    R_56FeOOH.a[R_56FeOOH.a<0.0]                         <- 0.0
    r_56FeOOH.a                                          <- Fe56_FeOOH.a$C/FeOOH.a$C
    r_56FeOOH.a[FeOOH.a$C<(C.lim*1e-9)]                  <- 0.0
    r_56FeOOH.a[Fe56_FeOOH.a$C<(C.lim*1e-9)]             <- 0.0
    r_56FeOOH.a[r_56FeOOH.a<0.0]                         <- 0.0
    
    R_56FeS                                              <- Fe56_FeS$C/(FeS$C-Fe56_FeS$C)
    R_56FeS[(FeS$C-Fe56_FeS$C)<(C.lim*1e-9)]             <- 0.0
    R_56FeS[(FeS$C)<(C.lim*1e-9)]                        <- 0.0
    R_56FeS[(Fe56_FeS$C)<(C.lim*1e-9)]                   <- 0.0
    R_56FeS[R_56FeS<0.0]                                 <- 0.0
    r_56FeS                                              <- Fe56_FeS$C/FeS$C
    r_56FeS[FeS$C<(C.lim*1e-9)]                          <- 0.0
    r_56FeS[Fe56_FeS$C<(C.lim*1e-9)]                     <- 0.0
    r_56FeS[r_56FeS<0.0]                                 <- 0.0
    
    R_56FeS2                                             <- Fe56_FeS2$C/(FeS2$C-Fe56_FeS2$C)
    R_56FeS2[(FeS2$C-Fe56_FeS2$C)<(C.lim*1e-9)]          <- 0.0
    R_56FeS2[(FeS2$C)<(C.lim*1e-9)]                      <- 0.0
    R_56FeS2[(Fe56_FeS2$C)<(C.lim*1e-9)]                 <- 0.0
    R_56FeS2[R_56FeS2<0.0]                               <- 0.0
    r_56FeS2                                             <- Fe56_FeS2$C/FeS2$C
    r_56FeS2[FeS2$C<(C.lim*1e-9)]                        <- 0.0
    r_56FeS2[Fe56_FeS2$C<(C.lim*1e-9)]                   <- 0.0
    r_56FeS2[r_56FeS2<0.0]                               <- 0.0
    
    R_56Fe                                               <- Fe56_Fe$C/(Fe$C-Fe56_Fe$C)
    R_56Fe[(Fe$C-Fe56_Fe$C)<(C.lim*1e-9)]                <- 0.0
    R_56Fe[(Fe$C)<(C.lim*1e-9)]                          <- 0.0
    R_56Fe[(Fe56_Fe$C)<(C.lim*1e-9)]                     <- 0.0
    R_56Fe[R_56Fe<0.0]                                   <- 0.0
    r_56Fe                                               <- Fe56_Fe$C/Fe$C
    r_56Fe[Fe$C<(C.lim*1e-9)]                            <- 0.0
    r_56Fe[Fe56_Fe$C<(C.lim*1e-9)]                       <- 0.0
    r_56Fe[r_56Fe<0.0]                                   <- 0.0
    
    R_56XFe                                              <- X.Fe56_Fe$C/(X.Fe$C-X.Fe56_Fe$C)
    R_56XFe[(X.Fe$C-X.Fe56_Fe$C)<(C.lim*1e-9)]           <- 0.0
    R_56XFe[(X.Fe$C)<(C.lim*1e-9)]                       <- 0.0
    R_56XFe[(X.Fe56_Fe$C)<(C.lim*1e-9)]                  <- 0.0
    R_56XFe[R_56XFe<0.0]                                 <- 0.0
    r_56XFe                                              <- X.Fe56_Fe$C/X.Fe$C
    r_56XFe[X.Fe$C<(C.lim*1e-9)]                         <- 0.0
    r_56XFe[X.Fe56_Fe$C<(C.lim*1e-9)]                    <- 0.0
    r_56XFe[r_56XFe<0.0]                                 <- 0.0
    
    #=========================================================================
    # TRANSPORT terms using tran.1D from ReacTran
    #=========================================================================
    
    if (summary.call) {ind <- c("tran","C.up","C.down","dif.flux","adv.flux","flux","flux.up","flux.down")} else {ind <- c("tran","flux.up","flux.down")}
    
    # Transport terms: Solids [umol cm-3 solid yr-1]
    
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
    
    FeOOH.f[ind]      <- tran.1D(C=FeOOH.f$C,     flux.up=F.FeOOH.f,     v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    Fe56_FeOOH.f[ind] <- tran.1D(C=Fe56_FeOOH.f$C,flux.up=F.Fe56_FeOOH.f,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeOOH.a[ind]      <- tran.1D(C=FeOOH.a$C,     flux.up=F.FeOOH.a,     v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    Fe56_FeOOH.a[ind] <- tran.1D(C=Fe56_FeOOH.a$C,flux.up=F.Fe56_FeOOH.a,v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeS[ind]          <- tran.1D(C=FeS$C,         flux.up=F.FeS,         v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    Fe56_FeS[ind]     <- tran.1D(C=Fe56_FeS$C,    flux.up=F.Fe56_FeS,    v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    FeS2[ind]         <- tran.1D(C=FeS2$C,        flux.up=F.FeS2,        v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    Fe56_FeS2[ind]    <- tran.1D(C=Fe56_FeS2$C,   flux.up=F.Fe56_FeS2,   v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    X.Fe[ind]         <- tran.1D(C=X.Fe$C,        flux.up=F.X.Fe,        v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    X.Fe56_Fe[ind]    <- tran.1D(C=X.Fe56_Fe$C,   flux.up=F.X.Fe56_Fe,   v=v.grid,D=Db.grid,VF=svf.grid,dx=grid,full.output=summary.call)
    
    # Transport terms: Solutes [umol cm-3 pore water yr-1]

    HCO3[ind]     <- tran.1D(C=HCO3$C,   C.up=HCO3$C.up,   C.down=HCO3$C.down,   v=u.grid,D=D.HCO3.grid,AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    O2[ind]       <- tran.1D(C=O2$C,     C.up=O2$C.up,     C.down=O2$C.down,     v=u.grid,D=D.O2.grid,  AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    SO4[ind]      <- tran.1D(C=SO4$C,    C.up=SO4$C.up,    C.down=SO4$C.down,    v=u.grid,D=D.SO4.grid, AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    Fe[ind]       <- tran.1D(C=Fe$C,     C.up=Fe$C.up,     C.down=Fe$C.down,     v=u.grid,D=D.Fe.grid,  AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    Fe56_Fe[ind]  <- tran.1D(C=Fe56_Fe$C,C.up=Fe56_Fe$C.up,C.down=Fe56_Fe$C.down,v=u.grid,D=D.Fe.grid,  AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    #Fe56_Fe[["tran"]] <- R_56Fe*Fe[["tran"]]; Fe56_Fe[["C.up"]] <- R_56Fe[1]*Fe[["C.up"]]; Fe56_Fe[["C.down"]] <- R_56Fe[N]*Fe[["C.down"]] 
    #Fe56_Fe[["flux.up"]] <- R_56Fe[1]*Fe[["flux.up"]]; Fe56_Fe[["flux.down"]] <- R_56Fe[1]*Fe[["flux.down"]]
    #if(summary.call){Fe56_Fe[["dif.flux"]] <- rep(NA,length(N+1)); Fe56_Fe[["adv.flux"]] <- rep(NA,length(N+1)); Fe56_Fe[["flux"]] <- rep(NA,length(N+1))}
    HS[ind]       <- tran.1D(C=HS$C,     C.up=HS$C.up,     C.down=HS$C.down,     v=u.grid,D=D.HS.grid,  AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    CH4[ind]      <- tran.1D(C=CH4$C,    C.up=CH4$C.up,    C.down=CH4$C.down,    v=u.grid,D=D.CH4.grid, AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    
    H2O[ind]      <- tran.1D(C=H2O$C,    C.up=H2O$C.up,    C.down=H2O$C.down,    v=u.grid,D=D.H2O.grid, AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
    H[ind]        <- tran.1D(C=H$C,      C.up=H$C.up,      C.down=H$C.down,      v=u.grid,D=D.H.grid,   AFDW=1,VF=por.grid,dx=grid,full.output=summary.call)
        
    #-----------------------------------------------------------------------------
    # Saturation state calculation
    #-----------------------------------------------------------------------------
    
    # Saturation state FeS 
    
    Omega.FeS$C      <- Fe$C*HS$C/(K_IS)#*FSAT(Fe$C,C.lim/100,5)*FSAT(HS$C,C.lim/100,5)
    Omega.FeS$C.up   <- Fe.ow*HS$C.up/K_IS
    Omega.FeS$C.down <- Fe$C.down*HS$C.down/K_IS
    
    #=========================================================================
    # Rates of Kinetic reactions 
    # Expressed per volume of bulk sediment [umol cm-3 yr-1]
    #=========================================================================
    
    # Mineralisation reactions
    
#    Cmin.f$R <- svf.grid$mid*k.f*OC.f$C*FSAT(OC.f$C,C.lim,5)
#    Cmin.s$R <- svf.grid$mid*k.s*OC.s$C*FSAT(OC.s$C,C.lim,5)

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
    
    Cmin$R <- Cmin.1$R + Cmin.2$R + Cmin.3$R + Cmin.4$R + Cmin.5$R + Cmin.6$R + Cmin.7$R + Cmin.8$R + Cmin.9$R + Cmin.10$R + Cmin.11$R + Cmin.12$R + Cmin.13$R + Cmin.14$R  
    
    O2.lim    <- O2$C/(O2$C+K_O2)*FSAT(O2$C,C.lim,5) #*(O2$C>sqrt(.Machine$double.eps))
    O2.inh    <- K_O2/(O2$C+K_O2)
    FeOOH.lim <- FeOOH.f$C/(FeOOH.f$C+K_FeOOH)*FSAT(FeOOH.f$C,C.lim,5)
    FeOOH.inh <- K_FeOOH/(FeOOH.f$C+K_FeOOH)
    SO4.lim   <- SO4$C/(SO4$C+K_SO4)*FSAT(SO4$C,C.lim,5)
    SO4.inh   <- K_SO4/(SO4$C+K_SO4)
    
    a.f <- O2.lim/(O2.lim + FeOOH.lim*O2.inh + SO4.lim*FeOOH.inh*O2.inh + SO4.inh*FeOOH.inh*O2.inh)
    f.f <- (FeOOH.lim*O2.inh)/(O2.lim + FeOOH.lim*O2.inh + SO4.lim*FeOOH.inh*O2.inh + SO4.inh*FeOOH.inh*O2.inh)
    s.f <- (SO4.lim*FeOOH.inh*O2.inh)/(O2.lim + FeOOH.lim*O2.inh + SO4.lim*FeOOH.inh*O2.inh + SO4.inh*FeOOH.inh*O2.inh)
    m.f <- (SO4.inh*FeOOH.inh*O2.inh)/(O2.lim + FeOOH.lim*O2.inh + SO4.lim*FeOOH.inh*O2.inh + SO4.inh*FeOOH.inh*O2.inh)
    
    AR$R  <- a.f*Cmin$R
    DIR$R <- f.f*Cmin$R 
    DIR_56Fe$R[FeOOH.f$C>=frac.limit] <- (a56Fe_DIR*R_56FeOOH.f[FeOOH.f$C>=frac.limit])/(1.0+a56Fe_DIR*R_56FeOOH.f[FeOOH.f$C>=frac.limit])*DIR$R[FeOOH.f$C>=frac.limit]# FeIII reduction
    DIR_56Fe$R[FeOOH.f$C<frac.limit]  <- r_56FeOOH.f[FeOOH.f$C<frac.limit]*DIR$R[FeOOH.f$C<frac.limit]# FeIII reduction
    SR$R  <- s.f*Cmin$R
    MG$R  <- m.f*Cmin$R 
    
    # Aerobic oxidation reactions
    
    #AAO$R  <- por.grid$mid*k.AAO*O2$C*FSAT(O2$C,C.lim/100,5)*NH4$C*FSAT(NH4$C,C.lim/100,5) # ammonium oxidation
    #SAO$R  <- svf.grid$mid*k.SAO*O2$C*FSAT(O2$C,C.lim/100,5)*X.NH4$C*FSAT(X.NH4$C,C.lim/100,5) # Sorbed ammonium oxidation
    
    FIO$R                         <- por.grid$mid*k.FIO*O2$C*FSAT(O2$C,C.lim,5)*Fe$C*FSAT(Fe$C,C.lim/100,5) # iron re-oxidation
    FIO_56Fe$R[Fe$C>=frac.limit]  <- (a56Fe_FIO*R_56Fe[Fe$C>=frac.limit])/(1+a56Fe_FIO*R_56Fe[Fe$C>=frac.limit])*FIO$R[Fe$C>=frac.limit]# FeIII reduction
    FIO_56Fe$R[Fe$C<frac.limit]   <- r_56Fe[Fe$C<frac.limit]*FIO$R[Fe$C<frac.limit]
    SIO$R                         <- svf.grid$mid*k.SIO*O2$C*FSAT(O2$C,C.lim,5)*X.Fe$C*FSAT(X.Fe$C,C.lim,5) # sorbed iron re-oxidation
    SIO_56Fe$R[X.Fe$C>=frac.limit] <- (a56Fe_FIO*R_56XFe[X.Fe$C>=frac.limit])/(1+a56Fe_FIO*R_56XFe[X.Fe$C>=frac.limit])*SIO$R[X.Fe$C>=frac.limit]# FeIII reduction
    SIO_56Fe$R[X.Fe$C<frac.limit]  <- r_56XFe[X.Fe$C<frac.limit]*SIO$R[X.Fe$C<frac.limit]
    
    CSO$R  <- por.grid$mid*k.CSO*O2$C*FSAT(O2$C,C.lim/100,5)*HS$C*FSAT(HS$C,C.lim/100,5) # Sulfide re-oxidation
    AMO$R  <- por.grid$mid*k.AMO*O2$C*FSAT(O2$C,C.lim/100,5)*CH4$C*FSAT(CH4$C,C.lim/100,5) # aerobic methane oxidation
 
    # Anaerobic oxidation reactions
    
    CSFO.f$R                             <- svf.grid$mid*k.CSFOf*FeOOH.f$C*FSAT(FeOOH.f$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim,5)# Sulfide re-oxidation by FeIII
    CSFO.f_56Fe$R[FeOOH.f$C>=frac.limit] <- (a56Fe_CSFO*R_56FeOOH.f[FeOOH.f$C>=frac.limit])/(1+a56Fe_CSFO*R_56FeOOH.f[FeOOH.f$C>=frac.limit])*CSFO.f$R[FeOOH.f$C>=frac.limit]# FeIII reduction
    CSFO.f_56Fe$R[FeOOH.f$C<frac.limit]  <- r_56FeOOH.f[FeOOH.f$C<frac.limit]*CSFO.f$R[FeOOH.f$C<frac.limit]
    
    CSFO.a$R                             <- svf.grid$mid*k.CSFOa*FeOOH.a$C*FSAT(FeOOH.a$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim,5)# Sulfide re-oxidation by FeIII
    CSFO.a_56Fe$R[FeOOH.a$C>=frac.limit] <- (a56Fe_CSFO*R_56FeOOH.a[FeOOH.a$C>=frac.limit])/(1+a56Fe_CSFO*R_56FeOOH.a[FeOOH.a$C>=frac.limit])*CSFO.a$R[FeOOH.a$C>=frac.limit]# FeIII reduction
    CSFO.a_56Fe$R[FeOOH.a$C<frac.limit]  <- r_56FeOOH.a[FeOOH.a$C<frac.limit]*CSFO.a$R[FeOOH.a$C<frac.limit]
    
    SMO$R  <- por.grid$mid*k.SMO*SO4$C*FSAT(SO4$C,C.lim/100,5)*CH4$C*FSAT(CH4$C,C.lim/100,5) # Sulfate methane oxidation
    
    # Reactions involving sulfide minerals
   
    ISO$R                         <- svf.grid$mid*k.ISO*O2$C*FSAT(O2$C,C.lim,5)*FeS$C*FSAT(FeS$C,C.lim,5)*(FeS$C>0)*(O2$C>0) # iron sulfide re-oxidation
    ISO_56Fe$R[FeS$C>=frac.limit] <- (a56Fe_ISO*R_56FeS[FeS$C>=frac.limit])/(1+a56Fe_ISO*R_56FeS[FeS$C>=frac.limit])*ISO$R[FeS$C>=frac.limit]# FeIII reduction
    ISO_56Fe$R[FeS$C<frac.limit]  <- r_56FeS[FeS$C<frac.limit]*ISO$R[FeS$C<frac.limit]
    
    PyO$R                          <- svf.grid$mid*k.PyO*O2$C*FSAT(O2$C,C.lim,5)*FeS2$C*FSAT(FeS2$C,C.lim,5)*(FeS2$C>0)*(O2$C>0) # Pyrite re-oxidation
    PyO_56Fe$R[FeS2$C>=frac.limit] <- (a56Fe_PyO*R_56FeS2[FeS2$C>=frac.limit])/(1+a56Fe_PyO*R_56FeS2[FeS2$C>=frac.limit])*PyO$R[FeS2$C>=frac.limit]# FeIII reduction
    PyO_56Fe$R[FeS2$C<frac.limit]  <- r_56FeS2[FeS2$C<frac.limit]*PyO$R[FeS2$C<frac.limit]
    
    PyP$R                         <- svf.grid$mid*k.PyP*FeS$C*FSAT(FeS$C,C.lim,5)*HS$C*FSAT(HS$C,C.lim,5)*FSAT(SO4$C,C.lim,5)*(FeS$C>0)*(HS$C>0)*(SO4$C>0) # sorbed iron re-oxidation
    PyP_56Fe$R[FeS$C>=frac.limit] <- (a56Fe_PyP*R_56FeS[FeS$C>=frac.limit])/(1+a56Fe_PyP*R_56FeS[FeS$C>=frac.limit])*PyP$R[FeS$C>=frac.limit]# FeIII reduction
    PyP_56Fe$R[FeS$C<frac.limit]  <- r_56FeS[FeS$C<frac.limit]*PyP$R[FeS$C<frac.limit]
    
    ISP$R                         <- svf.grid$mid*k.ISP*((Omega.FeS$C-1)^n.ISP)*(Omega.FeS$C>1)*FSAT(HS$C,C.lim,5)*FSAT(Fe$C,C.lim,5)*(Fe$C>0)*(HS$C>0)              #iron sulfide precipitation
    ISP_56Fe$R[Fe$C>=frac.limit]  <- (a56Fe_ISP*R_56Fe[Fe$C>=frac.limit])/(1+a56Fe_ISP*R_56Fe[Fe$C>=frac.limit])*ISP$R[Fe$C>=frac.limit]# FeIII reduction
    ISP_56Fe$R[Fe$C<frac.limit]   <- r_56Fe[Fe$C<frac.limit]*ISP$R[Fe$C<frac.limit]
    
    ISD$R                         <- svf.grid$mid*k.ISD*FeS$C*FSAT(FeS$C,C.lim,5)*(1-Omega.FeS$C)^n.ISD*(Omega.FeS$C<1)*(FeS$C>0) #iron sulfide dissolution
    ISD_56Fe$R[FeS$C>=frac.limit] <- (a56Fe_ISD*R_56FeS[FeS$C>=frac.limit])/(1+a56Fe_ISD*R_56FeS[FeS$C>=frac.limit])*ISD$R[FeS$C>=frac.limit]# FeIII reduction
    ISD_56Fe$R[FeS$C<frac.limit]  <- r_56FeS[FeS$C<frac.limit]*ISD$R[FeS$C<frac.limit]
    
    # Iron oxide ageing

    IOA$R                             <- svf.grid$mid*k.IOA*FeOOH.f$C*FSAT(FeOOH.f$C,C.lim,5) # Iron oxide ageing
    IOA_56Fe$R[FeOOH.f$C>=frac.limit] <- (a56Fe_IOA*R_56FeOOH.f[FeOOH.f$C>=frac.limit])/(1+a56Fe_IOA*R_56FeOOH.f[FeOOH.f$C>=frac.limit])*IOA$R[FeOOH.f$C>=frac.limit]# FeIII reduction
    IOA_56Fe$R[FeOOH.f$C<frac.limit]  <- r_56FeOOH.f[FeOOH.f$C<frac.limit]*IOA$R[FeOOH.f$C<frac.limit]
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
    
    FeOOH.f$kin.reac      <- (- 4*DIR$R      + 4*FIO$R      + 4*SIO$R      + ISO$R      - 8*CSFO.f$R      + PyO$R      - IOA$R     )/svf.grid$mid  
    Fe56_FeOOH.f$kin.reac <- (- 4*DIR_56Fe$R + 4*FIO_56Fe$R + 4*SIO_56Fe$R + ISO_56Fe$R - 8*CSFO.f_56Fe$R + PyO_56Fe$R - IOA_56Fe$R)/svf.grid$mid  
    FeOOH.a$kin.reac      <- (- 8*CSFO.a$R      + IOA$R     )/svf.grid$mid    
    Fe56_FeOOH.a$kin.reac <- (- 8*CSFO.a_56Fe$R + IOA_56Fe$R)/svf.grid$mid    
    FeS$kin.reac          <- (- ISO$R      - PyP$R     )/svf.grid$mid
    Fe56_FeS$kin.reac     <- (- ISO_56Fe$R - PyP_56Fe$R)/svf.grid$mid
    FeS2$kin.reac         <- (PyP$R      - PyO$R     )/svf.grid$mid
    Fe56_FeS2$kin.reac    <- (PyP_56Fe$R - PyO_56Fe$R)/svf.grid$mid
    X.Fe$kin.reac         <- (- 4*SIO$R     )/svf.grid$mid
    X.Fe56_Fe$kin.reac    <- (- 4*SIO_56Fe$R)/svf.grid$mid
    
    HCO3$kin.reac<- (AR$R + DIR$R + SR$R + 1/2*MG$R + AMO$R + SMO$R)/por.grid$mid
    O2$kin.reac  <- (- AR$R - FIO$R - SIO$R - 2*CSO$R - 9/4*ISO$R - 15/4*PyO$R - 2*AMO$R)/por.grid$mid
    SO4$kin.reac <- (- 0.5*SR$R + CSFO.f$R  + CSFO.a$R + CSO$R + ISO$R - 0.25*PyP$R + 2*PyO$R - SMO$R)/por.grid$mid
    
    Fe$kin.reac       <- (+ 4*DIR$R - 4*FIO$R      + 8*CSFO.f$R      + 8*CSFO.a$R     )/por.grid$mid
    Fe56_Fe$kin.reac  <- (+ 4*DIR_56Fe$R - 4*FIO_56Fe$R + 8*CSFO.f_56Fe$R + 8*CSFO.a_56Fe$R)/por.grid$mid
    HS$kin.reac       <- (+ 0.5*SR$R - CSFO.f$R - CSFO.a$R - CSO$R - 0.75*PyP$R + SMO$R)/por.grid$mid  
    CH4$kin.reac      <- (+ 1/2*MG$R - AMO$R - SMO$R)/por.grid$mid
    
    H2O$kin.reac <- (9*DIR$R - 1/4*MG$R - 6*FIO$R - 6*SIO$R + 12*CSFO.f$R + 12*CSFO.a$R - 3/2*ISO$R -5/4*PyP$R - 5/2*PyO$R + 2*AMO$R + 2*SMO$R)/por.grid$mid
    H$kin.reac   <- (+ AR$R - DIR$R + SR$R - MG$R + 8*FIO$R + 8*SIO$R + CSO$R - 7*CSFO.f$R - 7*CSFO.a$R + 2*ISO$R + PyP$R + 4*PyO$R - SMO$R)/por.grid$mid

    #=========================================================================
    # Equilibrium reactions
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================

    FeS$eq.reac      <- (ISP$R        - ISD$R     )/svf.grid$mid
    Fe56_FeS$eq.reac <- (ISP_56Fe$R   - ISD_56Fe$R)/svf.grid$mid
    Fe$eq.reac       <- (- ISP$R      + ISD$R     )/por.grid$mid
    Fe56_Fe$eq.reac  <- (- ISP_56Fe$R + ISD_56Fe$R)/por.grid$mid
    HS$eq.reac       <- (- ISP$R      + ISD$R)/por.grid$mid
    H$eq.reac        <- (  ISP$R      - ISD$R)/por.grid$mid

    #=========================================================================
    # Irrigation 
    # Solutes [umol cm-3 pore water yr-1] Solids [umol cm-3 solid yr-1] 
    #=========================================================================

    HCO3$irr    <- irr.grid$mid*irr.HCO3*(HCO3.ow - HCO3$C)
    O2$irr      <- irr.grid$mid*irr.O2*(O2.ow - O2$C)
    SO4$irr     <- irr.grid$mid*irr.SO4*(SO4.ow - SO4$C)

    Fe$irr      <- irr.grid$mid*irr.Fe*(Fe.ow - Fe$C)
    Fe56_Fe$irr <- irr.grid$mid*irr.Fe*(Fe56_Fe.ow - Fe56_Fe$C)
    HS$irr      <- irr.grid$mid*irr.HS*(HS.ow - HS$C)
    CH4$irr     <- irr.grid$mid*irr.CH4*(CH4.ow - CH4$C)
    H2O$irr     <- irr.grid$mid*(H2O.ow - H2O$C)
    
    #=========================================================================
    # Rate of changes composite species 
    #=========================================================================
    
    POC$tran     <- POC.1$tran     + POC.2$tran     + POC.3$tran     + POC.4$tran     + POC.5$tran     + POC.6$tran     + POC.7$tran     + POC.8$tran     + POC.9$tran     + POC.10$tran     + POC.11$tran     + POC.12$tran     + POC.13$tran     + POC.14$tran
    POC$kin.reac <- POC.1$kin.reac + POC.2$kin.reac + POC.3$kin.reac + POC.4$kin.reac + POC.5$kin.reac + POC.6$kin.reac + POC.7$kin.reac + POC.8$kin.reac + POC.9$kin.reac + POC.10$kin.reac + POC.11$kin.reac + POC.12$kin.reac + POC.13$kin.reac + POC.14$kin.reac
    POC$irr      <- POC.1$irr      + POC.2$irr      + POC.3$irr      + POC.4$irr      + POC.5$irr      + POC.6$irr      + POC.7$irr      + POC.8$irr      + POC.9$irr      + POC.10$irr      + POC.11$irr      + POC.12$irr      + POC.13$irr      + POC.14$irr
    POC$eq.reac  <- POC.1$eq.reac  + POC.2$eq.reac  + POC.3$eq.reac  + POC.4$eq.reac  + POC.5$eq.reac  + POC.6$eq.reac  + POC.7$eq.reac  + POC.8$eq.reac  + POC.9$eq.reac  + POC.10$eq.reac  + POC.11$eq.reac  + POC.12$eq.reac  + POC.13$eq.reac  + POC.14$eq.reac

    FeOOH$tran          <- FeOOH.f$tran          + FeOOH.a$tran
    Fe56_FeOOH$tran     <- Fe56_FeOOH.f$tran     + Fe56_FeOOH.a$tran
    FeOOH$kin.reac      <- FeOOH.f$kin.reac      + FeOOH.a$kin.reac
    Fe56_FeOOH$kin.reac <- Fe56_FeOOH.f$kin.reac + Fe56_FeOOH.a$kin.reac
    FeOOH$irr           <- FeOOH.f$irr           + FeOOH.a$irr
    Fe56_FeOOH$irr      <- Fe56_FeOOH.f$irr      + Fe56_FeOOH.a$irr
    FeOOH$eq.reac       <- FeOOH.f$eq.reac       + FeOOH.a$eq.reac
    Fe56_FeOOH$eq.reac  <- Fe56_FeOOH.f$eq.reac  + Fe56_FeOOH.a$eq.reac
    
    #=========================================================================
    # Rate of changes of all species
    #=========================================================================
    
  # Fe$irr <- alfa*(PL$Fe.ow - Fe$C)

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
    
    FeOOH.f$ddt      <- FeOOH.f$tran      + FeOOH.f$kin.reac      + FeOOH.f$eq.reac      + FeOOH.f$irr
    Fe56_FeOOH.f$ddt <- Fe56_FeOOH.f$tran + Fe56_FeOOH.f$kin.reac + Fe56_FeOOH.f$eq.reac + Fe56_FeOOH.f$irr
    FeOOH.a$ddt      <- FeOOH.a$tran      + FeOOH.a$kin.reac      + FeOOH.a$eq.reac      + FeOOH.a$irr
    Fe56_FeOOH.a$ddt <- Fe56_FeOOH.a$tran + Fe56_FeOOH.a$kin.reac + Fe56_FeOOH.a$eq.reac + Fe56_FeOOH.a$irr
    FeS$ddt          <- FeS$tran          + FeS$kin.reac          + FeS$eq.reac          + FeS$irr
    Fe56_FeS$ddt     <- Fe56_FeS$tran     + Fe56_FeS$kin.reac     + Fe56_FeS$eq.reac     + Fe56_FeS$irr
    FeS2$ddt         <- FeS2$tran         + FeS2$kin.reac         + FeS2$eq.reac         + FeS2$irr
    Fe56_FeS2$ddt    <- Fe56_FeS2$tran    + Fe56_FeS2$kin.reac    + Fe56_FeS2$eq.reac    + Fe56_FeS2$irr
    
    HCO3$ddt  <- HCO3$tran  + HCO3$kin.reac  + HCO3$eq.reac  + HCO3$irr
    O2$ddt    <- O2$tran    + O2$kin.reac    + O2$eq.reac    + O2$irr
    SO4$ddt   <- SO4$tran   + SO4$kin.reac   + SO4$eq.reac   + SO4$irr

    Fe$ddt       <- 1/(por.grid$mid+svf.grid$mid*K.FIS)*(por.grid$mid*(Fe$tran      + Fe$kin.reac      + Fe$eq.reac      + Fe$irr)      + svf.grid$mid*(X.Fe$tran      + X.Fe$kin.reac      + X.Fe$eq.reac      + X.Fe$irr))
    #Fe$ddt       <- Fe$tran      + Fe$kin.reac      + Fe$eq.reac      + Fe$irr
    Fe56_Fe$ddt  <- 1/(por.grid$mid+svf.grid$mid*K.FIS)*(por.grid$mid*(Fe56_Fe$tran + Fe56_Fe$kin.reac + Fe56_Fe$eq.reac + Fe56_Fe$irr) + svf.grid$mid*(X.Fe56_Fe$tran + X.Fe56_Fe$kin.reac + X.Fe56_Fe$eq.reac + X.Fe56_Fe$irr))
    #Fe56_Fe$ddt  <- Fe56_Fe$tran + Fe56_Fe$kin.reac + Fe56_Fe$eq.reac + Fe56_Fe$irr
    HS$ddt    <- HS$tran     + HS$kin.reac     + HS$eq.reac     + HS$irr
    CH4$ddt   <- CH4$tran    + CH4$kin.reac    + CH4$eq.reac    + CH4$irr
    H2O$ddt   <- H2O$tran    + H2O$kin.reac    + H2O$eq.reac    + H2O$irr
    H$ddt     <- H$tran      + H$kin.reac      + H$eq.reac      + H$irr
    
    POC$ddt        <- POC$tran        + POC$kin.reac        + POC$eq.reac        + POC$irr
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
  
  PL$var.names <- c(paste("POC.",1:14,sep=""),"FeOOH.f","Fe56_FeOOH.f","FeOOH.a","Fe56_FeOOH.a","FeS","Fe56_FeS","FeS2","Fe56_FeS2","X.Fe56_Fe","X.Fe",
                    "HCO3","O2","SO4","Fe","Fe56_Fe","HS","CH4","H2O","H")
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
    
    if ("FeOOH.f" %in% PL$var.names)      SV.guess[,"FeOOH.f"]      <- rep(0,length.out=N)
    if ("Fe56_FeOOH.f" %in% PL$var.names) SV.guess[,"Fe56_FeOOH.f"] <- rep(0,length.out=N)
    if ("FeOOH.a" %in% PL$var.names)      SV.guess[,"FeOOH.a"]      <- rep(0,length.out=N)
    if ("Fe56_FeOOH.a" %in% PL$var.names) SV.guess[,"Fe56_FeOOH.a"] <- rep(0,length.out=N)
    if ("FeS" %in% PL$var.names)          SV.guess[,"FeS"]          <- rep(0,length.out=N)
    if ("Fe56_FeS" %in% PL$var.names)     SV.guess[,"Fe56_FeS"]     <- rep(0,length.out=N)
    if ("FeS2" %in% PL$var.names)         SV.guess[,"FeS2"]         <- rep(0,length.out=N)
    if ("Fe56_FeS2" %in% PL$var.names)    SV.guess[,"Fe56_FeS2"]    <- rep(0,length.out=N)
    if ("X.Fe" %in% PL$var.names)         SV.guess[,"X.Fe"]         <- rep(0,length.out=N)
    if ("X.Fe56_Fe" %in% PL$var.names)    SV.guess[,"X.Fe56_Fe"]    <- rep(0,length.out=N)
    
    if ("HCO3" %in% PL$var.names)    SV.guess[,"HCO3"]    <- rep(HCO3.ow,length.out=N)
    if ("O2" %in% PL$var.names)      SV.guess[,"O2"]      <- rep(O2.ow,length.out=N)
    if ("SO4" %in% PL$var.names)     SV.guess[,"SO4"]     <- rep(SO4.ow,length.out=N)
    if ("Fe" %in% PL$var.names)      SV.guess[,"Fe"]      <- rep(Fe.ow,length.out=N)
    if ("Fe56_Fe" %in% PL$var.names) SV.guess[,"Fe56_Fe"] <- rep(Fe56_Fe.ow,length.out=N)
    if ("HS" %in% PL$var.names)      SV.guess[,"HS"]      <- rep(HS.ow,length.out=N)
    if ("CH4" %in% PL$var.names)     SV.guess[,"CH4"]     <- rep(CH4.ow,length.out=N)
    if ("H2O" %in% PL$var.names)     SV.guess[,"H2O"]     <- rep(H2O.ow,length.out=N)
    if ("H" %in% PL$var.names)       SV.guess[,"H"]       <- rep(H.ow,length.out=N)
    
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

