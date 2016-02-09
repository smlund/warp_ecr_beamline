# Parameters for FRIB front-end simulations using the x-y slice model 
# described in frib-front-xy.py 
# Notes:
#  * include inputs describing beam AND numerical (mesh, advance) parameters  

 
#########################
# Begin Inputs 
#########################

# Define species
#   Operate = Operating species targeted to be transported (e.g. U)
#   Support = Support   species used in ECR (e.g. O in U operation) 
# 
#     Syntax of later use:
#
#       for sp in species_name:
#         sp.ppzx()
#
#     Dictionaries also setup. 

# Setup U operation with O support gas 

# --- Define charge states 
Operate_charge_states = [33,34,25,26,27,28,29,30,31,32,35,36,37,38,39,40]
Support_charge_states = [1,2,3,4]

# --- Create species consistent with charge states 
Operate_species = [Species(type=Uranium,charge_state=i,name='U%d'%i) for i in Operate_charge_states]
Support_species = [Species(type=Oxygen, charge_state=i,name='O%d'%i) for i in Support_charge_states]

# --- Count species 
Operate_ns = len(Operate_species) 
Support_ns = len(Support_species) 

# --- Make abbreviated dictionary of species "sp" for cleaner use in later diagnostics 
#       sp.keys()
sp_Operate = {Operate_species[i].name:Operate_species[i] for i in range(Operate_ns)}
sp_Support = {Support_species[i].name:Support_species[i] for i in range(Support_ns)}

sp = {}
sp.update(sp_Operate)
sp.update(sp_Support)

# --- Assign species colors to help distinguish on plots: 
#        target U species: magenta/green
#        O        species: blue 
#        High   U species: cyan 
#        Low    U species: red 
sp['O1'].color = "blue"
sp['O2'].color = "blue"
sp['O3'].color = "blue"
sp['O4'].color = "blue"

sp['U25'].color = "red"
sp['U26'].color = "red"
sp['U27'].color = "red"
sp['U28'].color = "red"
sp['U29'].color = "red"
sp['U30'].color = "red"
sp['U31'].color = "red"
sp['U32'].color = "red"

sp['U33'].color = "magenta"
sp['U34'].color = "green"

sp['U35'].color = "cyan"
sp['U36'].color = "cyan"
sp['U37'].color = "cyan"
sp['U38'].color = "cyan"
sp['U39'].color = "cyan"
sp['U40'].color = "cyan"

# --- Target species tuple for later use in diagnostics 
sp_target = ['U33','U34']

# --- Set species properties for injection 
#        The beam kinetic energy (ekin) and axial velocity (vbeam) should not both
#        be set unless done so consistently.  If one is zero the code sets from
#        the other on generation:
#           vbeam = 0    => set vbeam from ekin
#           ekin  = 0    => set ekin  from vbeam
#

# --- --- unneutralized electric currents ... array elements corresponding to charge state arrays  
Operate_ibeam = array([0.210,0.205,0.035,0.051,0.068,0.088,0.115,0.150,0.175,0.192,0.178,0.142,0.110,0.072,0.043,0.031])*mA 
Support_ibeam = array([0.300,0.300,0.300,0.200])*mA 

# --- --- kinetic energy 
SourceBias = 35.*keV  # source voltage: set for Q_ref*SourceBias/A =>  4.9264706 keV/u for U 

Operate_ekin = array(Operate_charge_states)*SourceBias
Support_ekin = array(Support_charge_states)*SourceBias 

# --- --- initial ion thermal temp/phase-space area 
# 
#  Approach: Allow setting of thermal phase-space area by either a temperature specification 
#            OR a normalized rms emittance measure. Convert input values to a normalized rms 
#            edge emittance to work with the Warp loader.  
#
#  init_emit_spec = "emitn"    => set by init_emitn 
#                   or "temp"  => set by init_temp 
# 
#  init_emitn = Initial ion thermal rms normalized emittance [m-rad] ** NOT edge measure **
#  init_temp  = Initial ion temp [eV] 
#
#  Notes:
#  * Guilliaume:  Ions likely 2-3 eV and electrons few to 100s of keV.  
#                 Ions not equilibrated with electrons.
#  * Previous simulations used 0.4*mm*mr thermal phase-space area launched outside of magnetic fields.      
#  * Here we set only the thermal component of the emittance and the distribution is loaded without 
#    correlation terms. Later we adjust correlation terms to set the value of the canonical angular 
#    momentum as desired based on field values at launch and birth.    

init_emit_spec = "temp" 
init_temp  = 3.          # initial temp [eV] 
init_emitn = 0.4*mm*mr   # initial rms normalized emittance [m-rad]  

# --- --- initial beam size: elliptical cross-section beam 
#   
#   r_x  = 2*sqrt(<x^2>)           Initial x-plane beam rms edge radius [m] 
#   rp_x = 2*<x*x'>/sqrt(<x^2>)    Initial rms envelope angle   [rad] 
# 
#   r_y, rp_y                      Analogous y-plane measures 
#
# Note:
# * Now set r_x etc by fraction of ECR aperture size.  
# * Previous simulations used betatron functions with a specific emittance value to get 
#   desired beam size. 
# * Was set by betatron functions and an emittance specification.  Commented out code (left 
#   in for reference).   
#

# ####################### OLD ... set by betatron funcs ###############
#alpha_x = 0.
#beta_x  = 12.9*cm
#gamma_x = (1. + alpha_x**2)/beta_x 
#
#alpha_y = 0.
#beta_y  = 12.9*cm
#gamma_y = (1. + alpha_y**2)/beta_y
#
#emitn_edge = 0.4*mm*mr   # norm rms edge emittance used to set beam size in earlier simulations 
#
#v_z_ref   = sqrt(2.*jperev*Q_ref*SourceBias/m_ref)  # nonrel approx ref z-velocity 
#gamma_ref = 1./sqrt(1.-(v_z_ref/clight)**2)         # ref axial gamma factor (approx 1) 
#emit_edge = emitn_edge/(gamma_ref*v_z_ref/clight)   # unnormalized rms edge emittance 
#
#r_x  = sqrt(emit_edge*beta_x)             # envelope x-edge 
#r_y  = sqrt(emit_edge*beta_y)             # envelope y-edge 
#rp_x = -sqrt(emit_edge/beta_x)*alpha_x    # envelope x-angle 
#rp_y = -sqrt(emit_edge/beta_y)*alpha_y    # envelope y-angle 
# #######################

r_extractor = 4.*mm   # ECR extraction aperture radius [m]
                      #   Value both for Venus and Artemis A,B ECR Sources 

r_x = r_extractor  
r_y = r_extractor
rp_x = 0. 
rp_y = 0. 


## --- transverse thermal velocity and energy (eV) of nonrel ref particle from emittance 
#
#vt = v_z_ref*emit_edge/(2.*r_x) 
#Et = 0.5*m_ref*vt**2/jperev 
#
## --- intrinsic thermal emittance scale 
#Et_therm = 3.   # Guilliaume's estimated ion temp scale (eV) 
#vt_therm = sqrt(2.*(jperev*Et_therm)/m_ref)
#emit_therm  = 2.*r_x*vt_therm/v_z_ref
#emitn_therm = (gamma_ref*v_z_ref/clight)*emit_therm
#
## Ratio of thermal to edge emittance suggests value of P_theta contributing to effective emittance 
## emit_therm/emit_edge = 0.10  => most beam PS area from P_theta  

# --- Set properties of initial species load 

for i in range(Operate_ns):
  Osp = Operate_species[i]
  ekin_i  = Operate_ekin[i]
  betab_i = sqrt(2.*jperev*ekin_i/Osp.sm)/clight
  if init_emit_spec == "emitn":
    emitn_therm_edge = 4.*init_emitn 
  elif init_emit_spec == "temp": 
    emitn_therm_edge = betab_i*((2.*r_x)/sqrt(2.))*sqrt(init_temp/ekin_i)
  else:
    raise Exception("Error: init_emit_spec not set properly") 
  #
  Osp.ekin   = ekin_i              # kinetic energy of beam particle [eV]
  Osp.vbeam  = 0.                  # beam axial velocity [m/sec] (set from ekin if 0) 
  Osp.ibeam  = Operate_ibeam[i]    # beam current [amps] 
  Osp.emitnx = emitn_therm_edge    # beam x-emittance, rms edge [m-rad] 
  Osp.emitny = emitn_therm_edge    # beam y-emittance, rms edge [m-rad]
  Osp.vthz   = 0.                  # axial velocity spread [m/sec] 

for i in range(Support_ns):
  Ssp = Support_species[i]
  ekin_i  = Support_ekin[i]
  betab_i = sqrt(2.*jperev*ekin_i/Ssp.sm)/clight
  if init_emit_spec == "emitn":
    emitn_therm_edge = 4.*init_emitn 
  elif init_emit_spec == "temp": 
    emitn_therm_edge = betab_i*((2.*r_x)/sqrt(2.))*sqrt(init_temp/ekin_i)
  else:
    raise Exception("Error: init_emit_spec not set properly") 
  # 
  Ssp.ekin   = ekin_i  
  Ssp.vbeam  = 0.
  Ssp.ibeam  = Support_ibeam[i]
  Ssp.emitnx = emitn_therm_edge 
  Ssp.emitny = emitn_therm_edge
  Ssp.vthz   = 0.


#
# Beam centroid and rms envelope initial conditions at s=0      
#    
#   x0:   initial x-centroid xc = <x> [m]
#   y0:   initial y-centroid yc = <y> [m]
#   xp0:  initial x-centroid angle xc' = <x'> = d<x>/ds [rad]
#   yp0:  initial y-centroid angle yc' = <y'> = d<y>/ds [rad]
#
#   a0:   initial x-envelope edge a = 2*sqrt(<(x-xc)^2>) [m]
#   b0:   initial y-envelope edge b = 2*sqrt(<(y-yc)^2>) [m]
#   ap0:  initial x-envelope angle ap = a' = d a/ds [rad]
#   bp0:  initial y-envelope angle bp = b' = d b/ds [rad]

for i in range(Operate_ns):
  Osp = Operate_species[i]
  # --- centroid 
  Osp.x0  = 0.
  Osp.y0  = 0.
  Osp.xp0 = 0.
  Osp.yp0 = 0.
  # --- envelope 
  Osp.a0   = r_x
  Osp.b0   = r_y
  Osp.ap0  = rp_x
  Osp.bp0  = rp_y

for i in range(Support_ns):
  Ssp = Support_species[i]
  # --- centroid 
  Ssp.x0  = 0.   
  Ssp.y0  = 0.   
  Ssp.xp0 = 0.   
  Ssp.yp0 = 0.   
  # --- envelope 
  Ssp.a0   = r_x              
  Ssp.b0   = r_y           
  Ssp.ap0  = rp_x  
  Ssp.bp0  = rp_y  

#
# Define transverse simulation grid
#

# --- Symmetries.  Set for increased statistical efficiency.  These should
#     only be used in cases where lattice symmetries and initial beam
#     conditions and loads allow.  For no symmetry, set both options to false.
w3d.l2symtry = false     # 2-fold symmetry
w3d.l4symtry = false     # 4-fold symmetry

# -- Grid extent and increments
#     Chosen without consideration of symmetry. Consistently reset later 
#     if symmetry operations are set.  

r_grid = 8.*cm    # radial extent of grid 
n_grid = 400      # number grid cells (no symmetries) 


#
# Particle moving
#

z_launch  = 66.540938  # Axial position in lattice to launch  beam 
                       #   = ecr_z_extr defined in lattice file. 
z_adv     = 69.2       # Axial position in lattice to advance to 

ds = 2.*mm             # Axial advance step 

# 
# Particle loading 
#

nppg = 100.    # number of particles per grid cell


# Distribution loads type 
#  Comment: 
#   * See stptcl routine for how loading works. 
#   * Neutralization is NOT taken into account in the loading.    
#   * At present, all species are loaded with the same value of distrbtn.  
#   * rms equivalent beam loaded with the specified distribution from:  
#     KV => KV Distribution 
#
#     SG => semi-Gaussian distribution 
#             (KV density and local Gaussian angle spread about KV flutter) 
#
#     TE => Pseudoequilibrium with Thermal  Equilibrium form 
#     WB => Pseudoequilibrium with Waterbag Equilibrium form
#             The Pseudoequilibrium distributions use continuous focused 
#             equilibrium forms which are canoically transformed to AG 
#             symmetry of the lattice. 
#
#     For more info on loads, see review paper:
#       Lund, Kikuchi, and Davidson, PRSTAB 12, 114801 (2009) 


#w3d.distrbtn = "KV"          # initial KV distribution
#w3d.distrbtn = "TE"          # initial thermal distribution
w3d.distrbtn = "WB"           # initial waterbag distribution
#w3d.distrbtn = "SG"          # initial semi-Gaussian distribution 
