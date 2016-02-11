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

# --- Make charge state dictionary 
Operate_charge_states_dict = {'U%d'%i: i for i in Operate_charge_states} 
Support_charge_states_dict = {'O%d'%i: i for i in Support_charge_states} 

charge_states = {} 
charge_states.update(Operate_charge_states_dict)
charge_states.update(Support_charge_states_dict) 

# --- Count species 
Operate_ns = len(Operate_species) 
Support_ns = len(Support_species) 
nspecies = Operate_ns + Support_ns 

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

# --- Target species list for later use in diagnostics 
sp_target = ['U33','U34']

# --- Set species properties for injection 
#        The beam kinetic energy (ekin) and axial velocity (vbeam) should not both
#        be set unless done so consistently.  If one is zero the code sets from
#        the other on generation:
#           vbeam = 0    => set vbeam from ekin
#           ekin  = 0    => set ekin  from vbeam
#

# --- --- unneutralized electric currents
#           Set in a named dictionary ibeam. Data from ECR measurement data.    
Operate_ibeam = array([0.210,0.205,0.035,0.051,0.068,0.088,0.115,0.150,0.175,0.192,0.178,0.142,0.110,0.072,0.043,0.031])*mA 
Support_ibeam = array([0.300,0.300,0.300,0.200])*mA 

ibeam = {
'U25': 0.035*mA,
'U26': 0.051*mA,
'U27': 0.068*mA,
'U28': 0.088*mA,
'U39': 0.115*mA,
'U30': 0.150*mA,
'U31': 0.175*mA,
'U32': 0.192*mA,
'U33': 0.210*mA,  
'U34': 0.205*mA,
'U35': 0.178*mA,
'U36': 0.142*mA,
'U37': 0.110*mA,
'U38': 0.072*mA,
'U39': 0.043*mA,
'U40': 0.031*mA,
'O1': 0.300*mA,
'O2': 0.300*mA,
'O3': 0.300*mA,
'O4': 0.200*mA
         }

# --- --- kinetic energy
#          Set as Q*SourceBias in a dictionary ekin 
SourceBias = 35.*keV  # source voltage: set for Q_ref*SourceBias/A =>  4.9264706 keV/u for U 

Operate_ekin = array(Operate_charge_states)*SourceBias
Support_ekin = array(Support_charge_states)*SourceBias 

ekin = {key: SourceBias*charge_states[key] for key in sp.keys()} 

# --- --- initial ion thermal phase-space area 
#           Set in species-named dictionaries
#            * Can set individual values as in ibeam. 
# 
#  Approach: Allow setting of thermal phase-space area by either a temperature specification 
#            OR a normalized rms emittance measure. Convert input values to a normalized rms 
#            edge emittance to work with the Warp loader.  
#
#  init_ps_spec = "emitn"    => set by init_emitn 
#                 or "temp"  => set by init_temp 
# 
#  emitn = Initial ion thermal rms normalized emittance [m-rad] ** NOT edge measure **
#  temp  = Initial ion temp [eV] 
#
#  Notes:
#  * Guilliaume:  Ions likely 2-3 eV and electrons few to 100s of keV.  
#                 Ions not equilibrated with electrons.
#  * Previous simulations used 0.4*mm*mr thermal phase-space area launched outside of magnetic fields.      
#  * Here we set only the thermal component of the emittance and the distribution is loaded without 
#    correlation terms. Later we adjust correlation terms to set the value of the canonical angular 
#    momentum as desired based on field values at launch and birth. 
#  * Syntax is used to allow different beam sizes and temperatures for each species. Initial runs 
#    may employ a common value for all species    

init_ps_spec = "temp" 
temp_ion  = 3.          # initial ion temp [eV] 
emitn     = 0.4*mm*mr   # initial rms normalized emittance [m-rad]  ** NOT edge measure **

Operate_temp = temp_ion*ones(Operate_ns)  
Support_temp = temp_ion*ones(Support_ns) 

temp = {key: temp_ion for key in sp.keys()} 

Operate_emitnx = emitn*ones(Operate_ns)
Operate_emitny = emitn*ones(Operate_ns)

Support_emitnx = emitn*ones(Support_ns)
Support_emitny = emitn*ones(Support_ns)

emitnx = {key: emitn for key in sp.keys()}
emitny = {key: emitn for key in sp.keys()}

# --- --- initial beam centroid 
# 
#   xc  = <x>   Initial centroid coordinate [m]
#   xcp = <x'>  Initial centroid angle [rad] 
# 
#               Analogous y-plane measures 
#

Operate_xc  = zeros(Operate_ns) 
Operate_xcp = zeros(Operate_ns)

Operate_yc  = zeros(Operate_ns) 
Operate_ycp = zeros(Operate_ns)

Support_xc  = zeros(Operate_ns) 
Support_xcp = zeros(Operate_ns)

Support_yc  = zeros(Operate_ns) 
Support_ycp = zeros(Operate_ns)

xc  = {key: 0. for key in sp.keys()}
yc  = {key: 0. for key in sp.keys()}

xcp = {key: 0. for key in sp.keys()}
ycp = {key: 0. for key in sp.keys()}


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
#   desired beam size. Code used was (y-plane similar):   
#     alpha_x = 0.
#     beta_x  = 12.9*cm
#     gamma_x = (1. + alpha_x**2)/beta_x 
#
#     emitn_edge = 0.4*mm*mr   # norm rms edge emittance used to set beam size 
#
#     v_z_ref   = sqrt(2.*jperev*Q_ref*SourceBias/m_ref)  # nonrel approx ref z-velocity 
#     gamma_ref = 1./sqrt(1.-(v_z_ref/clight)**2)         # ref axial gamma factor (approx 1) 
#     emit_edge = emitn_edge/(gamma_ref*v_z_ref/clight)   # unnormalized rms edge emittance 
#
#     r_x  = sqrt(emit_edge*beta_x)             # envelope x-edge 
#     rp_x = -sqrt(emit_edge/beta_x)*alpha_x    # envelope x-angle 

r_extractor = 4.*mm   # ECR extraction aperture radius [m]
                      #   Value both for Venus and Artemis A,B ECR Sources 

Operate_r_x = r_extractor*ones(Operate_ns)  
Operate_r_y = r_extractor*ones(Operate_ns)

Operate_rp_x = 0.*ones(Operate_ns) 
Operate_rp_y = 0.*ones(Operate_ns) 

Support_r_x = r_extractor*ones(Operate_ns)  
Support_r_y = r_extractor*ones(Operate_ns)

Support_rp_x = 0.*ones(Operate_ns) 
Support_rp_y = 0.*ones(Operate_ns) 

rx  = {key: r_extractor for key in sp.keys()}
ry  = {key: r_extractor for key in sp.keys()}

rxp = {key: 0. for key in sp.keys()}
ryp = {key: 0. for key in sp.keys()}


# Specify ion birth method to set beam canonical angular momentum. 
# brith_mode = 0: Ions are born at launch point
#              1: Ions are born with a canonical angular momentum specified by the user
#              2: Ions are born in average field between the two peaks in B-field
#              3: Ions are born randomly between the two peaks in B-field
#
#  * Angular velocities are injected to the ions to conserve the 
#    statistical canonical angular momentum of the beam in method 2 & 3.  This only 
#    takes into account the solenoidal compoenent of the ECR field.  

birth_mode = 2

## Input for birth_mode == 1
rms_birth = 2.82
bz0_birth = 1.2

## Reference: rms_launch = 2.8 mm, bz0_launch = 2.15 T

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
