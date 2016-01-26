"""
Python input script for a Warp xy slice simulation of 
the FRIB front end.  For help and/or additional information contact:

     Steve Lund     lund@frib.msu.edu    (517) 908-7291
     Jonathan Wong  wong@nscl.msu.edu    (517) 908-7465 

For documentation see the Warp web-site:

     https://warp.lbl.gov

All code inputs are mks with the exception of particle kinetic energy (eV).   

To run this Warp script:
  Interactive (return to interpreter after executing):

    % python -i frib-front-xy.py 

  Non-interactive (exit interpreter after executing):

    % python frib-front-xy.py

This setup scripts for these simulations are large and have been broken 
apart into functional sections to aid maintainability.
 ** To be implemented **   

    frib-front-xy.py        :  Main setup/control for xy slice simulations 
    frib-front-xy-parms.py  :  Parameters for runs 

    frib-front-lat.py       :  Lattice setup/description (3D, xy etc) 
    frib-front-lat-diag.py  :  Diagnostics for lattice     

    frib-front-xy-load.py   :  Initial beam specification 
    frib-front-xy-diag.py   :  Diagnostics 

    frib-front-env.py       :  Axisymmetric Envelope solver
    frib-front-env-diag.py  :  Diagnostics for envelope solver  

"""

# Load Warp and various script packages 
from warp        import *               # Warp code 
from errorcheck  import checksymmetry   # Check for input errors
#from runcounters import *               # Counter for parametric runs 

# Set informational labels included on all output cgm plots.   
top.pline2   = "xy Slice Simulation: FRIB Front End" 
top.pline1   = " "   # Add more info, if desired.  

# Parameters 
mr = 0.001  # milli-rad unit 

# Invoke setup routine for graphics and output files (THIS IS MANDATORY)
setup()

# Set runmaker - included in informational labels on output plots
top.runmaker = "Wong and Lund"

# Beam parameters for simulation
#
#   Other than numerical parameters, the physical simulation is specified
#   by the numbers input immediately below.  


# --- Define species
# 
#     Syntax of later use:
#
#       for us in U_species:
#         us.ppzx()
#
#     Dictionaries also setup. 

# --- Define charge states 
U_charge_states = [33,34,25,26,27,28,29,30,31,32,35,36,37,38,39,40]
O_charge_states = [1,2,3,4]

# --- Create species consistent with charge states 
U_species = [Species(type=Uranium,charge_state=i,name='U%d'%i) for i in U_charge_states]
O_species = [Species(type=Oxygen, charge_state=i,name='O%d'%i) for i in O_charge_states]

# --- Count species 
U_ns = len(U_species) 
O_ns = len(O_species) 

# --- Make abbreviated dictionary of species "sp" for cleaner use in later diagnostics 
#       sp.keys()
sp_U = {U_species[i].name:U_species[i] for i in range(U_ns)}
sp_O = {O_species[i].name:O_species[i] for i in range(O_ns)}

sp = {}
sp.update(sp_U)
sp.update(sp_O)

# --- Setup for variable particle weights in slice code. 
#     Notes: 
#      * All particles will carry individual weight set with v-velocity for full slice 
#          consistency using pid arrays to store the initial axial velocity. This is 
#          setup by hand in slice code for multi-species due to code structure.  
#      * Way done will be consistent with axial velocity spread.   
#      * Will be put in weight adjustment each step using beforeloadrho() 
#      * pid array elements hold particle properties and are consistently 
#          mirrored when particles scraped  
#      * nextpid() gets next pid array index 
#  
# Via Dave Grote: need to setup scaling by hand for multi-species 
#  top.wpid = nextpid()     # setup/allocate pid array on generate 
#
#  species.getw()  = s.w   gets weights of species  
#                          [equivalent to species.w] (top.wpid must be set to work)  
#  species.getweights() returns product of species.sw*species.w 
# 
#  beam.w * beam.sw = total weight = beam.getweights() 
# 

top.wpid = nextpid()       # pid index for variable weights: initialize on generate  
uzp0pid  = nextpid() - 1   # pid index for initial uz to scale weights: initialize on generate  

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

# --- Target species tuple for later use 
sp_target = ['U33','U34']

# --- Reference species as average of target 
A_ref = 0. 
Q_ref = 0. 
for ii in sp_target:
  s = sp[ii]
  A_ref += s.mass/amu 
  Q_ref += s.charge/echarge 
A_ref = A_ref/len(sp_target) 
Q_ref = Q_ref/len(sp_target)  

m_ref = A_ref*amu 
q_ref = Q_ref*echarge 

# --- Set species properties for injection 
#        The beam kinetic energy (ekin) and axial velocity (vbeam) should not both
#        be set unless done so consistently.  If one is zero the code sets from
#        the other on generation:
#           vbeam = 0    => set vbeam from ekin
#           ekin  = 0    => set ekin  from vbeam
#

# --- --- unneutralized electric currents ... array elements corresponding to charge state arrays  
U_ibeam = array([0.210,0.205,0.035,0.051,0.068,0.088,0.115,0.150,0.175,0.192,0.178,0.142,0.110,0.072,0.043,0.031])*mA 
O_ibeam = array([0.300,0.300,0.300,0.200])*mA 

# --- --- kinetic energy 
SourceBias = 35.*keV  # source voltage: set for Q_ref*SourceBias/A =>  4.9264706 keV/u for U 

U_ekin = array(U_charge_states)*SourceBias
O_ekin = array(O_charge_states)*SourceBias 

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

for i in range(U_ns):
  Usp = U_species[i]
  ekin_i  = U_ekin[i]
  betab_i = sqrt(2.*jperev*ekin_i/Usp.sm)/clight
  if init_emit_spec == "emitn":
    emitn_therm_edge = 4.*init_emitn 
  elif init_emit_spec == "temp": 
    emitn_therm_edge = betab_i*((2.*r_x)/sqrt(2.))*sqrt(init_temp/ekin_i)
  else:
    raise Exception("Error: init_emit_spec not set properly") 
  #
  Usp.ekin   = ekin_i              # kinetic energy of beam particle [eV]
  Usp.vbeam  = 0.                  # beam axial velocity [m/sec] (set from ekin if 0) 
  Usp.ibeam  = U_ibeam[i]          # beam current [amps] 
  Usp.emitnx = emitn_therm_edge    # beam x-emittance, rms edge [m-rad] 
  Usp.emitny = emitn_therm_edge    # beam y-emittance, rms edge [m-rad]
  Usp.vthz   = 0.                  # axial velocity spread [m/sec] 

for i in range(O_ns):
  Osp = O_species[i]
  ekin_i  = O_ekin[i]
  betab_i = sqrt(2.*jperev*ekin_i/Osp.sm)/clight
  if init_emit_spec == "emitn":
    emitn_therm_edge = 4.*init_emitn 
  elif init_emit_spec == "temp": 
    emitn_therm_edge = betab_i*((2.*r_x)/sqrt(2.))*sqrt(init_temp/ekin_i)
  else:
    raise Exception("Error: init_emit_spec not set properly") 
  # 
  Osp.ekin   = ekin_i  
  Osp.vbeam  = 0.
  Osp.ibeam  = O_ibeam[i]
  Osp.emitnx = emitn_therm_edge 
  Osp.emitny = emitn_therm_edge
  Osp.vthz   = 0.

# Calculate vbeam and other species quantities for defined species 
derivqty()

# ---  Calculate and printout Q/M by species and store in a dictionary 
sp_qovm = {}
print("Species Charge to Mass Ratios:")
for ii in sort(sp.keys()):
  s = sp[ii]
  qovm = (s.charge/echarge)/(s.mass/amu)
  sp_qovm.update({ii:qovm})
  print("   Species: "+ii+" Q/A = %s"%qovm)

# ---  Calculate and printout injected rigidity [B rho] by species and store in a dictionary 
sp_brho = {}
print("Species Rigidity:")
for ii in sort(sp.keys()):
  s = sp[ii]
  gamma = 1./sqrt(1.-(s.vbeam/clight)**2)
  brho  = gamma*s.mass*s.vbeam/s.charge   # Define rigidity with mean axial velocity 
  sp_brho.update({ii:brho})
  print("   Species: "+ii+" [B rho] = %s T-m"%brho)


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

for i in range(U_ns):
  Usp = U_species[i]
  # --- centroid 
  Usp.x0  = 0.
  Usp.y0  = 0.
  Usp.xp0 = 0.
  Usp.yp0 = 0.
  # --- envelope 
  Usp.a0   = r_x
  Usp.b0   = r_y
  Usp.ap0  = rp_x
  Usp.bp0  = rp_y

for i in range(O_ns):
  Osp = O_species[i]
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


#
# Setup Lattice  
#
#  Read in lattice description from auxillary script file 

execfile("frib-front-lat.py") 


# Define transverse simulation grid

# --- Symmetries.  Set for increased statistical efficiency.  These should
#     only be used in cases where lattice symmetries and initial beam
#     conditions and loads allow.  For no symmetry, set both options to false.
w3d.l2symtry = false     # 2-fold symmetry
w3d.l4symtry = false     # 4-fold symmetry

# -- Grid increments
#      First choose number grid cells without symmetry and reset 
#      consistent with symmetry options

n_grid = 400 # 200 previous # number grid cells (no symmetries) 

sym_x = 1.
sym_y = 1.
if w3d.l4symtry:
  sym_x = 0.5
  sym_y = 0.5
elif w3d.l2symtry:
  sym_x = 0.5 

w3d.nx = int(sym_x*n_grid) 
w3d.ny = int(sym_y*n_grid)

# ---- Grid bounds 
#      Some bounds will be reset to zero by code on generation
#      if symmetry options are set.
l_grid = 2.*r_p               # length edge of simulation grid [m]      
w3d.xmmax =  l_grid/2.        # x-grid max limit [m] 
w3d.xmmin = -l_grid/2.        # x-grid min limit [m] 
w3d.ymmax =  l_grid/2.        # y-grid max limit [m] 
w3d.ymmin = -l_grid/2.        # y-grid min limit [m] 

# --- grid increments to use before code generation in setup
dx = l_grid/float(n_grid)


# Particle loading
#
# Set simulation macro-particle number (top.npmax) by specifying the
#  number of macro-particles to load per xy grid cell (nppg) using an
# rms equivalent uniform density beam measure.  This number is (re)set
# consistently with the symmetry options.  
#   sp['name'].np = particle #
#   top.nplive = total number live particles = len(sp)*top.npmax at time of load 
# 

nppg = 100.    # number of particles per grid cell
top.npmax = int(nppg*pi*(r_x*r_y)/dx**2*sym_x*sym_y) # max initial particles loaded (each species) 

# Distribution loads type 
#  Comment: 
#     See stptcl routine for how loading works.   
#     At present, all species are loaded with the same value of distrbtn.   zbeam can be 
#     set before the generate for the beam location. 

z_launch  = ecr_z_extr #+ 10.*cm 
top.zbeam = z_launch                 # present z of simulation, reset consistently 

# rms equivalent beam loaded with the specified distribution form 
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
w3d.distrbtn = "WB"          # initial waterbag distribution
#w3d.distrbtn = "SG"          # initial semi-Gaussian distribution 

# --- random number options to use in loading 
w3d.xrandom  = "pseudo" # "digitrev"    # load x,y,z  with digitreverse random numbers 
w3d.vtrandom = "pseudo" # "digitrev"    # load vx, vy with digitreverse random numbers
w3d.vzrandom = "pseudo" # "digitrev"    # load vz     with digitreverse random numbers 
w3d.cylinder = true          # load a cylinder

# Particle moving
#

z_adv = 69.2   # Axial position in lattice to advance to 

top.lrelativ   =  false    # turn off relativistic kinematics
top.relativity = 0         # turn off relativistic self-field correction
                           #   to account for approx diamagnetic B-field of beam

wxy.ds = 2.*mm             # ds for part adv [m] 
wxy.lvzchang = true        # Use iterative stepping, which is needed if
                           # the vz of the particles changes.
                           #  ... must change even in linear lattice 
                           #          for high-order energy conservation 
top.ibpush   = 2           # magnetic field particle push, 
                           #   0 - off, 1 - fast, 2 - accurate 


# Setup field solver using 2d multigrid field solver. 

w3d.boundxy = 0              # Neuman boundary conditions on edge of grid.
w3d.solvergeom = w3d.XYgeom  # fieldsolve type to 2d multigrid 

# --- Uncomment to turn off space-charge deposition for simulation of particles 
#     moving in applied field  
#top.depos = "none"


# Turn on x-window plots, if desired; use winkill() to close interactively.  
#winon()


################################
# Particle simulation
################################

# Generate the xy PIC code.  In the generate, particles are allocated and
# loaded consistent with initial conditions and load parameters
# set previously.  Particles are advanced with the step() command later
# after various diagnostics are setup.
package("wxy"); generate()

# Read in diagnostics for applied lattice fields 
execfile("frib-front-lat-diag.py") 

# Install conducting aperture on mesh
for i in aperture:
  installconductors(i,dfill=largepos) # will trigger "list has no attribute" error if aperture is used in installconductors

# Check that inputs are consistent with symmetries (errorcheck package function)
checksymmetry()

# ?? Why this here ??
# --- Weight setup 
for s in sp.values():       
  s.w   = 1.      # Need full charge: set relative weight to unity 

# Carry out an initial unneutralized field solve with conducting pipe after generate 
loadrho() 
fieldsolve() 

# Setup variable weight species needs for neutralization and acceleration 


# --- set initial weight factors consistent with neutralization factor 
#       w0 = initial weight (same for all particles in species) 
#       species.w = array of variable weight factors 
for s in sp.values():       
  s.w0  = 1.-neut_f1
  #s.w   = 1.-neut_f1   #?? why this commented out ?? 
  s.sw0    = s.sw       # save initial sw    (in case later changed)  
  s.vbeam0 = s.vbeam    # save initial vbeam (in case later changed)  

# --- save initial uzp for all species at once 
top.pgroup.pid[:,uzp0pid] = top.pgroup.uzp

# --- adjust weights  
@callfrombeforeloadrho
def adjustweights():
  for s in sp.values():
    s.w[:] = s.w0*s.pid[:,uzp0pid]/s.uzp

# Fix intitial history diagnostics to account for species weight changes
top.jhist = top.jhist-1   # needed to be minus 1 to reset save in right postion 
from getzmom import *
zmmnt() 
savehist(0.) 

# Read in diagnostics for slice run 
execfile("frib-front-xy-diag.py") 

diag_hist_hl()   # make sure initial diagnostics saved before any steps 

# Make plot of initial unneutralized beam potential profile 
         
diag_plt_phi_ax(label="Initial Unneutralized Beam Potential at y,x = 0 b,r") 
fma()

diag_plt_phi_ax(label="Initial Unneutralized Beam Potential at y,x = 0 b,r",xmax=1.5*r_x)
fma()

# Carry out explicit fieldsolve with adjusted rho consistent with neutralization 
loadrho() 
fieldsolve()


# Make plot of initial neutralized beam potential profile 
         
diag_plt_phi_ax(label="Initial f = %s Neutralized Beam Potential at y,x = 0 b,r"%(neut_f1))
fma()

diag_plt_phi_ax(label="Initial f = %s Neutralized Beam Potential at y,x = 0 b,r"%(neut_f1),xmax=1.5*r_x)
fma()

# Loading ion at launch point based on different assumptions on their birth in the ECR
# Modify intital distribution loaded on generate to include canonical angular momentum

execfile("frib-front-xy-load.py")

#raise Exception("to here")


# Make plot of initial Brho by species 
plt_diag_bro(label = "Initial Rigidity by Species") 

#raise exception("to here")


# Install diagnostic calls after simulation step
installafterstep(diag_calls)

# Step 0 diagnostics (if any) of the initial distribution loaded 
diag_calls() 

# Advance simulation specified steps 

#raise "to here"

# ---- to grated accel gap  
n_step = nint((neut_z1-z_launch)/wxy.ds) 
top.prwall = r_p_up    # consistent aperture 
step(n_step)

# --- reset species weights to turn off neutralization  
for s in sp.values():
   s.w0 = 1. 

loadrho()     # applies adjusted species weights  
fieldsolve()  # make field consistent with turned off neutralization 

# --- unneutralized advance in acceleration column  
n_step = nint((neut_z2-top.zbeam)/wxy.ds)
top.prwall = gag_rp   # consistent aperture 
step(n_step)

# --- reset species weights to turn on post accel gap neutralization  
for s in sp.values():
   s.w0 = 1.-neut_f2 

# --- neutralized advance after acceleration column to start of dipole   
n_step = nint((z_adv-top.zbeam)/wxy.ds) + 2   # add two extra steps in case of roundoff accumulation 
top.prwall = r_p_down    # consistent aperture 
step(n_step)



# Make additional history plots for final run if not already called 
#if not(top.it >= diag_hist_step.max()):
#  diag_hist()  # Add full arg list 

# Save restart dump of run.  By default the name of the dump file is
# top.runid (or script name if this is not set) with the step number (iii)
# and ".dump" appended to the name:
#       runidiii.pdb 
# To restart:
#   % python
#     >>> from warp import *
#     >>> restart("runidiii.dump") 
#
#dump() 

# Make plot of final Brho by species 
plt_diag_bro(label = "Final Rigidity by Species") 

# Output data to auxillary file
output_data = false 
output_data_file = "frib-front-xy_data.txt" 

if output_data:
  fout = open(output_data_file,"a")
  #
  #fout.write("Run Number %s :\n"%irun) 
  fout.write(" Solenoid Excitations:\n")
  fout.write("  S41 = %s \n"%s4p1_str)
  fout.write("  S41 = %s \n"%s4p2_str)
  for ii in sp_target:
    s  = sp[ii] 
    js = s.js 
    rmsx    = top.hxrms[0,top.jhist,js]
    drmsxds = top.hxxpbar[0,top.jhist,js]/top.hxrms[0,top.jhist,js]
    # Final rms envelope width and angles 
    fout.write("   "+ii+": \n")
    fout.write("     sqrt(<x^2>) = %s mm, d sqrt(<x^2>)/ds = %s mr \n"%(rmsx/mm,drmsxds/mr))     
  #
  fout.close() 


# Print out timing statistics of run 
printtimers() 



# Multi-Species Envelope Model 
#   Create overlay using axisymmetric, multi-species envelope model


CorrectionMode = 1 #set velocity correction method: 0 - no correction, 1 - dBdz only, 2 - dBdz + d2Edz2

integratewarp = 0 # integrate ode using real-time warp data; 0: no, 1: yes


execfile("frib-front-env.py")

#plotodeterms(0)
#plotwarpterms(0)
#termsodevswarp(0)
#termsdifference(0)

# Make sure that last plot is flushed from buffer
fma() 


