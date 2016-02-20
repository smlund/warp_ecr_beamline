"""
Script to generate field data for the Venus ECR source into WARP by reading in data from an ascii text file provided by the source group.  

Steve Lund 
lund@frib.msu.edu  
"""

from warp import *

setup() 


# Venus ECR solenoid on-axis field 

# Read in data file from magnet designer
#   Data has B_z(r=0,z) from the axial magnet midplane 
#   on a uniform mesh in z with B_z polarity negative.  
#   Use getdatafromtextfile() Warp function.
data = getdatafromtextfile("venus_bz0.txt",nskip=0,dims=[2,None],)

# Extract mesh and scaled B_z (normalized to peak value) from data 
#    d_ prefix for data,  _m suffix for mesh 
[z_m,bz0_m] = data[0:2]

# Plot of raw data 
plg(bz0_m,z_m)
limits('e','e',0.,'e')
ptitles('Venus ECR Raw Field Data: On-Axis B_z_','z [mm]','B_z_(r=0,z) [Tesla]',) 
fma() 

# Scale data for z_m in meters
z_m   = z_m*mm 

# Get max and extraction (z=0) field values 
bz_peak = bz0_m.max() 

iz_peak = 0 
while bz0_m[iz_peak] < bz_peak:
  iz_peak = iz_peak + 1
  
iz_extr = sum(where(z_m < 0.,1,0))-1 
bz_extr = bz0_m[iz_extr] 

# Calculate average bz in ecr region 
bz_avg = sum( bz0_m[iz_peak:iz_extr+1])/len(bz0_m[iz_peak:iz_extr+1])

# Calculate dz for uniform mesh  
dz_m = diff(unique(z_m))
dz = average(dz_m)


# Generate derivative of on-axis field from mesh 
bz0p_m = gradient(bz0_m)/dz 


## -- Extending the fringe using a Gaussian tail

# Input parameters

z_chop_over_sigma = 1.2 # should be >= 1 to ensure B'' > 0 in the tail
                        # the larger the value, the slower the tail decays
b_tail_end = 2e-6       # set the ending field of the tail

# Generating the tail

z_chop = z_m[-1]
bz0_chop = bz0_m[-1]
bz0p_chop = bz0p_m[-1]

z_centre = z_chop - abs(bz0_chop / bz0p_chop)*z_chop_over_sigma
a_Gaussian = bz0_chop*exp(z_chop_over_sigma/2)
sigma = abs(bz0_chop / bz0p_chop)*sqrt(z_chop_over_sigma)

zposition = z_chop + 0.001
fringefield = []
fringez = []
fringeprime = []
tailfield = 1.

while tailfield > 0.000002:
	tailfield = a_Gaussian*exp(-(zposition - z_centre)**2/2/sigma**2)
	taildBdz = tailfield*(z_centre-zposition)/sigma**2
	fringefield.append(tailfield)
	fringez.append(zposition)
	fringeprime.append(taildBdz)
	zposition += 0.001

import numpy as np

z_m   = np.append(z_m, fringez)
bz0_m   = np.append(bz0_m, fringefield)
bz0p_m  = np.append(bz0p_m, fringeprime)

fringe_extension = len(fringez)*0.001

print "field extended by %s meters using Gaussian tail" % fringe_extension


# Calculate mesh bounds and increments 
nz = nint((z_m.max()-z_m.min())/dz)
z_mmin = z_m.min()
z_mmax = z_m.max() 

# Output extraction plane value of B_z 
print("Extraction (z=0 plane) value of B_z = %s Tesla"%bz_extr)

# Plots of processed data 
# --- B_z(r=0,z) vs z 
plg(bz0_m,z_m)
limits('e','e',0,3.6)
ptitles('Venus ECR Solenoid: On-Axis B_z_','z [m]','B_z_(r=0,z) [Tesla]',) 
fma() 

# --- Linear optic B_r/r vs z 
plg(-.5*bz0p_m,z_m)
limits('e','e',-8.,8.)
ptitles('Venus ECR Solenoid: Linear Optics B_r_/r','z [m]','(B_r_/r)(r=0,z) [Tesla/m]',) 
fma() 

# Save linear field data to an external binary array 
fo = PWpickle.PW("ecr_venus.lin.20160218.pkl")
fo.ecr_venus_bz_peak = bz_peak 
fo.ecr_venus_bz_extr = bz_extr
fo.ecr_venus_dz = dz 
fo.ecr_venus_nz = nz 
fo.ecr_venus_z_m    = z_m 
fo.ecr_venus_z_extr = 0.      # z value of beam extraction 
fo.ecr_venus_bz0_m  = bz0_m
fo.ecr_venus_bz0p_m = bz0p_m 
