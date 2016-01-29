
# set z_begin and z_end, default equals to z_launch and z_adv in accordance with setup in frib-front-xy.py

z_begin = z_launch
z_end = z_adv

CorrectionMode = 1 #set velocity correction method: 0 - no correction, 1 - dBdz only, 2 - dBdz + d2Edz2

# set neutralization mode
# 0: same neutralization factor throughout
# 1: different neutralization factor in different regions
# default setup is neut_mode = 1 in accordance with frib-front-lat.py
# neut_f = 0 means no neutralization, neut_f = 1 means full neutralization

neut_mode = 1

# if neut_mode == 0, specify the neutralization factor

neut_f = 0.75 

# if neut_mode == 1, split into different regions and set the respective neut_f
# split the space into regions by specifying the region boundaries

# for example, if neut_f = 0.1 for 0 < z < 0.5 and neut_f = 0.9 for 0.5 < z < 1:
# neut_region_boundaries = [0., 0.5, 1.]
# neut_region_factors = [0.1, 0.9]

neut_region_boundaries = [z_begin, neut_z1, neut_z2, z_end]
neut_region_factors = [0.75, 0., 0.75]


# make sure that len(neut_region_boundaries) = len(neut_region_factors) + 1
if neut_mode == 1:
	if len(neut_region_boundaries) != len(neut_region_factors) + 1:
		raise exception("faulty neutralization region setup")






from scipy.integrate import odeint


# Make time array for solution

stepnum = int(round((z_end - z_begin)/0.005))

sss = linspace(z_begin, z_end, stepnum)

stepsize = (z_end - z_begin)/(stepnum - 1)




# Data needed in Env. Model

speciesq = append(U_charge_states, O_charge_states)*jperev

speciesI = append(U_ibeam, O_ibeam)

specieslist = append(U_species, O_species)



# state vector, for a total of n species
# ordering: same as rest of the code, U then O
# first lot of n elements are KE of the respective species
# second lot of n elements are the sigma
# third lot of elements are the dsigmadz

state_vector = [0]*3*top.ns

dct = {}

deltaz = stepsize/2.

termdata = []


# function that returns the neutralization factor at a given position

def get_neut(zz):
	if neut_mode == 0:
		neut_factor = 1 - neut_f
		
	if neut_mode == 1:
		for i in range(len(neut_region_factors)):
			if zz >= neut_region_boundaries[i]:
				if zz <= neut_region_boundaries[i+1]:
					neut_factor = 1 - neut_region_factors[i]
					break
		if zz < neut_region_boundaries[0]:
			neut_factor = 1 - neut_region_factors[0]
		if zz > neut_region_boundaries[len(neut_region_factors)]:
			neut_factor = 1 - neut_region_factors[len(neut_region_factors)-1]
	
	return neut_factor




## Set up the derivative of the state vector

def f(state_vector, rrr):
	
	efieldz = getappliedfields(0, 0, rrr)[2][0]
	bfieldz = getappliedfields(0, 0, rrr)[5][0]
	
	dEdz = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr - deltaz/2)[2][0])/deltaz
	dBdz = (getappliedfields(0, 0, rrr + deltaz/2)[5][0] - getappliedfields(0, 0, rrr - deltaz/2)[5][0])/deltaz
	
	d2Edz2 = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr)[2][0]*2 + getappliedfields(0, 0, rrr - deltaz/2)[2][0])/(deltaz**2)*4
	#dVdz = 0
	
	derivs = []	
	
	speciesbeta = []
	
	for i in range(top.ns):
		speciesbeta.append(sqrt((2*state_vector[i]*jperev)/(specieslist[i].mass*clight**2)))
		
	
	## build first lot in deriv output (i.e. dKEdz)
	
	for j in range(top.ns):
		
		if CorrectionMode == 0:
			derivs.append(speciesq[j]*(efieldz)/jperev )
			
		if CorrectionMode == 1:
			derivs.append(speciesq[j]/jperev*(efieldz + (speciesbeta[j]*clight*hl_pth[0,j]/2 - speciesq[j]*state_vector[j+top.ns]**2*bfieldz/4/specieslist[j].mass)*dBdz) )
		
		if CorrectionMode == 2:
			derivs.append(speciesq[j]/jperev*(efieldz - state_vector[j+top.ns]**2/4*d2Edz2 + (speciesbeta[j]*clight*hl_pth[0,j]/2 - speciesq[j]*state_vector[j+top.ns]**2*bfieldz/4/specieslist[j].mass)*dBdz) )

	for i in range(top.ns):
		dct["sigma%s" %(i)] = state_vector[i+top.ns]
		dct["dsigmadz%s" %(i)] = state_vector[i+2*top.ns]
		derivs.append(state_vector[i+2*top.ns]) #build second lot in deriv output
	
	
	## Set the neutralization factor
	
	neut_ode = get_neut(rrr)
	
	## build third lot in deriv output (i.e. sigma'' )
	
	for j in range(top.ns):
		
		scterm = 0
		
		for s in range(top.ns):
			QQQ = (speciesq[j]*speciesI[s])/(2*pi*eps0*specieslist[j].mass*speciesbeta[j]**2*speciesbeta[s]*clight**3)
			scterm += QQQ*neut_ode*state_vector[j+top.ns]/(state_vector[j+top.ns]**2 + state_vector[s+top.ns]**2)
		
		term1 = (speciesq[j]*-efieldz)/(2*state_vector[j]*jperev) * state_vector[j+2*top.ns]
		
		term2 = (speciesq[j]*-dEdz)/(4*state_vector[j]*jperev) * state_vector[j+top.ns]
		
		term3 = ((speciesq[j]*bfieldz)/(2*specieslist[j].mass*speciesbeta[j]*clight))**2*state_vector[j+top.ns]

		emitterm = ((top.hepsr[0,0,j]/2)**2 + (hl_pthn[0,j] /speciesbeta[j])**2) / state_vector[j+top.ns]**3
		
		# top.hepsr equals two times the rms-r-emittance
		
		d2sigmadz2 = term1 + term2 - term3 + scterm + emitterm
		
		derivs.append(d2sigmadz2)
		
		termdata.append([j, rrr, term1, term2, -term3, scterm, emitterm, state_vector[j]])
		
	return derivs




# Set up initial states ( dimensions = 3 * number of species, 1st lot KE, 2nd lot sigma-r, 3rd lot sigma-r-prime)

initialstates = []

## Kinetic energy

for i in range(U_ns):
	initialstates.append(U_ekin[i])
	
for i in range(O_ns):
	initialstates.append(O_ekin[i])

## Initial rms-radius

for i in range(top.ns):
	initialstates.append(sqrt((specieslist[i].a0/2)**2 + (specieslist[i].b0/2)**2))

## Initial envelope angle

for i in range(top.ns):
	initialstates.append((specieslist[i].a0*specieslist[i].ap0 + specieslist[i].b0*specieslist[i].bp0) / initialstates[i+top.ns])




## Numerical Solution of the Env. Model

psoln = odeint (f, initialstates, sss, hmax = stepsize, mxstep=5000)


# Load the diagnostic associated with the Env. Model

execfile("frib-front-env-diag.py")



