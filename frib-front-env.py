
from scipy.integrate import odeint


# Make time array for solution

stepnum = int(round((z_adv - z_launch)/0.005))

sss = linspace(z_launch, z_adv, stepnum)

stepsize = (z_adv - z_launch)/(stepnum - 1)




# Data needed in ODE

speciesq = U_charge_states
speciesq.extend(O_charge_states)
speciesq[:] = [x*jperev for x in speciesq]


speciesI = append(U_ibeam, O_ibeam)

#speciesI = U_ibeam
#speciesI.extend(O_ibeam)

speciesCombined = U_species
speciesCombined.extend(O_species)

# state vector, for a total of n species
# ordering: same as rest of the code, U then O
# first lot of n elements are KE of the respective species
# second lot of n elements are the sigma
# third lot of elements are the dsigmadz

xyz = [0]*3*top.ns

dct = {}

deltaz = stepsize/2.

termdata = []

## Set up the derivative of the state vector

def f(xyz, rrr):
	
	efieldz = getappliedfields(0, 0, rrr)[2][0]
	bfieldz = getappliedfields(0, 0, rrr)[5][0]
	
	dEdz = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr - deltaz/2)[2][0])/deltaz
	dBdz = (getappliedfields(0, 0, rrr + deltaz/2)[5][0] - getappliedfields(0, 0, rrr - deltaz/2)[5][0])/deltaz
	
	d2Edz2 = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr)[2][0]*2 + getappliedfields(0, 0, rrr - deltaz/2)[2][0])/(deltaz**2)*4
	#dVdz = 0
	
	derivs = []	
	
	speciesbeta = []
	
	for i in range(top.ns):
		speciesbeta.append(sqrt((2*xyz[i]*jperev)/(speciesCombined[i].mass*clight**2)))
		
	
	## build first lot in deriv output (i.e. dKEdz)
	
	for j in range(top.ns):
		
		if CorrectionMode == 0:
			derivs.append(speciesq[j]*(efieldz)/jperev )
			
		if CorrectionMode == 1:
			derivs.append(speciesq[j]/jperev*(efieldz + (speciesbeta[j]*clight*hl_pth[0,j]/2 - speciesq[j]*xyz[j+top.ns]**2*bfieldz/4/speciesCombined[j].mass)*dBdz) )
		
		if CorrectionMode == 2:
			derivs.append(speciesq[j]/jperev*(efieldz - xyz[j+top.ns]**2/4*d2Edz2 + (speciesbeta[j]*clight*hl_pth[0,j]/2 - speciesq[j]*xyz[j+top.ns]**2*bfieldz/4/speciesCombined[j].mass)*dBdz) )

	for i in range(top.ns):
		dct["sigma%s" %(i)] = xyz[i+top.ns]
		dct["dsigmadz%s" %(i)] = xyz[i+2*top.ns]
		derivs.append(xyz[i+2*top.ns]) #build second lot in deriv output
	
	
	## Turn off neutralization in the ES Gap
	
	if neut_z1 < rrr < neut_z2:
		neut_ode = 1 - 0
	else:
		neut_ode = 1 - neut_f1
	
	
	## build third lot in deriv output (i.e. sigma'' )
	
	for j in range(top.ns):
		
		scterm = 0
		
		for s in range(top.ns):
			QQQ = (speciesq[j]*speciesI[s])/(2*pi*eps0*speciesCombined[j].mass*speciesbeta[j]**2*speciesbeta[s]*clight**3)
			scterm += QQQ*neut_ode*xyz[j+top.ns]/(xyz[j+top.ns]**2 + xyz[s+top.ns]**2)
		
		term1 = (speciesq[j]*-efieldz)/(2*xyz[j]*jperev) * xyz[j+2*top.ns]
		
		term2 = (speciesq[j]*-dEdz)/(4*xyz[j]*jperev) * xyz[j+top.ns]
		
		term3 = ((speciesq[j]*bfieldz)/(2*speciesCombined[j].mass*speciesbeta[j]*clight))**2*xyz[j+top.ns]

		emitterm = ((top.hepsr[0,0,j]/2)**2 + (hl_pthn[0,j] /speciesbeta[j])**2) / xyz[j+top.ns]**3
		
		# top.hepsr equals two times the rms-r-emittance
		
		d2sigmadz2 = term1 + term2 - term3 + scterm + emitterm
		
		derivs.append(d2sigmadz2)
		
		termdata.append([j, rrr, term1, term2, -term3, scterm, emitterm, xyz[j]])
		
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
	initialstates.append(sqrt((speciesCombined[i].a0/2)**2 + (speciesCombined[i].b0/2)**2))

## Initial envelope angle

for i in range(top.ns):
	initialstates.append((speciesCombined[i].a0*speciesCombined[i].ap0 + speciesCombined[i].b0*speciesCombined[i].bp0) / initialstates[i+top.ns])




## Numerical Solution of the ODE

psoln = odeint (f, initialstates, sss, hmax = stepsize, mxstep=5000)










## Compute the terms in the ODE using real-time WARP data

termdatawarp = []

for iii in range(0, top.jhist+1):
	
	xyz = 0
	
	zzz = top.hzbeam[iii]
	deltaz = diff(top.hzbeam)[0]/2
	
	efieldz = getappliedfields(0, 0, zzz)[2][0]
	bfieldz = getappliedfields(0, 0, zzz)[5][0]
	dEdz = (getappliedfields(0, 0, zzz + deltaz/2)[2][0] - getappliedfields(0, 0, zzz - deltaz/2)[2][0])/deltaz
	
	kelist = array([hl_ekin[iii,j] for j in range(top.ns)])
	speciesbeta = array([sqrt((2*kelist[j]*jperev)/(speciesCombined[j].mass*clight**2)) for j in range(top.ns)])
	anglelist = array([top.hxxpbar[0,iii,j]/(top.hxrms[0,iii,j])*sqrt(2) for j in range(top.ns)])
	envelopelist = array([top.hxrms[0,iii,j]*sqrt(2) for j in range(top.ns)])
	pthetalist = array([hl_pthn[iii,j] for j in range(top.ns)])
	
	xyz = concatenate((kelist, envelopelist, anglelist))
	
	if neut_z1 < zzz < neut_z2:
		neut_ode = 1 - 0
	else:
		neut_ode = 1 - neut_f1
	
	for j in range(top.ns):
				
		scterm = 0
		
		for s in range(top.ns):
			QQQ = (speciesq[j]*speciesI[s])/(2*pi*eps0*speciesCombined[j].mass*speciesbeta[j]**2*speciesbeta[s]*clight**3)
			scterm += QQQ*neut_ode*xyz[j+top.ns]/(xyz[j+top.ns]**2 + xyz[s+top.ns]**2)
		
		term1 = (speciesq[j]*-efieldz)/(2*xyz[j]*jperev) * xyz[j+2*top.ns]
		
		term2 = (speciesq[j]*-dEdz)/(4*xyz[j]*jperev) * xyz[j+top.ns]
		
		term3 = ((speciesq[j]*bfieldz)/(2*speciesCombined[j].mass*speciesbeta[j]*clight))**2*xyz[j+top.ns]

		emitterm = ((top.hepsr[0,iii,j]/2)**2 + (pthetalist[j] /speciesbeta[j])**2) / xyz[j+top.ns]**3
		
		#d2sigmadz2 = term1 + term2 - term3 + scterm + emitterm
		
		#derivs.append(d2sigmadz2)
		
		termdatawarp.append([j, zzz, term1, term2, -term3, scterm, emitterm, xyz[j]])









from scipy import interpolate

xyzxyz = [0]*3*top.ns

dct = {}

deltaz = stepsize/2.



## Set up function used to integrate the ODE using real-time WARP data

zlist = []
		
zlist = array([top.hzbeam[kkk] for kkk in range(0, top.jhist+1)])

def fwarp(xyzxyz, rrr):
	
	efieldz = getappliedfields(0, 0, rrr)[2][0]
	bfieldz = getappliedfields(0, 0, rrr)[5][0]
	
	dEdz = (getappliedfields(0, 0, rrr + deltaz/2)[2][0] - getappliedfields(0, 0, rrr - deltaz/2)[2][0])/deltaz
	#dVdz = 0
	
	derivs = []
	
	for i in range(U_ns):
		derivs.append(U_species[i].charge/jperev*efieldz) #build first lot in deriv output
	
	for i in range(O_ns):
		derivs.append(O_species[i].charge/jperev*efieldz) #build first lot in deriv output
	
	for i in range(top.ns):
		derivs.append(xyzxyz[i+2*top.ns]) #build second lot in deriv output
	
	speciesbeta = []
	
	for i in range(top.ns):
		speciesbeta.append(sqrt((2*xyzxyz[i]*jperev)/(speciesCombined[i].mass*clight**2)))

	if neut_z1 < rrr < neut_z2:
		neut_ode = 1 - 0
	else:
		neut_ode = 1 - neut_f1

	for j in range(top.ns):
		
		scterm = 0
		
		
		
		emittancelist = []
		
		#for kkk in range(len(top.hepsr)):
			#emittancelist = emittancelist.append(top.hepsr[0,kkk,j])
		
		emittancelist = array([top.hepsr[0,kkk,j]/2 for kkk in range(0, top.jhist+1)])
		
		pthetaLIST = []
		
		pthetaLIST = array([hl_pthn[kkk,j] for kkk in range(0, top.jhist+1)])
		
		kineticenergylist = []
		
		kineticenergylist = array([hl_ekin[kkk,j] for kkk in range(0, top.jhist+1)])
		

		emitinter = interpolate.interp1d(zlist, emittancelist, kind='slinear')
		
		pthetainter = interpolate.interp1d(zlist, pthetaLIST, kind='slinear')
		
		keinter = interpolate.interp1d(zlist, kineticenergylist, kind='slinear')
		
		for s in range(top.ns):
			QQQ = (speciesq[j]*speciesI[s])/(2*pi*eps0*speciesCombined[j].mass*speciesbeta[j]**2*speciesbeta[s]*clight**3)
			scterm += QQQ*neut_ode*xyzxyz[j+top.ns]/(xyzxyz[j+top.ns]**2 + xyzxyz[s+top.ns]**2)
		
		term1 = (speciesq[j]*-efieldz)/(2*xyzxyz[j]*jperev) * xyzxyz[j+2*top.ns]
		
		term2 = (speciesq[j]*-dEdz)/(4*xyzxyz[j]*jperev) * xyzxyz[j+top.ns]
		
		term3 = ((speciesq[j]*bfieldz)/(2*speciesCombined[j].mass*speciesbeta[j]*clight))**2*xyzxyz[j+top.ns]

		emitterm = ((emitinter(rrr))**2 + (pthetainter(rrr) /speciesbeta[j])**2) / xyzxyz[j+top.ns]**3
		
		d2sigmadz2 = term1 + term2 - term3 + scterm + emitterm
		
		derivs.append(d2sigmadz2)
		
	return derivs
		
		
if integratewarp == 1:
	psolnwarp = odeint (fwarp, initialstates, sss, hmax = stepsize, mxstep=5000)













termdata2 = []
termdatawarp2=[]

for i in range(0, len(termdata), 20):
	termdata2.append(termdata[i])

termdata2 = array(termdata2)

termdatawarp2 = array([termdatawarp[i] for i in  range(0, len(termdatawarp), 20)])




# --- Target Species, x-plane 
hpenvx(titles=false,yscale=1./(2.*mm),lhzbeam=true)
lab = ""    
for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      lab+= ii + "("+co+"), "
      hpenvx(js=js,color=co,titles=false,yscale=1./(2.*mm),lhzbeam=true)    
ptitles("ODE (dashed), Warp (solid)","z [m]","RMS x-Envelope [mm]", )

plg(psoln[:,20]*1000/sqrt(2),sss, color="magenta", type="dash")
plg(psoln[:,21]*1000/sqrt(2),sss, color="green", type="dash")
#limits(plotxmin,plotxmax, plotymin, plotymax)


combinedenvelop = [0]*stepnum

totalweight=0.
for j in speciesCombined:
	totalweight += j.sw

for i in range(stepnum):
	
	SumofR2overN = 0.
	
	for j in range(top.ns):
		#combinedenvelop[i] += psoln[i, j+top.ns]*1000/sqrt(2)*speciesI[j]/sum(speciesI)
		SumofR2overN += psoln[i, j+top.ns]**2*(speciesCombined[j].sw/totalweight)
	
	combinedenvelop[i] = SumofR2overN**0.5*1000/sqrt(2)

plg(combinedenvelop,sss, color="black", type="dash")

fma()	







SpeciesDictionary = {'U33':0, 'U34':1}

for kkk in range(2,10):
	abc = 'U%d' % (kkk+23)
	DummyDict = {abc:kkk}
	SpeciesDictionary.update(DummyDict)

for kkk in range(1,7):
	abc = 'U%d' % (kkk+34)
	DummyDict = {abc:kkk+9}
	SpeciesDictionary.update(DummyDict)
	
for kkk in range(1,5):
	abc = 'O%d' % kkk
	DummyDict = {abc:kkk+15}
	SpeciesDictionary.update(DummyDict)


def SpeciesEnvelope(iii):
    s = sp[iii]
    js = s.js
    co = s.color
    hpenvx(js=js,color=co,titles=false,yscale=1./(2.*mm),lhzbeam=true)

    
    
    RefNumber = SpeciesDictionary[iii]
    
    plg(psoln[:,RefNumber+20]*1000/sqrt(2),sss, color=co, type="dash")
    ptitles("ODE (dashed), Warp (solid), %s" % iii,"z [m]","RMS x-Envelope [mm]", )
    #limits(plotxmin,plotxmax, plotymin, plotymax)
    
    fma()


s = sp['O1']
js = s.js
hpenvx(js=js,color='black',titles=false,yscale=1./(2.*mm),lhzbeam=true)
RefNumber = SpeciesDictionary['O1']
plg(psoln[:,RefNumber+20]*1000/sqrt(2),sss, color='black', type="dash")

s = sp['O4']
js = s.js
hpenvx(js=js,color='red',titles=false,yscale=1./(2.*mm),lhzbeam=true)
RefNumber = SpeciesDictionary['O4']
plg(psoln[:,RefNumber+20]*1000/sqrt(2),sss, color='red', type="dash")

ptitles("ODE (dashed), Warp (solid), O1 black, O4 red","z [m]","RMS x-Envelope [mm]", )
#limits(plotxmin,plotxmax, plotymin, plotymax)

fma()


for jjj in sp:
	SpeciesEnvelope(jjj)










for ii in sp_target:
      s = sp[ii]
      js = s.js
      co = s.color
      A  = s.mass/amu
      #hpekin(js=js,color=co,titles=false,yscale=1./A,lhzbeam=true)
      plg(hl_ekin[0:top.jhist+1,js]/(A*kV),top.hzbeam[0:top.jhist+1],color=co)        
ptitles("Kinetic Energy, ODE (dashed), Warp (solid)","z [m]","KeV/u", )

plg(psoln[:,0]/1000/238.02891,sss, color="magenta", type="dash")
plg(psoln[:,1]/1000/238.02891,sss, color="green", type="dash")

fma()



# --- Target Species, x-plane 
#lab = ""    
#for ii in sp_target:
      #s = sp[ii]
      #js = s.js
      #co = s.color
      #lab+= ii + "("+co+"), "
      #plg(top.hepsr[0, 0:top.jhist+1,js]*CorrectionFactor,top.hzbeam[0:top.jhist+1],color=co)
      #plg(sqrt(4*top.hepsx[0, 0:top.jhist+1,js]**2 - hl_lz_bar2[0:top.jhist+1,js]**2),top.hzbeam[0:top.jhist+1],color=co, type="dash")    
#ptitles("History: "+lab+"RMS x-Envelope","z [m]","RMS Width [mm]", )

#fma()




## Plot the terms in the Envelope Equation as calculated by the ODE

def plotodeterms(jj):
	
	termdata4 = []
	termdatawarp4=[]

	
	termdata4 = array([termdata[i] for i in range(jj, len(termdata), top.ns)])

	termdatawarp4 = array([termdatawarp[i] for i in  range(jj, len(termdatawarp), top.ns)])
	
	plg(termdata4[:,2], termdata4[:,1], color="black", type="dash")
	plg(termdata4[:,3], termdata4[:,1], color="green")
	plg(termdata4[:,4], termdata4[:,1], color="black")
	plg(termdata4[:,5], termdata4[:,1], color="red")
	plg(termdata4[:,6], termdata4[:,1], color="blue")
	
	ptitles("ODE Terms, black-dash:E green:dEdz black:B red:SC blue:Emit","z [m]","sigma'' [rad/m]",)
	
	fma()


plotodeterms(0)




## Plot the terms in the Envelope Equation as calculated by the WARP

def plotwarpterms(jj):
	
	termdata4 = []
	termdatawarp4=[]

	
	termdata4 = array([termdata[i] for i in range(jj, len(termdata), top.ns)])

	termdatawarp4 = array([termdatawarp[i] for i in  range(jj, len(termdatawarp), top.ns)])
	
	plg(termdatawarp4[:,2], termdatawarp4[:,1], color="black", type="dash")
	plg(termdatawarp4[:,3], termdatawarp4[:,1], color="green")
	plg(termdatawarp4[:,4], termdatawarp4[:,1], color="black")
	plg(termdatawarp4[:,5], termdatawarp4[:,1], color="red")
	plg(termdatawarp4[:,6], termdatawarp4[:,1], color="blue")
	
	ptitles("Terms with WARP input, black-dash:E green:dEdz black:B red:SC blue:Emit","z [m]","sigma'' [rad/m]",)
	
	fma()


plotwarpterms(0)





## Compare the value of each term in the Envelope Equation from WARP and ODE

def termsodevswarp(jj):
	
	termdata3 = []
	termdatawarp3=[]

	for i in range(jj, len(termdata), top.ns):
		termdata3.append(termdata[i])

	termdata3 = array(termdata3)

	termdatawarp3 = array([termdatawarp[i] for i in  range(jj, len(termdatawarp), top.ns)])


	plg(termdata3[:,2],termdata3[:,1])
	plg(termdatawarp3[:,2],termdatawarp3[:,1], type="dash")
	ptitles("Term 1 (E-field), ODEinput (solid), WARPinput (dashed)","z [m]","sigma'' [rad/m]",)
	fma()
	
	plg(termdata3[:,3],termdata3[:,1])
	plg(termdatawarp3[:,3],termdatawarp3[:,1], type="dash")
	ptitles("Term 2 (dEdz), ODEinput (solid), WARPinput (dashed)","z [m]","sigma'' [rad/m]",)
	fma()
	
	plg(termdata3[:,4],termdata3[:,1])
	plg(termdatawarp3[:,4],termdatawarp3[:,1], type="dash")
	ptitles("Term 3 (B-field), ODEinput (solid), WARPinput (dashed)","z [m]","sigma'' [rad/m]",)
	fma()
	
	plg(termdata3[:,5],termdata3[:,1])
	plg(termdatawarp3[:,5],termdatawarp3[:,1], type="dash")
	ptitles("Space Charge Term, ODEinput (solid), WARPinput (dashed)","z [m]","sigma'' [rad/m]",)
	fma()
	
	plg(termdata3[:,6],termdata3[:,1])
	plg(termdatawarp3[:,6],termdatawarp3[:,1], type="dash")
	ptitles("Emittance Term, ODEinput (solid), WARPinput (dashed)","z [m]","sigma'' [rad/m]",)
	fma()



termsodevswarp(0)









## Plot the difference between the ODE terms and the WARP terms

def termsdifference(jj):
	
	termdata4 = []
	termdatawarp4=[]

	
	termdata4 = array([termdata[i] for i in range(jj, len(termdata), top.ns)])

	termdatawarp4 = array([termdatawarp[i] for i in  range(jj, len(termdatawarp), top.ns)])
	
	warpterm1inter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,2], bounds_error = False, fill_value = 0., kind='linear')
	warpterm2inter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,3], bounds_error = False, fill_value = 0., kind='linear')
	warpterm3inter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,4], bounds_error = False, fill_value = 0., kind='linear')
	warpscterminter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,5], bounds_error = False, fill_value = 0., kind='linear')
	warpemitterminter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,6], bounds_error = False, fill_value = 0., kind='linear')
	
	warpterm1 = array([ warpterm1inter(kkk) for kkk in termdata4[:,1] ])
	warpterm2 = array([ warpterm2inter(kkk) for kkk in termdata4[:,1] ])
	warpterm3 = array([ warpterm3inter(kkk) for kkk in termdata4[:,1] ])
	warpscterm = array([ warpscterminter(kkk) for kkk in termdata4[:,1] ])
	warpemitterm = array([ warpemitterminter(kkk) for kkk in termdata4[:,1] ])
	
	if len(termdata4[:,2]) != len(warpterm1):
		print "termdata4[:,2] not same length as warpterm1"
	
	term1diff = termdata4[:,2] - warpterm1
	term2diff = termdata4[:,3] - warpterm2
	term3diff = termdata4[:,4] - warpterm3
	sctermdiff = termdata4[:,5] - warpscterm
	emittermdiff = termdata4[:,6] - warpemitterm
	
	plg(term1diff, termdata4[:,1], color="black", type="dash")
	plg(term2diff, termdata4[:,1], color="green")
	plg(term3diff, termdata4[:,1], color="black")
	plg(sctermdiff, termdata4[:,1], color="red")
	plg(emittermdiff, termdata4[:,1], color="blue")
	
	ptitles("ODEinput minus WARPinput for different terms","z [m]","sigma'' [rad/m]",)
	
	fma()

termsdifference(0)




## Plot the KE of a species from WARP and ODE

def kedifference(jj):
	
	termdata4 = []
	termdatawarp4=[]

	
	termdata4 = array([termdata[i] for i in range(jj, len(termdata), top.ns)])

	termdatawarp4 = array([termdatawarp[i] for i in  range(jj, len(termdatawarp), top.ns)])
	
	warpkeinter = interpolate.interp1d(termdatawarp4[:,1], termdatawarp4[:,7], bounds_error = False, fill_value = 0., kind='linear')
	
	warpke = array([ warpkeinter(kkk) for kkk in termdata4[:,1] ])
	
	if len(termdata4[:,2]) != len(warpke):
		print "termdata4[:,2] not same length as warpterm1"
	
	kediff = termdata4[:,7] - warpke
	
	kediffcut = array(kediff[4:-4])
	warpkecut = array(warpke[4:-4])
	
	kepercentagediff = kediffcut/warpkecut*100
	
	if jj == 0:
		plg(kepercentagediff, termdata4[:,1][4:-4], color="magenta")
	elif jj == 1:
		plg(kepercentagediff, termdata4[:,1][4:-4], color="green")
	else:
		plg(kepercentagediff, termdata4[:,1][4:-4], color="black")
	

kedifference(0)
kedifference(1)

ptitles("KE Difference between ODE and WARP","z [m]","% difference",)
	
fma()




if integratewarp == 1:
	
	
	
	## ODE with WARP input versus WARP Simulations

	# --- Target Species, x-plane 
	hpenvx(titles=false,yscale=1./(2.*mm),lhzbeam=true)
	lab = ""    
	for ii in sp_target:
	      s = sp[ii]
	      js = s.js
	      co = s.color
	      lab+= ii + "("+co+"), "
	      hpenvx(js=js,color=co,titles=false,yscale=1./(2.*mm),lhzbeam=true)    
	ptitles("ODE with WARP input (dashed), WARP Simulations (solid)","z [m]","RMS x-Envelope [mm]", )
	
	plg(psolnwarp[:,20]*1000/sqrt(2),sss, color="magenta", type="dash")
	plg(psolnwarp[:,21]*1000/sqrt(2),sss, color="green", type="dash")
	
	
	combinedenvelop = [0]*stepnum
	
	TotalWeight=0.
	for j in speciesCombined:
		TotalWeight += j.sw
	
	for i in range(stepnum):
		
		SumofR2overN = 0.
		
		for j in range(top.ns):
			#combinedenvelop[i] += psoln[i, j+top.ns]*1000/sqrt(2)*speciesI[j]/sum(speciesI)
			SumofR2overN += psolnwarp[i, j+top.ns]**2*(speciesCombined[j].sw/TotalWeight)
		
		combinedenvelop[i] = SumofR2overN**0.5*1000/sqrt(2)
	
	plg(combinedenvelop,sss, color="black", type="dash")	
	
	fma()
	
	
	
	## Compare ODE with ODE using WARP input
	
	# --- Target Species, x-plane 
	
	ptitles("ODE (solid), ODE with WARP input (dashed)","z [m]","RMS x-Envelope [mm]", )
	
	plg(psoln[:,20]*1000/sqrt(2),sss, color="magenta")
	plg(psoln[:,21]*1000/sqrt(2),sss, color="green")
	plg(psolnwarp[:,20]*1000/sqrt(2),sss, color="magenta", type="dash")
	plg(psolnwarp[:,21]*1000/sqrt(2),sss, color="green", type="dash")
	
	combinedenvelop = [0]*stepnum
	
	TotalWeight=0.
	for j in speciesCombined:
		TotalWeight += j.sw
	
	for i in range(stepnum):
		
		SumofR2overN = 0.
		
		for j in range(top.ns):
			#combinedenvelop[i] += psoln[i, j+top.ns]*1000/sqrt(2)*speciesI[j]/sum(speciesI)
			SumofR2overN += psolnwarp[i, j+top.ns]**2*(speciesCombined[j].sw/TotalWeight)
		
		combinedenvelop[i] = SumofR2overN**0.5*1000/sqrt(2)
	
	plg(combinedenvelop,sss, color="black", type="dash")	
	
	
	
	combinedenvelop = [0]*stepnum
	
	TotalWeight=0.
	for j in speciesCombined:
		TotalWeight += j.sw
	
	for i in range(stepnum):
		
		SumofR2overN = 0.
		
		for j in range(top.ns):
			#combinedenvelop[i] += psoln[i, j+top.ns]*1000/sqrt(2)*speciesI[j]/sum(speciesI)
			SumofR2overN += psoln[i, j+top.ns]**2*(speciesCombined[j].sw/TotalWeight)
		
		combinedenvelop[i] = SumofR2overN**0.5*1000/sqrt(2)
	
	plg(combinedenvelop,sss, color="black")	
	
	fma()


ptitles("Fields: %s, f_neut: %s, Correction: %s, CentralDifference" %(abc, neut_f1, CorrectionMode))

fma()
