# Diagnostics for lattice field 
# Steve Lund 
# FRIB Front End 
# 2014 

# Grated Accel Gap 
if gag_typ == "nl":
  plotegrd(ie=gag[0],component='z',ix=0,iy=0,titles=false)
  ptitles("Grated Accel Gap: on-axis Ez vs z","z [m]","Ez [V/m]",)
  fma()

  plotegrd(ie=gag[0],component='x',ix=gag_nr/2,iy=0,titles=false)
  ptitles("Grated Accel Gap: Er vs z at r = %s mm"%(gag_r_m[gag_nr/2]/mm),"z [m]","Er [V/m]",)
  fma()


# On Axis Bz and Ez  
nz = 1000  
z_mmin = 65.5 
z_mmax = 69.
dz = (z_mmax-z_mmin)/nz   
x_mesh = zeros(nz+1)
y_mesh = zeros(nz+1)
z_mesh = z_mmin + arange(nz+1)*dz  

(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x=x_mesh,y=y_mesh,z=z_mesh)

plg(bz_mesh,z_mesh)
ptitles("Axial Magnetic Field at r=0","z [m]","Bz [Tesla]",)
fma() 

plg(ez_mesh,z_mesh) 
ptitles("Axial Electric Field at r=0","z [m]","Ez [V/m]",) 
fma() 


# Br and Er at fraction of aperture within applied field grid    

r_f = r_p/2.  
x_mesh = r_f*ones(nz+1) 

(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x=x_mesh,y=y_mesh,z=z_mesh)

plg(ex_mesh,z_mesh)
ptitles("Radial Electric Field at r=%s mm"%(r_f/mm),"z [m]","Er [V/m]",)
fma()

# Symmetry check of Br and Er that signs are correct along principal axes  

x_mesh = r_f*ones(nz+1) 

(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x= x_mesh,y=y_mesh,z=z_mesh)
plg(bx_mesh,z_mesh,color="black")
(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x=-x_mesh,y=y_mesh,z=z_mesh)
plg(bx_mesh,z_mesh,color="red"  )
ptitles("Bx Magnetic Field at y=0, x = +-%s mm B/R"%(r_f/mm),"z [m]","Bx [Tesla]",)
fma()

(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x= x_mesh,y=y_mesh,z=z_mesh)
plg(ex_mesh,z_mesh,color="black")
(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x=-x_mesh,y=y_mesh,z=z_mesh)
plg(ex_mesh,z_mesh,color="red"  )
ptitles("Ex Electric Field at y=0, x = +-%s mm B/R"%(r_f/mm),"z [m]","Ex [V/m]",)
fma()

x_mesh =    zeros(nz+1)
y_mesh = r_f*ones(nz+1) 

(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x=x_mesh,y= y_mesh,z=z_mesh)
plg(by_mesh,z_mesh,color="black")
(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x=x_mesh,y=-y_mesh,z=z_mesh)
plg(by_mesh,z_mesh,color="red"  )
ptitles("By Magnetic Field at x=0, y = +-%s mm B/R"%(r_f/mm),"z [m]","By [Tesla]",)
fma()

(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x=x_mesh,y= y_mesh,z=z_mesh)
plg(ey_mesh,z_mesh,color="black")
(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x=x_mesh,y=-y_mesh,z=z_mesh)
plg(ey_mesh,z_mesh,color="red"  )
ptitles("Ey Electric Field at x=0, y = +-%s mm B/R"%(r_f/mm),"z [m]","Ey [V/m]",)
fma()


# On axis Bz element by element 

def diag_field_ele_Bz(type):
  dict_id  = {'s4p1':s4p1[0],  's4p2':s4p2[0], 'ecr':ecr[0]}
  dict_typ = {'s4p1':s4p1_typ, 's4p2':s4p2_typ,'ecr':ecr_typ}
  id  = dict_id[type]
  typ = dict_typ[type]
  if typ == "lin":
    plotmmlt(im=id,titles=false)
  elif typ == "nl":
    plotbgrd(ib=id,component='z',ix=0,iy=0,titles=false)
  else:
    return

diag_field_ele_Bz('ecr')
diag_field_ele_Bz('s4p1')
diag_field_ele_Bz('s4p2')
ptitles("On Axis Bz of Elements","z [m]","Bz [Tesla]",)
fma() 


# Plot of Bz and Ez on Axis on a similar scale to help understand field structure

nz = 1000  
z_mesh = linspace(ecr_z_extr,z_adv,nz+1)

x_mesh = zeros(nz+1)
y_mesh = zeros(nz+1) 

(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x=x_mesh,y= y_mesh,z=z_mesh)

ez_max = maxnd( abs(ez_mesh) ) 
bz_max = maxnd( abs(bz_mesh) ) 

plg([1.,0.],[z_launch,z_launch], type = "dash")
plg(bz_mesh/bz_max, z_mesh)
plg(ez_mesh/ez_max, z_mesh, color = "red")  
ptitles("Scaled Axial Field of Elements (B: Bz) and (R: Ez)","z [m]","Bz and Ez",) 
fma() 


z_mesh = linspace(ecr_z_extr - 1.1,z_adv,nz+1)
(ex_mesh,ey_mesh,ez_mesh,bx_mesh,by_mesh,bz_mesh) = getappliedfields(x=x_mesh,y= y_mesh,z=z_mesh)

ez_max = maxnd( abs(ez_mesh) ) 
bz_max = maxnd( abs(bz_mesh) ) 

plg([1.,0.],[z_launch,z_launch], type = "dash")
plg(bz_mesh/bz_max, z_mesh)
plg(ez_mesh/ez_max, z_mesh, color = "red")  
ptitles("Scaled Axial Field of Elements (B: Bz) and (R: Ez)","z [m]","Bz and Ez",) 
fma() 


diag_field_ele_Bz('ecr')
diag_field_ele_Bz('s4p1')
diag_field_ele_Bz('s4p2')
plg(2.*ez_mesh/ez_max, z_mesh, color = "red")  
plg([3.5,0.],[z_launch,z_launch], type = "dash")
limits('e','e',0.,3.5)
ptitles("On Axis Bz (B) of Elements and Scaled Ez (R)","z [m]","Bz [Tesla]",)
fma() 
  
