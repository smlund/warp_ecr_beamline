warp_ion_frontend
Git repository of Warp simulation scripts/tools to model front-ends for 
heavy ion accelerators. Setup by default to be directly applicable to 
the FRIB front-end with ECR type sources.  However, tools are more broadly 
applicable/adaptable to simulation of multi-species heavy-ion beams.  

Professor Steven M. Lund
Physics and Astronomy Department
Facility for Rare Isotope Beams
Michigan State University
lund@frib.msu.edu
517-908-7291

Chun Yan (Jonathan) Wong 
Physics and Astronomy Department
Michigan State University
wong@nscl.msu.edu
517-908-7465

Kei Fukushima
Facility for Rare Isotope Beams
Michigan State University
fukushim@frib.msu.edu
517-908-7254


To initialize the repository, 

   % git clone https://github.com/smlund/warp_ion_frontend

This will create a directory, ./warp_ion_frontend where command was 
run with the archive files.   

To get the latest version descend into directory warp_ion_beamline and run:

  % git pull 

When modifying the repository (for those with edit privilege, please contact 
me if you want to contribute and I will add you) 

  ... edit files etc then checkin (use readme.txt as example here) using: 
  % git add readme.txt 
  % git commit -m "SML: updated readme.txt file" 
  % git push 

To remove a file from git control (to not include in future pulls), use 

  % git rm file
  % git commit -m "SML: removing file from repo"
  % git push 

The local copy of the file removed can be retained on disk by subsituting


  % git rm --cached file 

in the above.  Another way to remove a file from the master node repo is to
go to the github web side, click on a file, click the delete button and then
confirm at the bottom of the page.   

Field files needed for the lattices described by the simulation are stored
in the Dropbox file sharing system. Code users do NOT need dropbox access
to retrieve the needed field description files for lattice elements.
A unix shell script "frib-lat-fields-dropbox-fetch" is provided that
employs wget with dropbox links to download the needed field element
description files in various lat_element_name subdirectories
     e.g.  lat_ecr_venus   Venus source 
	   lat_s4          s4 solenoid 
	   lat_gag         grated acceleration gap
	   ...
Run this executable script using:
 % ./frib-lat-fields-dropbox-fetch

This script and the linked Dropbox account must be maintained consistently
as lattice element descriptions change. Code users may want access to the
Dropbox account which contains much info/input on the generation of lattice
element field files, code input (Poisson, CST Studio, ...), element plots,
etc. Contact Steve Lund for access.


