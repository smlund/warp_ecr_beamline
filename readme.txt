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

To initialize the repository, 
 https form
   % git clone https://github.com/smlund/warp_ion_frontend.git 
 ssh form (can have permission problems depending on setup) 
   % git clone git@github.com:smlund/warp_ion_frontend.git 

This will create a directory, ./warp_ion_frontend where command was 
run with the archive files.   

To get the latest version descend into directory warp_ion_frontend and run:

  % git pull 

When modifying the repository (for those with edit privilege, please contact 
me if you want to contribute and I will add you) 

  ... edit files etc then checkin (use readme.txt as example here) using: 
  % git add readme.txt 
  % git commit -m "SML: updated readme.txt file" 
  % git push 

You may be need to enter your git username and password in this step.  

To remove a file from git control (to not include in future pulls), use 

  % git rm file
  % git commit -m "SML: removing file from repo"
  % git push 

The local copy of the file removed can be retained on disk by subsituting

  % git rm --cached file 

in the above.  Another way to remove a file from the master node repo is to
go to the github web side, click on a file, click the delete button and then
confirm at the bottom of the page.

The file removal procedures above  will NOT remove the file from the
history within git. This is preferable in most cases.  However, when a git
repo is cloned, if a large file (say binary) was contained at one point and
deleted, it will be downloaded (and then deleted) in progression used
to generate the repo on the local machine. This can be a problem. To remove
all traces of such a file "file.pkl" from the repo, you can run:

  % git filter-branch -f --index-filter 'git rm --cached --ignore-unmatch file.pkl'

Another command that can help:

  % git filter-branch -f --tree-filter "rm -rf *.pkl" --prune-empty -- --all

This should ONLY be used in cases where it is really needed (large files that
should not be there).  

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
Dropbox account which contains much info/input on the generation of the 
lattice. Contact Steve Lund for access if you have a dropbox account and this 
account can be "shared" with you. Note that this directories is large due to the 
size of field element arrays and shared account contribute to a dropbox users 
quota. So a larger (paid) dropbox account would likely be needed for full sharing. 
element field files, code input (Poisson, CST Studio, ...), element plots,
etc.


