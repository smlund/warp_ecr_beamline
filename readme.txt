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

To initialize the repository, 

   % git clone  git@github.com:smlund/warp_ion_frontend.git

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

The local copy of the file removed can be retained on disk by using 

  # git rm --cached file 

