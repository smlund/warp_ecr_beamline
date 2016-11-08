# Run as stand alone script after scanning data is generated
# Find the point closest to target values among a set of alpha-beta points

import numpy as np

noSClist = np.loadtxt('ScanSC60Zoom.out')

len1 = len(noSClist)

min_diff = 10.
min_pos_1 = 0
min_pos_2 = 0

alpha_target = 0.306
beta_target = 4.383

for ii in range(len1):
	penalty = abs(noSClist[ii,2] - beta_target) + 2*abs(noSClist[ii,3] - alpha_target)
	
	if penalty < min_diff:
		min_diff = penalty
		min_pos_1 = ii

print noSClist[min_pos_1, :]
print min_diff
		
