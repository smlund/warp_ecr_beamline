# output Warp and envelope solver results for the target species


## Warp results
# in sequence: z-position, centroid, rms envelope size, emittance, ptheta

last_nonzero = amax(nonzero(hl_epsx[:,0]))

outfile = open("ThermalnoSCwarp.out", "a")

savetxt(outfile, hl_zbeam[0:last_nonzero][newaxis])
savetxt(outfile, hl_xcen[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_xcenp[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_ycen[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_ycenp[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_xrms[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_yrms[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_rrms[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_epsx[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_epsxn[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_epsy[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_epsyn[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_epsr[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_epsrn[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_pth[0:last_nonzero,1][newaxis])
savetxt(outfile, hl_pthn[0:last_nonzero,1][newaxis])

savetxt(outfile, hl_xcen[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_xcenp[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_ycen[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_ycenp[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_xrms[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_yrms[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_rrms[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_epsx[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_epsxn[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_epsy[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_epsyn[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_epsr[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_epsrn[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_pth[0:last_nonzero,0][newaxis])
savetxt(outfile, hl_pthn[0:last_nonzero,0][newaxis])

outfile.close()


## Envelope solver results

outfile2 = open("ThermalnoSCenv.out", "a")

savetxt(outfile2, sss[newaxis])
savetxt(outfile2, psoln[:,21][newaxis])
savetxt(outfile2, psoln[:,41][newaxis])

outfile2.close()
