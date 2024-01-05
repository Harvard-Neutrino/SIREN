import numpy as np
fluxTotal = np.loadtxt('mbflux.txt')
Emid = 0.5*(fluxTotal[:,0]+fluxTotal[:,1])
bnb_flux = {}
for i,key in enumerate(['numu','numubar','nue','nuebar']):
	bnb_flux[key] = np.zeros((len(Emid),2))
	bnb_flux[key][:,0] = Emid
	bnb_flux[key][:,1] = fluxTotal[:,i+2]*(1000./50.)
	np.savetxt('BNB_%s_flux.txt'%key,bnb_flux[key])
