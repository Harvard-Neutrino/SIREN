all: inject_HNLs
	
inject_HNLs:
	g++ inject_HNLs.cpp -o inject_dipole -I/usr/include/hdf5/serial/ -I${PREFIX}/include/ -I${PWD}/extern/rk/include/ -L${PREFIX}/lib/ -L${PREFIX}/lib64/ -L${PREFIX}/lib/ -L${PREFIX}/lib/pkgconfig/ -lcfitsio -lLeptonInjector -g
