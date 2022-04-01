all: inject_HNLs inject_HNLs_MB
	
inject_HNLs:
	g++ inject_HNLs.cpp -o inject_dipole -I/usr/include/hdf5/serial/ -I${PREFIX}/include/ -I${PWD}/extern/rk/include/ -L${PREFIX}/lib/ -L${PREFIX}/lib64/ -L${PREFIX}/lib/ -L${PREFIX}/lib/pkgconfig/ -lcfitsio -lLeptonInjector -g

inject_HNLs_MB:
	g++ inject_HNLs_MB.cpp -o inject_dipole_MB -I/usr/include/hdf5/serial/ -I${PREFIX}/include/ -I${PWD}/extern/rk/include/ -L${PREFIX}/lib/ -L${PREFIX}/lib64/ -L${PREFIX}/lib/ -L${PREFIX}/lib/pkgconfig/ -lcfitsio -lLeptonInjector -g
