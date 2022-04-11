all: inject_HNLs_MINERvA inject_HNLs_MB
	
inject_HNLs_MINERvA:
	g++ inject_HNLs_MINERvA.cpp -o inject_dipole_MINERvA -I/usr/include/hdf5/serial/ -I${PREFIX}/include/ -I${PWD}/extern/rk/include/ -L${PREFIX}/lib/ -L${PREFIX}/lib64/ -L${PREFIX}/lib/ -L${PREFIX}/lib/pkgconfig/ -lcfitsio -lLeptonInjector -g

inject_HNLs_MB:
	g++ inject_HNLs_MB.cpp -o inject_dipole_MB -I/usr/include/hdf5/serial/ -I${PREFIX}/include/ -I${PWD}/extern/rk/include/ -L${PREFIX}/lib/ -L${PREFIX}/lib64/ -L${PREFIX}/lib/ -L${PREFIX}/lib/pkgconfig/ -lcfitsio -lLeptonInjector -g
