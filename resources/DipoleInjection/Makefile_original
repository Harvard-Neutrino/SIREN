all: inject_HNLs_MINERvA inject_HNLs_MB inject_HNLs_CCM

args=-I/usr/include/hdf5/serial/ -I${PREFIX}/include/ -I${PWD}/extern/rk/include/ -I/opt/homebrew/opt/cfitsio/include/ -L/opt/homebrew/lib/ -L${PREFIX}/lib/ -L${PREFIX}/lib/ -L${PREFIX}/lib/pkgconfig/ -lcfitsio -lLeptonInjector -g -std=c++14

inject_HNLs_MINERvA:
	g++ inject_HNLs_MINERvA.cpp inject_HNLs_MINERvA ${args}

inject_HNLs_MB:
	g++ inject_HNLs_MB.cpp -o inject_HNLs_MB ${args}

inject_HNLs_CCM:
	g++ inject_HNLs_CCM.cpp -o inject_HNLs_CCM ${args}
