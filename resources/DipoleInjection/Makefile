all: inject_HNLs_CCM

args=-L/n/home01/awen/LeptonInjector/local/lib -I/n/home01/awen/LeptonInjector/local/include/ -I/usr/include/hdf5/serial/ -I/n/home01/awen/LeptonInjector/sources/LeptonInjector/vendor/cereal/include/ -I${PREFIX}/include/ -I${PWD}/extern/rk/include/ -I/n/home01/awen/homebrew/opt/cfitsio/include -L/n/home01/awen/homebrew/lib -L${PREFIX}/lib/ -L${PREFIX}/lib/pkgconfig/ -lcfitsio -g -std=c++14

inject_HNLs_CCM:
	g++ inject_HNLs_CCM.cpp -o inject_HNLs_CCM ${args}
