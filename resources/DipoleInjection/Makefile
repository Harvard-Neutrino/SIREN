all: inject_HNLs_CCM

#args=-L/n/home01/awen/LeptonInjector/local/lib/ -L/n/home01/awen/LeptonInjector/local/lib64/ -I/n/home01/awen/LeptonInjector/local/include/ -I/n/home01/awen/miniconda3/envs/LI_env/include/ -I/n/home01/awen/homebrew/opt/hdf5/include/ -I/n/home01/awen/LeptonInjector/sources/LeptonInjector/vendor/cereal/include/ -I/n/home01/awen/LeptonInjector/sources/LeptonInjector/vendor/rk/rk/ -L/n/home01/awen/LeptonInjector/local/lib/pkgconfig/ -L/n/home01/awen/LeptonInjector/sources/LeptonInjector/ -g -std=c++14 -lm -ldl

#args=-L/n/home01/awen/LeptonInjector/local/lib/ -I/n/home01/awen/LeptonInjector/local/include/ -I/n/home01/awen/homebrew/opt/hdf5/include/ -I/usr/include/hdf5/serial/ -I/n/home01/awen/LeptonInjector/sources/LeptonInjector/vendor/cereal/include/ -I/n/home01/awen/LeptonInjector/sources/LeptonInjector/vendor/rk/rk/ -I/n/home01/awen/homebrew/opt/cfitsio/include/ -L/n/home01/awen/homebrew/lib -L/n/home01/awen/LeptonInjector/local/lib/pkgconfig/ -lcfitsio -L/n/home01/awen/LeptonInjector/sources/LeptonInjector/ -g -std=c++14 -lm -ldl

#args=-I/usr/include/hdf5/serial/ -I${PREFIX}/include/ -I${PWD}/extern/rk/include/ -I/opt/homebrew/opt/cfitsio/include/ -L/opt/homebrew/lib/ -L${PREFIX}/lib/ -L${PREFIX}/lib/ -L${PREFIX}/lib/pkgconfig/ -lcfitsio -lLeptonInjector -g -std=c++14

args=-I/usr/include/hdf5/serial/ -I${PREFIX}/include/ -I${PWD}/extern/rk/include/ -I/opt/homebrew/opt/cfitsio/include/ -L/opt/homebrew/lib/ -L${PREFIX}/lib/ -L${PREFIX}/lib64/ -L${PREFIX}/lib/pkgconfig/ -L${SROOT}/lib/ -lcfitsio -lLeptonInjector -g -std=c++14 -L/usr/local/lib64/

inject_HNLs_CCM:
	g++ inject_HNLs_CCM.cpp -o inject_HNLs_CCM ${args}
