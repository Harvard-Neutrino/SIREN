all: inject_muons inject_HNLs
	
inject_muons:
	g++ inject_muons.cpp -o inject_muons -I/usr/include/hdf5/serial/ -I${PREFIX}/include/ -L${PREFIX}/lib/ -L${PREFIX}/lib64/ -L${PREFIX}/lib/pkgconfig/ -lLeptonInjector -g
	
inject_HNLs:
	g++ inject_HNLs.cpp -o inject_HNLs -I/usr/include/hdf5/serial/ -I${PREFIX}/include/ -L${PREFIX}/lib/ -L${PREFIX}/lib64/ -L${PREFIX}/lib/pkgconfig/ -lLeptonInjector -g
