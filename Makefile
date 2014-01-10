all: mktest

hworld: hWorld.cpp 
	g++ -O3 -g hWorld.cpp -o hworld -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl
mktest: mktest.cpp mktest.hpp
	icc -O3 -g mktest.cpp -o BooTQT -mkl -I/home/nwk2/include -I/opt/intel/composerexe-2011.5.220/mkl/include -L/home/nwk2/lib -L/opt/intel/composerexe-2011.5.220/mkl/lib/intel64/ -lhdf5 -lhdf5_hl 
multest: multest.cpp mktest.hpp
	icc -O3 -g multest.cpp -fp-model fast=2 -ipo -o multest -mkl=parallel -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl
test:	unit_tests.cpp
	icc  -O0 -g unit_tests.cpp -debug -o unit_tests -mkl -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl
import: snpexpfiles.cpp
	icc -O2 -g snpexpfiles.cpp -o import_snpgenes -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl
stest: stest.cpp mktest.hpp
	icc -O2 -g stest.cpp -o stest -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl -mkl=parallel
maptest: maptest.cpp
	icc -O2 -g maptest.cpp -o maptest -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl -mkl=parallel
sortest: sortest.cpp mktest.hpp
	icc -O0 -g sortest.cpp  -o sortest -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl -mkl=parallel
mattest: mattest.cpp mat.o BooT.o
	icc -O0 -g mattest.cpp mat.o BooT.o -o mattest -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl -mkl=parallel
BooT.o: BooT.cpp BooT.hpp strideiter.hpp mat.o
	icc -O0 -g BooT.cpp mat.o -c -o BooT.o -I/home/nwk2/include -L/home/nwk2/lib -mkl=parallel
mat.o: arrayclass.cpp arrayclass.hpp
	icc -O0 -g arrayclass.cpp -c -o mat.o -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl -mkl=parallel


clean:
	rm -f mktest unit_tests import_snpgenes
