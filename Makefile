all: mktest



anno_import: import_anno.cpp 
	g++ -O3 -g import_anno.cpp -o import_anno -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl
BeQTL: mktest.cpp mktest.hpp
	icc -O3 -g mktest.cpp -o BeQTL -mkl -I/home/nwk2/include -I/opt/intel/composerexe-2011.5.220/mkl/include -L/home/nwk2/lib -L/opt/intel/composerexe-2011.5.220/mkl/lib/intel64/ -lhdf5 -lhdf5_hl 

import: snpexpfiles.cpp
	icc -O2 -g snpexpfiles.cpp -o import_snpgenes -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl

BooT.o: BooT.cpp BooT.hpp strideiter.hpp mat.o
	icc -O0 -g BooT.cpp mat.o -c -o BooT.o -I/home/nwk2/include -L/home/nwk2/lib -mkl=parallel
mat.o: arrayclass.cpp arrayclass.hpp
	icc -O0 -g arrayclass.cpp -c -o mat.o -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl -mkl=parallel


clean:
	rm -f mktest unit_tests import_snpgenes
