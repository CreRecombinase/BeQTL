all: mktest

mktest: mktest.cpp
	icc -O3 -g mktest.cpp -o mktest.o -mkl=parallel -I/home/nwk2/include -L/home/nwk2/lib -lhdf5 -lhdf5_hl

clean:
	rm -f mktest.o