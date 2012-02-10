CC=/usr/bin/g++ -g -Wall
CXX=/curc/tools/nonfree/redhat_5_x86_64/ics_2012.0.032/composer_xe_2011_sp1.6.233/bin/intel64/icpc

HDF5_INCLUDE=-I/curc/tools/free/redhat_5_x86_64/hdf5-1.8.7/include
HDF5_LIB_DIR=-L/curc/tools/free/redhat_5_x86_64/hdf5-1.8.7/lib
HDF5_LIBS=-lhdf5
BLOCK=-DUSE_BLOCK=1
INTEL=-DUSE_INTEL=1

#Intel Stuff:
MKL_INCLUDE=-I/curc/tools/nonfree/redhat_5_x86_64/ics_2012.0.032/composer_xe_2011_sp1.6.233/mkl/include
MKL_LIB_DIR=-L/curc/tools/nonfree/redhat_5_x86_64/ics_2012.0.032/composer_xe_2011_sp1.6.233/mkl/lib/intel64
MKL_LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core


INCLUDE=$(HDF5_INCLUDE)
LIB_DIR=$(HDF5_LIB_DIR)
LIBS=$(HDF5_LIBS)


all:
	make clean
	make naive
	make block
	make intel
	make read_write

naive:
	$(CC) lab2.cpp $(INCLUDE) $(LIB_DIR) $(LIBS) -o naive_multiply

block:
	$(CC) $(BLOCK) lab2.cpp $(INCLUDE) $(LIB_DIR) $(LIBS) -o block_multiply 

intel:
	$(CXX) $(INTEL) lab2.cpp $(INCLUDE) $(LIB_DIR) $(LIBS) $(MKL_INCLUDE) $(MKL_LIB_DIR) $(MKL_LIBS) -o intel_multiply 

read_write:
	$(CC) lab2.cpp $(INCLUDE) $(LIB_DIR) $(LIBS) -o read_write_matrix

clean:
	rm -f naive_multiply block_multiply intel_multiply read_write_matrix



