CC=/usr/bin/g++

HDF5_INCLUDE=-I/curc/tools/free/redhat_5_x86_64/hdf5-1.8.7/include
HDF5_LIB_DIR=-L/curc/tools/free/redhat_5_x86_64/hdf5-1.8.7/lib
HDF5_LIBS=-lhdf5

INCLUDE=$(HDF5_INCLUDE)
LIB_DIR=$(HDF5_LIB_DIR)
LIBS=$(HDF5_LIBS)

all:
	rm -f naive_multiply
	make naive_multiply

naive_multiply:
	$(CC) naive_multiply.cpp $(INCLUDE) $(LIB_DIR) $(LIBS) -o naive_multiply

clean:
	rm naive_multiply data-1d.hdf5



