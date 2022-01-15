GCC				= gcc
CFLAGS =  -ggdb -fopenmp -fPIC -O3 -Wall -Werror -Wno-error=cpp -fno-strict-aliasing
CFLAGS += -Wno-error=unused-function
CFLAGS += -Wno-error=unused-variable

HDF5DIR = /usr/lib/x86_64-linux-gnu/hdf5/serial
TOMLCDIR = /opt/mnt

CINCFLAGS = -I"$(HDF5DIR)/include" -I"$(TOMLCDIR)/include"
CLDFLAGS = -L"$(HDF5DIR)/lib" -lhdf5
CLDFLAGS += -L"$(TOMLCDIR)/lib" -ltoml

UVH5_HEADERS = uvh5.h uvh5_bool_t.h h5_dataspace.h
UVH5_TOML_HEADERS = uvh5_toml.h
UVH5_TOML_SOURCES = uvh5_toml.c
UVH5_TOML_FILES = $(UVH5_TOML_HEADERS) $(UVH5_TOML_SOURCES)

anew: clean run clean_build
	rm -f *.h.gch *.o

main.o: main.c $(UVH5_HEADERS) $(UVH5_TOML_FILES)
	$(GCC) $(CINCFLAGS) $(CFLAGS) -c $? $(CLDFLAGS)

main: main.o uvh5_toml.o
	$(GCC) $(CINCFLAGS) $(CFLAGS) $? -o $@ $(CLDFLAGS)

run: main
	./main

clean_build:
	rm -f *.h.gch *.o

clean: clean_build
	rm -f main