TARGET := run
SRC := $(shell ls *.cc) $(shell ls cg/*.cc)
HEADERS := $(shell ls *.h)

target := gcc
#target := darwin
#target := mic
target := intel

MPI := $(HOME)/Programs/install/mpich2/bin/
MPI :=

CXX :=    $(MPI)mpicxx
CC :=     $(MPI)mpicc
CXXFLAGS := -fPIC -O0 -g -Wall -Wextra
CPPFLAGS := -I. -std=c++11 -mkl
LIBDIR :=   
PREFIX :=   
LDFLAGS :=  -Wl,-rpath,$(PREFIX)/lib

ifeq ($(target),gcc)
CXXFLAGS := $(CXXFLAGS) -fopenmp
endif

ifeq ($(target),darwin)
LDFLAGS := $(LDFLAGS) -stdlib=libc++
CPPFLAGS := $(CPPFLAGS) -Dno_affinity
endif

ifeq ($(target),intel)
CXX := $(CXX) -cxx=icpc
CC := $(CC) -cc=icc
CXXFLAGS := $(CXXFLAGS) -openmp
LDFLAGS := $(LDFLAGS) -lmkl_rt
endif

ifeq ($(target),mic)
CXXFLAGS := $(CXXFLAGS) -mmic -openmp
CXX := $(CXX) -cxx=icpc
CC := $(CC) -cc=icc
LDFLAGS := $(LDFLAGS) -L$(MKLROOT)/lib/mic -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core
endif

OBJ := $(SRC:.cc=.o) 
OBJ := $(OBJ:.cpp=.o)
OBJ := $(OBJ:.c=.o)

.PHONY: clean install 

all: $(TARGET)

$(TARGET): $(OBJ) 
	$(CXX) -o $@ $+ $(CXXFLAGS) $(LDFLAGS) $(LIBS) 

%.o: %.cc  $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean: 
	rm -f $(TARGET) $(OBJ) 

install: $(TARGET)
	cp $< $(PREFIX)/bin

