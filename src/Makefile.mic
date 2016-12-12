TARGET := run
#SRC := $(shell ls *.cc) $(shell ls cg/*.cc) cholesky/cholesky.cc cholesky/choleskyprofile.cc
SRC := $(shell ls *.cc) $(shell ls cg/*.cc) cholesky/cholesky.cc cholesky/choleskyprofile.cc lu/lu.cc lu/luprofile.cc
HEADERS := $(shell ls *.h) cholesky/cholesky.h lu/lu.h

CXX :=    icpc
CXXFLAGS := -fPIC -O2 -g -mkl -w3 -diag-disable remark -Wunused-variable
CPPFLAGS := -I. -std=gnu++11 -qoffload -pthread -qopenmp
LDFLAGS :=  -Wl,-rpath,$(PREFIX)/lib
LIBS := -lmkl_rt -lmicmgmt -lrt

OBJ := $(SRC:.cc=.o) 

.PHONY: clean install 

all: $(TARGET)

$(TARGET): $(OBJ) 
	$(CXX) -o $@ $+ $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(LIBS)

%.o: %.cc  $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean: 
	rm -f $(TARGET) $(OBJ)

