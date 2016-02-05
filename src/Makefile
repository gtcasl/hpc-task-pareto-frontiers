TARGET := spmv
SRC := main.cc setup.cc mult.cc

# the sstmacro-config script must be found in PATH
CXX :=   g++
CC :=    gcc
CXXFLAGS := -O3 -fopenmp
CPPFLAGS := -I. -I/home/eric/.local/stow/lwperf-host/include -DLWPERF_USE_CSV
LIBDIR :=   
PREFIX :=   
LDFLAGS := -L/home/eric/.local/stow/lwperf-host/lib
LIBS=-llwperf

OBJ := $(SRC:.cc=.o) 

.PHONY: clean install 

all: $(TARGET)

$(TARGET): $(OBJ) 
	$(CXX) -o $@ $+ $(LDFLAGS) $(CXXFLAGS) $(LIBS) 

%.o: %.cc index.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


clean: 
	rm -f $(TARGET) $(OBJ) 


