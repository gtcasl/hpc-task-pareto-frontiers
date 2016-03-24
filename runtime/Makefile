TARGET := run
SRC := $(shell ls *.cc) $(shell ls cg/*.cc)
HEADERS := $(shell ls *.h)

#target := gcc
target := darwin
#target := mic
#target := intel

CXX :=    mpicxx
CC :=     mpicc
CXXFLAGS := -fPIC
CPPFLAGS := -I. -std=c++11 
LIBDIR :=   
PREFIX :=   
LDFLAGS :=  -Wl,-rpath,$(PREFIX)/lib 

ifeq ($(target),darwin)
LDFLAGS := $(LDFLAGS) -stdlib=libc++
CPPFLAGS := $(CPPFLAGS) -Dno_affinity
endif

ifeq ($(target),intel)
CXX := $(CXX) -cxx=icpc
CC := $(CC) -cc=icc
endif

ifeq ($(target),mic)
CXXFLAGS := $(CXXFLAGS) -mmic
CXX := $(CXX) -cxx=icpc
CC := $(CC) -cc=icc
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

