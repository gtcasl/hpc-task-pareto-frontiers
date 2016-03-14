TARGET := run
SRC := $(shell ls *.cc) 
HEADERS := $(shell ls *.h)

target := darwin

CXX :=    mpicxx
CC :=     mpicc
CXXFLAGS := -fPIC
CPPFLAGS := -I. -std=c++11
LIBDIR :=   
PREFIX :=   
LDFLAGS :=  -Wl,-rpath,$(PREFIX)/lib -stdlib=libc++

ifeq ($(target),darwin)
LDFLAGS := $(LDFLAGS) -stdlib=libc++
CPPFLAGS := $(CPPFLAGS) -Dno_affinity
endif

OBJ := $(SRC:.cc=.o) 
OBJ := $(OBJ:.cpp=.o)
OBJ := $(OBJ:.c=.o)

.PHONY: clean install 

all: $(TARGET)

$(TARGET): $(OBJ) 
	$(CXX) -o $@ $+ $(LDFLAGS) $(LIBS) 

%.o: %.cc  $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean: 
	rm -f $(TARGET) $(OBJ) 

install: $(TARGET)
	cp $< $(PREFIX)/bin

