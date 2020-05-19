src = $(wildcard src/*.c) 
ccsrc = $(wildcard src/*.cpp)
obj = $(src:.c=.o) $(ccsrc:.cpp=.o)
objDebug = $(src:.c=.oDebug) $(ccsrc:.cpp=.oDebug)

CXXFLAGS =-O3 -D NDEBUG -fopenmp -std=gnu++11 -march=native
LDFLAGS =-O3 -fopenmp -std=gnu++11 -lgsl -lgslcblas -march=native

CXXDEBUGFLAGS =-g -O0 -std=gnu++11 -march=native 
LDDEBUGFLAGS =-g -O0 -std=gnu++11 -lgsl -lgslcblas -march=native

INC=-I/home/ivancicr/gsl/include
LIB=-L/home/ivancicr/gsl/lib

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LIB) $(INC) -o $@ -c $<

%.o: %.c
	$(CC) $(CXXFLAGS) $(LIB) $(INC) -o $@ -c $<

%.oDebug: %.c
	$(CC) $(CXXDEBUGFLAGS) $(LIB) $(INC) -o $@ -c $<

%.oDebug: %.cpp
	$(CXX) $(CXXDEBUGFLAGS) $(LIB) $(INC) -o $@ -c $<

a.out: $(obj)
	$(CXX) $(LIB) $(INC) -o $@ $^ $(LDFLAGS)
all: a.out

debug.out: $(objDebug)
	$(CXX) $(LIB) $(INC) -o $@ $^ $(LDDEBUGFLAGS)
debug: debug.out


.PHONY: clean
clean:
	rm -f $(obj) $(objDebug) debug.out a.out

