src = $(wildcard src/*.c) 
ccsrc = $(wildcard src/*.cpp)
obj = $(src:.c=.o) $(ccsrc:.cpp=.o)
objDebug = $(src:.c=.oDebug) $(ccsrc:.cpp=.oDebug)

CXXFLAGS =-O3 -D NDEBUG -fopenmp -std=gnu++11 -march=native
LDFLAGS =-O3 -fopenmp -std=gnu++11 -lmgl -lgsl -lgslcblas -march=skylake-avx512

CXXDEBUGFLAGS =-g -O0 -std=gnu++11 -march=native 
LDDEBUGFLAGS =-g -O0 -std=gnu++11 -lmgl -lgsl -lgslcblas -march=skylake-avx512

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

%.o: %.c
	$(CC) $(CXXFLAGS) -o $@ -c $<

%.oDebug: %.c
	$(CC) $(CXXDEBUGFLAGS) -o $@ -c $<

%.oDebug: %.cpp
	$(CXX) $(CXXDEBUGFLAGS) -o $@ -c $<

a.out: $(obj)
	$(CXX) -o $@ $^ $(LDFLAGS)
all: a.out

debug.out: $(objDebug)
	$(CXX) -o $@ $^ $(LDDEBUGFLAGS)
debug: debug.out


.PHONY: clean
clean:
	rm -f $(obj) $(objDebug) debug.out a.out

