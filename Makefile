src = $(wildcard src/*.c) 
ccsrc = $(wildcard src/*.cpp)
cusrc = $(wildcard src/*.cu)

obj = $(src:.c=.o) $(ccsrc:.cpp=.o) $(cusrc:.cu=.o)

CXXFLAGS =-O3 -D NDEBUG -fopenmp -std=gnu++11 -march=native
LDFLAGS =-O3 -fopenmp -std=gnu++11 -lmgl -lgsl -lgslcblas -lnlopt -march=native

CUDA_ROOT_DIR=/usr/local/cuda
NVCC=nvcc
NVCCFLAGS=--cudart static --relocatable-device-code=false -gencode arch=compute_61,code=compute_61 -gencode arch=compute_61,code=sm_61 
NVCC_LIBS=

# CUDA library directory:
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include
# CUDA linking libraries:
CUDA_LINK_LIBS= -lcudart


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

%.o: %.c _kiss_fft_guts.h
	$(CC) $(CXXFLAGS) -o $@ -c $<

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

GridSim_cuda: $(obj)
	$(NVCC) $(NVCCFLAGS) $(obj) -link -o $@ $(CUDA_INC_DIR) $(CUDA_LIB_DIR) $(CUDA_LINK_LIBS) -lmgl -lgsl -lgslcblas -lnlopt -lgomp
all: GridSim_cuda

.PHONY: clean
clean:
	rm -f $(obj) GridSim_cuda

