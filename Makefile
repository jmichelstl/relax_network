CC=g++
FLAGS=-std=c++11
CUFLAG=--gpu-architecture=sm_60 

netrelax:
	$(CC) $(FLAGS) -o netrelax src/relax_network.cpp -pthread
curelax:
	nvcc src/cu_relax.cu -ccbin $(CC) $(FLAGS) $(CUFLAG) -I/usr/local/cuda-10.0/include -L/usr/local/lib -lcublas -o curelax
relclean:
	rm netrelax
crclean:
	rm curelax
