#Parameters
#NPROC=$(shell sed -n 1p Data/parameters.txt | tr -d -c 0-9.-)
#N=$(shell sed -n 2p Data/parameters.txt | tr -d -c 0-9.-)
#D=$(shell sed -n 3p Data/parameters.txt | tr -d -c 0-9.-)
#SIGMA=$(shell sed -n 4p Data/parameters.txt | tr -d -c 0-9.-) 
#EPSILON=$(shell sed -n 5p Data/parameters.txt | tr -d -c 0-9.-)
#M=$(shell sed -n 6p Data/parameters.txt | tr -d -c 0-9.-)
#CUTT=$(shell sed -n 7p Data/parameters.txt | tr -d -c 0-9.-)
#T=$(shell sed -n 8p Data/parameters.txt | tr -d -c 0-9.-)
#DAMP=$(shell sed -n 9p Data/parameters.txt | tr -d -c 0-9.-)
#ITER=$(shell sed -n 10p Data/parameters.txt | tr -d -c 0-9.-)
#DT=$(shell sed -n 11p Data/parameters.txt | tr -d -c 0-9.-)
#LX=$(shell sed -n 12p Data/parameters.txt | tr -d -c 0-9.-)
#LY=$(shell sed -n 13p Data/parameters.txt | tr -d -c 0-9.-)
#LZ=$(shell sed -n 14p Data/parameters.txt | tr -d -c 0-9.-)
#SAMPLE=$(shell sed -n 15p Data/parameters.txt | tr -d -c 0-9.-)
#BINS=$(shell sed -n 16p Data/parameters.txt | tr -d -c 0-9.-)
#SEED=$(shell sed -n 17p Data/parameters.txt | tr -d -c 0-9.-)
#Thr=$(shell sed -n 18p Data/parameters.txt | tr -d -c 0-9.-)

#Compiling parameters
CXX=g++
FLAGS= -std=c++17
#RUN=mpirun --use-hwthread-cpus -np $(PROCCESORS) ./
RUN=./
SOURCES=$(wildcard code/*.cpp)
DEPENDENCIES=$(SOURCES:Code/%.cpp=objects/%.o)

.PHONY: all clean oclean graph movie

all: main.x
	@echo -e '\n'
	@echo -e 'Running program ... \n'
	@$(RUN)$< 
	@echo -e '\n'
	@echo -e '\nDone!'

main.x: $(DEPENDENCIES)
	@echo -e '\n'
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

objects/%.o: code/%.cpp
	@echo -e '\n'
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'
graph:
	@echo -e '\n'
	@echo -e 'Making graphs...\n'
	@echo -e 'Done!\n'

clean:
	@rm -rf *.x results/* code/*~ data/*~ scripts/*~ main.x

oclean:
	@rm -rf objects/*.o
