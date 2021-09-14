#Parameters
T=$(shell sed -n 1p data/parameters.txt | tr -d -c 0-9.-)
DT=$(shell sed -n 2p data/parameters.txt | tr -d -c 0-9.-)
FRAGMENTS=$(shell sed -n 3p data/parameters.txt | tr -d -c 0-9.-)
R=$(shell sed -n 4p data/parameters.txt | tr -d -c 0-9.-)
FRAMES=$(shell sed -n 5p data/parameters.txt | tr -d -c 0-9.-)
MASSLOSS=$(shell sed -n 6p data/parameters.txt | tr -d -c 0-9.-)
ENERGY=$(shell sed -n 7p data/parameters.txt | tr -d -c 0-9.-)
TE=$(shell sed -n 8p data/parameters.txt | tr -d -c 0-9.-)
G=$(shell sed -n 9p data/parameters.txt | tr -d -c 0-9.-)

#Compiling parameters
CXX=g++
FLAGS= -std=c++17 -g
LIBS=
#RUN=mpirun --use-hwthread-cpus -np $(PROCCESORS) ./
RUN=./
SOURCES=$(wildcard code/*.cpp)
DEPENDENCIES=$(SOURCES:code/%.cpp=objects/%.o)

.PHONY: all clean oclean graph movie

all: main.x
	@echo -e '\n'
	@echo -e 'Running program ... \n'
	@$(RUN)$< $(T) $(DT) $(FRAGMENTS) $(R) $(FRAMES) $(MASSLOSS) $(ENERGY) $(TE) $(G)
	@echo -e '\n'
	@echo -e '\nDone!'

main.x: $(DEPENDENCIES)
	@echo -e '\n'
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ $(LIBS) -o $@
	@echo -e 'Done!\n'

objects/%.o: code/%.cpp
	@echo -e '\n'
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(LIBS) -o $@
	@echo -e 'Done!\n'
graph:
	@echo -e '\n'
	@echo -e 'Making graphs...\n'
	@gnuplot scripts/Baricenter.gp
	@gnuplot scripts/Mass.gp
	@gnuplot scripts/Senergy.gp
	@gnuplot scripts/Tenergy.gp
	@python3 scripts/coes.py
	@echo -e 'Done!\n'

clean:
	@rm -rf *.x results/* code/*~ data/*~ scripts/*~ main.x

oclean:
	@rm -rf objects/*.o
