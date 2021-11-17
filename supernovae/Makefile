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
YES=$(shell sed -n 10p data/parameters.txt | cut -d '#' -f 1)

#Compiling parameters
CXX=g++
XXFLAGS=-std=c++17 -g -O3
SANITFLAGS =-fsanitize=address -fsanitize=leak -fsanitize=undefined
LIBS=
RUN=./
SOURCES=$(wildcard code/*.cpp)
DEPENDENCIES=$(SOURCES:code/%.cpp=objects/%.o)
FLAGS = $(XXFLAGS)

.PHONY: all clean oclean graph movie

all: main.x
	@echo -e '\n'
	@echo -e 'Running program ... \n'
	@$(RUN)$< $(T) $(DT) $(FRAGMENTS) $(R) $(FRAMES) $(MASSLOSS) $(ENERGY) $(TE) $(G) $(YES)
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
	@gnuplot scripts/Ederivative.gp
	@gnuplot scripts/CMvel.gp
	@echo -e 'Done!\n'

Gif:
	ffmpeg   -framerate 60   -pattern_type glob   -i 'Animation/Supernovae.*.png'   -r 15   -vf scale=512:-1   Supernovae.gif ;

clean:
	@rm -rf *.x results/* code/*~ data/*~ scripts/*~ main.x

oclean:
	@rm -rf objects/*.o
	
