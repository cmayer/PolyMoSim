# Makefile for PolyMoSim
# SPDX-License-Identifier: GPL-3.0-only
#
#   make               Build the program PolyMoSim
#   make install       Copy PolyMoSim to $(DESTDIR)$(PREFIX)/bin (default: /usr/local/bin)
#   make test          Build PolyMoSim and run the tests in tests/
#   make clean         Remove the build products
#
# CXX, CXXFLAGS and LDFLAGS can be set in the environment or on the command line, e.g.
#   make CXX=clang++ CXXFLAGS="-O2 -g"
#   make install PREFIX=$HOME/.local

CXX      ?= g++
CXXFLAGS ?= -O3 -Wall -Wextra # -static ## Static linking is recommended e.g. for mingw,

# Older build scripts pass the compiler as CPP=... on the command line. This still works.
ifeq ($(origin CPP),command line)
CXX := $(CPP)
endif

PREFIX   ?= /usr/local
BINDIR   ?= $(PREFIX)/bin

BIN      = PolyMoSim

LINKOBJ  = tree_admin.o model_admin.o mymodel.o BasicNode.o BasicTree.o CRandom.o global-types-and-parameters.o split_admin.o BSplit.o $(RES) discrete_gamma.o models-aa.o
OBJ      = $(LINKOBJ) PolyMoSim.o

RM       = rm -f

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o $(BIN) $(LIBS)

install: $(BIN)
	mkdir -p $(DESTDIR)$(BINDIR)
	cp $(BIN) $(DESTDIR)$(BINDIR)/$(BIN)
	chmod 755 $(DESTDIR)$(BINDIR)/$(BIN)

test: $(BIN)
	sh tests/run_tests.sh ./$(BIN)

clean:
	${RM} $(OBJ) $(BIN)

tree_admin.o: tree_admin.cpp tree_admin.h BasicTree.h PolyMoSim.h
	$(CXX) -c tree_admin.cpp -o tree_admin.o $(CXXFLAGS)

model_admin.o: model_admin.cpp model_admin.h mymodel.h PolyMoSim.h
	$(CXX) -c model_admin.cpp -o model_admin.o $(CXXFLAGS)

mymodel.o: mymodel.cpp mymodel.h PolyMoSim.h # gamma_random.h
	$(CXX) -c mymodel.cpp -o mymodel.o $(CXXFLAGS)

models-aa.o: models-aa.cpp mymodel.h PolyMoSim.h # gamma_random.h
	$(CXX) -c models-aa.cpp -o models-aa.o $(CXXFLAGS)

BasicNode.o: BasicNode.cpp BasicNode.h mymodel.h
	$(CXX) -c BasicNode.cpp -o BasicNode.o $(CXXFLAGS)

BasicTree.o: BasicTree.cpp BasicTree.h BasicNode.h model_admin.h PolyMoSim.h
	$(CXX) -c BasicTree.cpp -o BasicTree.o $(CXXFLAGS)

discrete_gamma.o: discrete_gamma.cpp discrete_gamma.hpp
	$(CXX) -c discrete_gamma.cpp -o discrete_gamma.o $(CXXFLAGS)

CRandom.o: CRandom.cpp CRandom.h
	$(CXX) -c CRandom.cpp -o CRandom.o $(CXXFLAGS)

split_admin.o: split_admin.cpp split_admin.h fast-dynamic-bitset.h BSplit.h
	$(CXX) -c split_admin.cpp -o split_admin.o $(CXXFLAGS)

BSplit.o: BSplit.cpp BSplit.h fast-dynamic-bitset.h
	$(CXX) -c BSplit.cpp -o BSplit.o $(CXXFLAGS)

global-types-and-parameters.o: global-types-and-parameters.cpp global-types-and-parameters.hpp
	$(CXX) -I . -c global-types-and-parameters.cpp -o global-types-and-parameters.o $(CXXFLAGS)

PolyMoSim.o: PolyMoSim.cpp PolyMoSim.h symbol-cartesian-product.h
	$(CXX) -c PolyMoSim.cpp -o PolyMoSim.o $(CXXFLAGS)

.PHONY: all install test clean
