CPP      = g++ -Wall -Wextra
# CPP      = g++-4 -Wall -Wextra


# WINDRES = windres.exe
RES      = 
#OBJ      = tree_admin.o sequence_gen.o pvartree.o model_admin.o mymodel.o BasicNode.o BasicTree.o CRandom.o global-types-and-parameters.o $(RES) PolyMoSim.o pvarmain.o
#LINKOBJ  = tree_admin.o sequence_gen.o pvartree.o model_admin.o mymodel.o BasicNode.o BasicTree.o CRandom.o global-types-and-parameters.o $(RES)

OBJ      = tree_admin.o model_admin.o mymodel.o BasicNode.o BasicTree.o CRandom.o global-types-and-parameters.o split_admin.o BSplit.o $(RES) discrete_gamma.o PolyMoSim.o 
LINKOBJ  = tree_admin.o model_admin.o mymodel.o BasicNode.o BasicTree.o CRandom.o global-types-and-parameters.o split_admin.o BSplit.o $(RES) discrete_gamma.o

LIBS     =   
## INCS     = -I ../../../Klassen/
## CXXINCS  = -I ../../../Klassen/ 
BIN      =   PolyMoSim-v1.1.4

CXXFLAGS = $(CXXINCS)  -g3 -ggdb # -static ## Static linking is recommended e.g. for mingw, 
CFLAGS   = $(INCS)  -g3 -ggdb    # -static ## Static linking is recommended e.g. for mingw, 

#CXXFLAGS = $(CXXINCS)  -O2 # -static ## Static linking is recommended e.g. for mingw, 
#CFLAGS   = $(INCS)     -O2 # -static ## Static linking is recommended e.g. for mingw, 

RM       = rm -f

all: PolyMoSim-v1.1.4


clean:
	${RM} $(OBJ) $(BIN)

# $(BIN): $(OBJ)
# 	$(CPP) $(LINKOBJ) -o $(BIN) $(LIBS)

PolyMoSim-v1.1.4: $(OBJ) PolyMoSim.o
	$(CPP) $(LINKOBJ) PolyMoSim.o -o PolyMoSim-v1.1.4 $(LIBS)

tree_admin.o: tree_admin.cpp tree_admin.h BasicTree.h PolyMoSim.h
	$(CPP) -c tree_admin.cpp -o tree_admin.o $(CXXFLAGS)

model_admin.o: model_admin.cpp model_admin.h mymodel.h PolyMoSim.h
	$(CPP) -c model_admin.cpp -o model_admin.o $(CXXFLAGS)

mymodel.o: mymodel.cpp mymodel.h PolyMoSim.h # gamma_random.h
	$(CPP) -c mymodel.cpp -o mymodel.o $(CXXFLAGS)

BasicNode.o: BasicNode.cpp BasicNode.h mymodel.h
	$(CPP) -c BasicNode.cpp -o BasicNode.o $(CXXFLAGS)

BasicTree.o: BasicTree.cpp BasicTree.h BasicNode.h model_admin.h PolyMoSim.h
	$(CPP) -c BasicTree.cpp -o BasicTree.o $(CXXFLAGS)

discrete_gamma.o: discrete_gamma.c discrete_gamma.h
	$(CPP) -c discrete_gamma.c -o discrete_gamma.o $(CXXFLAGS)

CRandom.o: CRandom.cpp CRandom.h
	$(CPP) -c CRandom.cpp -o CRandom.o $(CXXFLAGS)

split_admin.o: split_admin.cpp split_admin.h fast-dynamic-bitset.h BSplit.h
	$(CPP) -c split_admin.cpp -o split_admin.o $(CXXFLAGS)

BSplit.o: BSplit.cpp BSplit.h fast-dynamic-bitset.h
	$(CPP) -c BSplit.cpp -o BSplit.o $(CXXFLAGS)

global-types-and-parameters.o: global-types-and-parameters.cpp global-types-and-parameters.h
	$(CPP)	 -I . -c global-types-and-parameters.cpp -o global-types-and-parameters.o $(CXXFLAGS)

PolyMoSim.o: PolyMoSim.cpp PolyMoSim.h symbol-cartesian-product.h
	$(CPP) -c PolyMoSim.cpp -o PolyMoSim.o $(CXXFLAGS)

