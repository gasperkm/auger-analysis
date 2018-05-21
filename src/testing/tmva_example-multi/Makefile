ROOTCONF=$(shell root-config --cflags --glibs)
ROOTINC=$(shell root-config --incdir)
ROOTLIB=$(shell root-config --libs)
INC=-I. -I$(ROOTINC) -I$(TMVASYS)/inc
LIBS=-L. -L$(TMVASYS)/lib -lTMVA $(ROOTLIB) -lMLP -lMinuit -lstdc++

all: tmva_simple tmvagui

tmva_simple: tmva_simple.cpp
	g++ $(INC) -o tmva_simple tmva_simple.cpp separate_functions.cpp mvaefficiency.cpp $(LIBS)

tmvagui: TMVAGui.C
	g++ -o tmvagui TMVAGui.C $(ROOTCONF)

tarball:
	mkdir tmva_example
	cp *.C example.root Makefile *.cpp tmva_example/
	tar -zcf tmva_example.tar.gz tmva_example/
	rm -fr tmva_example/

clean:
	rm -fr tmva_simple tmvagui *C.d *C.so
