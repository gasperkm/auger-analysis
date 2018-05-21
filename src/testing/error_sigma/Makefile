ROOTCONF=$(shell root-config --cflags --glibs)
ROOTINC=$(shell root-config --incdir)
ROOTLIB=$(shell root-config --libs)
INC=-I. -I$(ROOTINC)
LIBS=-L. $(ROOTLIB) -lMLP -lMinuit -lstdc++

all: error_sigma

error_sigma: error_sigma.cpp
	$(CXX) $(INC) -o error_sigma error_sigma.cpp $(LIBS)

clean:
	rm -fr error_sigma
