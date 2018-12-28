 
# Make variables ----------------------------------------------------

# Setup for wxWidgets and ROOT
ROOTCONF=$(shell root-config --cflags --libs)
ROOTINC=$(shell root-config --incdir)
ROOTLIB=$(shell root-config --libs)

# Code position directories
SRC=../../../src
IDIR=../../../include

# Compiler options
COMPOPT=-fPIC -g -std=c++11 -Wno-deprecated
#COMPOPT=-fPIC -std=c++11

# Include directories
INC=-I. -I$(IDIR) -I$(ADSTROOT)/include/adst -I$(ROOTINC)

# Used libraries
OFFLINELIBS=-lRecEventKG -lAnalysisKG
LIBS=-L. -L$(ADSTROOT)/lib $(OFFLINELIBS) $(ROOTLIB) -lMLP -lMinuit -lstdc++
# -------------------------------------------------------------------

# Base rules --------------------------------------------------------

# General make options
all: adst_quick_plot

adst_quick_plot: adst_quick_plot.cpp $(SRC)/separate_functions.cpp $(SRC)/root_style.cpp $(SRC)/mva_result_read.cpp $(SRC)/mva_methods.cpp $(SRC)/primary_type.cpp
	$(CXX) $(COMPOPT) $(INC) adst_quick_plot.cpp $(SRC)/separate_functions.cpp $(SRC)/root_style.cpp $(SRC)/mva_result_read.cpp $(SRC)/mva_methods.cpp $(SRC)/primary_type.cpp -o $@ $(LIBS)

# Rules for cleaning the installation
clean:
	rm -fr adst_quick_plot
# -------------------------------------------------------------------
