 
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
# -------------------------------------------------------------------

# Base rules --------------------------------------------------------

# General make options
all: multiplot_fractions multiplot_fractions_pub

multiplot_fractions: multiplot_fractions.cpp $(SRC)/separate_functions.cpp $(SRC)/root_style.cpp $(SRC)/mva_result_read.cpp $(SRC)/mva_methods.cpp $(SRC)/primary_type.cpp
	$(CXX) $(COMPOPT) -I. -I$(IDIR) -I$(ROOTINC) multiplot_fractions.cpp $(SRC)/separate_functions.cpp $(SRC)/root_style.cpp $(SRC)/mva_result_read.cpp $(SRC)/mva_methods.cpp $(SRC)/primary_type.cpp -o $@ -L. $(ROOTCONF) -lMLP -lMinuit -lstdc++

multiplot_fractions_pub: multiplot_fractions_pub.cpp $(SRC)/separate_functions.cpp $(SRC)/root_style.cpp $(SRC)/mva_result_read.cpp $(SRC)/mva_methods.cpp $(SRC)/primary_type.cpp
	$(CXX) $(COMPOPT) -I. -I$(IDIR) -I$(ROOTINC) multiplot_fractions_pub.cpp $(SRC)/separate_functions.cpp $(SRC)/root_style.cpp $(SRC)/mva_result_read.cpp $(SRC)/mva_methods.cpp $(SRC)/primary_type.cpp -o $@ -L. $(ROOTCONF) -lMLP -lMinuit -lstdc++

# Rules for cleaning the installation
clean:
	rm -fr multiplot_fractions multiplot_fractions_pub
# -------------------------------------------------------------------
