######################################################################
# Project: TMVA - Toolkit for Multivariate Data Analysis             #
# Code   : Example usage for TMVA analyses                           #
###################################################################### 

MAKEFLAGS = --no-print-directory -r -s

INCLUDE = -I$(TMVASYS)/ $(shell root-config --cflags)
LIBS    = -L$(TMVASYS)/lib -lTMVA.1 $(shell root-config --libs) -lMLP -lTreePlayer -lMinuit

BINS = TMVAClassification \
	TMVAClassificationCategory \
	TMVAClassificationApplication \
	TMVAClassificationCategoryApplication \
	TMVARegression \
	TMVARegressionApplication \
	TMVAMulticlass \
	TMVAMulticlassApplication \
	TMVAMultipleBackgroundExample

TMVALIB = $(TMVASYS)/lib/libTMVA.1.so

PUBLISH = $(TMVASYS)/www


all: $(TMVALIB) $(BINS)

run: TMVA.root 

apply: TMVApp.root


clean:
	rm -f $(BINS)

$(TMVALIB): $(TMVASYS)/src/*.cxx $(TMVASYS)/inc/*.h
	$(MAKE) -C ..

install:
	@mkdir -p $(PUBLISH)
	@cp -u $(TMVASYS)/test/pad.inc $(PUBLISH)
	@cp -u $(TMVASYS)/test/tabView.js $(PUBLISH)
	@cp -u $(TMVASYS)/test/index.php $(PUBLISH)
	@cp -u $(TMVASYS)/macros/tmva_logo.gif $(PUBLISH)


TMVA.root: $(TMVALIB) TMVAClassification 
	./TMVAClassification -b $(TESTCLASSIFIERS)

TMVApp.root: TMVAClassificationApplication TMVA.root
	./TMVAClassificationApplication $(TESTCLASSIFIERS)


$(BINS): % : %.cxx $(TMVASYS)/inc/*.h
	@echo -n "Building $@ ... "
	$(CXX) $(CCFLAGS) $< $(INCLUDE) $(LIBS) -g -o $@
	@echo "Done"

.PHONY : all run apply web clean install ut

