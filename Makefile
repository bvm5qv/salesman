ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
CXXFLAGS  += $(ROOTCFLAGS) -I$(ODELIB) -Wall -O3
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS) -Wl,-rpath,$(ODELIB) -L$(ODELIB) -lode
GSLFLAGS   = -lgsl -lgslcblas
GXX	   = g++ $(CXXFLAGS)

EXES = datareader sales

all: $(EXES)

datareader: datareader.cpp
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)

sales: sales.cpp
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)

clean:
	rm -f datareader *~ *png
