ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)
INCLUDES    := -I../include

CXX       := g++
#CXXFLAGS  += -std=c++11 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS)
CXXFLAGS  += -std=c++11 -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

OBJECTS2   := Subtraction.o run_genie_analysis.o Fiducial.o genie_analysis.o
OBJECTS3   := Subtraction.o run_systematics.o Fiducial.o systematics.o
OBJECTS4   := gst.o run_gst.o

all: genie_analysis systematics

genie_analysis: $(OBJECTS2)
		$(CXX) -o genie_analysis $(OBJECTS2) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

systematics: $(OBJECTS3)
		$(CXX) -o systematics $(OBJECTS3) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)
genie_true: $(OBJECTS4)
		$(CXX) -o genie_true $(OBJECTS4) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o run_genie_analysis run_systematics *~

%.o: %.C
	$(CXX) -c $< -O2 $(CXXFLAGS) $(INCLUDES)
