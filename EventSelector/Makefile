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

OBJECTS1   := event_selector.o run_event_selector.o

all: event_selector

event_selector: $(OBJECTS1)
		$(CXX) -o event_selector $(OBJECTS1) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o

%.o: %.C
	$(CXX) -c $< -O2 $(CXXFLAGS) $(INCLUDES)
