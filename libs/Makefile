#
# mock_galaxy_library/lib
#

LIBRARY   = libmockgal.a libmockgal.dylib

all: $(LIBRARY)

CPPFLAGS := -O2
LIBS     := -lm -lgsl -lgslcblas -lhdf5

# add non standark library path if necessary
DIR_PATH := # e.g. /opt/gsl

CFLAGS   += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS     += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)
LIBS     += 

OBJS     := msg.o power.o mf.o

libmockgal.a: $(OBJS)
	ar r $@ $(OBJS)

libmockgal.dylib: $(OBJS)
	$(CXX) -dynamiclib $(LIBS) $(OBJS) -o $@

.PHONY: clean dependence
clean:
	rm -f $(LIBRARY) $(OBJS)

dependence:
	g++ -MM -MG *.cpp
