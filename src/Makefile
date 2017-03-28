CXX = g++
CC = gcc
RR = R CMD SHLIB
CFLAGS = -std=c11 -fPIC -I. -I/usr/include/R/
CXXFLAGS = -std=c++11 -fopenmp -fPIC -I. -I${mkEigenInc} -I/usr/include/R/ -fpermissive
RFLAGS = -fPIC -lstdc++ -fopenmp -Wl,-z,defs
LDFLAGS = -lm


CSOURCES = FEMtriangulation.c triangle.c

CXXSOURCES = integration.cpp mesh_objects.cpp FEMeval.cpp fdaPDE.cpp
	
TARGET_LIB = fdaPDE.so

COBJS = $(CSOURCES:.c=.o)
CXXOBJS = $(CXXSOURCES:.cpp=.o)

OBJS = $(COBJS) $(CXXOBJS)


$(TARGET_LIB) : $(OBJS) 
	$(RR) $(RFLAGS) -o $@ $^ $(LDFLAGS)

$(CXXOBJS) : %.o : %.cpp 
	$(CXX) $(CXXFLAGS) -c $< 
	
$(COBJS) : %.o : %.c 
	$(CC) $(CFLAGS) -c $< 

.PHONY: clean distclean

clean :
	$(RM) $(OBJS)

distclean : clean
	$(RM) $(TARGET_LIB)
