
DIMENSION=2

OFLAGS  = -g -O0 -Wall
#OFLAGS  = -O3 -Wall

# Set this to wherever your HDF5 installation resides
# Make sure that this is compatible with your compiler below.
HDFBASE = /usr/lib/x86_64-linux-gnu/hdf5/mpich

# Your C++ compiler or MPI C++ compiler wrapper
# CXX     = g++
CXX     = mpic++

CXXFLAGS = $(OFLAGS)

INCLUDE = -I/usr/local/include -I$(HDFBASE)/include

TARGET=vellamo$(DIMENSION)d

SOURCES = src/boundary.cpp \
  src/diagnostic.cpp \
  src/euler_solver.cpp \
  src/hydro_fields.cpp \
  src/vellamo.cpp


OBJECTS = $(SOURCES:.cpp=.o)

LDFLAGS = -L$(HDFBASE)/lib -Wl,-rpath,$(HDFBASE)/lib

LOADLIBS = -lhdf5 -lschnek -lm
BINDIR = bin
OBJDIR = obj

FULLTARGET = $(BINDIR)/$(TARGET)

all: $(FULLTARGET)

$(FULLTARGET): $(OBJECTS) 
	@mkdir -p $(BINDIR)
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)


%.o: %.cpp
	$(CXX) -DENV_DIMENSION=$(DIMENSION) -o $@ -c $(CXXFLAGS) $(INCLUDE) $<


clean:
	-rm -f $(OBJECTS) core $(FULLTARGET)


