
TARGET=vellamo

#OFLAGS  = -g -O0 -Wall
OFLAGS  = -O3 -Wall

INCLUDE = -I/usr/local/include
#CXX     = $(X_CXX)
CXX     = mpic++

CXXFLAGS = $(OFLAGS)

SOURCES = src/boundary.cpp \
  src/diagnostic.cpp \
  src/euler_solver.cpp \
  src/vellamo.cpp


OBJECTS = $(SOURCES:.cpp=.o)

LDFLAGS = 

LOADLIBS = -lhdf5 -lschnek -lm
BINDIR = bin
OBJDIR = obj

FULLTARGET = $(BINDIR)/$(TARGET)

all: $(FULLTARGET)

$(FULLTARGET): $(OBJECTS) 
	@mkdir -p $(BINDIR)
	$(CXX) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)


%.o: %.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) $(INCLUDE) $<


clean:
	-rm -f $(OBJECTS) core $(FULLTARGET)


