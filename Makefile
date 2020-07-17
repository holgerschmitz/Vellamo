
TARGET_BASE = vellamo
OFLAGS  = -g -O0 -Wall
#OFLAGS  = -O3 -Wall


DIMENSIONS = 1 2 3

SOURCES = $(wildcard src/*.cpp) \
  huerto/hydrodynamics/hydro_fields.cpp\
  huerto/hydrodynamics/euler/adiabatic_knp.cpp


# Set this to wherever your HDF5 installation resides
# Make sure that this is compatible with your compiler below.
HDFBASE = /usr/lib/x86_64-linux-gnu/hdf5/openmpi

# Your C++ compiler or MPI C++ compiler wrapper
# CXX     = g++
CXX		 = mpic++
LINK	 = mpic++
CXXFLAGS = $(OFLAGS)

INCLUDE = -I/usr/local/include -I$(HDFBASE)/include

BUILD_DIR = build
BIN_DIR = bin
OBJECTS = $(addprefix $(BUILD_DIR)/,$(patsubst %.cpp,%.o,$(SOURCES)))

LDFLAGS = -L$(HDFBASE)/lib -Wl,-rpath,$(HDFBASE)/lib

LOADLIBS = -lhdf5 -lschnek -lm

# FULLTARGET = $(BIN_DIR)/$(TARGET)

DIM1_FLAGS = -DHUERTO_ONE_DIM
DIM2_FLAGS = -DHUERTO_TWO_DIM
DIM3_FLAGS = -DHUERTO_THREE_DIM

define PROGRAM_template =
 TARGET$(1)D_OBJS = $(addprefix $(BUILD_DIR)/$(1)d/,$(patsubst %.cpp,%.o,$(SOURCES)))
 $(BIN_DIR)/$(TARGET_BASE)$(1)d: $$(TARGET$(1)D_OBJS)
	$(LINK) $^ -o $@ $(OFLAGS) $(LDFLAGS) $(LOADLIBS)
 $$(TARGET$(1)D_OBJS): $(BUILD_DIR)/$(1)d/%.o: %.cpp
	$(CXX) -o $$@ -c $(CXXFLAGS) $(INCLUDE) $$(DIM$(1)_FLAGS) $$<
 FULLTARGET += $(BIN_DIR)/$(TARGET_BASE)$(1)d
 ALL_OBJS   += $$(TARGET_BASE$(1)D_OBJS)
endef


$(foreach dimension,$(DIMENSIONS),$(eval $(call PROGRAM_template,$(dimension))))

all: $(FULLTARGET)



clean:
	-rm -f $(ALL_OBJS) core $(FULLTARGET)


