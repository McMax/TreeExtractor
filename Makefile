CC	= g++-4.9
LD	= g++-4.9

CCFLAGS = -O2 `root-config --cflags` -Wall -I./inc -I$(PEV_INC)
LDFLAGS = -O2 `root-config --libs` -Wall -L./lib -I$(PEV_LIB)

TOPDIR = .
SRC_DIR = $(TOPDIR)/src
OBJ_DIR = $(TOPDIR)/lib
INC_DIR = $(TOPDIR)/inc

PEV_DIR = $(TOPDIR)/../Particle_Event_Clusters
PEV_LIB = $(PEV_DIR)/lib
PEV_SRC = $(PEV_DIR)/src
PEV_INC = $(PEV_DIR)/inc

SOURCES := $(shell find $(SRC_DIR) -type f -name "*.cpp")
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))
PEV_OBJECTS = $(PEV_LIB)/Particle.o $(PEV_LIB)/Event.o $(PEV_LIB)/ParticleTree.o $(PEV_LIB)/Dict.o

all: extractor merger

extractor: $(OBJ_DIR)/Extractor.o $(OBJ_DIR)/Prefifi.o $(OBJ_DIR)/RootWriter.o $(PEV_OBJECTS)
	$(LD) -o $@ $^ $(LDFLAGS)

merger: $(OBJ_DIR)/Merger.o $(PEV_OBJECTS)
	$(LD) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(OBJ_DIR)
	$(CC) -c $(CCFLAGS) $< -o $@ 

$(PEV_OBJECTS):
	@echo "No base libs. Create them"

clean:
	@rm -rf merger extractor ./lib
