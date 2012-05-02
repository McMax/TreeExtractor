CC	= g++
LD	= g++

CCFLAGS = -g -O0 `root-config --cflags` -lEG -Wall -I./inc
LDFLAGS = -g -O0 `root-config --libs` -lEG -Wall -L./lib

TOPDIR = .
SRC_DIR = $(TOPDIR)/src
OBJ_DIR = $(TOPDIR)/lib
INC_DIR = $(TOPDIR)/inc

PROGRAM = extractor

SOURCES := $(shell find $(SRC_DIR) -type f -name "*.cpp" | grep -v "Prefifi.cpp")
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))
#INCLUDES := $(shell find $(INC_DIR) -type f -name "*.h" | grep -v "linkdef.h" | grep -v "Prefifi.h")
#INCLUDES += $(INC_DIR)/linkdef.h
INCLUDES = inc/Particle.h inc/Event.h inc/linkdef.h

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS) Dict.o
	$(LD) $(LDFLAGS) $^ -o $@ 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp Dict.cpp
	mkdir -p $(OBJ_DIR)
	$(CC) -c $(CCFLAGS) $< -o $@ 

Dict.cpp: $(INCLUDES)
	@echo "Generating dictionary..."
	rootcint -f Dict.cpp -c -P -I$(ROOTSYS) -I/usr/local/include $(INCLUDES)
	g++ -c $(CCFLAGS) Dict.cpp -o Dict.o	

clean:
	@rm $(PROGRAM) -r ./lib

realclean:
	@rm $(PROGRAM)
	@rm -r ./lib
	@rm Dict.*
