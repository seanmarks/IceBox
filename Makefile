# NOTES:
# -g	Extra symbolic debugging info for use with the gdb debugger
# -o	Specify output file name
# -c	Compile into object file
# -Wall	Print all warning messages
#
# $@	Name of the target
# $<	First dependency of the target
# $^	All dependencies of the target
#
# Linking occurs when the object files (.o) are put together into an executable (.exe)

INSTALL_DIR=$(HOME)/programs/IceBox

# Compilers
CC=gcc
CXX=g++

# Output binary's name
PROJECT=IceBox


### !!! Shouldn't have to change anything below this line !!! ###


#############
### Flags ###
#############

# Compiler flags
# - Misc: -march=native -Wno-comment -Wno-sign-compare -DPLUMED_MODE
CXXFLAGS += -g -std=c++11 -DCPLUSPLUS -O3 -Wall

# Linking
LDFLAGS += -g
LIBS    += -lm


#############################
### Directories and Files ###
#############################

START_DIR := $(PWD)

SRC_DIR         := src
BUILD_DIR       := build
BUILD_BIN_DIR   := $(BUILD_DIR)/bin
INSTALL_BIN_DIR := $(INSTALL_DIR)/bin

# Get source and target objects
SRC_FILES = $(shell find $(SRC_DIR) -name '*.cpp')
SRC_DIRS  = $(shell find $(SRC_DIR) -type d | sed 's/$(SRC_DIR)/./g' )
OBJECTS   = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_FILES))


############
### Make ###
############

all : buildrepo $(PROJECT)

.PHONY : install
install :
	@{ \
	if [[ ! -d $(INSTALL_BIN_DIR) ]]; then \
		echo "Creating $(INSTALL_BIN_DIR) ..." ;\
		mkdir -p $(INSTALL_BIN_DIR) ;\
	fi ;\
	echo "Installing at $(INSTALL_DIR) ..." ;\
	cp $(BUILD_BIN_DIR)/* $(INSTALL_BIN_DIR) ;\
	echo "Done" ;\
	}

# Link project
$(PROJECT) : $(OBJECTS)
	$(CXX) -o $(BUILD_BIN_DIR)/$@ $(LDFLAGS) $(LIBS) $(OBJECTS)

# Compile
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) -o $@ $(CXXFLAGS) -c $<

.PHONY : clean	
clean :
	rm -rf $(BUILD_DIR)/*

.PHONY : clean_install
clean_install :
	rm -rf $(INSTALL_DIR)/*

# Make directories 
buildrepo :
	@{ \
	for dir in $(BUILD_DIR) $(BUILD_BIN_DIR) $(INSTALL_DIR)/bin ; do \
		if [[ ! -d $$dir ]]; then \
			echo "Making directory $$dir ..." ;\
			mkdir -p $$dir ;\
		fi ;\
	done ;\
	}


###### !!!

# NOTE: old contents are below
#  buildrepo:
#    @$(call make-repo)

# Create obj directory structure
define make-repo
	mkdir -p $(BUILD_DIR)
	for dir in $(SRC_DIRS); \
	do \
		mkdir -p $(BUILD_DIR)/$$dir; \
	done
endef

