# Compiler selection (supported options are "gnu", "intel", and "clang")
COMPILER = "gnu"
# COMPILER = "intel"
# COMPILER = "clang"

# List of Gems that should be compiled
GEM = GRHayL_Core Con2Prim EOS Induction Reconstruction

# Set build directory name
BUILD_DIR = build

# HDF5 configuration
HDF5_DIR     = /usr/lib/x86_64-linux-gnu/hdf5/serial
HDF5_INC_DIR = $(HDF5_DIR)/include
HDF5_LIB_DIR = $(HDF5_DIR)/lib

# Get directory structures for Gems
DIR = $(shell find $(GEM) -type d)

# Find all source files
SRC = $(wildcard $(addsuffix /*.c, $(DIR)))

# Find all header files
INC = $(wildcard $(addsuffix /*.h, $(DIR)))

# Set all object files
OBJ = $(addprefix $(BUILD_DIR)/, $(SRC:.c=.o))

# Different compiler options
ifeq ($(COMPILER), "gnu")
	CC       = gcc
	CFLAGS   = -Wall -O2 -march=native -std=c99 -fopenmp
	LD_FLAGS = -lm
endif

ifeq ($(COMPILER), "intel")
	CC       = icc
	CFLAGS   = -Wall -O2 -march=native -std=c99 -fopenmp
	LD_FLAGS = -lm
endif

ifeq ($(COMPILER), "clang")
	CC       = clang
	CFLAGS   = -Wall -O2 -march=native -std=c99 -fopenmp
	LD_FLAGS = -lm
endif

# Now adjust CFLAGS and LD_FLAGS
CFLAGS   += -I./include -I$(HDF5_INC_DIR) -shared -fPIC
LD_FLAGS += -L$(HDF5_LIB_DIR) -lhdf5

BUILD_DIRS = lib/ $(addprefix $(BUILD_DIR)/, /$(addsuffix /, $(DIR)))

all: $(BUILD_DIRS) lib/libgrhayl.so

lib/libgrhayl.so: $(OBJ)
	@echo "Linking GRHayL object files"
	@$(CC) $(CFLAGS) $< -o $@ $(LD_FLAGS)
	@echo "All done!"

$(BUILD_DIRS):
	@echo "Creating directory $@"
	@mkdir -p $@

$(OBJ): build/%.o : %.c $(INC)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -c $< -o $@

clean:
	@rm -f $(OBJ) lib/libgrhayl.a
	@echo "Removing objects and library file"

veryclean: clean
	@rm -rf $(BUILD_DIRS)
	@echo "Removing build and lib directories"
