# GRHayL Makefile
# Author(s): Leo Werneck and Sam Cupp

# Compiler selection (supported options are "gnu", "intel", and "clang")
COMPILER = "gnu"
# COMPILER = "intel"
# COMPILER = "clang"

# List of Gems that should be compiled
GEM = GRHayL_Core Con2Prim EOS Induction Reconstruction Neutrinos

# Set build directory name
BUILD_DIR = build

# HDF5 configuration
HDF5_DIR     = /usr/lib/x86_64-linux-gnu/hdf5/serial
HDF5_INC_DIR = $(HDF5_DIR)/include
HDF5_LIB_DIR = $(HDF5_DIR)/lib

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
	CFLAGS   = -Wall -O2 -march=native -std=c99 -fopenmp=libgomp
	LD_FLAGS = -lm
endif

# Now adjust CFLAGS and LD_FLAGS
CFLAGS   += -I./include -I$(HDF5_INC_DIR)
LD_FLAGS += -L$(HDF5_LIB_DIR) -lhdf5

# Get directory structures for Gems
DIR = $(shell find $(GEM) -type d)

# Find all source files
SRC = $(wildcard $(addsuffix /*.c, $(DIR)))

# Find all header files
INC = $(wildcard $(addsuffix /*.h, $(DIR)))

# Set all object files
OBJ = $(addprefix $(BUILD_DIR)/, $(SRC:.c=.o))

# Get all sources and headers in the unit test directory
UNIT_TEST_DIR = Unit_Tests
UNIT_TEST_SRC = $(wildcard $(UNIT_TEST_DIR)/*.c)
UNIT_TEST_OBJ = $(addprefix $(BUILD_DIR)/, $(UNIT_TEST_SRC:.c=.o))
UNIT_TEST_INC = $(wildcard $(UNIT_TEST_DIR)/*.h)

# Now we must separate the tests from auxiliary files
UNIT_TEST_MAIN_SRC = $(wildcard $(UNIT_TEST_DIR)/unit_test_*.c)
UNIT_TEST_AUXS_SRC = $(filter-out $(UNIT_TEST_MAIN_SRC), $(UNIT_TEST_SRC))

# Finally, get the object files
UNIT_TEST_MAIN_OBJ = $(addprefix $(BUILD_DIR)/, $(UNIT_TEST_MAIN_SRC:.c=.o))
UNIT_TEST_AUXS_OBJ = $(addprefix $(BUILD_DIR)/, $(UNIT_TEST_AUXS_SRC:.c=.o))

# Every MAIN src/obj will result in a unit test executable
UNIT_TEST_MAIN_EXE = $(addprefix exe/, $(notdir $(UNIT_TEST_MAIN_SRC:.c=)))

# Build directories (these are created if they do not already exist)
BUILD_DIRS = lib/ $(addprefix $(BUILD_DIR)/, $(addsuffix /, $(DIR))) exe/ $(addprefix $(BUILD_DIR)/, $(UNIT_TEST_DIR))

# .---------------.
# | Build Targets |
# .---------------.

# Build all of GRHayL (standard library + unit tests)
all: $(BUILD_DIRS) lib/libgrhayl.so $(UNIT_TEST_MAIN_EXE)
	@echo "All done!"

# Compile the GRHayL standard library
lib/libgrhayl.so: $(OBJ)
	@echo "Linking GRHayL object files"
	@$(CC) $(CFLAGS) -shared -fPIC $(OBJ) -o $@ $(LD_FLAGS)

# Create all directories in BUILD_DIRS if they do not already exist
$(BUILD_DIRS):
	@echo "Creating directory $@"
	@mkdir -p $@

# Compile all object files for all gems in GEM
$(OBJ): $(BUILD_DIR)/%.o : %.c $(INC)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -fPIC -c $< -o $@

# Compile all auxiliary object files for the unit tests
$(UNIT_TEST_AUXS_OBJ): $(BUILD_DIR)/%.o : %.c $(INC)  $(UNIT_TEST_INC)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -I$(UNIT_TEST_DIR) -c $< -o $@

# Compile all main object files for the unit tests
$(UNIT_TEST_MAIN_OBJ): $(BUILD_DIR)/%.o : %.c $(INC)  $(UNIT_TEST_INC)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -I$(UNIT_TEST_DIR) -c $< -o $@

# Create the unit test executables
$(UNIT_TEST_MAIN_EXE): exe/% : $(BUILD_DIR)/$(UNIT_TEST_DIR)/%.o $(UNIT_TEST_AUXS_OBJ) $(INC) $(UNIT_TEST_INC) lib/libgrhayl.so
	@echo "Creating unit test $@"
	@$(CC) $(CFLAGS) -I$(UNIT_TEST_DIR) $(UNIT_TEST_AUXS_OBJ) $< -o $@ $(LD_FLAGS) -L./lib -lgrhayl

# Clean: remove libraries, all objects, and all executables
clean:
	@rm -f $(OBJ) lib/libgrhayl.a $(UNIT_TEST_OBJ) $(UNIT_TEST_MAIN_EXE)
	@echo "Removing objects and library file"

# Thorough clean: remove all files in clean and all build, library, and executable directories
veryclean: clean
	@rm -rf $(BUILD_DIRS) lib/ exe/
	@echo "Removing build directories"
