# GRHayL Makefile

# List of Gems that should be compiled
GEM = GRHayL_Core Con2Prim EOS Induction Reconstruction Neutrinos Flux_Source

# Set build directory name
BUILD_DIR = build

# Compiler selection (supported options are gnu, intel, and clang)
ifeq ($(COMPILER),)
  COMPILER = gnu
endif

UNAME_S := $(shell uname -s)

# Different compiler options
CFLAGS  += -Wall -O2 -march=native -std=c99
LDFLAGS += -lm

ifeq ($(COMPILER), gnu)
  ifeq ($(UNAME_S),Linux)
    CC = gcc
  endif
  ifeq ($(UNAME_S),Darwin)
    ifeq ($(CC),cc)
      CC := $(shell find `brew --prefix`/opt/gcc/bin | egrep 'opt/gcc/bin/gcc-[0-9]{2}$$')
    endif
    ifeq ($(CC),gcc)
      CC := $(shell find `brew --prefix`/opt/gcc/bin | egrep 'opt/gcc/bin/gcc-[0-9]{2}$$')
    endif
  endif
  CFLAGS += -fopenmp
endif
ifeq ($(COMPILER), intel)
  CC       = icx
  CFLAGS += -fopenmp
endif
ifeq ($(COMPILER), clang)
  ifeq ($(UNAME_S),Linux)
    CC      = clang-14
    CFLAGS += -fopenmp=libgomp
  endif
  ifeq ($(UNAME_S),Darwin)
    ifeq ($(CC),cc)
      CC := $(shell find `brew --prefix`/opt/llvm/bin | egrep 'opt/llvm/bin/clang-[0-9]{2}$$')
    endif
    ifeq ($(CC),gcc)
      CC := $(shell find `brew --prefix`/opt/llvm/bin | egrep 'opt/llvm/bin/clang-[0-9]{2}$$')
    endif
    CFLAGS += -fopenmp=libomp
  endif
endif

# HDF5 configuration
ifeq ($(HDF5_INC_DIR),)
  ifeq ($(HDF5_DIR),)
    ifeq ($(UNAME_S),Linux)
      HDF5_DIR = /usr/lib/x86_64-linux-gnu/hdf5/serial
    else
      ifeq ($(UNAME_S),Darwin)
        HDF5_DIR := `brew --prefix`/opt/hdf5
      endif
    endif
  endif
  HDF5_INC_DIR = $(HDF5_DIR)/include
endif
ifeq ($(HDF5_LIB_DIR),)
  HDF5_LIB_DIR = $(HDF5_DIR)/lib
endif

# Now adjust CFLAGS and LDFLAGS
CFLAGS  += -I./GRHayL/include -I$(HDF5_INC_DIR)
LDFLAGS += -L$(HDF5_LIB_DIR) -lhdf5

# Get directory structures for Gems
DIR := $(shell find $(addprefix GRHayL/, $(GEM) include) -type d)

# Find all source files
SRC := $(wildcard $(addsuffix /*.c, $(DIR)))

# Find all header files
INC := $(wildcard $(addsuffix /*.h, $(DIR)))

# Set all object files
OBJ := $(addprefix $(BUILD_DIR)/, $(SRC:.c=.o))

# Get all sources and headers in the unit test directory
UNIT_TEST_DIR = Unit_Tests
UNIT_TEST_SRC := $(wildcard $(UNIT_TEST_DIR)/*.c)
UNIT_TEST_OBJ := $(addprefix $(BUILD_DIR)/, $(UNIT_TEST_SRC:.c=.o))
UNIT_TEST_INC := $(wildcard $(UNIT_TEST_DIR)/*.h)

# Now we must separate the tests from auxiliary files
UNIT_TEST_MAIN_SRC := $(wildcard $(UNIT_TEST_DIR)/unit_test_*.c)
UNIT_TEST_AUXS_SRC := $(filter-out $(UNIT_TEST_MAIN_SRC), $(UNIT_TEST_SRC))

# Finally, get the object files
UNIT_TEST_MAIN_OBJ := $(addprefix $(BUILD_DIR)/, $(UNIT_TEST_MAIN_SRC:.c=.o))
UNIT_TEST_AUXS_OBJ := $(addprefix $(BUILD_DIR)/, $(UNIT_TEST_AUXS_SRC:.c=.o))

# Every MAIN src/obj will result in a unit test executable
UNIT_TEST_MAIN_EXE := $(addprefix exe/, $(notdir $(UNIT_TEST_MAIN_SRC:.c=)))

# Build directories (these are created if they do not already exist)
BUILD_DIRS := lib/ $(addprefix $(BUILD_DIR)/, $(addsuffix /, $(DIR))) exe/ $(addprefix $(BUILD_DIR)/, $(UNIT_TEST_DIR))

# .---------------.
# | Build Targets |
# .---------------.

# Build all of GRHayL (standard library + unit tests)
all: $(BUILD_DIRS) lib/libgrhayl.so $(UNIT_TEST_MAIN_EXE)
	@echo "All done!"

# Compile the GRHayL standard library
lib/libgrhayl.so: $(OBJ)
	@echo "Linking GRHayL object files"
	@$(CC) $(CFLAGS) -shared -fPIC $(OBJ) -o $@ $(LDFLAGS)

# Create all directories in BUILD_DIRS if they do not already exist
$(BUILD_DIRS):
	@echo "Creating directory $@"
	@mkdir -p $@

# Compile all object files for all gems in GEM
$(OBJ): $(BUILD_DIR)/%.o : %.c $(INC)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -fPIC -c $< -o $@

# Compile all auxiliary object files for the unit tests
$(UNIT_TEST_AUXS_OBJ): $(BUILD_DIR)/%.o : %.c $(INC) $(UNIT_TEST_INC)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -I$(UNIT_TEST_DIR) -c $< -o $@

# Compile all main object files for the unit tests
$(UNIT_TEST_MAIN_OBJ): $(BUILD_DIR)/%.o : %.c $(INC) $(UNIT_TEST_INC)
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -I$(UNIT_TEST_DIR) -c $< -o $@

# Create the unit test executables
$(UNIT_TEST_MAIN_EXE): exe/% : $(BUILD_DIR)/$(UNIT_TEST_DIR)/%.o $(UNIT_TEST_AUXS_OBJ) $(INC) $(UNIT_TEST_INC) lib/libgrhayl.so
	@echo "Creating unit test $@"
	@$(CC) $(CFLAGS) -I$(UNIT_TEST_DIR) $(UNIT_TEST_AUXS_OBJ) $< -o $@ $(LDFLAGS) -L./lib -lgrhayl

# Clean: remove libraries, all objects, and all executables
clean:
	@rm -f $(OBJ) lib/libgrhayl.a $(UNIT_TEST_OBJ) $(UNIT_TEST_MAIN_EXE)
	@echo "Removing objects and library file"

# Thorough clean: remove all files in clean and all build, library, and executable directories
realclean: clean
	@rm -rf $(BUILD_DIRS) build/ lib/ exe/
	@echo "Removing build directories"
