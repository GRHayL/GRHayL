# Compiler selection (supported options are "gnu", "intel", and "clang")
COMPILER = "gnu"
# COMPILER = "intel"
# COMPILER = "clang"

# Set name of the build directory
BUILD_DIR = build/

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

# HDF5 configuration
HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial
HDF5_INC_DIR=$(HDF5_DIR)/include
HDF5_LIB_DIR=$(HDF5_DIR)/lib

# Now adjust CFLAGS and LD_FLAGS
CFLAGS   += -I./include -I$(HDF5_INC_DIR) -shared -fPIC
LD_FLAGS += -L$(HDF5_LIB_DIR) -lhdf5

# Source files
SRC := $(wildcard GRHayL_Core/*.c Con2Prim/*.c Con2Prim/C2P_Routines/Font_Fix/*.c Con2Prim/C2P_Routines/Hybrid/*.c EOS/Hybrid/*.c EOS/Tabulated/*.c EOS/Tabulated/interpolators/*.c Induction/*.c Reconstruction/*.c)

# Object files
OBJ = $(patsubst %.c,$(BUILD_DIR)%.o,$(SRC))

# Header files
INC = ./EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h ./include/induction.h ./include/con2prim.h ./include/unit_tests.h ./include/GRHayL.h ./include/NRPyEOS_Hybrid.h ./include/NRPyEOS_Tabulated.h ./Con2Prim/C2P_Routines/harm_u2p_util.h

BUILD_DIRS = lib/ build/Con2Prim/ build/Con2Prim/C2P_Routines/Font_Fix/ build/Con2Prim/C2P_Routines/Hybrid/ build/EOS/Hybrid/ build/EOS/Tabulated/ build/EOS/Tabulated/interpolators/ build/GRHayL_Core/ build/Induction/ build/Reconstruction/ build/Unit_Tests/

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
	@rm -rf build/ lib/
	@echo "Removing build and lib directories"
