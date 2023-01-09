# generate_makefile.sh - Generates a Makefile for GRHayL
# Author: Leo Werneck
#
# Basic usage (may require manually editing the generated Makefile):
#
#  ./generate_makefile.sh
#
# Usage with HDF5 directory specification:
#
#  ./generate_makefile.sh /usr/local/opt/hdf5

if [[ $# -gt 1 ]]; then
    printf "(GRHayL) Error: correct usage is $0 [hdf5_dir (optional)]\n"
    exit 1
fi

mcd_locs=`find . | grep make.code.defn`

awk -v hdf5_dir=$1 '
BEGIN { print "(GRHayL) Beginning Makefile automatic generation..." }
/SRCS/,/[^\\] *$/ {

  if( match(FILENAME, /ET\//) )
    next

  if( FNR != pFNR ) {
    n = split(FILENAME, name, "/")
    dir = name[1]"/"
    for(i=2;i<n;i++)
      dir = dir""name[i]"/"
    pFNR = FNR
  }

  sub(/SRCS *= */, "", $0)
  sub(/\\/, "", $0)
  sub(/^ */, "", $0)
  sub(/ *$/, "", $0)

  n = split($0, current_srcs, " ")
  for(key in current_srcs) {
    src = src" "dir""current_srcs[key]
    sub(/\.c/, ".o", current_srcs[key])
    obj = obj" obj/"current_srcs[key]
  }
}
END {
  n = split(src, srcs, " ")
  m = split(obj, objs, " ")
  if( n != m ) {
    printf("Error: number of source files (%d) does not match number of object files (%d).\n", n, m)
    exit
  }
  printf("(GRHayL) Found %d source and object files.\n", n)
  printf("(GRHayL) Writing Makefile...\n")

  print "# This Makefile was automatically generated using the GRHayL script generate_makefile.sh" > "Makefile"
  print "# Compiler options" >> "Makefile"
  print "CC       = gcc" >> "Makefile"
  print "CFLAGS   = -Wall -O2 -march=native -std=c99 -fopenmp -I./include" >> "Makefile"
  print "LD_FLAGS = -lm\n" >> "Makefile"
  print "# HDF5 configuration" >> "Makefile"
  print "HDF5_DIR="hdf5_dir >> "Makefile"
  print "HDF5_INC_DIR=$(HDF5_DIR)/include\nHDF5_LIB_DIR=$(HDF5_DIR)/lib" >> "Makefile"
  print "CFLAGS += -I$(HDF5_INC_DIR)" >> "Makefile"
  print "LD_FLAGS += -L$(HDF5_LIB_DIR) -lhdf5\n" >> "Makefile"
  print "# Source files\nSRC="src"\n" >> "Makefile"
  print "# Object files\nOBJ="obj"\n" >> "Makefile"
  print "all: objdir $(OBJ)\n" >> "Makefile"
  print "objdir:\n\tmkdir -p obj\n" >> "Makefile"
  for(i=1;i<=n;i++) {
    print objs[i]": "srcs[i]" $(INC)" >> "Makefile"
    print "\t$(CC) $(CFLAGS) -c "srcs[i]" -o "objs[i]"\n" >> "Makefile"
  }
  print "clean:\n\trm -f $(OBJ)\n" >> "Makefile"
  print "veryclean: clean\n\trm -rf obj/" >> "Makefile"

  print "(GRHayL) Finished writing Makefile."
  if( length(hdf5_dir) > 0 )
    print "(GRHayL) Set HDF5 directory to "hdf5_dir
  else
    print "(GRHayL) Please do not forget to add the path to your HDF5 installation to the Makefile."
  print "(GRHayL) All done!"
}
' $mcd_locs
