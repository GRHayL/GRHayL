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

mcd_locs=`find . | grep make.code.defn | sort`
inc_locs=`find . | egrep '\.h$'`

awk -v hdf5_dir=$1 -v inc="`echo $inc_locs`" '

FNR==1 { fileid++ }

BEGIN { print "(GRHayL) Beginning Makefile automatic generation..."}
FILENAME !~ /ET\// && /SRCS/,/[^\\] *$/ {

  if( fileid != prev_fileid ) {
    n = split(FILENAME, name, "/")
    srcdir = name[1]"/"
    for(i=2;i<n;i++)
      srcdir = srcdir""name[i]"/"
    prev_fileid = fileid
    objdir = srcdir
    sub(/\.*\//, "", objdir)
    src_dirs[++nsrcdirs] = srcdir
    if( srcdir != "./" && srcdir != "../" )
      obj_dirs[++nobjdirs] = "build/"objdir
  }

  sub(/SRCS *= */, "", $0)
  sub(/\\/, "", $0)
  sub(/^ */, "", $0)
  sub(/ *$/, "", $0)

  n = split($0, current_srcs, " ")
  for(key in current_srcs) {
    src = src" "srcdir""current_srcs[key]
    sub(/\.c/, ".o", current_srcs[key])
    src_filecount[srcdir]++

    objdir = srcdir
    sub(/\.\//, "", objdir)
    obj = obj" build/"objdir""current_srcs[key]
  }
}
END {
  nsrc = split(src, srcs, " ")
  nobj = split(obj, objs, " ")
  if( nsrc != nobj ) {
    printf("Error: number of source files (%d) does not match number of object files (%d).\n", nsrc, nobj)
    exit
  }

  # Find the include directories
  n = split(inc_files, incs, " ")
  for(i=1;i<=n;i++) {
    m = split(incs[i], dirs, "/")
    dir = dirs[1]
    for(j=2;j<m;j++)
      dir = dir"/"dirs[j]
    inc_filecount[dir]++
    print nincdirs, inc_dirs[nincdirs], dir, inc_filecount[dir]
    if( dir != inc_dirs[nincdirs] )
      inc_dirs[++nincdirs] = dir
  }

  bubble_sort(nincdirs, inc_dirs)
  bubble_sort(nsrcdirs, src_dirs)
  bubble_sort(nobjdirs, obj_dirs)

  printf("(GRHayL) File summary:\n", nsrc)
  printf("(GRHayL)   Source files:\n")
  for(i=1;i<=nsrcdirs;i++)
    if(src_filecount[src_dirs[i]]>0)
      printf("(GRHayL)     Found %02d file(s) in directory %s\n", src_filecount[src_dirs[i]], src_dirs[i])
  printf("(GRHayL)   Header files:\n")
   for(i=1;i<=nincdirs;i++)
    printf("(GRHayL)     Found %02d file(s) in directory %s\n", inc_filecount[inc_dirs[i]], inc_dirs[i])

  printf("(GRHayL) Writing Makefile...\n")
  print "# This Makefile was automatically generated using the GRHayL script generate_makefile.sh" > "Makefile"
  print "# Compiler selection (supported options are \"gnu\", \"intel\", and \"clang\"" >> "Makefile"
  print "COMPILER = \"gnu\"" >> "Makefile"
  print "# COMPILER = \"intel\"" >> "Makefile"
  print "# COMPILER = \"clang\"" >> "Makefile"

  print "# GNU C compiler" >> "Makefile"
  print "GNU_CC       = gcc" >> "Makefile"
  print "GNU_CFLAGS   = -Wall -O2 -march=native -std=c99 -fopenmp" >> "Makefile"
  print "GNU_LD_FLAGS = -lm\n" >> "Makefile"

  print "# Intel C compiler" >> "Makefile"
  print "INTEL_CC       = icc" >> "Makefile"
  print "INTEL_CFLAGS   = -Wall -O2 -march=native -std=c99 -qopenmp" >> "Makefile"
  print "INTEL_LD_FLAGS = -lm\n" >> "Makefile"

  print "# Clang C compiler" >> "Makefile"
  print "CLANG_CC       = clang" >> "Makefile"
  print "CLANG_CFLAGS   = -Wall -O2 -march=native -std=c99 -fopenmp" >> "Makefile"
  print "CLANG_LD_FLAGS = -lm\n" >> "Makefile"

  print "ifeq ($(COMPILER), \"gnu\")" >> "Makefile"
  print "\tCC=$(GNU_CC)\n\tCFLAGS=$(GNU_CFLAGS)\n\tLD_FLAGS=$(GNU_LD_FLAGS)" >> "Makefile"
  print "endif\n" >> "Makefile"

  print "ifeq ($(COMPILER), \"intel\")" >> "Makefile"
  print "\tCC=$(INTEL_CC)\n\tCFLAGS=$(INTEL_CFLAGS)\n\tLD_FLAGS=$(INTEL_LD_FLAGS)" >> "Makefile"
  print "endif\n" >> "Makefile"

  print "ifeq ($(COMPILER), \"clang\")" >> "Makefile"
  print "\tCC=$(CLANG_CC)\n\tCFLAGS=$(CLANG_CFLAGS)\n\tLD_FLAGS=$(CLANG_LD_FLAGS)" >> "Makefile"
  print "endif\n" >> "Makefile"

  print "# HDF5 configuration" >> "Makefile"
  print "HDF5_DIR="hdf5_dir >> "Makefile"
  print "HDF5_INC_DIR=$(HDF5_DIR)/include\nHDF5_LIB_DIR=$(HDF5_DIR)/lib\n" >> "Makefile"
  print "# Now adjust CFLAGS and LD_FLAGS" >> "Makefile"
  print "CFLAGS += -I./include -I$(HDF5_INC_DIR)" >> "Makefile"
  print "LD_FLAGS += -L$(HDF5_LIB_DIR) -lhdf5\n" >> "Makefile"
  print "# Source files\nSRC ="src"\n" >> "Makefile"
  print "# Object files\nOBJ ="obj"\n" >> "Makefile"
  print "# Header files\nINC = "inc"\n" >> "Makefile"
  printf("all: libdir") >> "Makefile"
  for(i=1;i<=nobjdirs;i++)
    printf(" %s", obj_dirs[i]) >> "Makefile"
  print " libgrhayl\n" >> "Makefile"
  print "libgrhayl: $(OBJ)\n\t@echo \"Linking GRHayL object files\"\n\t@ar rcs lib/libgrhayl.a $(OBJ)\n\t@echo \"All done!\"\n" >> "Makefile"
  print "libdir:\n\t@mkdir -p lib\n\t@echo \"Creating lib directory\"\n" >> "Makefile"
  for(i=1;i<=nobjdirs;i++)
    print obj_dirs[i]":\n\t@mkdir -p "obj_dirs[i]"\n\t@echo \"Creating build directory "obj_dirs[i]"\"\n" >> "Makefile"
  for(i=1;i<=nsrc;i++) {
    print objs[i]": "srcs[i]" $(INC)" >> "Makefile"
    print "\t@$(CC) $(CFLAGS) -c "srcs[i]" -o "objs[i] >> "Makefile"
    print "\t@echo \"Compiling source file "srcs[i]"\"\n" >> "Makefile"
  }
  print "clean:\n\t@rm -f $(OBJ) lib/libgrhayl.a\n\t@echo \"Removing objects and library file\"\n" >> "Makefile"
  print "veryclean: clean\n\t@rm -rf build/ lib/\n\t@echo \"Removing build and lib directories\"" >> "Makefile"

  print "(GRHayL) Finished writing Makefile."
  if( length(hdf5_dir) > 0 )
    print "(GRHayL) Set HDF5 directory to "hdf5_dir
  else
    print "(GRHayL) Please do not forget to add the path to your HDF5 installation to the Makefile."
  print "(GRHayL) All done!"
}

function bubble_sort(n, y) {
  for(i=1;i<n;i++)
    for(j=i+1;j<=n;j++)
      if( y[i] > y[j] ) {
        tmp = y[i]
        y[i] = y[j]
        y[j] = tmp
      }
}
' $mcd_locs
