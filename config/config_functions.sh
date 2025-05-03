#!/bin/sh

normalize_spaces() {
      echo "$1" | awk '{$1=$1; print}'
}

append() {
  printf "%s\n" "$1" >> Makefile
}

tab_append() {
  printf "\t%s\n" "$1" >> Makefile
}

gen_rule() {
  target="$1"
  deps="$2"
  recipe="$3"
  echo "$target: $deps" >> Makefile
  echo "  $recipe" >> Makefile
  echo "" >> Makefile
}

red()
{
    printf "\033[1m\033[31m$1\033[0m"
}

green()
{
    printf "\033[1m\033[32m$1\033[0m"
}

blue()
{
    printf "\033[1m\033[34m$1\033[0m"
}

YES()
{
    green "YES"
    if [ $# -eq 1 ]; then
        printf " ($1)\n"
    else
        printf "\n"
    fi
}

NO()
{
    red "NO"
    if [ $# -eq 1 ]; then
        printf " ($1)\n"
    else
        printf "\n"
    fi
}

error()
{
    red "ERROR: "
    printf "$1\n"
    exit 1
}

check_command()
{
    printf "Checking if command $1 exists... "
    tmp_=$(command -v $1)
    if [ "$tmp_" = "" ]; then
        NO
        eval have_$(printf $1 | tr - _)=0
        return 0
    else
        YES $tmp_
        eval have_$(printf $1 | tr - _)=1
        return 1
    fi
}

check_compiler()
{
    printf "Checking if compiler $1 exists... "
    tmp_=$(command -v $1)
    if [ "$tmp_" = "" ]; then
        NO
    else
        CC=$1
        YES "$tmp_"
    fi
}

c_compiles()
{
    printf "$5"
    printf "$1" | $2 $3 -x c - $4 >/dev/null 2>&1
    if [ $? -eq 0 ]; then
        YES "$2$3<objs>$4"
        rm -f a.out
        return 0
    else
        NO "$2$3<objs>$4"
        rm -f a.out
        return 1
    fi
}

Realpath()
{
    echo $(cd $(dirname $1); pwd)/$(basename $1) | \
        awk '{gsub(/\/\.\./, "", $0); gsub(/\/\./, "", $0); print}'
}

NL='
'
probe_file()
{
    # Source: https://serverfault.com/questions/225798/
    tmp_=$(find $1 -name "$2" 2>/dev/null | grep .)
    return $?
}

set_hdf5()
{
    for dir in $1 $2; do
        if [ ! -d $dir ]; then
            error "Directory $dir not found"
        fi
    done
    printf "Checking if hdf5.h can be found in $1... "
    probe_file $1/ hdf5.h
    err=$?
    if [ $err -ne 0 ]; then
        case $err in
            1) NO;;
            2) printf "\nFound multiple hdf5.h files:\n"
               printf "$tmp_\n"
               printf "Please specify the desired one using the --hdf5dir, --hdf5inc,\n"
               printf "and/or --hdf5lib options. Run ./configure -h for further details.\n";;
        esac
        return $err
    else
        YES "$tmp_"
        cflags="$cflags -I$(dirname $tmp_)"
    fi
    for f in $hdf5libfiles; do
        printf "Checking if $f can be found in $2... "
        probe_file $2/ $f
        if [ $? -ne 0 ]; then
            NO;
            if [ "$f" = "libhdf5.a" ]; then
                return 1;
            fi
        else
            YES "$tmp_"
            clibs="$clibs -L$(dirname $tmp_) -lhdf5"
            break
        fi
    done
    return 0
}

parse_mcd() {
    local indent="$1"
    local path="$2"
    local depth="$3"
    local name=$(basename "$path")
    case "$depth" in
        0) printf "$indent- %s\n" $(blue "$name") ;;
        1) mkdir -p "$builddir/${path#GRHayL/}"
           printf "$indent- %s\n" $(green $name) ;;
        *) mkdir -p "$builddir/${path#GRHayL/}"
           printf "$indent- $name\n" ;;
    esac
    file="$2/make.code.defn"
    if [ ! -f $file ]; then
        error "File $file not found"
    fi
    for sd in $(./scripts/parser awk $file "subdirs"); do
        local sdpath="$path/$sd"
        sdsrcs=$(./scripts/parser awk "$sdpath/make.code.defn" "sources")
        sdincs=$(./scripts/parser awk "$sdpath/make.code.defn" "headers")
        sdihds=$(./scripts/parser awk "$sdpath/make.code.defn" "install_headers")
        srcs="$srcs $sdsrcs"
        incs="$incs $sdincs"
        ihds="$ihds $sdihds"
        parse_mcd "$indent  " $sdpath $(($depth+1))
    done
}
