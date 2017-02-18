#!/bin/bash


# NAME=`echo "$FILE" | cut -d'.' -f1`
# EXTENSION=`echo "$FILE" | cut -d'.' -f2`

SUFF="x"
include_flag="-I./include/"
source_dir="source"
program_dir="programs"
compiler="mpiCC ${include_flag}"

rm -f ./obj/*.o
rm -f ./lib/main/*.o
echo "==========================================================="
echo " Compiling script for HST Code "
echo "==========================================================="
echo 
echo "Compiling sources..."
mkdir -p bin
mkdir -p obj
rm -f bin/*.$SUFF
for src_name in `ls $source_dir/*.cpp | cut -d"/" -f2 |  cut -d"." -f1`; do
	src="$source_dir/$src_name.cpp"
  	obj="obj/$src_file.o"

	echo "compiling $src"
	$compiler -c $src 
	mv $src_name.o obj/ 
done


#================================================================
# compiling codes 
#================================================================
echo
echo "Compiling programs..."

for prog_name in `ls $program_dir/*.cpp | cut -d"/" -f2 |  cut -d"." -f1`; do
	program="$program_dir/$prog_name.cpp"
	mkdir -p lib/$prog_name


	echo "  compiling $program"
	$compiler  -c $program 
	mv $prog_name.o lib/$prog_name/ 
done

#================================================================
# linking programs
#================================================================
echo
echo "Linking programs..."

for prog_name in `ls $program_dir/*.cpp | cut -d"/" -f2 |  cut -d"." -f1`; do
    echo "  linking $prog_name"
    $compiler -ldl -o ./bin/$prog_name.$SUFF `ls obj/*` `ls lib/$prog_name/$prog_name*.o`
done
#================================================================



# 2> /dev/null &