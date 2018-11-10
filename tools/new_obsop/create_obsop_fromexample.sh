#!/bin/sh
example="Example"
example_dir="example"

echo "Usage: create_obsop_fromexample.sh <ObsOperatorName> <path/to/obsoperator>"
echo "Example: create_obsop_fromexample.sh Radiance atmosphere/radiance2"
echo "Please use CamelCase for the ObsOperatorName and path relative to src/ufo/"
 
if [ "$#" -lt 2 ]; then
    echo "Error: must call with 2 parameters"
    exit 1
fi

generate=$1
generate_dir=../../src/ufo/$2
generate_path=$2

example_lc=`echo ${example} | perl -ne 'print lc'`
generate_lc=`echo ${generate} | perl -ne 'print lc'`

example_cpp_define=`echo ${example_dir}_OBS${example} | perl -ne 'print uc'`
generate_cpp_define=`echo ${generate_path}_OBS${generate} | perl -ne 'print uc' | perl -pe "s/\//_/g"`

mkdir -p ${generate_dir}
cp ${example_dir}/CMakeLists.txt ${generate_dir}
cp ${example_dir}/Obs${example}.cc     ${generate_dir}/Obs${generate}.cc
cp ${example_dir}/Obs${example}.h      ${generate_dir}/Obs${generate}.h
cp ${example_dir}/Obs${example}TLAD.cc ${generate_dir}/Obs${generate}TLAD.cc
cp ${example_dir}/Obs${example}TLAD.h  ${generate_dir}/Obs${generate}TLAD.h
cp ${example_dir}/ufo_${example_lc}_interface.F90  ${generate_dir}/ufo_${generate_lc}_interface.F90
cp ${example_dir}/ufo_${example_lc}_mod.F90        ${generate_dir}/ufo_${generate_lc}_mod.F90
cp ${example_dir}/ufo_${example_lc}_tlad_interface.F90  ${generate_dir}/ufo_${generate_lc}_tlad_interface.F90
cp ${example_dir}/ufo_${example_lc}_tlad_mod.F90        ${generate_dir}/ufo_${generate_lc}_tlad_mod.F90

# replace the defines in *h files
perl -p -i -e "s/${example_cpp_define}/${generate_cpp_define}/g" ${generate_dir}/Obs*.h
# replace Example class name with new name
perl -p -i -e "s/${example}/${generate}/g" ${generate_dir}/*
# replace example struct and routine names in the Fortran calls
perl -p -i -e "s/ufo_${example_lc}/ufo_${generate_lc}/g" ${generate_dir}/*
# replace include headers in *cc files
perl -p -i -e "s#${example_dir}#${generate_path}#g" ${generate_dir}/Obs*.cc
# replace example in the rest of the files
perl -p -i -e "s/${example_lc}/${generate_lc}/g" ${generate_dir}/*
