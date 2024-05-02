#!/bin/sh
example="Example"
example_path="tools/new_obsop/example"

echo "Usage: create_obsop_fromexample.sh <ObsOperatorName> <path/to/obsoperator>"
echo "Example: create_obsop_fromexample.sh CRTM crtm2"
echo "Please use CamelCase for the ObsOperatorName and path relative to src/ufo/"
 
if [ "$#" -lt 2 ]; then
    echo "Error: must call with 2 parameters"
    exit 1
fi

generate=$1
generate_dir=$2

example_lc=`echo ${example} | perl -ne 'print lc'`
generate_lc=`echo ${generate} | perl -ne 'print lc'`

example_cpp_define=`echo ${example_path}_OBS${example} | perl -ne 'print uc' | perl -pe "s/\//_/g"`
generate_cpp_define=`echo ${generate_dir}_OBS${generate} | perl -ne 'print uc' | perl -pe "s/\//_/g"`

mkdir -p ${generate_dir}
cp example/CMakeLists.txt ${generate_dir}
cp example/Obs${example}.cc     ${generate_dir}/Obs${generate}.cc
cp example/Obs${example}.h      ${generate_dir}/Obs${generate}.h
cp example/Obs${example}TLAD.cc ${generate_dir}/Obs${generate}TLAD.cc
cp example/Obs${example}TLAD.h  ${generate_dir}/Obs${generate}TLAD.h
cp example/Obs${example}.interface.F90  ${generate_dir}/Obs${generate}.interface.F90
cp example/Obs${example}.interface.h    ${generate_dir}/Obs${generate}.interface.h
cp example/Obs${example}TLAD.interface.F90  ${generate_dir}/Obs${generate}TLAD.interface.F90
cp example/Obs${example}TLAD.interface.h    ${generate_dir}/Obs${generate}TLAD.interface.h
cp example/ufo_${example_lc}_mod.F90        ${generate_dir}/ufo_${generate_lc}_mod.F90
cp example/ufo_${example_lc}_tlad_mod.F90   ${generate_dir}/ufo_${generate_lc}_tlad_mod.F90
cp example/Obs${example}Parameters.h ${generate_dir}/Obs${generate}Parameters.h

# replace Example class name with new name
perl -p -i -e "s/${example}/${generate}/g" ${generate_dir}/*
# replace example struct and routine names in the Fortran calls
perl -p -i -e "s/ufo_${example_lc}/ufo_${generate_lc}/g" ${generate_dir}/*
# replace include headers in *cc and *h files
perl -p -i -e "s#${example_path}#${generate_dir}#g" ${generate_dir}/Obs*
# replace example in the rest of the files
perl -p -i -e "s/${example_lc}/${generate_lc}/g" ${generate_dir}/*

echo "Directory created"
