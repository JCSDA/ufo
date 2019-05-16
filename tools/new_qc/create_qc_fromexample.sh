#!/bin/sh
example="Example"
example_dir="example"

echo "Usage: create_qc_fromexample.sh <QcName>"
echo "Example: create_qc_fromexample.sh Background2"
echo "Please use CamelCase for the QcName"
 
if [ "$#" -lt 1 ]; then
    echo "Error: must call with 1 parameter"
    exit 1
fi

generate=$1
generate_dir=../../src/ufo/filters/

example_lc=`echo ${example} | perl -ne 'print lc'`
generate_lc=`echo ${generate} | perl -ne 'print lc'`

example_cpp_define=`echo ${example} | perl -ne 'print uc'`
generate_cpp_define=`echo ${generate} | perl -ne 'print uc'`

cp ${example_dir}/${example}.cc     ${generate_dir}/${generate}.cc
cp ${example_dir}/${example}.h      ${generate_dir}/${generate}.h
cp ${example_dir}/${example}.interface.F90 ${generate_dir}/${generate}.interface.F90
cp ${example_dir}/${example}.interface.h   ${generate_dir}/${generate}.interface.h
cp ${example_dir}/ufo_${example_lc}_mod.F90  ${generate_dir}/ufo_${generate_lc}_mod.F90

# replace the defines in *h files
perl -p -i -e "s/${example_cpp_define}/${generate_cpp_define}/g" ${generate_dir}/${generate}.*h
# replace Example class name with new name
perl -p -i -e "s/${example}/${generate}/g" ${generate_dir}/${generate}.*
# replace example struct and routine names in the Fortran calls
perl -p -i -e "s/ufo_${example_lc}/ufo_${generate_lc}/g" ${generate_dir}/ufo_${generate_lc}_mod.F90
perl -p -i -e "s/ufo_${example_lc}/ufo_${generate_lc}/g" ${generate_dir}/${generate}.*
# replace include headers in *cc files
perl -p -i -e "s#${example_dir}#${generate_path}#g" ${generate_dir}/${generate}.cc
# replace example in the rest of the files
perl -p -i -e "s/${example_lc}/${generate_lc}/g" ${generate_dir}/ufo_${generate_lc}_mod.F90

