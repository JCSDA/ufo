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
generate_dir=../../src/ufo/

example_lc=`echo ${example} | perl -ne 'print lc'`
generate_lc=`echo ${generate} | perl -ne 'print lc'`

example_cpp_define=`echo ${example} | perl -ne 'print uc'`
generate_cpp_define=`echo ${generate} | perl -ne 'print uc'`

cp ${example_dir}/${example}Check.cc     ${generate_dir}/${generate}Check.cc
cp ${example_dir}/${example}Check.h      ${generate_dir}/${generate}Check.h
cp ${example_dir}/${example}Check.interface.F90 ${generate_dir}/${generate}Check.interface.F90
cp ${example_dir}/${example}Check.interface.h   ${generate_dir}/${generate}Check.interface.h
cp ${example_dir}/ufo_${example_lc}check_mod.F90  ${generate_dir}/ufo_${generate_lc}check_mod.F90

# replace the defines in *h files
perl -p -i -e "s/${example_cpp_define}/${generate_cpp_define}/g" ${generate_dir}/${generate}Check.*h
# replace Example class name with new name
perl -p -i -e "s/${example}/${generate}/g" ${generate_dir}/${generate}Check.*
# replace example struct and routine names in the Fortran calls
perl -p -i -e "s/ufo_${example_lc}/ufo_${generate_lc}/g" ${generate_dir}/ufo_${generate_lc}check_mod.F90
perl -p -i -e "s/ufo_${example_lc}/ufo_${generate_lc}/g" ${generate_dir}/${generate}Check.*
# replace include headers in *cc files
perl -p -i -e "s#${example_dir}#${generate_path}#g" ${generate_dir}/${generate}Check.cc
# replace example in the rest of the files
perl -p -i -e "s/${example_lc}/${generate_lc}/g" ${generate_dir}/ufo_${generate_lc}check_mod.F90

