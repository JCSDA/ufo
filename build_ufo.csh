#!/bin/csh -f

#User defined stuff
setenv BASE "/home/vagrant/jedi/"

#All codes and builds
setenv ALLSRC ${BASE}"code/"
setenv ALLBUILDS ${BASE}"build/"

#Sources
setenv UFO_SRC ${ALLSRC}"ufo/"
setenv OOPS_SRC ${ALLSRC}"oops/"

#Builds 
setenv UFO_BUILD ${ALLBUILDS}"ufo/"
setenv OOPS_BUILD ${ALLBUILDS}"oops/" 

#Point Env. var. to CRTM build ...
setenv CRTM_PATH ${BASE}/code/crtm/libsrc/
setenv CRTM_LIBRARIES ${BASE}/code/crtm/libsrc/libcrtm.a
setenv CRTM_INCLUDE ${BASE}/code/crtm/libsrc/

#CLEAN BUILD
rm -rf ${UFO_BUILD}
mkdir ${UFO_BUILD}
cd ${UFO_BUILD}

ecbuild -DOOPS_PATH=${OOPS_BUILD} \
        -DLAPACK_LIBRARIES=$LAPACK_LIBRARIES \
	-DCRTM_LIBRARIES=$CRTM_LIBRARIES \
        -DCRTM_INCLUDE=$CRTM_INCLUDE \
        --build=release \
	${UFO_SRC}
make -j4
exit 0
