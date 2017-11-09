#!/bin/csh -f
setenv BASE "/home/gvernier/Sandboxes/codesprint-ufo/jedi/"
setenv ALLSRC ${BASE}"code/"
setenv ALLBUILDS ${BASE}"build/"

#Sources
setenv UFO_SRC ${ALLSRC}"ufo/"
setenv OOPS_SRC ${ALLSRC}"oops/"

#Builds 
setenv UFO_BUILD ${ALLBUILDS}"ufo/"
setenv OOPS_BUILD ${ALLBUILDS}"oops/" 

#CRTM
setenv CRTM_PATH /home/gvernier/Sandboxes/codesprint-ufo/jedi/code/crtm-release/libsrc/
#${ALLBUILDS}"gsi/lib/" 
setenv CRTM_LIBRARIES /home/gvernier/Sandboxes/codesprint-ufo/jedi/code/crtm-release/libsrc/libcrtm.a
#"$CRTM_PATH/libcrtm_v2.2.3.a"
setenv CRTM_INCLUDE /home/gvernier/Sandboxes/codesprint-ufo/jedi/code/crtm-release/libsrc/
#${ALLBUILDS}"gsi/include/" 

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
#exit 0
