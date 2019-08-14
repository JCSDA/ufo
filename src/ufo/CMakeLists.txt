# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND ufo_src_files
    Fortran.h
    GeoVaLs.cc
    GeoVaLs.h
    GeoVaLs.interface.F90
    GeoVaLs.interface.h
    LinearObsOperator.cc
    LinearObsOperator.h
    LinearObsOperatorBase.cc
    LinearObsOperatorBase.h
    Locations.cc
    Locations.h
    Locations.interface.F90
    Locations.interface.h
    ObsBias.h
    ObsBiasCovariance.h
    ObsBiasIncrement.h
    ObsOperator.cc
    ObsOperator.h
    ObsOperatorBase.cc
    ObsOperatorBase.h
    UfoTrait.h
    ufo_geovals_mod.F90
    ufo_locs_mod.F90
    ufo_variables_mod.F90
    ufo_constants_mod.F90
)

# macro to prepend a prefix with relative path
# can this be added to ecbuild for use elsewhere?
function(PREPEND var prefix )
    set ( listVar "" )
    foreach (f ${ARGN})
        list (APPEND listVar "${prefix}/${f}")
    endforeach(f)
    set ( ${var} "${listVar}" PARENT_SCOPE )
endfunction(PREPEND)

add_subdirectory( utils )
add_subdirectory( basis )
add_subdirectory( filters )
add_subdirectory( identity )
add_subdirectory( atmvertinterp )
add_subdirectory( atmvertinterplay )
add_subdirectory( atmsfcinterp )
add_subdirectory( crtm )
if( ${RTTOV_FOUND} )
    add_subdirectory( rttov )
endif( ${RTTOV_FOUND} )
add_subdirectory( gnssro )
add_subdirectory( surface )
add_subdirectory( ncdiag )
add_subdirectory( marine )
list( APPEND ufo_src_files
    ${utils_src_files}
    ${basis_src_files}
    ${filters_src_files}
    ${identity_src_files}
    ${atmvertinterp_src_files}
    ${atmvertinterplay_src_files}
    ${atmsfcinterp_src_files}
    ${crtm_src_files}
    ${rttov_src_files}
    ${gnssro_src_files}
    ${surface_src_files}
    ${ncdiag_src_files}
)
list( APPEND ufo_src_files ${marine_src_files} )
list( APPEND UFO_LIBS_DEP oops ioda crtm ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} )
if( ${GSW_FOUND} )
    list( APPEND UFO_LIBS_DEP gsw )
endif( ${GSW_FOUND} )
if( ${RTTOV_FOUND} )
    list( APPEND UFO_LIBS_DEP rttov )
endif( ${RTTOV_FOUND} )

if( ${ROPP-UFO_FOUND} )
    list( APPEND UFO_LIBS_DEP ropp-ufo )
endif( ${ROPP-UFO_FOUND} )

ecbuild_add_library( TARGET   ufo
                     SOURCES  ${ufo_src_files}
                     LIBS     ${UFO_LIBS_DEP}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )
