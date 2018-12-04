# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND ufo_src_files
    BackgroundCheck.cc
    BackgroundCheck.h
    BackgroundCheck.interface.F90
    BackgroundCheck.interface.h
    Fortran.h
    Locations.cc
    Locations.h
    Locations.interface.F90
    Locations.interface.h
    FortranGeoVals.h
    GeoVaLs.cc
    GeoVaLs.h
    LinearObsOperator.cc
    LinearObsOperator.h
    LinearObsOperatorBase.cc
    LinearObsOperatorBase.h
    ObsBias.h
    ObsBiasCovariance.h
    ObsBiasIncrement.h
    ObsOperator.cc
    ObsOperator.h
    ObsOperatorBase.cc
    ObsOperatorBase.h
    UfoTrait.h
    ufo_bgcheck_mod.F90
    ufo_locs_mod.F90
    ufo_geovals_interface.F90
    ufo_geovals_mod.F90
    ufo_variables_mod.F90
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
add_subdirectory( atmosphere )
add_subdirectory( constituents )
if( ${GSW_FOUND} )
    add_subdirectory( marine )
endif( ${GSW_FOUND} )
list( APPEND ufo_src_files
    ${utils_src_files}
    ${basis_src_files}
    ${atmosphere_src_files}
    ${constituents_src_files}
)
if( ${GSW_FOUND} )
    list( APPEND ufo_src_files ${marine_src_files} )
endif( ${GSW_FOUND} )
list( APPEND UFO_LIBS_DEP oops ioda crtm ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} )
if( ${GSW_FOUND} )
    list( APPEND UFO_LIBS_DEP gsw )
endif( ${GSW_FOUND} )

if( ${ROPP-UFO_FOUND} )
    list( APPEND UFO_LIBS_DEP ropp-ufo )
endif( ${ROPP-UFO_FOUND} )

message(STATUS ${ufo_src_files})

ecbuild_add_library( TARGET   ufo
                     SOURCES  ${ufo_src_files}
                     LIBS     ${UFO_LIBS_DEP}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )
