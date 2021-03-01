# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND ufo_src_files
    Fortran.h
    AnalyticInit.cc
    AnalyticInit.h
    GeoVaLs.cc
    GeoVaLs.h
    GeoVaLs.interface.F90
    GeoVaLs.interface.h
    instantiateObsFilterFactory.h
    LinearObsBiasOperator.cc
    LinearObsBiasOperator.h
    LinearObsOperator.cc
    LinearObsOperator.h
    LinearObsOperatorBase.cc
    LinearObsOperatorBase.h
    Locations.cc
    Locations.h
    ObsBias.cc
    ObsBias.h
    ObsBiasCovariance.cc
    ObsBiasCovariance.h
    ObsBiasIncrement.cc
    ObsBiasIncrement.h
    ObsBiasOperator.cc
    ObsBiasOperator.h
    ObsDiagnostics.cc
    ObsDiagnostics.h
    ObsOperator.cc
    ObsOperator.h
    ObsOperatorBase.cc
    ObsOperatorBase.h
    ObsTraits.h
    locations_f.cc
    locations_f.h
    ufo_geovals_mod.F90
    ufo_locations_mod.F90
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
add_subdirectory( predictors )
add_subdirectory( filters )
add_subdirectory( identity )
add_subdirectory( atmvertinterp )
add_subdirectory( atmvertinterplay )
add_subdirectory( atmsfcinterp )
add_subdirectory( avgkernel )
if( ${CRTM_FOUND} )
    add_subdirectory( crtm )
endif( ${CRTM_FOUND} )
add_subdirectory( radarreflectivity )
add_subdirectory( radarradialvelocity )
if( ${RTTOV_FOUND} )
    add_subdirectory( rttov )
    add_subdirectory( rttovcpp )
endif( ${RTTOV_FOUND} )
add_subdirectory( gnssro )
add_subdirectory( sfcpcorrected )
add_subdirectory( marine )
add_subdirectory( aerosols )
add_subdirectory( timeoper )
add_subdirectory( profile )
add_subdirectory( groundgnss )
add_subdirectory( calculate )
add_subdirectory( scatwind )
add_subdirectory( compositeoper )
add_subdirectory( backgrounderrorvertinterp )
add_subdirectory( backgrounderroridentity )

list( APPEND ufo_src_files
    ${utils_src_files}
    ${basis_src_files}
    ${predictor_src_files}
    ${filters_src_files}
    ${identity_src_files}
    ${atmvertinterp_src_files}
    ${atmvertinterplay_src_files}
    ${atmsfcinterp_src_files}
    ${avgkernel_src_files}
    ${crtm_src_files}
    ${rttov_src_files}
    ${rttovcpp_src_files}
    ${gnssro_src_files}
    ${aerosols_src_files}
    ${sfcpcorrected_src_files}
    ${radarreflectivity_src_files}
    ${radarradialvelocity_src_files}
    ${timeoper_src_files}
    ${profile_src_files}
    ${groundgnss_src_files}
    ${calculate_src_files}
    ${scatwind_src_files}
    ${compositeoper_src_files}
    ${backgrounderrorvertinterp_src_files}
    ${backgrounderroridentity_src_files}
)
list( APPEND ufo_src_files ${marine_src_files} )
list( APPEND UFO_LIBS_DEP oops ioda ioda_engines ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} )
if( ${CRTM_FOUND} )
    list( APPEND UFO_LIBS_DEP crtm )
endif( ${CRTM_FOUND} )
if( ${GSW_FOUND} )
    list( APPEND UFO_LIBS_DEP gsw )
endif( ${GSW_FOUND} )
if( ${RTTOV_FOUND} )
    list( APPEND UFO_LIBS_DEP rttov )
endif( ${RTTOV_FOUND} )

if( ${ROPP-UFO_FOUND} )
    list( APPEND UFO_LIBS_DEP ropp-ufo )
endif( ${ROPP-UFO_FOUND} )

if( ${GEOS-AERO_FOUND} )
    list( APPEND UFO_LIBS_DEP geos-aero )
endif( ${GEOS-AERO_FOUND} )


ecbuild_add_library( TARGET   ufo
                     SOURCES  ${ufo_src_files}
                     LIBS     ${UFO_LIBS_DEP}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${INSTALL_INCLUDE_DIR}/ufo
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

if (${RTTOV_FOUND})
  target_compile_definitions(ufo PUBLIC RTTOV_FOUND)
endif()
