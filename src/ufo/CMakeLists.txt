# (C) Copyright 2017-2021 UCAR.
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
    instantiateObsErrorFactory.h
    instantiateObsFilterFactory.h
    instantiateObsLocFactory.h
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
    ObsBiasParameters.h
    ObsDiagnostics.cc
    ObsDiagnostics.h
    ObsOperator.cc
    ObsOperator.h
    ObsOperatorBase.cc
    ObsOperatorBase.h
    ObsOperatorParametersBase.h
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
add_subdirectory( errors )
add_subdirectory( filters )
add_subdirectory( obslocalization )
add_subdirectory( identity )
add_subdirectory( atmvertinterp )
add_subdirectory( atmvertinterplay )
add_subdirectory( atmsfcinterp )
add_subdirectory( avgkernel )
if( crtm_FOUND )
    add_subdirectory( crtm )
endif()
add_subdirectory( radarradialvelocity )
if( rttov_FOUND )
    add_subdirectory( rttov )
    add_subdirectory( rttovcpp )
endif()
add_subdirectory( gnssro )
add_subdirectory( sfcpcorrected )
add_subdirectory( marine )
add_subdirectory( aerosols )
add_subdirectory( timeoper )
add_subdirectory( profile )
add_subdirectory( groundgnss )
add_subdirectory( variabletransforms )
add_subdirectory( scatwind )
add_subdirectory( sattcwv )
add_subdirectory( compositeoper )
add_subdirectory( backgrounderrorvertinterp )
add_subdirectory( backgrounderroridentity )
add_subdirectory( categoricaloper )

list( APPEND ufo_src_files
    ${utils_src_files}
    ${basis_src_files}
    ${predictor_src_files}
    ${errors_src_files}
    ${filters_src_files}
    ${obslocalization_src_files}
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
    ${radarradialvelocity_src_files}
    ${timeoper_src_files}
    ${profile_src_files}
    ${groundgnss_src_files}
    ${variabletransforms_src_files}
    ${scatwind_src_files}
    ${sattcwv_src_files}
    ${compositeoper_src_files}
    ${backgrounderrorvertinterp_src_files}
    ${backgrounderroridentity_src_files}
    ${categoricaloper_src_files}
)
list( APPEND ufo_src_files ${marine_src_files} )

## Create Library target
ecbuild_add_library( TARGET   ufo
                     SOURCES  ${ufo_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${INSTALL_INCLUDE_DIR}/${PROJECT_NAME}
                    )

## Link dependencies
target_link_libraries(ufo PUBLIC MPI::MPI_C MPI::MPI_CXX MPI::MPI_Fortran)
target_link_libraries(ufo PUBLIC Boost::boost)
target_link_libraries(ufo PUBLIC NetCDF::NetCDF_C NetCDF::NetCDF_Fortran)
target_link_libraries(ufo PUBLIC Eigen3::Eigen)
target_link_libraries(ufo PUBLIC gsl::gsl-lite)
target_link_libraries(ufo PUBLIC eckit)
target_link_libraries(ufo PUBLIC fckit)
target_link_libraries(ufo PUBLIC ioda)
target_link_libraries(ufo PUBLIC oops)

# Optional dependencies
if(crtm_FOUND)
    target_link_libraries(ufo PUBLIC crtm)
endif()

if(gsw_FOUND)
    target_link_libraries(ufo PUBLIC gsw)
endif()

if(rttov_FOUND)
    target_link_libraries(ufo PUBLIC rttov)
    target_compile_definitions(ufo PUBLIC RTTOV_FOUND)
endif()

if(ropp-ufo_FOUND)
    target_link_libraries(ufo PUBLIC ropp-ufo)
endif()

if(geos-aero_FOUND)
    target_link_libraries(ufo PUBLIC geos-aero)
endif()


## Include paths
target_include_directories(ufo PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
                                      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)


## Fortran modules
set(MODULE_DIR ${PROJECT_NAME}/module)
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${CMAKE_INSTALL_LIBDIR}/${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}/${MODULE_DIR}>)

