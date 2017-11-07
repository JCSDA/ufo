list( APPEND ufo_src_files
Fortran.h
GeoVaLs.cc
GeoVaLs.h
Locations.h
ObsBias.h
ObsSpace.cc
ObsSpace.h
ObsVector.cc
ObsVector.h
ObsWSpeed.cc
ObsWSpeed.h
UfoTrait.h
Variables.h
)

ecbuild_add_library( TARGET   ufo
                     SOURCES  ${ufo_src_files}
                     LIBS     oops ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

