file( GLOB gsw_src_files ./marine/gsw/*.f90 )
list( APPEND ufo_src_files ${gsw_src_files} )

list( APPEND ufo_src_files
Fortran.h
GeoVaLs.cc
GeoVaLs.h
ObsBias.h
ObsBiasCovariance.h
ObsBiasIncrement.h
ObsCheck.cc
ObsCheck.h
UfoTrait.h
instantiateObsOperatorFactory.h
instantiateLinearObsOpFactory.h
ufo_obscheck_mod.F90
GeoVaLs.interface.F90
ufo_geovals_mod.F90
ufo_variables_mod.F90
ufo_obs_data_basis_mod.F90
ufo_obs_data.F90
ufo_obs_data_mod.F90
# Atmos: utils
atmosphere/utils/vert_interp.F90
# Atmos: Radiosondes
atmosphere/radiosonde/ObsRadiosonde.h
atmosphere/radiosonde/ObsRadiosondeTLAD.h
atmosphere/radiosonde/ObsRadiosonde.interface.F90
atmosphere/radiosonde/ObsRadiosondeTLAD.interface.F90
atmosphere/radiosonde/ufo_radiosonde_mod.F90
atmosphere/radiosonde/ufo_radiosonde_tlad_mod.F90
# Atmos: Radiances
atmosphere/radiance/ObsRadiance.h
atmosphere/radiance/ObsRadianceTLAD.h
atmosphere/radiance/ObsRadiance.interface.F90
atmosphere/radiance/ObsRadianceTLAD.interface.F90
atmosphere/radiance/ufo_radiance_mod.F90
atmosphere/radiance/ufo_radiance_tlad_mod.F90
# Marine: Sea-ice fraction
./marine/seaicefraction/ObsSeaIceFraction.h
./marine/seaicefraction/ObsSeaIceFractionTLAD.h
./marine/seaicefraction/ObsSeaIceFraction.interface.F90
./marine/seaicefraction/ObsSeaIceFractionTLAD.interface.F90
./marine/seaicefraction/ufo_seaicefrac_mod.F90
./marine/seaicefraction/ufo_seaicefrac_tlad_mod.F90
# Marine: Sea-ice thickness
./marine/seaicethickness/ObsSeaIceThickness.h
./marine/seaicethickness/ObsSeaIceThicknessTLAD.h
./marine/seaicethickness/ObsSeaIceThickness.interface.F90
./marine/seaicethickness/ObsSeaIceThicknessTLAD.interface.F90
./marine/seaicethickness/ufo_seaicethick_mod.F90
./marine/seaicethickness/ufo_seaicethick_tlad_mod.F90
# Marine: Steric height
./marine/stericheight/ObsStericHeight.h
./marine/stericheight/ObsStericHeightTLAD.h
./marine/stericheight/ObsStericHeight.interface.F90
./marine/stericheight/ObsStericHeightTLAD.interface.F90
./marine/stericheight/ufo_stericheight_mod.F90
./marine/stericheight/ufo_stericheight_tlad_mod.F90
./marine/obsop/ufo_steric.F90
# Constituents: AOD
constituents/aod/ObsAod.h
constituents/aod/ObsAodTLAD.h
constituents/aod/ObsAod.interface.F90
constituents/aod/ObsAodTLAD.interface.F90
constituents/aod/ufo_aod_mod.F90
constituents/aod/ufo_aod_tlad_mod.F90
constituents/aod/ufo_aod_misc.F90
)

message(STATUS ${ufo_src_files})
ecbuild_add_library( TARGET   ufo
                     SOURCES  ${ufo_src_files}
                     LIBS     oops ioda ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${CRTM_LIBRARIES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

