/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_RTTOV_OBSRADIANCERTTOVPARAMETERS_H_
#define UFO_RTTOV_OBSRADIANCERTTOVPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

class RTTOVObsOptionsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(RTTOVObsOptionsParameters, Parameters)

 public:
  /// This is a string used to create the coefficient file name for example "noaa_20_atms" to make
  /// RTTOV use rtcoef_noaa_20_atms.dat
  oops::RequiredParameter<std::string> sensorID{"Sensor_ID", this};

  /// The path to the coefficient files.
  oops::RequiredParameter<std::string> coefficientPath{"CoefficientPath", this};

  /// Should RTTOV convert from mixing ratio to ppmv for rttov processing
  /// RTTOV can handle either but the output BTs obtained can be slightly different
  oops::Parameter<bool> RTTOVGasUnitConv{"RTTOV_GasUnitConv", false, this};

  /// The default option to setup the code for a particular version of RTTOV
  oops::Parameter<std::string> RTTOVDefaultOpts{"RTTOV_default_opts", "RTTOV", this};

  /// Set a series of options in the radiative transfer interface so that it matches with
  /// the OPS SatRad setup.
  oops::Parameter<bool> satRadCompatibility{"SatRad_compatibility", true, this};

  /// Use liquid water in the q saturation calculations, only used when
  /// SatRad_compatibility is true.
  oops::Parameter<bool> useRHwaterForQC{"UseRHwaterForQC", true, this};

  /// Use a check on cold surfaces to update the surface temperature and skin temperature
  /// This is legacy code from the Met Office OPS and should generally not be used.
  /// It is included at the moment for backwards compatibility testing.  Only used when
  /// SatRad_compatibility is true.
  oops::Parameter<bool> useColdSurfaceCheck{"UseColdSurfaceCheck", false, this};

  /// Maximum number of profiles to be processed by RTTOV per pass.  If this is true then
  /// one profile is run at a time.  This will be slow but useful for debugging.  The normal number
  /// of profiles per pass is nchan_max_sim / nchan_inst.
  oops::Parameter<bool> profByProf{"prof_by_prof", false, this};

  /// This allows the user to specify the maximum number of channels to simulate per batch
  /// (nchan_max_sim)
  oops::Parameter<int> maxChanPerBatch{"max_channels_per_batch", 10000, this};

  /// If performing the qsplit calculation for qtotal within the interface this option
  /// allows for the inclusion of rain in the calculation.  Although the default is false
  /// it is turned to true in the interface if simulating a microwave instrument with
  /// SatRad_compatibility set to true.
  oops::OptionalParameter<bool> qTSplitRain{"QtSplitRain", this};

  /// Check the rttov profile using the rttov check profile routine and flag if the check is failed.
  /// This check makes sure the values read from the geovals are within certain bounds.
  oops::Parameter<bool> RTTOVProfileCheckInput{"RTTOV_profile_checkinput", false, this};

  /// Specify the profiles where further diagnostics are needed
  /// for debugging e.g. [1, 2, 3]
  /// Note this is numbering from 1 as it is Fortran that will provide output.
  oops::OptionalParameter<std::vector<int>> inspectProfileNumber{"InspectProfileNumber", this};

  /// -----------------------------------------------------------------------------------
  /// RTTOV options all of these are loaded using set_options_rttov and are not required
  /// because there are defaults within the RTTOV code.  The values are taken from
  /// rttov/src/rttov/main/rttov_types.F90
  /// -------------------------------------------------
  /// General radiative transfer options: opts % rt_all
  /// -------------------------------------------------
  /// If true RTTOV calculations account for atmospheric refraction
  /// default is false
  oops::OptionalParameter<bool> RTTOVAddRefrac{"RTTOV_addrefrac", this};

  /// Determines input perturbation for AD/K routines: if true
  /// radiance_ad/k % bt is used for channels with wavelengths >
  /// 3µm, otherwise radiance_ad/k % total is used.
  /// default is false
  oops::OptionalParameter<bool> RTTOVSwitchRad{"RTTOV_switchrad", this};

  /// If true activate use of surface humidity
  /// default is true
  oops::OptionalParameter<bool> RTTOVUseq2m{"RTTOV_use_q2m", this};

  /// If true activate treatment of surface as Lambertian instead of
  /// specular reflector for downwelling emitted radiance.
  /// The specularity variable determines how specular or
  /// Lambertian the surface is.
  /// default is false
  oops::OptionalParameter<bool> RTTOVDoLambertian{"RTTOV_do_lambertian", this};

  /// If true the Lambertian downwelling radiance is computed using a
  /// fixed effective angle. If false the angle is derived from a
  /// parameterisation computed from the total atmospheric optical
  /// depth. Only relevant if do_lambertian is true.
  /// default is true
  oops::OptionalParameter<bool> RTTOVLambertianFixedAngle{"RTTOV_lambertian_fixed_angle", this};

  /// If true treat atmosphere as strictly plane-parallel (no curvature or
  /// refraction); automatically activated for all simulated radiances if
  /// the DOM scattering solver is used.
  /// default is false
  oops::OptionalParameter<bool> RTTOVPlaneParallel{"RTTOV_plane_parallel", this};

  /// If true the downwelling atmospheric emission is computed using
  /// the linear-in-tau approximation for the Planck source term. If
  /// false a simpler layer average source term is used which is slightly
  /// faster and has negligible impact on radiances (RTTOV 12 default = true,
  /// recommended = false). Note that the Lambertian downwelling
  /// radiances always use the layer average source term.
  oops::OptionalParameter<bool> RTTOVRadDownLinTau{"RTTOV_rad_down_lin_tau", this};

  /// This switch enables/disables a check on the delta-transmittance in
  /// the near-surface layer when considering this layer contribution in
  /// the integration of the RTE. By setting to false (omitting the
  /// check) we reduce discontinuities in the direct model with very
  /// small (mK) impact on radiances. Option deprecated.
  /// default is true
  oops::OptionalParameter<bool> RTTOVDtauTest{"RTTOV_dtau_test", this};

  /// -------------------------------------------------
  /// MW-only radiative transfer options: opts % rt_mw
  /// -------------------------------------------------
  /// Select the MW sea surface emissivity model to use. Valid values
  /// are 0-6; values 1-6 specify FASTEM-1 to FASTEM-6; value 0
  /// specifies TESSEM2.
  /// default is 6
  oops::OptionalParameter<int> RTTOVFastemVersion{"RTTOV_fastem_version", this};

  /// If true use the foam fraction (0-1) value specified in
  /// profiles(:)%skin%foam_fraction in FASTEM which is loaded from the geovals.
  /// Does not apply to TESSEM2.
  /// default is false.
  oops::OptionalParameter<bool> RTTOVSupplyFoamFraction{"RTTOV_supply_foam_fraction", this};

  /// If true user is supplying cloud liquid water profiles.
  /// This applies to “clear-sky” simulations only: the cloud
  /// is treated as a purely absorbing medium. For full scattering use
  /// RTTOV-SCATT instead.
  /// default is false
  oops::OptionalParameter<bool> RTTOVCLWData{"RTTOV_clw_data", this};

  /// Choose CLW permittivity parameterisation: 1=>Liebe (1989),
  /// 2=>Rosekranz (2015), 3=>Turner,Kneifel,Cadeddu(2016)
  /// default is 1
  oops::OptionalParameter<int> RTTOVCLWScheme{"RTTOV_clw_scheme", this};

  /// Apply MW CLW calculations on coef/user levels (true/false respectively)
  /// default is true
  oops::OptionalParameter<bool> RTTOVCLWCalcOnCoefLev{"RTTOV_clw_calc_on_coef_lev", this};

  /// MW-only. Lower pressure limit for MW CLW calculations (hPa): any CLW
  /// in levels with pressures smaller than this value is ignored.
  /// default in RTTOV12 is 322.0 hPa
  oops::OptionalParameter<float> RTTOVCLWCloudTop{"RTTOV_clw_cloud_top", this};

  /// Apply band-correction for Planck radiance and BT calculations.
  /// default is true
  oops::OptionalParameter<bool> RTTOVApplyBandCorrection{"RTTOV_apply_band_correction", this};

  /// ---------------------------------------------------------------------------
  /// Options related to interpolation on the vertical grid: opts % interpolation
  /// ---------------------------------------------------------------------------
  /// If true input profiles may be supplied on user-defined levels, and
  /// internal interpolation is used.
  /// default is false
  oops::OptionalParameter<bool> RTTOVAddInterp{"RTTOV_addinterp", this};

  /// Set the interpolation mode. See table 7 of the RTTOV12 userguide for more details.
  /// Valid values are 1-5.
  /// default is 1
  oops::OptionalParameter<int> RTTOVInterpMode{"RTTOV_interp_mode", this};

  /// Allow TL/AD/K of user pressure levels ; only applies if addinterp is true.
  /// default is false
  oops::OptionalParameter<bool> RTTOVLGradP{"RTTOV_lgradp", this};

  /// If true treat user's model top as space boundary (default = true,
  /// it is not recommended to set this false). Option deprecated.
  oops::OptionalParameter<bool> RTTOVSpaceTop{"RTTOV_spacetop", this};

  /// Extrapolate input profiles up to top coefficient level maintaining relative
  /// values with respect to regression limits - see section 7.3 of RTTOV 12 user
  /// guide (default = false); only applies if addinterp is true. Option deprecated.
  oops::OptionalParameter<bool> RTTOVRegLimitExtrap{"RTTOV_reg_limit_extrap", this};

  /// ----------------------------------------------
  /// General configurations options: opts % config
  /// ----------------------------------------------
  /// If true input profiles outside the limits specified in the coefficient
  /// files are reset to the min/max. If false such profiles will generate
  /// warning messages unless verbose is false. (default = false)
  oops::OptionalParameter<bool> RTTOVApplyRegLimits{"RTTOV_apply_reg_limits", this};

  /// If false, only messages for fatal errors are output
  /// (default = true)
  oops::OptionalParameter<bool> RTTOVVerbose{"RTTOV_verbose", this};

  /// If true checks whether input profiles are within both absolute and
  /// regression limits. If false no check is performed. (default = true)
  oops::OptionalParameter<bool> RTTOVDoCheckInput{"RTTOV_do_checkinput", this};

  /// If false the input surface elevation is assigned to the pressure
  /// level at or immediately below the input 2m (surface) pressure. If
  /// true the surface elevation is assigned to the specified surface
  /// pressure (default = false). Option deprecated.
  oops::OptionalParameter<bool> RTTOVFixHgpl{"RTTOV_fix_hgpl", this};

  /// --------------------------------------------------------
  /// Visible/IR-only radiative transfer options: opts % rt_ir
  /// --------------------------------------------------------
  /// Select the solar sea BRDF model: 1=>JONSWAP, 2=>Elfouhaily et al.
  /// default is 1
  oops::OptionalParameter<int> RTTOVSolarSeaBrdfModel{"RTTOV_solar_sea_brdf_model", this};

  /// Select the IR sea surface emissivity model to use.
  /// Valid values: 1 => ISEM; 2 => IREMIS.
  /// default is 2
  oops::OptionalParameter<int> RTTOVIRSeaEmisModel{"RTTOV_ir_sea_emis_model", this};

  /// If true enable solar calculations for solar-affected channels
  /// default is false
  oops::OptionalParameter<bool> RTTOVAddSolar{"RTTOV_addsolar", this};

  /// Vis/IR options. If false disables the Rayleigh single-scattering calculation for
  /// visible channels when addsolar is true (not
  /// recommended to set this to false without good reason).
  /// default is true
  oops::OptionalParameter<bool> RTTOVRayleighSingleScatt{"RTTOV_rayleigh_single_scatt", this};

  /// Vis/IR options. If true includes non-LTE bias correction for hi-res sounders
  /// This is independent of addsolar.
  /// default is false
  oops::OptionalParameter<bool> RTTOVDoNLTECorrection{"RTTOV_do_nlte_correction", this};

  /// Vis/IR options. If true account for scattering due to aerosols.
  /// default in RTTOV12 is false
  oops::OptionalParameter<bool> RTTOVAddAerosl{"RTTOV_addaerosl", this};

  /// Vis/IR options. If true and addaerosl is true the user specifies the aerosol
  /// scattering optical parameters instead of supplying number density
  /// profiles for pre-defined particle types.
  /// default is false
  oops::OptionalParameter<bool> RTTOVUserAerOptParam{"RTTOV_user_aer_opt_param", this};

  /// Vis/IR options. If true account for scattering due to clouds.
  /// default is false
  oops::OptionalParameter<bool> RTTOVAddClouds{"RTTOV_addclouds", this};

  /// Vis/IR options. If true and addaerosl is true the user specifies the aerosol
  /// scattering optical parameters instead of supplying number density
  /// profiles for pre-defined particle types.
  /// default is false
  oops::OptionalParameter<bool> RTTOVUserCldOptParam{"RTTOV_user_cld_opt_param", this};

  /// If true input cloud concentrations should represent grid box
  /// average values (i.e. not accounting for cloud fraction). If false the
  /// cloud concentration should represent the value for the cloudy
  /// fraction of the layer (i.e. the grid box average divided by the layer
  /// cloud fraction). Only applies when addclouds is true and
  /// user_cld_opt_param is false.
  /// default is false
  oops::OptionalParameter<bool> RTTOVGridBoxAvgCloud{"RTTOV_grid_box_avg_cloud", this};

  /// Ignore cloud streams with weights lower than this.
  /// default is -1.0
  oops::OptionalParameter<float> RTTOVCldStrThreshold{"RTTOV_cldstr_threshold", this};

  /// Switch for simplified cloud stream option - USE WITH CAUTION
  /// default is false
  oops::OptionalParameter<bool> RTTOVCldStrSimple{"RTTOV_cldstr_simple", this};

  /// Upper pressure limit for cldstr_simple_option (hPa)
  /// default is 750.0 hPa
  oops::OptionalParameter<float> RTTOVCldStrLowCloudTop{"RTTOV_cldstr_low_cloud_top", this};

  /// Scattering model to use for thermal emission source term: 1 =>
  /// DOM; 2 => Chou-scaling; only applies when
  /// addclouds or addaerosl is true
  /// default is 2
  oops::OptionalParameter<int> RTTOVIRScattModel{"RTTOV_ir_scatt_model", this};

  /// Scattering model to use for solar source term:
  /// 1 => DOM; 2 => single-scattering; only applies
  /// when addclouds or addaerosl is true and addsolar is true.
  /// default is 1
  oops::OptionalParameter<int> RTTOVVisScattModel{"RTTOV_vis_scatt_model", this};

  /// Number of streams (discrete ordinates) to use with DOM
  /// scattering solver. Must be >=2 and even; only
  /// applies when addclouds or addaerosl is true and DOM is
  /// selected as a scattering solver.
  /// default is 8
  oops::OptionalParameter<int> RTTOVDomNstreams{"RTTOV_dom_nstreams", this};

  /// Parameter to determine convergence criterion for DOM
  /// azimuthal loop for solar scattering simulations. If zero or less the
  /// loop never exits early; only applies when
  /// addclouds or addaerosl is true and DOM is selected as a
  /// scattering solver
  /// default is 0.0
  oops::OptionalParameter<float> RTTOVDomAccuracy{"RTTOV_dom_accuracy", this};

  /// DOM ignores layers below this total level-to-space absorption
  /// optical depth. If zero or less all layers are treated in solver;
  /// only applies when addclouds or addaerosl is true
  /// and DOM is selected as a scattering solver.
  /// default is 0.0
  oops::OptionalParameter<float> RTTOVDomOpdepThreshold{"RTTOV_dom_opdep_threshold", this};

  /// If true user is supplying ozone profiles. Relevant
  /// for all VIS/IR sensors, and a limited number of MW sensors
  /// default is false
  oops::OptionalParameter<bool> RTTOVOzoneData{"RTTOV_ozone_data", this};

  /// If true user is supplying CO2 profiles. Currently
  /// only relevant for VIS/IR sensors
  /// default is false
  oops::OptionalParameter<bool> RTTOVCO2Data{"RTTOV_co2_data", this};

  /// If true user is supplying N2O profiles. Currently
  /// only relevant for VIS/IR sensors.
  /// default is false
  oops::OptionalParameter<bool> RTTOVN2OData{"RTTOV_n2o_data", this};

  /// If true user is supplying CO profiles. Currently
  /// only relevant for VIS/IR sensors.
  /// default is false
  oops::OptionalParameter<bool> RTTOVCOData{"RTTOV_co_data", this};

  /// If true user is supplying CH4 profiles. Currently
  /// only relevant for VIS/IR sensors.
  /// default is false
  oops::OptionalParameter<bool> RTTOVCH4Data{"RTTOV_ch4_data", this};

  /// If true user is supplying SO2 profiles. Currently
  /// only relevant for VIS/IR sensors
  /// default is false
  oops::OptionalParameter<bool> RTTOVSO2Data{"RTTOV_so2_data", this};

  /// ----------------------------------------------------------------------
  /// Principal Components-only radiative transfer options: ops % rt_ir % pc
  /// This is being deliberately omitted because the interface is currently
  /// not setup to deal with PC's
  /// ----------------------------------------------------------------------

  /// ----------------------------------------------------------------------
  /// Options related to HTFRTC: ops % htfrtc_opts
  /// This is being deliberately omitted because the interface is currently
  /// not setup to use htfrtc
  /// ----------------------------------------------------------------------

  /// ---------------------------------------------------------------------------
  /// Options related RTTOV-SCATT: opts_scatt
  /// This is being deliberately omitted because the interface is currently
  /// not setup to use RTTOV-SCATT
  /// ---------------------------------------------------------------------------
};

class ObsRadianceRTTOVParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsRadianceRTTOVParameters, ObsOperatorParametersBase)

 public:
  /// The options used by the radiative transfer model.  The rttov settings and some
  /// interface specific settings.
  oops::RequiredParameter<RTTOVObsOptionsParameters> obsOptions{"obs options", this};

  /// Write out detailed information to help with diagnosing problems with the interface and/or
  /// RTTOV
  oops::Parameter<bool> debug{"Debug", false, this};

  /// The type of geoval that is expected by the interface.  Hopefully this can be removed
  /// in the future but for now the oprions are: MetO, SatRad and CRTM.  This is to allow
  /// for Met Office specific setups and comparison with CRTM.
  oops::OptionalParameter<std::string> geovalType{"GeoVal_type", this};

  /// The absorbers needed by RTTOV for this observation type and RTTOV model.
  oops::OptionalParameter<std::vector<std::string>> absorbers{"Absorbers", this};

  /// Options specific to the Linear observation operator
  oops::OptionalParameter<std::vector<std::string>> linearModelAbsorbers{
      "linear model absorbers", this};
};

}  // namespace ufo

#endif  // UFO_RTTOV_OBSRADIANCERTTOVPARAMETERS_H_
