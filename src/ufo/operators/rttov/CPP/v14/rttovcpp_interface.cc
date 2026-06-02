/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/rttov/CPP/v14/rttovcpp_interface.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"

#include "rttov/wrapper/RttovProfile.h"
#include "rttov/wrapper/RttovSafe.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;
}

namespace ufo {

// -----------------------------------------------------------------------------
/*! \brief Set one single RttovSafe object for a single sensor
*
* \details **rttovcpp_interface()** loads rttov coef files, set profiles and
* surface emissivity/reflectance, call Jacobian function, and perform some quality
* control. It is called by both simulateObs() and setTrajectory().
*
* \param[in] geovals reference to the input full model state at observation locations
* \param[in] odb_ reference to the input observational information
* \param[in] CoefFileName rttov coef file to be loaded
* \param[in] channels_ indexes for a subset of channels for a single sensor
* \param[out] nlevels number of model vertical levels
* \param[out] aRttov_ reference to the output rttov object
* \param[out] skip_profile logical to determine if a profile is used or not
*
* \author Zhiquan (Jake) Liu (NCAR/MMM), initial version for clear-sky DA 
*
* \date Feb. 2021, initial version
*
* \date July 2025, updated from rttov12 to rttov14
* \author Zhiquan (Jake) Liu and Lipeng Jiang (NCAR/MMM)
*
*/
void rttovcpp_interface(const GeoVaLs & geovals, const ioda::ObsSpace & odb_,
                        rttov::RttovSafe & aRttov_, const std::string CoefFileName,
                        const std::vector<int> channels_, std::size_t & nlevels,
                        std::vector<bool> & skip_profile) {
  // 1. Set options for a RttovSafe instance:
  //-----------------------------------------------
  // 1.1 general setting for all sensors: clear-sky

  aRttov_.setFileCoef(CoefFileName);
  aRttov_.options.setVerboseWrapper(true);  // more output info
  aRttov_.options.setVerbose(true);
  aRttov_.options.setCO2Data(false);
  aRttov_.options.setStoreRad(true);
  aRttov_.options.setStoreRad2(true);
  aRttov_.options.setStoreEmisTerms(true);
  // Rttov_.options.setStoreEmisRefl(true);   // Enable emissivity/reflectance retrieval
  // aRttov_.options.setStoreTrans(true);      // Enable Transmittance retrieval
  // aRttov_.options.setStoreDiagOutput(true); // Enable diag retrieval
  aRttov_.options.setCheckProfiles(false);
  aRttov_.options.setHydrometeors(false);
  aRttov_.options.setStoreEmisRefl(true);
  aRttov_.options.setStoreTrans(true);

  // 1.2 for Microwave Sensors
  aRttov_.options.setMwSeaEmisModel(2);  // 2=FASTEM6, 3=SURFEM-Ocean
  aRttov_.options.setUseFoamFraction(false);

  // 1.3 Load coef for subset channels of an instrument
  //-------------------------------------------------
  try {
      aRttov_.loadInst(channels_);
  }
  catch (std::exception& e) {
      oops::Log::error() << "Error loading instrument " << e.what() << std::endl;
  }

  oops::Log::info() << "ObsRadianceRTTOVCPP channels: " << channels_ << std::endl;

  // 2. Allocate profiles
  //---------------------------------------------------------------------------------
  nlevels   = geovals.nlevs(oops::Variable{"air_temperature"});   // set private data member
  std::size_t nprofiles = odb_.nlocs();
  std::size_t nchannels = aRttov_.getNchannels();
  std::size_t nsurfaces = 1;

  std::vector <rttov::Profile> profiles;   // RTTOV Profile object
  for (std::size_t p = 0; p < nprofiles; p++) {
      rttov::Profile aProfile(nlevels+1, nsurfaces);  // rttov14 needs half-level P
      profiles.push_back(aProfile);
  }

  // global scope for it accessable by reflectance* cos(sunzen)
  std::vector<double> sunzen(nprofiles, 0.0);

  // 3. Populate the profiles object
  //---------------------------------------------------------------------------------
  try {
      std::vector<double>  tmpvar1d(nlevels, 0.0);    // one single vertical profile
      std::vector<double>  tmpvar2d(nprofiles, 0.0);  // one single level field

  // 3.1 Common 3D fields needed
  //----------------------------------------------------
    // 3.1.1 Retrieve pressure in hPa
      std::vector<std::vector<double>> tmpvar3d;    // [nlevels][nprofiles]
      std::vector<double> phalf1d(nlevels+1, 0.0);  // one single vertical profile
      std::vector<double> ps(nprofiles, 0.0);

      for (std::size_t i = 0; i < nlevels; ++i) {
         geovals.getAtLevel(tmpvar2d, oops::Variable{"air_pressure_levels"}, i);  // get 1 P_half
         tmpvar3d.push_back(tmpvar2d);  // push one level P into 3D P
      }

      geovals.get(ps, oops::Variable{"air_pressure_at_surface"});  // get one level Ps
      tmpvar3d.push_back(ps);

      for (std::size_t i = 0; i < nprofiles; ++i) {
          for (std::size_t k = 0; k < nlevels+1; ++k) {
            // get one vertical profile, rttov level index is from top to bottom
            phalf1d[k]  = tmpvar3d[k][i]*0.01;
          }
          profiles[i].setPHalf(phalf1d);
      }
     // release memory of tmpvar3d variable
      std::vector<std::vector<double>>().swap(tmpvar3d);

    // 3.1.2 Retrieve temperature in K
      std::vector<std::vector<double>> tmpvar3d_T;  // [nlevels][nprofiles]
      for (std::size_t i = 0; i < nlevels; ++i) {
         geovals.getAtLevel(tmpvar2d, oops::Variable{"air_temperature"}, i);  // get one level T
         tmpvar3d_T.push_back(tmpvar2d);   // push one level T into 3D T
      }
      for (std::size_t i = 0; i < nprofiles; ++i) {
          for (std::size_t k = 0; k < nlevels; ++k) {
              tmpvar1d[k] = tmpvar3d_T[k][i];
          }
          profiles[i].setT(tmpvar1d);
      }
      std::vector<std::vector<double>>().swap(tmpvar3d_T);

    // 3.1.3 Retrieve specific humidity in kg/kg
      std::vector<std::vector<double>> tmpvar3d_Q;  // [nlevels][nprofiles]
      for (std::size_t i = 0; i < nlevels; ++i) {
         geovals.getAtLevel(tmpvar2d, oops::Variable{
             "water_vapor_mixing_ratio_wrt_moist_air"}, i);
         tmpvar3d_Q.push_back(tmpvar2d);
      }
      for (std::size_t i = 0; i < nprofiles; ++i) {
          profiles[i].setGasUnits(rttov::kg_per_kg);
          for (std::size_t k = 0; k < nlevels; ++k) {
              tmpvar1d[k] = tmpvar3d_Q[k][i];
          }
          profiles[i].setQ(tmpvar1d);
      }
      std::vector<std::vector<double>>().swap(tmpvar3d_Q);

    // 3.2 2D surface fields at obs locations
    //-------------------------------------------
      std::vector<double> t2m(nprofiles, 0.0);
      std::vector<double> q2m(nprofiles, 0.0);
      std::vector<double> u10(nprofiles, 0.0);
      std::vector<double> v10(nprofiles, 0.0);
      std::vector<double> tskin(nprofiles, 0.0);
      std::vector<int>    landmask(nprofiles);  // 1: land, 0:ocean
      std::vector<double> seaice_frac(nprofiles, 0.0);
      std::vector<double> elev(nprofiles, 0.0);

    // Retrieve surface variables
      geovals.get(elev, oops::Variable{"geopotential_height_at_surface"});  // in m
      geovals.get(ps, oops::Variable{"air_pressure_at_surface"});  // in Pa, get one level Ps
      geovals.get(t2m, oops::Variable{"air_temperature_at_2m"});  // Kelvin
      geovals.get(q2m, oops::Variable{
          "water_vapor_mixing_ratio_wrt_moist_air_at_2m"});  // kg/kg
      geovals.get(u10, oops::Variable{"eastward_wind_at_10m"});
      geovals.get(v10, oops::Variable{"northward_wind_at_10m"});
      geovals.get(tskin, oops::Variable{"skin_temperature_at_surface"});  // Kelvin
      geovals.get(landmask, oops::Variable{"landmask"});  // 1: land, 0:ocean, var_sfc_landmask
      geovals.get(seaice_frac, oops::Variable{"seaice_fraction"});

    // 3.3 Obs metadata
    //-----------------------------------------------
      std::vector<double> satzen(nprofiles, 0.0);  // always needed
      std::vector<double> satazi(nprofiles, 0.0);  // not always needed
      std::vector<double> sunazi(nprofiles, 0.0);  // not always needed
      std::vector<double> lat(nprofiles, 0.0);
      std::vector<double> lon(nprofiles, 0.0);
      std::vector<util::DateTime> times(nprofiles);

      odb_.get_db("MetaData", "sensorZenithAngle",  satzen);  // in degree
      odb_.get_db("MetaData", "sensorAzimuthAngle", satazi);  // in degree
      odb_.get_db("MetaData", "solarZenithAngle",   sunzen);  // in degree
      odb_.get_db("MetaData", "solarAzimuthAngle",  sunazi);  // in degree
      odb_.get_db("MetaData", "latitude",  lat);  // -90~90 in degree
      odb_.get_db("MetaData", "longitude", lon);  // 0~360 in degree
      odb_.get_db("MetaData", "dateTime", times);

  // 4. Call rttov set functions
  //---------------------------------------------------------------------------------
      util::DateTime time1;
      int year, month, day, hour, minute, second;
      int surftype;
      int isurf = 0;

      for (std::size_t i = 0; i < nprofiles; i++) {
         profiles[i].setGasUnits(rttov::kg_per_kg);

         // convert mpas landmask/xice to rttov surface type
         // may need to make this more generic for different models
         if ( landmask[i] == 0 )      surftype=1;  // sea
         if ( landmask[i] == 1 )      surftype=0;  // land
         if ( seaice_frac[i] >= 0.5 ) surftype=2;  // sea-ice
         profiles[i].setSurfGeom(lat[i], lon[i], 0.001*elev[i]);

         time1 = times[i];
         time1.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);
         profiles[i].setDateTimes(year, month, day, hour, minute, second);

         // 0:land, 1:sea, 2:sea-ice, (sea, fresh water) temporary
         profiles[i].setSurfType(isurf, surftype, 0);

         // isurf, t2m (K), q2m (kg/kg), u10/v10 (m/s), wind fetch
         profiles[i].setNearSurface(isurf, t2m[i], q2m[i], u10[i], v10[i], 100000.);

         // isurf, tskin (k), salinity (35), snow_fraction, foam_fraction, fastem_coef_1-5, spec.
         // over sea/land
         profiles[i].setSkin(isurf, tskin[i], 35., 0., 0., 3.0, 5.0, 15.0, 0.1, 0.3);
         if ( surftype == 2 )  // over seaice, newice(no snow)
           profiles[i].setSkin(isurf, tskin[i], 35., 0., 0., 2.9, 3.4, 27.0, 0.0, 0.0);

         profiles[i].setAngles(std::abs(satzen[i]), satazi[i], sunzen[i], sunazi[i]);
      }
  }  // end try
  catch (std::exception& e) {
      oops::Log::error() << "Error defining the profile data " << e.what() << std::endl;
  }

  // 4.1 Associate the profiles with each RttovSafe instance: the profiles undergo
  //    some checks so use a try block to catch any errors that are thrown up
  try {
      aRttov_.setTheProfiles(profiles);
  }
  catch (std::exception& e) {
      oops::Log::error() << "Error setting the profiles " << e.what() << std::endl;
  }

  // 5. Set the surface emissivity/reflectance arrays
  //    and associate with the Rttov objects
  //--------------------------------------------------
  double surfemisrefl[5][nprofiles][nsurfaces][nchannels];

  aRttov_.setSurfEmisRefl(reinterpret_cast<double *>(surfemisrefl));

// Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
// Negative values will cause RTTOV to supply emissivity/BRDF values (i.e. equivalent to
// calcemis/calcrefl TRUE - see RTTOV user guide)
  for (int j = 0; j < 5; j++) {
      for (int p = 0; p < nprofiles; p++) {
          for (int s = 0; s < nsurfaces; s++) {
              for (int c = 0; c < nchannels; c++) surfemisrefl[j][p][s][c] = -1.;
          }
      }
  }

// 6. Call the RTTOV K model for one instrument for all profiles:
// no arguments are supplied so all 'loaded' channels are simulated
//----------------------------------------------------------------------
  try {
      aRttov_.runK();
  }
  catch (std::exception& e) {
      oops::Log::error() << "Error running RTTOV K model " << e.what() << std::endl;
  }

// 7. Check if Jacobian or Bt has any NaN and set to skip bad profiles
//----------------------------------------------------------------------
  std::vector<double> var_k(nlevels, 0.0);
  std::vector <double> bt;

  for (size_t p = 0; p < nprofiles; p++) skip_profile.push_back(false);

  for (size_t p = 0; p < nprofiles; p++) {
    for (size_t c = 0; c < nchannels; c++) {
      var_k = aRttov_.getTK(p, c);               // T Jacobian for a single profile/channel
      for (size_t l = 0; l < nlevels; ++l) if (std::isnan(var_k[l])) {skip_profile[p] = true;}

      var_k = aRttov_.getItemK(rttov::Q, p, c);  // Q Jacobian for a single profile/channel
      for (size_t l = 0; l < nlevels; ++l) {
          if (std::isnan(var_k[l]) || std::abs(var_k[l]) > 200000.0) {
              skip_profile[p] = true;
          }
      }
    }
  }

  for (size_t p = 0; p < nprofiles; p++) {
    bt = aRttov_.getBtRefl(p);
    for (size_t c = 0; c < nchannels; c++) if (std::isnan(bt[c])) {skip_profile[p] = true;}
  }

  oops::Log::trace() << "rttovcpp_interface done" << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief Save RTTOV diagnostic outputs (YDiag) into ObsDiagnostics
 *
 * \details **rttovcpp_setYdiag()** extracts RTTOV forward and Jacobian
 * diagnostics (e.g., surface emissivity, brightness temperature Jacobians,
 * transmittances) from a single `RttovSafe` object and stores them in the
 * `ObsDiagnostics` container. Variables are allocated only if requested by
 * JEDI yaml configuration, and values are written per channel and per profile.
 * Both 2D diagnostics (surface emissivity, BT clear-sky, Jacobians wrt surface
 * parameters) and 3D diagnostics (layer transmittances, T and Q Jacobian) are supported.
 *
 * \param[in]  geovals   Reference to input model state (`GeoVaLs`) at obs locations
 * \param[in]  aRttov_   Reference to the RTTOV wrapper object already run forward/Jacobian
 * \param[inout] d       ObsDiagnostics object in which YDiag fields are stored
 * \param[in]  channels_ Vector of channel indices for the sensor
 * \param[in]  nlevels   Number of vertical levels in the model/RTTOV setup
 *
 */

  void rttovcpp_setYdiag(const GeoVaLs & geovals, rttov::RttovSafe & aRttov_,
                        ObsDiagnostics & d,
                        const std::vector<int> channels_,
                        std::size_t & nlevels) {
    std::size_t nprofiles = geovals.nlocs();
    std::size_t nchannels = aRttov_.getNchannels();
    // handle ydiag
    // default 2d variables
    const std::vector<std::string> ydiag_varnames_2d{
      "surface_emissivity",
      "brightness_temperature_jacobian_surface_emissivity",
      "brightness_temperature_jacobian_skin_temperature_at_surface",
      "brightness_temperature_assuming_clear_sky",
    };
    // default 3d variables
    const std::vector<std::string> ydiag_varnames_3d{
      "transmittances_of_atmosphere_layer",
      "brightness_temperature_jacobian_air_temperature",
      "brightness_temperature_jacobian_water_vapor_mixing_ratio_wrt_dry_air",
    };
    // combine 2d and 3d variables to ydiag_varnames
    std::vector<std::string> ydiag_varnames = ydiag_varnames_2d;
    ydiag_varnames.insert(ydiag_varnames.end(),
                        ydiag_varnames_3d.begin(),
                        ydiag_varnames_3d.end());
    // get ydiag object
    ufo::GeoVaLs & ydiag = d.geovals();

    // get variables in ydiag; just for debug
    const oops::Variables & vars = ydiag.getVars();
    for (const oops::Variable var : vars) {
      oops::Log::trace() << "ydiag name: " << var << std::endl;
    }

    size_t iSurface = 0;  // support only one surface type for now.
    for (const std::string &varname : ydiag_varnames) {
      // check varname dimension
      const bool is3d = std::find(ydiag_varnames_3d.begin(), ydiag_varnames_3d.end(), varname)
                      != ydiag_varnames_3d.end();
      const bool is2d = std::find(ydiag_varnames_2d.begin(), ydiag_varnames_2d.end(), varname)
                      != ydiag_varnames_2d.end();
      int this_nlevels = 0;
      if (is3d) {
        this_nlevels = nlevels;
      } else if (is2d) {
        this_nlevels = 1;
      }

      bool save_ydiag = true;
      // allocate ydiag for all channels
      for (std::size_t ich = 0; ich < nchannels; ich++) {
        int ch = channels_[ich];
        oops::Variables this_var({varname+"_"+std::to_string(ch)});
        if (ydiag.has(this_var[0])) {
          ydiag.allocate(this_nlevels, this_var);
        } else {
          // this ydiag is not needed according to JEDI Yaml settings.
          // no ydiag related QC or BC procedures; no ydiag output request.
          save_ydiag = false;
        }
      }
      // save RTTOV output to ydiag
      if (save_ydiag) {
        std::vector<double> res;
        if (is2d) {
          std::vector<double> prof_1(1, 0.0);
          for (size_t p = 0; p < nprofiles; p++) {
            if (varname.find("brightness_temperature_assuming_clear_sky") != std::string::npos) {
              res = aRttov_.getBtClear(p);
            } else if (varname.find("brightness_temperature_jacobian_surface_emissivity") !=
                       std::string::npos) {
              res = aRttov_.getSurfEmisK(p, iSurface);
              oops::Log::trace() << "emisK= " << res   << std::endl;
            } else if (varname.find("brightness_temperature_jacobian_skin_temperature_at_surface")
                       != std::string::npos) {
              res = aRttov_.getTskinEffK(p, iSurface);
              oops::Log::trace() << "TskinEffK= " << res   << std::endl;
            } else if (varname.find("surface_emissivity") != std::string::npos) {
              res = aRttov_.getSurfEmis(p, iSurface);
              oops::Log::trace() << "SurfEmis= " << res   << std::endl;
            }
            // set values for each channel
            for (std::size_t ich = 0; ich < nchannels; ich++) {
              const int ch = channels_[ich];
              oops::Variables this_var({varname+"_"+std::to_string(ch)});
              prof_1[0] = res[ich];
              ydiag.putProfile(prof_1, this_var[0], p);
            }
          }
        } else if (is3d) {
          std::vector<double> prof_n(nlevels, 0.0);
          for (int p = 0; p < nprofiles; p++) {
            for (int ich = 0; ich < nchannels; ich++) {
              const int ch = channels_[ich];
              oops::Variables this_var({varname+"_"+std::to_string(ch)});
              if (varname.find("transmittances_of_atmosphere_layer") != std::string::npos) {
                res = aRttov_.getTauLevels(p, ich);
                std::vector<double> subprof(res.begin(), res.end() - 1);
                ydiag.putProfile(subprof, this_var[0], p);
              } else if (varname.find("brightness_temperature_jacobian_air_temperature") !=
                         std::string::npos) {
                res = aRttov_.getTK(p, ich);
                std::vector<double> subprof(res.begin(), res.end());
                oops::Log::trace() << "BtJacbAirTemp= " << res   << std::endl;
                ydiag.putProfile(subprof, this_var[0], p);
              } else if (varname.find("brightness_temperature_jacobian_"
                                     "water_vapor_mixing_ratio_wrt_dry_air") != std::string::npos) {
                  res = aRttov_.getItemK(rttov::Q, p, ich);
                  std::vector<double> subprof(res.begin(), res.end());
                  ydiag.putProfile(subprof, this_var[0], p);
              }
            }
          }
        }
      }
    }
  }
  // -----------------------------------------------------------------------------

}  // namespace ufo
