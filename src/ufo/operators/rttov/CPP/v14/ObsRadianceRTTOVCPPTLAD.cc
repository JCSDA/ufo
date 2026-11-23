/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/rttov/CPP/ObsRadianceRTTOVCPPTLAD.h"
#include "ufo/operators/rttov/CPP/v14/rttovcpp_interface.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadianceRTTOVCPPTLAD> makerRTTOVCPPTL_("RTTOVCPP");
// -----------------------------------------------------------------------------

ObsRadianceRTTOVCPPTLAD::ObsRadianceRTTOVCPPTLAD(const ioda::ObsSpace & odb,
                                           const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), varin_()
{
  useHydrometeors_ = parameters.useHydrometeors.value();
  hydrometeorVars_ = parameters.hydrometeors.value();

  // 1. Always push the base clear-sky variables
  const std::vector<std::string> vv_clr{
     "air_temperature",
     "water_vapor_mixing_ratio_wrt_moist_air"
  };

  for (const auto& var : vv_clr) {
    varin_.push_back(var);
  }

  // 2. Conditionally push the extra variables if all-sky is enabled
  //    (cloud fraction is an RTTOV input only, not an increment variable)
  if (useHydrometeors_) {
    for (const auto& var : hydrometeorVars_) {
      if (var == "cloud_area_fraction_in_atmosphere_layer") continue;
      varin_.push_back(var);
    }
  }

  // get channels from observations
  const oops::ObsVariables & observed = odb.assimvariables();
  channels_ = observed.channels();  // set private data member channels_

  // get optical depth coef file name from yaml
  const std::string CoefPath = parameters.CoefPath;
  const std::string SensorID = parameters.SensorID;
  const std::string sep = (!CoefPath.empty() && CoefPath.back() == '/') ? "" : "/";
  CoefFileName = CoefPath + sep + "rtcoef_" + SensorID + ".dat";
  HydrotableName = CoefPath + sep + "rttov_hydrotable_" + SensorID + ".dat";

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceRTTOVCPPTLAD::~ObsRadianceRTTOVCPPTLAD() {
  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPPTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics & d,
                                            const QCFlags_t & qc_flags) {
  ufo::rttovcpp_interface(geovals, obsspace(), aRttov_, useHydrometeors_, hydrometeorVars_,
                          CoefFileName, HydrotableName, channels_,
                          nlevels, skip_profile);
  ufo::rttovcpp_setYdiag(geovals, aRttov_, d, channels_, nlevels, useHydrometeors_,
                          hydrometeorVars_);

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPPTLAD::simulateObsTL(const GeoVaLs & dx, ioda::ObsVector & dy) const {
//
  std::size_t nprofiles = dy.nlocs();
  std::size_t nchannels = aRttov_.getNchannels();

  std::vector<double>  tmpvar2d(nprofiles, 0.0);  // one single level field
  std::vector<std::vector<double>> dT;  // [nlevels][nprofiles]
  std::vector<std::vector<double>> dQ;  // [nlevels][nprofiles]

  // Retrieve temperature increment in K
  for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{"air_temperature"}, i);  // get one level T
      dT.push_back(tmpvar2d);  // push one level T into 3D T
  }

  // Retrieve specific humidity increment in kg/kg
  for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{
          "water_vapor_mixing_ratio_wrt_moist_air"}, i);  // get one level Q
      dQ.push_back(tmpvar2d);  // push one level Q into 3D Q
  }

  // Retrieve hydrometeor increments in kg/kg
  const std::vector<RttovHydrometeor> hydros =
      rttovcpp_activeHydrometeors(useHydrometeors_, hydrometeorVars_);
  std::vector<std::vector<std::vector<double>>> dHydro;  // [nhydro][nlevels][nprofiles]
  for (const RttovHydrometeor & hydro : hydros) {
    std::vector<std::vector<double>> singleVarData;
    for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{hydro.varName}, i);
      singleVarData.push_back(tmpvar2d);
    }
    dHydro.push_back(singleVarData);
  }

//-------------------------------------------
  ASSERT(dx.nlocs() == dy.nlocs());
  ASSERT(nchannels == dy.nvars());
  dy.zero();

  std::vector <double> var_k(nlevels, 0.0);

  for (size_t p = 0; p < nprofiles; p++) {
    if (skip_profile[p]) continue;
    for (size_t c = 0; c < nchannels; c++) {
      var_k = aRttov_.getTK(p, c);              // T Jacobian for a single profile/channel
      for (size_t l = 0; l < nlevels; l++)
          dy[p*nchannels+c] += var_k[l]*dT[l][p];

      var_k = aRttov_.getItemK(rttov::Q, p, c);   // Q Jacobian
      for (size_t l = 0; l < nlevels; l++)
          dy[p*nchannels+c] += var_k[l]*dQ[l][p];

      for (size_t v = 0; v < hydros.size(); v++) {   // Hydrometeor Jacobians
        var_k = aRttov_.getItemK(hydros[v].itemId, p, c);
        for (size_t l = 0; l < nlevels; l++)
            dy[p*nchannels+c] += var_k[l]*dHydro[v][l][p];
      }
    }
  }

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

  void ObsRadianceRTTOVCPPTLAD::simulateObsAD(GeoVaLs & dx, const ioda::ObsVector & dy) const {
  std::size_t nprofiles = dy.nlocs();
  std::size_t nchannels = aRttov_.getNchannels();

  std::vector<double>  tmpvar2d(nprofiles, 0.0);  // one single level field
  std::vector<std::vector<double>> dT;  // [nlevels][nprofiles]
  std::vector<std::vector<double>> dQ;  // [nlevels][nprofiles]

  // Retrieve temperature increment in K
  for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{"air_temperature"}, i);  // get one level T
      dT.push_back(tmpvar2d);  // push one level T into 3D T
  }

  // Retrieve specific humidity increment in kg/kg
  for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{
          "water_vapor_mixing_ratio_wrt_moist_air"}, i);  // get one level Q
      dQ.push_back(tmpvar2d);  // push one level Q into 3D Q
  }

  // Hydrometeor adjoint accumulators [nhydro][nlevels][nprofiles]
  const std::vector<RttovHydrometeor> hydros =
      rttovcpp_activeHydrometeors(useHydrometeors_, hydrometeorVars_);
  std::vector<std::vector<std::vector<double>>> dHydro(hydros.size(),
      std::vector<std::vector<double>>(nlevels, std::vector<double>(nprofiles, 0.0)));

//-------------------------------------------
  ASSERT(dx.nlocs() == dy.nlocs());

  const double missing = util::missingValue<double>();

  std::vector <double> var_k(nlevels, 0.0);

  for (size_t p = 0; p < nprofiles; p++) {
    if (skip_profile[p]) continue;
    for (size_t c = 0; c < nchannels; c++) {
      double dy_val = dy[p*nchannels+c];

      if (dy_val == missing) continue;

      var_k = aRttov_.getTK(p, c);               // T Jacobian, nlevels
      for (size_t l = 0; l < nlevels; ++l) dT[l][p] += var_k[l] * dy_val;

      var_k = aRttov_.getItemK(rttov::Q, p, c);  // Q Jacobian, nlevels
      for (size_t l = 0; l < nlevels; l++) dQ[l][p] += var_k[l] * dy_val;

      for (size_t v = 0; v < hydros.size(); v++) {  // Hydrometeor Jacobians, nlevels
        var_k = aRttov_.getItemK(hydros[v].itemId, p, c);
        for (size_t l = 0; l < nlevels; l++) dHydro[v][l][p] += var_k[l] * dy_val;
      }
    }
  }

  // Put temperature increment in K
  for (std::size_t l = 0; l < nlevels; ++l) {
      for (size_t p = 0; p < nprofiles; p++) {
          tmpvar2d[p] = dT[l][p];
      }
      dx.putAtLevel(tmpvar2d, oops::Variable{"air_temperature"}, l);  // put one level T
  }

  // Put specific humidity increment in kg/kg
  for (std::size_t l = 0; l < nlevels; ++l) {
      for (size_t p = 0; p < nprofiles; p++) {
          tmpvar2d[p] = dQ[l][p];
      }
      dx.putAtLevel(tmpvar2d, oops::Variable{
          "water_vapor_mixing_ratio_wrt_moist_air"}, l);  // put one level Q
  }

  // Put hydrometeor increments in kg/kg
  for (std::size_t v = 0; v < hydros.size(); ++v) {
    for (std::size_t l = 0; l < nlevels; ++l) {
      dx.getAtLevel(tmpvar2d, oops::Variable{hydros[v].varName}, l);
      for (std::size_t p = 0; p < nprofiles; ++p) tmpvar2d[p] += dHydro[v][l][p];
      dx.putAtLevel(tmpvar2d, oops::Variable{hydros[v].varName}, l);
    }
  }

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPPTLAD::print(std::ostream & os) const {
  os << "ObsRadianceRTTOVCPPTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
