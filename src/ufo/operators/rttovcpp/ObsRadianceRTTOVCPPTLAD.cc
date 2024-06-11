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
#include "ufo/operators/rttovcpp/ObsRadianceRTTOVCPPTLAD.h"
#include "ufo/operators/rttovcpp/rttovcpp_interface.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadianceRTTOVCPPTLAD> makerRTTOVCPPTL_("RTTOVCPP");
// -----------------------------------------------------------------------------

ObsRadianceRTTOVCPPTLAD::ObsRadianceRTTOVCPPTLAD(const ioda::ObsSpace & odb,
                                           const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), varin_()
{
  // Increment fields to be requested from getvalues and stored in geovals
  const std::vector<std::string> vv{
     "air_temperature",
     "specific_humidity"
  };

  for (size_t jvar = 0; jvar < vv.size(); ++jvar) {
     varin_.push_back(vv[jvar]);  // set private data member varin_
  }

  // get channels from observations
  const oops::ObsVariables & observed = odb.assimvariables();
  channels_ = observed.channels();  // set private data member channels_

  // get optical depth coef file name from yaml
  const std::string CoefPath = parameters.CoefPath;
  const std::string SensorID = parameters.SensorID;
  CoefFileName = CoefPath + "rtcoef_" + SensorID + ".dat";

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceRTTOVCPPTLAD::~ObsRadianceRTTOVCPPTLAD() {
  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPPTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                            const QCFlags_t & qc_flags) {
//
  ufo::rttovcpp_interface(geovals, obsspace(), aRttov_, CoefFileName, channels_,
                          nlevels, skip_profile);

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPPTLAD::simulateObsTL(const GeoVaLs & dx, ioda::ObsVector & dy,
                                            const QCFlags_t & qc_flags) const {
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
      dx.getAtLevel(tmpvar2d, oops::Variable{"specific_humidity"}, i);  // get one level Q
      dQ.push_back(tmpvar2d);  // push one level Q into 3D Q
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
    }
  }

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPPTLAD::simulateObsAD(GeoVaLs & dx, const ioda::ObsVector & dy,
                                            const QCFlags_t & qc_flags) const {
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
      dx.getAtLevel(tmpvar2d, oops::Variable{"specific_humidity"}, i);  // get one level Q
      dQ.push_back(tmpvar2d);  // push one level Q into 3D Q
  }

//-------------------------------------------
  ASSERT(dx.nlocs() == dy.nlocs());

  const double missing = util::missingValue<double>();

  std::vector <double> var_k(nlevels, 0.0);

  for (size_t p = 0; p < nprofiles; p++) {
    if (skip_profile[p]) continue;
    for (size_t c = 0; c < nchannels; c++) {
      var_k = aRttov_.getTK(p, c);               // T Jacobian, nlevels
      for (size_t l = 0; l < nlevels; ++l) {
        if (dy[p*nchannels+c] != missing) {
            dT[l][p] += dy[p*nchannels+c] * var_k[l];
        }
      }

      var_k = aRttov_.getItemK(rttov::Q, p, c);  // Q Jacobian, nlevels
      for (size_t l = 0; l < nlevels; l++) {
        if (dy[p*nchannels+c] != missing) {
            dQ[l][p] += dy[p*nchannels+c] * var_k[l];
        }
      }
    }
  }

  // Put temperature increment in kg/kg
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
      dx.putAtLevel(tmpvar2d, oops::Variable{"specific_humidity"}, l);  // put one level Q
  }

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPPTLAD::print(std::ostream & os) const {
  os << "ObsRadianceRTTOVCPPTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
