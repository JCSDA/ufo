/*
 * (C) Copyright 2017-2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/rttovcpp/ObsRadianceRTTOVCPP.h"
#include "ufo/operators/rttovcpp/rttovcpp_interface.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceRTTOVCPP> makerRTTOVCPP_("RTTOVCPP");

// -----------------------------------------------------------------------------

ObsRadianceRTTOVCPP::ObsRadianceRTTOVCPP(const ioda::ObsSpace & odb,
                                         const Parameters_ & parameters)
  : ObsOperatorBase(odb), odb_(odb), varin_()
{
  // Fields to be requested from getvalues and stored in geovals
  // need to be consistent with those defined in ufo_variables_mod.F90
  //-----------------------------------------------------------------------------
  const std::vector<std::string> vv{
    "air_pressure",
    "air_temperature",
    "specific_humidity",
    "surface_pressure",
    "surface_temperature",   // this is actually var_sfc_t2m
    "specific_humidity_at_two_meters_above_surface",
    "uwind_at_10m",
    "vwind_at_10m",
    "skin_temperature",
    "seaice_fraction",
    "landmask",
    "surface_geopotential_height"
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
  oops::Log::info() << CoefFileName << std::endl;

  oops::Log::trace() << "ObsRadianceRTTOVCPP created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceRTTOVCPP::~ObsRadianceRTTOVCPP() {
  oops::Log::trace() << "ObsRadianceRTTOVCPP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPP::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & hofx,
                                   ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
//
  std::vector<bool>  skip_profile;
  ufo::rttovcpp_interface(geovals, odb_, aRttov_, CoefFileName, channels_,
                          nlevels, skip_profile);

// ------------------------------------------------------------------------
// Obtain calculated brightness temperature for all profiles/channels
// ------------------------------------------------------------------------
  std::size_t nprofiles = geovals.nlocs();
  std::size_t nchannels = aRttov_.getNchannels();

  ASSERT(geovals.nlocs() == hofx.nlocs());
  hofx.zero();  // this may not be necessary

  const double missing = util::missingValue<double>();

  std::vector <double> bt;

  for (size_t p = 0; p < nprofiles; p++) {
      for (size_t c = 0; c < nchannels; c++) hofx[p*nchannels+c] = missing;
      if (skip_profile[p]) continue;
      bt = aRttov_.getBtRefl(p);
      for (size_t c = 0; c < nchannels; c++) hofx[p*nchannels+c] = bt[c];
  }

  oops::Log::trace() << "ObsRadianceRTTOVCPP::simulateObs done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPP::print(std::ostream & os) const {
  os << "ObsRadianceRTTOVCPP::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
