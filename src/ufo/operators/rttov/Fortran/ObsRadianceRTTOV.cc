/*
 * (C) Copyright 2026 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/rttov/Fortran/ObsRadianceRTTOV.h"

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceRTTOV> makerRTTOV_("RTTOV");

// -----------------------------------------------------------------------------

ObsRadianceRTTOV::ObsRadianceRTTOV(const ioda::ObsSpace & odb,
                                   const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOperRadianceRTTOV_(0), odb_(odb), varin_()
{
  oops::Log::trace() << "ObsRadianceRTTOV constructor start." << std::endl;
  // parse channels from the config and create variable names
  const oops::ObsVariables & observed = odb.assimvariables();
  std::vector<int> channels_list = observed.channels();

  // values for passed qc
  std::vector<int> qc_passed = {QCflags::pass, QCflags::passive};

  // call Fortran setup routine
  ufo_radiancerttov_setup_f90(keyOperRadianceRTTOV_, parameters.toConfiguration(),
                              channels_list.size(), channels_list[0], varin_,
                              qc_passed.size(), qc_passed[0]);

  // Remove ozone from varin_ if calculate from ref is switched on
  if (parameters.obsOptions.value().RTTOVScaleRefOzone.value())
      varin_ -= oops::Variable{"mole_fraction_of_ozone_in_air"};

  // Sanity check the provided vector if inspectProfileLatLonBox is used
  // latlonbox = [latmin, latmax, lonmin, lonmax]
  if (parameters.obsOptions.value().inspectProfileLatLonBox.value().has_value()) {
      std::vector<float> latlonbox =
          parameters.obsOptions.value().inspectProfileLatLonBox.value().get();
      assert(size(latlonbox) == 4);
      assert(latlonbox[0] < latlonbox[1]);
      assert(latlonbox[2] < latlonbox[3]);
  }

  oops::Log::info() << "ObsRadianceRTTOV channels: " << channels_list << std::endl;
  oops::Log::trace() << "ObsRadianceRTTOV constructor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceRTTOV::~ObsRadianceRTTOV() {
  ufo_radiancerttov_delete_f90(keyOperRadianceRTTOV_);
  oops::Log::trace() << "ObsRadianceRTTOV destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOV::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                  ObsDiagnostics & dvec, const QCFlags_t& qc_flags) const {
  ufo_radiancerttov_simobs_f90(keyOperRadianceRTTOV_, gom.toFortran(), odb_,
                               ovec.nvars(), ovec.nlocs(), ovec.toFortran(),
                               dvec.toFortran(), qc_flags);
  oops::Log::trace() << "ObsRadianceRTTOV::simulateObs done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOV::print(std::ostream & os) const {
  os << "ObsRadianceRTTOV::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
