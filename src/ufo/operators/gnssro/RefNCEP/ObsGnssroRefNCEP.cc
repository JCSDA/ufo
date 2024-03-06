/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/gnssro/RefNCEP/ObsGnssroRefNCEP.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroRefNCEP> makerGnssroRefNCEP_("GnssroRefNCEP");
// -----------------------------------------------------------------------------

ObsGnssroRefNCEP::ObsGnssroRefNCEP(const ioda::ObsSpace & odb, const Parameters_ & params)
  : ObsOperatorBase(odb), keyOperGnssroRefNCEP_(0), odb_(odb), varin_()
{
  const std::vector<std::string> vv{"air_temperature", "specific_humidity", "air_pressure",
                                    "geopotential_height"};
  varin_.reset(new oops::Variables(vv));

  ufo_gnssro_refncep_setup_f90(keyOperGnssroRefNCEP_,
                               params.options.value().toConfiguration());
  oops::Log::trace() << "ObsGnssroRefNCEP created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroRefNCEP::~ObsGnssroRefNCEP() {
  ufo_gnssro_refncep_delete_f90(keyOperGnssroRefNCEP_);
  oops::Log::trace() << "ObsGnssroRefNCEP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNCEP::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                               ObsDiagnostics &, const QCFlags_t & qc_flags) const {
  ufo_gnssro_refncep_simobs_f90(keyOperGnssroRefNCEP_, gom.toFortran(), odb_,
                            ovec.size(), ovec.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNCEP::print(std::ostream & os) const {
  os << "ObsGnssroRefNCEP::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
