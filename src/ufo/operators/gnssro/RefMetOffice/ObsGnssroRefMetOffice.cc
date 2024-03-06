/*
 * (C) British Crown Copyright 2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/gnssro/RefMetOffice/ObsGnssroRefMetOffice.h"
#include "ufo/operators/gnssro/RefMetOffice/ObsGnssroRefMetOfficeParameters.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroRefMetOffice> makerGnssroRefMetOffice_("GnssroRefMetOffice");
// -----------------------------------------------------------------------------

ObsGnssroRefMetOffice::ObsGnssroRefMetOffice(const ioda::ObsSpace & odb,
                                             const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOperGnssroRefMetOffice_(0), odb_(odb), varin_(),
    parameters_(parameters)
{
  ObsGnssroRefMetOfficeOptions obsOptions = parameters_.obsOptions.value();

  const std::vector<std::string> vv{"air_pressure_levels", "specific_humidity",
                                    "geopotential_height", "geopotential_height_levels"};
  varin_.reset(new oops::Variables(vv));

  ufo_gnssro_refmetoffice_setup_f90(keyOperGnssroRefMetOffice_,
                                    obsOptions.vertInterpOPS,
                                    obsOptions.pseudoLevels,
                                    obsOptions.minTempGrad);

  oops::Log::trace() << "ObsGnssroRefMetOffice created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroRefMetOffice::~ObsGnssroRefMetOffice() {
  ufo_gnssro_refmetoffice_delete_f90(keyOperGnssroRefMetOffice_);
  oops::Log::trace() << "ObsGnssroRefMetOffice destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRefMetOffice::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                     ObsDiagnostics & ydiags, const QCFlags_t & qc_flags) const {
  ufo_gnssro_refmetoffice_simobs_f90(keyOperGnssroRefMetOffice_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran(), ydiags.toFortran());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefMetOffice::print(std::ostream & os) const {
  os << "ObsGnssroRefMetOffice: config = " <<  parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
