/*
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/scatwind/NeutralMetOffice/ObsScatwindNeutralMetOffice.h"

#include <ostream>

#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsScatwindNeutralMetOffice>
  makerScatwindNeutralMetOffice_("ScatwindNeutralMetOffice");
// -----------------------------------------------------------------------------

ObsScatwindNeutralMetOffice::ObsScatwindNeutralMetOffice(const ioda::ObsSpace & odb,
                                                         const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOperScatwindNeutralMetOffice_(0),
    odb_(odb), varin_()
{
  ufo_scatwind_neutralmetoffice_setup_f90(keyOperScatwindNeutralMetOffice_,
                                          odb.assimvariables(), varin_);

  oops::Log::trace() << "ObsScatwindNeutralMetOffice created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsScatwindNeutralMetOffice::~ObsScatwindNeutralMetOffice() {
  ufo_scatwind_neutralmetoffice_delete_f90(keyOperScatwindNeutralMetOffice_);
  oops::Log::trace() << "ObsScatwindNeutralMetOffice destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsScatwindNeutralMetOffice::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                              ObsDiagnostics &) const {
    oops::Log::trace() << "ObsScatwindNeutralMetOffice::simulateObs entered" << std::endl;

    ufo_scatwind_neutralmetoffice_simobs_f90(keyOperScatwindNeutralMetOffice_,
                                             gom.toFortran(), odb_,
                                             ovec.nvars(), ovec.nlocs(), ovec.toFortran());

    oops::Log::trace() << "ObsScatwindNeutralMetOffice::simulateObs exit" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsScatwindNeutralMetOffice::print(std::ostream & os) const {
  os << "ObsScatwindNeutralMetOffice::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
