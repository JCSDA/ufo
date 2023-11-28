/*
 * (C) Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

/* 1-D var qc
 *   J(x) = (x-xb)T B-1 (x-xb) + (y-H(x))T R-1 (y-H(x))
 *   Code adapted from Met Office OPS System
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <vector>

#include "ufo/filters/gnssroonedvarcheck/GNSSROOneDVarCheck.h"
#include "ufo/filters/gnssroonedvarcheck/GNSSROOneDVarCheck.interface.h"
#include "ufo/GeoVaLs.h"

#include "oops/util/IntSetParser.h"

namespace ufo {

// -----------------------------------------------------------------------------

GNSSROOneDVarCheck::GNSSROOneDVarCheck(ioda::ObsSpace & obsdb,
                                       const Parameters_ & parameters,
                                       std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                       std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::debug() << "GNSSROOneDVarCheck contructor starting" << std::endl;

  std::set<int> channelset = oops::parseIntSet(parameters_.channelList);
  std::vector<int> channels;
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels));

  // Setup fortran object
  ufo_gnssroonedvarcheck_create_f90(key_,
                                    obsdb,
                                    parameters_.bmatrix_filename.value().size(),
                                    parameters_.bmatrix_filename.value().c_str(),
                                    parameters_.capsupersat.value(),
                                    parameters_.cost_funct_test.value(),
                                    parameters_.Delta_ct2.value(),
                                    parameters_.Delta_factor.value(),
                                    parameters_.min_temp_grad.value(),
                                    parameters_.n_iteration_test.value(),
                                    parameters_.OB_test.value(),
                                    parameters_.pseudo_ops.value(),
                                    parameters_.vert_interp_ops.value(),
                                    parameters_.y_test.value(),
                                    GNSSROOneDVarCheck::qcFlag(),
                                    channels.size(),
                                    channels[0],
                                    parameters_.noSuperCheck.value());

  oops::Log::debug() << "GNSSROOneDVarCheck contructor complete. " << std::endl;
}

// -----------------------------------------------------------------------------

GNSSROOneDVarCheck::~GNSSROOneDVarCheck() {
  ufo_gnssroonedvarcheck_delete_f90(key_);
  oops::Log::trace() << "GNSSROOneDVarCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void GNSSROOneDVarCheck::applyFilter(const std::vector<bool> & apply,
                               const Variables & filtervars,
                               std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "GNSSROOneDVarCheck Filter starting" << std::endl;

  // Get GeoVaLs
  const ufo::GeoVaLs * gvals = data_.getGeoVaLs();

  // Convert apply to char for passing to fortran
  std::vector<char> apply_char(apply.size(), 'F');
  for (size_t i = 0; i < apply_char.size(); i++) {
    if (apply[i]) {apply_char[i]='T';}
  }

  // Save qc flags to database for retrieval in fortran - needed for channel selection
  flags_->save("FortranQC");      // temporary measure as per ROobserror qc
  obserr_->save("FortranERR");    // Pass latest errors to 1DVar

  // Pass it all to fortran
  ufo_gnssroonedvarcheck_apply_f90(key_,
                                  gvals->toFortran(),
                                  apply_char.size(), apply_char[0]);

  // Read qc flags from database
  flags_->read("FortranQC");    // temporary measure as per ROobserror qc

  oops::Log::trace() << "GNSSROOneDVarCheck Filter complete" << std::endl;
}

// -----------------------------------------------------------------------------

void GNSSROOneDVarCheck::print(std::ostream & os) const {
  os << "GNSSROOneDVarCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
