/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/PreQC.h"

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"

namespace ufo {

PreQC::PreQC(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
             std::shared_ptr<ioda::ObsDataVector<int> > flags,
             std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::debug() << "PreQC: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

void PreQC::applyFilter(const std::vector<bool> & apply,
                        const Variables & filtervars,
                        std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "PreQC applyFilter starting " << std::endl;

  const int missing = util::missingValue<int>();

  // Read QC flags from pre-processing
  ioda::ObsDataVector<int> preqc(obsdb_,
                                 filtervars.toOopsObsVariables(),
                                 parameters_.inputQC);
  oops::Log::debug() << "PreQC::PreQC preqc: " << preqc;

  // Get min and max values and reject outside range
  const int qcmin = parameters_.minvalue;
  const int qcmax = parameters_.maxvalue;

  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    const ioda::ObsDataRow<int> &currentPreQC = preqc[jv];
    std::vector<bool> &currentFlagged = flagged[jv];
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (apply[jobs] &&
          (currentPreQC[jobs] == missing ||
           currentPreQC[jobs] > qcmax ||
           currentPreQC[jobs] < qcmin)) {
        currentFlagged[jobs] = true;
      }
    }
  }

  oops::Log::trace() << "PreQC applyFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

void PreQC::print(std::ostream & os) const {
  os << "PreQC: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
