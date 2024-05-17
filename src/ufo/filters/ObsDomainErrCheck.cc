/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/ObsDomainErrCheck.h"

#include <algorithm>
#include <string>
#include <type_traits>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/ObsFunctionScattering.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsDomainErrCheck::ObsDomainErrCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                               std::shared_ptr<ioda::ObsDataVector<int> > flags,
                               std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr),
    parameters_(parameters)
{
  oops::Log::debug() << "ObsDomainErrCheck: config = " << parameters_ << std::endl;
  ASSERT(obserr);
}

// -----------------------------------------------------------------------------

ObsDomainErrCheck::~ObsDomainErrCheck() {}

// -----------------------------------------------------------------------------

void ObsDomainErrCheck::applyFilter(const std::vector<bool> & inside,
                                    const Variables & filtervars,
                                    std::vector<std::vector<bool>> & flagged) const {
  const oops::ObsVariables observed = obsdb_.obsvariables();

  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsObsVariables(), "ObsValue");
  size_t nlocs = obsdb_.nlocs();

// compute function
  std::vector<float> values(nlocs);
  ioda::ObsDataVector<float> vals(obsdb_, "Scattering");
  ObsFunctionScattering obsdiag;
  obsdiag.compute(data_, vals);
  for (size_t jj = 0; jj < nlocs; ++jj) {
    values[jj] = vals["Scattering"][jj];
  }

  size_t count = 0;
  const float parameter = parameters_.infltparameter;
  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    size_t iv = observed.find(filtervars.variable(jv).variable());
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (!inside[jobs]) {
        flagged[jv][jobs] = true;
      } else {
        // style note: use missingValue<decltype(foo)> because obserr_ is declared in another file
        // and its underlying type could in principle change without it being obvious here.
        ASSERT((*obserr_)[iv][jobs] != util::missingValue<std::remove_reference_t<
                                                          decltype((*obserr_)[iv][jobs])>>());
        ASSERT(obs[jv][jobs] != util::missingValue<float>());
        float bound = 2.5 * (*obserr_)[iv][jobs];
        float obserrinc = parameter * std::max((values[jobs]-9.0), 0.0) * (*obserr_)[iv][jobs];
        obserrinc = std::max((*obserr_)[iv][jobs], bound);
        (*obserr_)[iv][jobs] = sqrt(pow((*obserr_)[iv][jobs], 2) + pow(obserrinc, 2));
        ++count;
        }
      }
  }
  oops::Log::info() << "count=" << count << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDomainErrCheck::print(std::ostream & os) const {
  os << "ObsDomainErrCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
