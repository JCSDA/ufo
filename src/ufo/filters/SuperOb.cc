/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/missingValues.h"

#include "ufo/filters/SuperOb.h"
#include "ufo/superob/SuperObBase.h"

namespace ufo {

// -----------------------------------------------------------------------------

SuperOb::SuperOb(ioda::ObsSpace & obsdb,
                 const Parameters_ & parameters,
                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr),
    options_(parameters) {
  oops::Log::trace() << "SuperOb constructor" << std::endl;

  // The ObsSpace must have been divided into records.
  if (obsdb.obs_group_vars().empty()) {
    throw eckit::BadParameter("group variables configuration is empty.", Here());
  }

  allvars_ += Variables(filtervars_, "HofX");

  oops::Log::debug() << "SuperOb constructor finished" << std::endl;
}

SuperOb::~SuperOb() {
  oops::Log::trace() << "SuperOb destructor" << std::endl;
}

void SuperOb::applyFilter(const std::vector<bool> & apply,
                          const Variables & filtervars,
                          std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "SuperOb applyFilter" << std::endl;

  const SuperObParametersWrapper & params = options_.algorithmParameters.value();

  std::unique_ptr<SuperObBase> superOb =
    SuperObFactory::create(params.superObName,
                           data_, apply, filtervars, *flags_, flagged);

  superOb->runAlgorithm();
}

void SuperOb::print(std::ostream & os) const {
  os << "SuperOb: config = " << options_ << std::endl;
}

}  // namespace ufo
