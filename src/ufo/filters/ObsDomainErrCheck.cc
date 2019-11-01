/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/ObsDomainErrCheck.h"

#include <algorithm>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsDomainErrCheck::ObsDomainErrCheck(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                               boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                               boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr),
    parameter_(0.0)
{
  oops::Log::debug() << "ObsDomainErrCheck: config = " << config_ << std::endl;
  ASSERT(obserr);

  const float missing = util::missingValue(missing);
  float parameter_ = config.getFloat("infltparameter", missing);
  ASSERT(parameter_ != missing);
}

// -----------------------------------------------------------------------------

ObsDomainErrCheck::~ObsDomainErrCheck() {}

// -----------------------------------------------------------------------------

void ObsDomainErrCheck::applyFilter(const std::vector<bool> & inside,
                                    const oops::Variables & filtervars,
                                    std::vector<std::vector<bool>> & flagged) const {
  const oops::Variables observed = obsdb_.obsvariables();

  ioda::ObsDataVector<float> obs(obsdb_, filtervars, "ObsValue");
  size_t nlocs = obsdb_.nlocs();

// compute function
  std::vector<eckit::LocalConfiguration> masks;
  config_.get("where", masks);
  std::vector<float> values(nlocs);
  for (size_t jm = 0; jm < masks.size(); ++jm) {
//  Get variable and group
    const std::string vargrp(masks[jm].getString("variable"));
    std::string fvar, grp;
    std::string obgrp = "MetaData";
    splitVarGroup(vargrp, fvar, grp);
    if (fvar == "Scattering" && grp == "ObsFunction") {
      ioda::ObsDataVector<float> vals(obsdb_, fvar);
      ObsFunction obsdiag(fvar);
      obsdiag.compute(data_, vals);
      for (size_t jj = 0; jj < nlocs; ++jj) {
        values[jj] = vals[fvar][jj];
      }
    }
  }

  size_t count = 0;
  for (size_t jv = 0; jv < filtervars.size(); ++jv) {
    size_t iv = observed.find(filtervars[jv]);
    for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
      if (!inside[jobs]) {
        flagged[jv][jobs] = true;
      } else {
        ASSERT(obserr_[iv][jobs] != util::missingValue(obserr_[iv][jobs]));
        ASSERT(obs[jv][jobs] != util::missingValue(obs[jv][jobs]));
        float bound = 2.5*obserr_[iv][jobs];
        float obserrinc = parameter_ * std::max((values[jobs]-9.0), 0.0) * obserr_[iv][jobs];
        obserrinc = std::max(obserr_[iv][jobs], bound);
        obserr_[iv][jobs] = sqrt(pow(obserr_[iv][jobs], 2) + pow(obserrinc, 2));
        ++count;
        }
      }
  }
  oops::Log::info() << "count=" << count << std::endl;
}

// -----------------------------------------------------------------------------

void ObsDomainErrCheck::print(std::ostream & os) const {
  os << "ObsDomainErrCheck: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
