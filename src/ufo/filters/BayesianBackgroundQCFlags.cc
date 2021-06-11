/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/BayesianBackgroundQCFlags.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"

namespace ufo {

BayesianBackgroundQCFlags::BayesianBackgroundQCFlags
(ioda::ObsSpace & obsdb,
 const Parameters_ & parameters,
 std::shared_ptr<ioda::ObsDataVector<int> > flags,
 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "BayesianBackgroundQCFlags constructor" << std::endl;
}

BayesianBackgroundQCFlags::~BayesianBackgroundQCFlags() {
  oops::Log::trace() << "BayesianBackgroundQCFlags destructor" << std::endl;
}

std::string BayesianBackgroundQCFlags::getPGEsubstituteName(const std::string &varname) const {
  const auto PGEsubstituteNames = parameters_.PGEsubstituteNames.value();
  const auto it = PGEsubstituteNames.find(varname);
  if (it != PGEsubstituteNames.end())
    return it->second;
  else
    return varname;
}

void BayesianBackgroundQCFlags::setFlags(const std::string& varname,
                                         const std::vector<bool>& apply,
                                         std::vector<bool>& flagged) const {
  const float missingValueFloat = util::missingValue(missingValueFloat);
  const int missingValueInt = util::missingValue(missingValueInt);
  const size_t nlocs = obsdb_.nlocs();
  // PGE multiplication factor used to store PGE values for later use.
  const float PGEMult = 1000.0;
  // Missing data indicator for packed PGEs (for compatibility with OPS).
  const float PGEMDI = 1.111;
  // PGE rejection limit.
  const float PGECrit = parameters_.PGEParameters.PGE_PGECrit.value();
  // Packed PGE rejection limit.
  const float packedPGECrit = PGEMult * PGECrit;

  // Sometimes the PGE of one variable is used to set the QC flags of another;
  // this happens for (e.g.) wind u and v components.
  // By default varPGE = varname.
  const std::string varPGE = getPGEsubstituteName(varname);

  // Get QC flags.
  std::vector <int> QCflags(nlocs, missingValueInt);
  if (obsdb_.has("QCFlags", varname))
    obsdb_.get_db("QCFlags", varname, QCflags);
  else
    throw eckit::BadValue(varname + "@QCFlags not present", Here());

  // Get the PGE values that each observation had after
  // the Bayesian background and buddy checks were applied.
  // If the checks were not applied, default to using
  // the existing PGE values.
  if (obsdb_.has("GrossErrorProbabilityBuddyCheck", varPGE)) {
    // Buddy check was applied.
    std::vector <float> buddyCheckPGEs(nlocs, missingValueFloat);
    obsdb_.get_db("GrossErrorProbabilityBuddyCheck", varPGE, buddyCheckPGEs);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (!apply[iloc]) continue;
      if (buddyCheckPGEs[iloc] == PGEMDI) {
        QCflags[iloc] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
        flagged[iloc] = true;
      } else if (buddyCheckPGEs[iloc] >= PGECrit) {
        QCflags[iloc] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
        QCflags[iloc] |= ufo::MetOfficeQCFlags::Elem::BuddyRejectFlag;
        flagged[iloc] = true;
      }
    }
  } else if (obsdb_.has("GrossErrorProbability", varPGE)) {
    // Get the PGE values that each observation had before
    // the Bayesian background and buddy checks were applied.
    std::vector <float> PGEs(nlocs, missingValueFloat);
    obsdb_.get_db("GrossErrorProbability", varPGE, PGEs);
    // Buddy check was not applied.
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (apply[iloc] && PGEs[iloc] >= packedPGECrit) {
        QCflags[iloc] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
        flagged[iloc] = true;
      }
    }
  } else {
    std::stringstream errormessage;
    errormessage << "At least one of "
                 << varname + "@GrossErrorProbability or "
                 << varname + "@GrossErrorProbabilityBuddyCheck must be present"
                 << std::endl;
    throw eckit::BadValue(errormessage.str(), Here());
  }

  // Save modified flags.
  obsdb_.put_db("QCFlags", varname, QCflags);
}

void BayesianBackgroundQCFlags::applyFilter(const std::vector<bool> & apply,
                                            const Variables & filtervars,
                                            std::vector<std::vector<bool>> & flagged) const {
  print(oops::Log::trace());

  for (size_t ivar = 0; ivar < filtervars.nvars(); ++ivar) {
    const std::string varname = filtervars.variable(ivar).variable();
    setFlags(varname, apply, flagged[ivar]);
  }
}

void BayesianBackgroundQCFlags::print(std::ostream & os) const {
  os << "BayesianBackgroundQCFlags: config = " << parameters_ << std::endl;
}

}  // namespace ufo
