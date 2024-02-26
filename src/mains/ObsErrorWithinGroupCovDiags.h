/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/TimeWindow.h"

#include "ufo/errors/ObsErrorWithinGroupCov.h"

namespace ufo {

/// \brief Options used to configure which diagnostics to save and where
class DiagParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DiagParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> filename{"output filename",
                         "filename where the diagnostics will be saved", this};
  oops::RequiredParameter<size_t> recnum{"record number", this};
};

/// \brief Options used to configure an application saving diagnostics for the
/// ObsErrorWithinGroupCov correlations.
class ObsErrorWithinGroupCovDiagsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorWithinGroupCovDiagsParameters, Parameters)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> timeWindow{"time window",
        "options used to configure the assimilation time window", this};
  oops::RequiredParameter<ioda::ObsTopLevelParameters> obsSpace{"obs space",
        "options used to configure the observation space", this};
  oops::RequiredParameter<ObsErrorWithinGroupCovParameters> obsError{"obs error",
        "options configuring ObsErrorWithinGroupCov", this};
  oops::RequiredParameter<DiagParameters> diags{"obs error diagnostics",
        "options configuring which diagnostics to save and where", this};
};

// -----------------------------------------------------------------------------

class ObsErrorWithinGroupCovDiags : public oops::Application {
 public:
  explicit ObsErrorWithinGroupCovDiags(const eckit::mpi::Comm & comm = oops::mpi::world()):
    Application(comm) {}

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig, bool validate) const {
    ObsErrorWithinGroupCovDiagsParameters params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    ioda::ObsSpace obsdb(params.obsSpace.value().toConfiguration(), this->getComm(),
                         util::TimeWindow(fullConfig.getSubConfiguration("time window")),
                         this->getComm());
    ioda::ObsVector randomVec(obsdb);
    randomVec.random();

    ObsErrorWithinGroupCov obserr(params.obsError, obsdb, this->getComm());
    obserr.saveCorrelations(params.diags.value().filename, params.diags.value().recnum,
                            randomVec);
    return 0;
  }

// -----------------------------------------------------------------------------

 private:
  std::string appname() const { return "ufo::ObsErrorWithinGroupCovDiags"; }
};

}  // namespace ufo
