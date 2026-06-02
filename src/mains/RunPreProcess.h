/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
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

#include "ufo/ObsFilters.h"

namespace ufo {

/// \brief Options used to configure an application saving pre-processed files
class RunPreProcessParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(RunPreProcessParameters, Parameters)

  typedef typename std::vector<ObsFilterParametersWrapper> FilterParams_;

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> timeWindow{"time window",
        "options used to configure the assimilation time window", this};
  oops::RequiredParameter<ioda::ObsTopLevelParameters> obsSpace{"obs space",
        "options used to configure the observation space", this};
  ObsFiltersParameters filtersParams{this};
};

// -----------------------------------------------------------------------------

class RunPreProcess : public oops::Application {
 public:
  explicit RunPreProcess(const eckit::mpi::Comm & comm = oops::mpi::world()):
    Application(comm) {}

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const {
    RunPreProcessParameters params;
    params.deserialize(fullConfig);

    ioda::ObsSpace obsdb(params.obsSpace.value().toConfiguration(), this->getComm(),
                         util::TimeWindow(params.timeWindow.value()),
                         this->getComm());

    std::shared_ptr<ioda::ObsDataVector<float>> obserr_ =
      std::make_shared<ioda::ObsDataVector<float>>(obsdb, obsdb.obsvariables(), "ObsError");
    std::shared_ptr<ioda::ObsDataVector<int>> qcflags_ =
      std::make_shared<ioda::ObsDataVector<int>>(obsdb, obsdb.obsvariables());

    ufo::ObsFilters filters(obsdb,
                            params.filtersParams.toConfiguration(),
                            qcflags_, obserr_);

    filters.preProcess();

    obsdb.save();

    return 0;
  }

// -----------------------------------------------------------------------------

 private:
  std::string appname() const { return "ufo::RunPreProcess"; }
};

}  // namespace ufo
