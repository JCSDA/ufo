/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/categoricaloper/ObsCategorical.h"

#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/categoricaloper/ObsCategoricalParameters.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsCategorical> obsCategoricalMaker_("Categorical");
// -----------------------------------------------------------------------------

ObsCategorical::ObsCategorical(const ioda::ObsSpace & odb,
                               const Parameters_ & params)
  : ObsOperatorBase(odb), odb_(odb)
{
  oops::Log::trace() << "ObsCategorical constructor starting" << std::endl;

  data_.configure(odb, params);

  oops::Log::trace() << "ObsCategorical constructor finished" << std::endl;
}

// -----------------------------------------------------------------------------

ObsCategorical::~ObsCategorical() {
  oops::Log::trace() << "ObsCategorical destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCategorical::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                 ObsDiagnostics & ydiags, const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsCategorical: simulateObs entered" << std::endl;

  oops::Log::debug() << "Running operators" << std::endl;

  // Container of ObsVectors produced by each operator.
  std::map <std::string, ioda::ObsVector> ovecs;
  // Run each operator and store output in ovecs.
  for (const auto& component : data_.components()) {
    ioda::ObsVector ovecTemp(ovec);
    component.second->simulateObs(gv, ovecTemp, ydiags, qc_flags);
    ovecs.insert({component.first, ovecTemp});
  }

  oops::Log::debug() << "Producing final ObsVector" << std::endl;

  data_.fillHofX(ovecs, ovec);

  oops::Log::trace() << "ObsCategorical: simulateObs finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsCategorical::print(std::ostream & os) const {
  data_.print(os);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
