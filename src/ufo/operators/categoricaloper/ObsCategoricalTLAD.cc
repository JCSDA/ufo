/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/categoricaloper/ObsCategoricalTLAD.h"

#include <algorithm>
#include <map>
#include <ostream>
#include <utility>

#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/categoricaloper/ObsCategoricalParameters.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsCategoricalTLAD> makerCategoricalTL_("Categorical");
// -----------------------------------------------------------------------------

ObsCategoricalTLAD::ObsCategoricalTLAD(const ioda::ObsSpace & odb,
                                       const Parameters_ & params)
  : LinearObsOperatorBase(odb, VariableNameMap(params.AliasFile.value())), odb_(odb)
{
  oops::Log::trace() << "ObsCategoricalTLAD constructor starting" << std::endl;

  data_.configure(odb, params);

  oops::Log::trace() << "ObsCategoricalTLAD created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsCategoricalTLAD::~ObsCategoricalTLAD() {
  oops::Log::trace() << "ObsCategoricalTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsCategoricalTLAD::setTrajectory(const GeoVaLs & geovals,
                                       ObsDiagnostics & ydiags,
                                       const QCFlags_t & qc_flags) {
  oops::Log::trace() << "ObsCategoricalTLAD: setTrajectory entered" << std::endl;

  // Set trajectory for each operator.
  for (const auto& component : data_.components())
    component.second->setTrajectory(geovals, ydiags, qc_flags);

  oops::Log::trace() << "ObsCategoricalTLAD: setTrajectory finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsCategoricalTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                       const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsCategoricalTLAD: simulateObsTL entered" << std::endl;

  oops::Log::debug() << "Running TL operators" << std::endl;

  // Container of ObsVectors produced by each TL operator.
  std::map <std::string, ioda::ObsVector> ovecs;
  // Run each TL operator and store output in ovecs.
  for (const auto& component : data_.components()) {
    ioda::ObsVector ovecTemp(ovec);
    component.second->simulateObsTL(geovals, ovecTemp, qc_flags);
    ovecs.insert({component.first, ovecTemp});
  }

  oops::Log::debug() << "Producing final TL" << std::endl;

  data_.fillHofX(ovecs, ovec);

  oops::Log::trace() << "ObsCategoricalTLAD: simulateObsTL finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsCategoricalTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                       const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsCategoricalTLAD: simulateObsAD entered" << std::endl;

  oops::Log::debug() << "Running AD operators" << std::endl;

  // Container of GeoVaLs produced by each AD operator.
  std::map <std::string, GeoVaLs> gvals;
  // Run each AD operator and store output in gvals.
  for (const auto& component : data_.components()) {
    GeoVaLs gvalTemp(geovals);
    component.second->simulateObsAD(gvalTemp, ovec, qc_flags);
    gvals.insert({component.first, gvalTemp});
  }

  oops::Log::debug() << "Producing final AD" << std::endl;

  // Insert values into geovals according to the categorical variable.
  // Use the fallback operator when necessary.
  const std::vector<std::string> &varnames = ovec.varnames().variables();
  std::vector <double> vecgv;
  for (size_t jloc = 0; jloc < ovec.nlocs(); ++jloc) {
    const auto &operName = data_.locOperNames()[jloc];
    const auto &gvaloper = gvals.at(operName);
    // Loop over each variable at this location.
    for (const auto& varname : varnames) {
      vecgv.resize(gvaloper.nlevs(nameMap_.convertName(varname)));
      gvaloper.getAtLocation(vecgv, nameMap_.convertName(varname), jloc);
      geovals.putAtLocation(vecgv, nameMap_.convertName(varname), jloc);
    }
  }

  oops::Log::trace() << "ObsCategoricalTLAD: simulateObsAD finished" <<  std::endl;
}

// -----------------------------------------------------------------------------

void ObsCategoricalTLAD::print(std::ostream & os) const {
  data_.print(os);
}

// -----------------------------------------------------------------------------

}  // namespace ufo
