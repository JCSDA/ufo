/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsOperatorBase.h"

#include <vector>

#include "ioda/ObsSpace.h"
#include "oops/base/Locations.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsOperatorBase::Locations_ ObsOperatorBase::locations() const {
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  std::vector<float> lons(odb_.nlocs());
  std::vector<float> lats(odb_.nlocs());
  std::vector<util::DateTime> times(odb_.nlocs());
  odb_.get_db("MetaData", "latitude", lats);
  odb_.get_db("MetaData", "longitude", lons);
  odb_.get_db("MetaData", "dateTime", times);
  return SampledLocations_(
        std::make_unique<SampledLocations>(lons, lats, times, odb_.distribution()));
}

// -----------------------------------------------------------------------------

oops::ObsVariables ObsOperatorBase::simulatedVars() const {
  return odb_.assimvariables();
}

// -----------------------------------------------------------------------------

void ObsOperatorBase::computeReducedVars(const oops::Variables & /*vars*/, GeoVaLs & gvals) const {
  // Ask the GeoVaLs object to store a reduced representation of all variables for which this
  // representation has not yet been allocated (the Doxygen comment explains why the `vars` argument
  // is ignored)
  const oops::Variables allVars = gvals.getVars();
  oops::Variables missingVars = allVars;
  missingVars -= gvals.getReducedVars();
  std::vector<size_t> sizes(missingVars.size());
  for (size_t i = 0; i < missingVars.size(); ++i) {
    // Variable `missingVars[i]` must already exist in the sampled format; otherwise the program
    // will be aborted.
    sizes[i] = gvals.nlevs(missingVars[i], GeoVaLFormat::SAMPLED);
  }
  gvals.addReducedVars(missingVars, sizes);

  // Throw an exception if the reduced version of any variable is not aliased with (stored in the
  // same piece of memory as) the sampled version -- it means we're expected to fill it, but we
  // don't know how. (If reduction of this variable is really needed, the obs operator should
  // override this method appropriately.)
  for (const auto & var : allVars) {
    if (!gvals.areReducedAndSampledFormatsAliased(var))
      throw eckit::NotImplemented("Unable to reduce variable " + var.name(), Here());
  }
}

// -----------------------------------------------------------------------------

ObsOperatorFactory::ObsOperatorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in ufo::ObsOperatorFactory." << std::endl;
    ABORT("Element already registered in ufo::ObsOperatorFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

ObsOperatorBase * ObsOperatorFactory::create(const ioda::ObsSpace & odb,
                                             const ObsOperatorParametersBase & params) {
  oops::Log::trace() << "ObsOperatorBase::create starting" << std::endl;
  const std::string &id = params.name.value().value();
  typename std::map<std::string, ObsOperatorFactory*>::iterator jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in ufo::ObsOperatorFactory." << std::endl;
    ABORT("Element does not exist in ufo::ObsOperatorFactory.");
  }
  ObsOperatorBase * ptr = jloc->second->make(odb, params);
  oops::Log::trace() << "ObsOperatorBase::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

std::unique_ptr<ObsOperatorParametersBase>
ObsOperatorFactory::createParameters(const std::string &name) {
  typename std::map<std::string, ObsOperatorFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in ufo::ObsOperatorFactory");
  }
  return it->second->makeParameters();
}

// -----------------------------------------------------------------------------

}  // namespace ufo
