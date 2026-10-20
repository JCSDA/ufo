/*
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/sfccorrected/SurfaceOperatorBase.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"

namespace ufo {

SurfaceOperatorBase::SurfaceOperatorBase(const std::string & name, const Parameters_ & params)
  : params_(params), name_(name)
{}

void SurfaceOperatorBase::settraj(const GeoVaLs &, const ioda::ObsSpace &) {
  throw eckit::NotImplemented(name_ + "::settraj not yet implemented", Here());
}

void SurfaceOperatorBase::TL(const GeoVaLs &, std::vector<double> &) const {
  throw eckit::NotImplemented(name_ + "::TL not yet implemented", Here());
}

void SurfaceOperatorBase::AD(GeoVaLs &, std::vector<double> &) const {
  throw eckit::NotImplemented(name_ + "::AD not yet implemented", Here());
}

SurfaceOperatorFactory::SurfaceOperatorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end())
    throw eckit::BadParameter(name + " already registered in ufo::SurfaceOperatorFactory.", Here());
  getMakers()[name] = this;
}

std::unique_ptr<SurfaceOperatorBase>
SurfaceOperatorFactory::create(const std::string & name,
                               const Parameters_ & params) {
  oops::Log::trace() << "SurfaceOperatorBase::create start" << std::endl;
  typename std::map<std::string, SurfaceOperatorFactory*>::iterator jloc = getMakers().find(name);
  if (jloc == getMakers().end()) {
    std::string makerNameList;
    for (const auto & makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
    throw eckit::BadParameter(name + " does not exist in ufo::SurfaceOperatorFactory. "
                              "Possible values:" + makerNameList, Here());
  }
  std::unique_ptr<SurfaceOperatorBase> ptr = jloc->second->make(name, params);
  oops::Log::trace() << "SurfaceOperatorBase::create done" << std::endl;
  return ptr;
}

std::vector<std::unique_ptr<SurfaceOperatorBase>>
SurfaceOperatorFactory::createForVariables(const oops::ObsVariables & operatorVars,
                                           const Parameters_ & params,
                                           oops::Variables & requiredVars) {
  // Get surface correction type as a string
  std::string methodname;
  for (util::NamedEnumerator<SfcCorrectionType> namedValue :
       SfcCorrectionTypeParameterTraitsHelper::namedValues) {
    if (namedValue.value == params.correctionType.value()) {
      methodname = namedValue.name;
    }
  }

  // Create operators for each variable / surface correction type combination
  std::vector<std::unique_ptr<SurfaceOperatorBase>> operators;
  for (const std::string & var : operatorVars.variables()) {
    std::unique_ptr<SurfaceOperatorBase> oper = create(var + "_" + methodname, params);
    requiredVars += oper->requiredVars();
    operators.push_back(std::move(oper));
  }
  return operators;
}

}  // namespace ufo
