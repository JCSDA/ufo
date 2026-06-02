/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsFilter.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"

#include "ufo/ObsFilterBase.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsFilter::ObsFilter(ioda::ObsSpace & os,
                     const ObsFilterParametersBase & parameters,
                     ObsDataPtr_<int> flags, ObsDataPtr_<float> obserr)
  : ofilt_(), filterName_("ObsFilter::"+parameters.filter.value().value())
{
  oops::Log::trace() << "ObsFilter::ObsFilter starting" << std::endl;
  util::Timer timer(classname(), "ObsFilter");
  util::Timer timef(filterName_, "ObsFilter");
  ofilt_ = FilterFactory::create(os, parameters, flags, obserr);
  oops::Log::trace() << "ObsFilter::ObsFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsFilter::~ObsFilter() {
  oops::Log::trace() << "ObsFilter::~ObsFilter starting" << std::endl;
  util::Timer timer(classname(), "~ObsFilter");
  ofilt_.reset();
  oops::Log::trace() << "ObsFilter::~ObsFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsFilter::preProcess() {
  oops::Log::trace() << "ObsFilter::preProcess starting" << std::endl;
  util::Timer timer(classname(), "preProcess");
  util::Timer timef(filterName_, "preProcess");
  ofilt_->preProcess();
  oops::Log::trace() << "ObsFilter::preProcess done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsFilter::priorFilter(const GeoVaLs & gv) {
  oops::Log::trace() << "ObsFilter::priorFilter starting" << std::endl;
  util::Timer timer(classname(), "priorFilter");
  util::Timer timef(filterName_, "priorFilter");
  ofilt_->priorFilter(gv);
  oops::Log::trace() << "ObsFilter::priorFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsFilter::postFilter(const GeoVaLs & gv,
                           const ioda::ObsVector & ov,
                           const ioda::ObsVector & bv,
                           const ObsDiagnostics & dv) {
  oops::Log::trace() << "ObsFilter::postFilter starting" << std::endl;
  util::Timer timer(classname(), "postFilter");
  util::Timer timef(filterName_, "postFilter");
  ofilt_->postFilter(gv, ov, bv, dv);
  oops::Log::trace() << "ObsFilter::postFilter done" << std::endl;
}

// -----------------------------------------------------------------------------

oops::Variables ObsFilter::requiredVars() const {
  oops::Log::trace() << "ObsFilter::requiredVars" << std::endl;
  return ofilt_->requiredVars();
}

// -----------------------------------------------------------------------------

oops::ObsVariables ObsFilter::requiredHdiagnostics() const {
  oops::Log::trace() << "ObsFilter::requiredHdiagnostics" << std::endl;
  return ofilt_->requiredHdiagnostics();
}

// -----------------------------------------------------------------------------

void ObsFilter::print(std::ostream & os) const {
  oops::Log::trace() << "ObsFilter::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *ofilt_;
  oops::Log::trace() << "ObsFilter::print done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsFilter::checkFilterData(const FilterStage filterStage) {
  oops::Log::trace() << "ObsFilter::checkFilterData starting" << std::endl;
  util::Timer timer(classname(), "checkFilterData");
  ofilt_->checkFilterData(filterStage);
  oops::Log::trace() << "ObsFilter::checkFilterData done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
