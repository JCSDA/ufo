/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ObsProcessorBase.h"

#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"

#include "ufo/filters/actions/FilterAction.h"
#include "ufo/filters/GenericFilterParameters.h"
#include "ufo/filters/processWhere.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsProcessorBase::ObsProcessorBase(ioda::ObsSpace & os, bool deferToPost,
                                   std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                   std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : obsdb_(os),
    flags_(flags), obserr_(obserr),
    data_(obsdb_), prior_(false), post_(false),
    deferToPost_(deferToPost)
{
  oops::Log::trace() << "ObsProcessorBase constructor" << std::endl;
  ASSERT(flags);
  ASSERT(obserr);
  data_.associate(*flags_, "QCflagsData");
  data_.associate(*obserr_, "ObsErrorData");
}

// -----------------------------------------------------------------------------

ObsProcessorBase::~ObsProcessorBase() {
  oops::Log::trace() << "ObsProcessorBase destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProcessorBase::preProcess() {
  oops::Log::trace() << "ObsProcessorBase preProcess begin" << std::endl;
// Cannot determine earlier when to apply filter because subclass
// constructors add to allvars
  if (allvars_.hasGroup("HofX") || allvars_.hasGroup("ObsDiag") || deferToPost_) {
    post_ = true;
  } else {
    if (allvars_.hasGroup("GeoVaLs")) {
      prior_ = true;
    } else {
      this->doFilter();
    }
  }
  oops::Log::trace() << "ObsProcessorBase preProcess end" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProcessorBase::priorFilter(const GeoVaLs & gv) {
  oops::Log::trace() << "ObsProcessorBase priorFilter begin" << std::endl;
  if (prior_ || post_) data_.associate(gv);
  if (prior_) this->doFilter();
  oops::Log::trace() << "ObsProcessorBase priorFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProcessorBase::postFilter(const ioda::ObsVector & hofx, const ObsDiagnostics & diags) {
  oops::Log::trace() << "ObsProcessorBase postFilter begin" << std::endl;
  if (post_) {
    data_.associate(hofx, "HofX");
    data_.associate(diags);
    this->doFilter();
  }
  oops::Log::trace() << "ObsProcessorBase postFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
