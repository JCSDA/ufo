/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ObsProcessorBase.h"

#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/filters/actions/FilterAction.h"
#include "ufo/filters/GenericFilterParameters.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ScopedDefaultGeoVaLFormatChange.h"

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
  if (allvars_.hasGroup("HofX") || allvars_.hasGroup("ObsDiag") ||
      allvars_.hasGroup("ObsBiasData") || deferToPost_) {
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
  ScopedDefaultGeoVaLFormatChange change(gv, GeoVaLFormat::REDUCED);
  if (prior_ || post_) data_.associate(gv);
  if (prior_) this->doFilter();
  oops::Log::trace() << "ObsProcessorBase priorFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProcessorBase::postFilter(const GeoVaLs & gv,
                                  const ioda::ObsVector & hofx,
                                  const ioda::ObsVector & bias,
                                  const ObsDiagnostics & diags) {
  oops::Log::trace() << "ObsProcessorBase postFilter begin" << std::endl;
  if (post_) {
    ScopedDefaultGeoVaLFormatChange change(gv, GeoVaLFormat::REDUCED);
    data_.associate(gv);
    data_.associate(hofx, "HofX");
    data_.associate(bias, "ObsBiasData");
    data_.associate(diags);
    this->doFilter();
  }
  oops::Log::trace() << "ObsProcessorBase postFilter end" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsProcessorBase::checkFilterData(const oops::FilterStage filterStage) {
  // Return if filters have been automatically designated as pre, prior or post.
  if (filterStage == oops::FilterStage::AUTO)
    return;

  // Pre filters (run before GetValues) cannot request quantities in the
  // GeoVaLs, HofX, ObsDiag or ObsBiasData groups.
  if (filterStage == oops::FilterStage::PRE &&
      (allvars_.hasGroup("GeoVaLs") ||
       allvars_.hasGroup("HofX") ||
       allvars_.hasGroup("ObsDiag") ||
       allvars_.hasGroup("ObsBiasData"))) {
    throw eckit::UserError("Invalid pre filter requested", Here());
  }

  // Prior filters (run after GetValues and before observation operator) cannot request
  // quantities in the HofX, ObsDiag or ObsBiasData groups.
  if (filterStage == oops::FilterStage::PRIOR &&
      (allvars_.hasGroup("HofX") ||
       allvars_.hasGroup("ObsDiag") ||
       allvars_.hasGroup("ObsBiasData"))) {
    throw eckit::UserError("Invalid prior filter requested", Here());
  }

  // There are no requirements on post filters (run after observation operator).

  // Depending on filter stage set prior_ and/or post_ to true.
  // This ensures priorFilter and postFilter will run doFilter().
  if (filterStage == oops::FilterStage::POST) {
    prior_ = true;
    post_ = true;
  }
  if (filterStage == oops::FilterStage::PRIOR) {
    prior_ = true;
  }
}

}  // namespace ufo
