/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsFilters.h"

#include <memory>
#include <set>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/AssociativeContainers.h"
#include "oops/util/IntSetParser.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsFilters::ObsFilters(ioda::ObsSpace & os,
                       const eckit::Configuration & config,
                       ObsDataPtr_<int> qcflags, ObsDataPtr_<float> obserr,
                       const int iteration)
  : obsspace_(os), autoFilters_(), preFilters_(), priorFilters_(), postFilters_(),
    geovars_(), diagvars_(), qcflags_(qcflags), obserrfilter_(obserr),
    iteration_(iteration) {
  oops::Log::trace() << "ObsFilters::ObsFilters starting";

  ObsFiltersParameters params;
  params.deserialize(config);

  const FilterParams_ & autoFiltersParams = params.obsFilters;
  const FilterParams_ & preFiltersParams = params.obsPreFilters;
  const FilterParams_ & priorFiltersParams = params.obsPriorFilters;
  const FilterParams_ & postFiltersParams = params.obsPostFilters;

  const bool atLeastOneAutoFilterConfigured = autoFiltersParams.size() > 0;
  const bool atLeastOnePrePriorPostFilterConfigured =
    preFiltersParams.size() +
    priorFiltersParams.size() +
    postFiltersParams.size() > 0;

  if (atLeastOneAutoFilterConfigured && atLeastOnePrePriorPostFilterConfigured) {
    throw eckit::UserError("It is not possible to use both an `obs filters` option "
                           "and one or more of the `obs pre filters`, "
                           "`obs prior filters` and `obs post filters` options", Here());
  }

  // If at least one filter has been configured, prepend the QC manager and
  // append the Final Check to the list of filters.
  const bool atLeastOneFilterConfigured =
    atLeastOneAutoFilterConfigured || atLeastOnePrePriorPostFilterConfigured;

// Prepare QC handling and statistics if any filters are present
  if (atLeastOneFilterConfigured) {
    eckit::LocalConfiguration conf;
    conf.set("filter", "QCmanager");
    ObsFilterParametersWrapper filterParams;
    filterParams.deserialize(conf);
    if (atLeastOneAutoFilterConfigured)
      autoFilters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrfilter_);
    else
      preFilters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrfilter_);
  }

// Create the filters, only at 0-th iteration, or at iterations specified in "apply at iterations"
  if (atLeastOneAutoFilterConfigured)
    appendToFiltersList(autoFiltersParams, autoFilters_);
  if (atLeastOnePrePriorPostFilterConfigured) {
    appendToFiltersList(preFiltersParams, preFilters_);
    appendToFiltersList(priorFiltersParams, priorFilters_);
    appendToFiltersList(postFiltersParams, postFilters_);
  }

// Create the final filter run at the end of the pipeline
  if (atLeastOneFilterConfigured) {
    eckit::LocalConfiguration conf;
    conf.set("filter", "Final Check");
    ObsFilterParametersWrapper filterParams;
    filterParams.deserialize(conf);
    if (atLeastOneAutoFilterConfigured)
      autoFilters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrfilter_);
    else
      postFilters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrfilter_);
  }

  oops::Log::trace() << "ObsFilters::ObsFilters done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsFilters::appendToFiltersList(const FilterParams_ & filtersParams,
                                     std::vector<ObsFilter> & filters) {
  for (const ObsFilterParametersWrapper & filterParams : filtersParams) {
    // Only create filters for the 0-th iteration by default
    bool apply = (iteration_ == 0);
    // If "apply at iterations" is set, check if this is the right iteration
    if (filterParams.applyAtIterations.value() != boost::none) {
      const std::set<int> iters = oops::parseIntSet(*filterParams.applyAtIterations.value());
      apply = oops::contains(iters, iteration_);
    }
    if (apply) {
      filters.emplace_back(obsspace_, filterParams.filterParameters, qcflags_, obserrfilter_);
      geovars_ += filters.back().requiredVars();
      diagvars_ += filters.back().requiredHdiagnostics();
    }
  }
}

// -----------------------------------------------------------------------------

void ObsFilters::preProcess() {
  if (autoFilters_.size() > 0) {
    for (ObsFilter & filter : autoFilters_) {
      filter.checkFilterData(FilterStage::AUTO);
      filter.preProcess();
    }
  } else {
    for (ObsFilter & filter : preFilters_) {
      filter.checkFilterData(FilterStage::PRE);
      filter.preProcess();
    }
  }
}

// -----------------------------------------------------------------------------

void ObsFilters::priorFilter(const GeoVaLs & gv) {
  if (autoFilters_.size() > 0) {
    for (ObsFilter & filter : autoFilters_) {
      filter.checkFilterData(FilterStage::AUTO);
      filter.priorFilter(gv);
    }
  } else {
    for (ObsFilter & filter : priorFilters_) {
      filter.checkFilterData(FilterStage::PRIOR);
      filter.priorFilter(gv);
    }
  }
}

// -----------------------------------------------------------------------------

void ObsFilters::postFilter(const GeoVaLs & gv,
                            const ioda::ObsVector & hofx,
                            const ioda::ObsVector & bias,
                            const ObsDiagnostics & diags) {
  if (autoFilters_.size() > 0) {
    for (ObsFilter & filter : autoFilters_) {
      filter.checkFilterData(FilterStage::AUTO);
      filter.postFilter(gv, hofx, bias, diags);
    }
  } else {
    for (ObsFilter & filter : postFilters_) {
      filter.checkFilterData(FilterStage::POST);
      filter.postFilter(gv, hofx, bias, diags);
    }
  }
}

// -----------------------------------------------------------------------------

void ObsFilters::print(std::ostream & os) const {
  if (autoFilters_.size() > 0) {
    os << "ObsFilters: " << autoFilters_.size() << " elements:" << std::endl;
    for (const ObsFilter & filter : autoFilters_)
      os << filter << std::endl;
  } else {
    os << "Pre filters: " << preFilters_.size() << " elements:" << std::endl;
    for (const ObsFilter & preFilter : preFilters_)
      os << preFilter << std::endl;
    os << "Prior filters: " << priorFilters_.size() << " elements:" << std::endl;
    for (const ObsFilter & priorFilter : priorFilters_)
      os << priorFilter << std::endl;
    os << "Post filters: " << postFilters_.size() << " elements:" << std::endl;
    for (const ObsFilter & postFilter : postFilters_)
      os << postFilter << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
