/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

#include "ufo/ObsFilter.h"

namespace ioda {
class ObsSpace;
class ObsVector;
template <typename DATA> class ObsDataVector;
}

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

/// Configuration options for the ObsFilters class.
class ObsFiltersParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsFiltersParameters, Parameters)

  typedef typename std::vector<ObsFilterParametersWrapper> FilterParams_;

 public:
  /// Options used to configure observation filters whose stage of operation
  /// (pre, prior or post) is determined automatically.
  oops::Parameter<FilterParams_> obsFilters{"obs filters", {}, this};

  /// Options used to configure observation pre filters.
  /// These filters are called before GetValues.
  /// Both GeoVaLs and H(x) are unavailable.
  oops::Parameter<FilterParams_> obsPreFilters{"obs pre filters", {}, this};

  /// Options used to configure observation filters.
  /// These filters are called after GetValues and before the observation operator.
  /// GeoVaLs are available and H(x) is unavailable.
  oops::Parameter<FilterParams_> obsPriorFilters{"obs prior filters", {}, this};

  /// Options used to configure observation filters.
  /// These filters are called after the observation operator.
  ///  Both GeoVaLs and H(x) are available.
  oops::Parameter<FilterParams_> obsPostFilters{"obs post filters", {}, this};
};

/// Holds observation filters (usually QC) for one observation type

// -----------------------------------------------------------------------------

class ObsFilters : public util::Printable,
                   private boost::noncopyable {
  typedef std::vector<ObsFilterParametersWrapper> FilterParams_;
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ioda::ObsDataVector<DATA> >;

 public:
  /// Initialize all filters for \p obspace, from parameters, using
  /// \p qcflags and \p obserr (observation error variances)
  /// \p iteration argument indicates outer loop iteration in the variational
  /// assimilation.
  ObsFilters(ioda::ObsSpace &,
             const eckit::Configuration &,
             ObsDataPtr_<int> qcflags, ObsDataPtr_<float> obserr,
             const int iteration = 0);

  void preProcess();
  void priorFilter(const GeoVaLs &);
  void postFilter(const GeoVaLs &,
                  const ioda::ObsVector &,
                  const ioda::ObsVector &,
                  const ObsDiagnostics &);

  oops::Variables requiredVars() const {return geovars_;}
  oops::ObsVariables requiredHdiagnostics() const {return diagvars_;}

 private:
  void print(std::ostream &) const override;

  /// Configure a filter and append it to a vector of filters.
  void appendToFiltersList(const FilterParams_ & filtersParams,
                           std::vector<ObsFilter> & filters);

  ioda::ObsSpace & obsspace_;
  // List of filters for which the stage (pre/prior/post) will be determined automatically.
  std::vector<ObsFilter> autoFilters_;
  // List of filters which have been designated to run at the pre stage.
  std::vector<ObsFilter> preFilters_;
  // List of filters which have been designated to run at the prior stage.
  std::vector<ObsFilter> priorFilters_;
  // List of filters which have been designated to run at the post stage .
  std::vector<ObsFilter> postFilters_;
  oops::Variables geovars_;
  oops::ObsVariables diagvars_;
  ObsDataPtr_<int> qcflags_;
  ObsDataPtr_<float> & obserrfilter_;
  const int iteration_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
