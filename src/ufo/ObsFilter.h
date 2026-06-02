/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/ObsFilterBase.h"

namespace ioda {
class ObsSpace;
class ObsVector;
template <typename DATA> class ObsDataVector;
}

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

// -----------------------------------------------------------------------------

/// \brief A filter processing observations.
class ObsFilter : public util::Printable,
                  private util::ObjectCounter<ObsFilter > {
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ioda::ObsDataVector<DATA> >;

 public:
  static const std::string classname() {return "ufo::ObsFilter";}

  /// \brief Create a new observation filter.
  ///
  /// \param obsspace
  ///   Space containing the observations to process.
  /// \param params
  ///   The filter's configuration parameters.
  /// \param qcflags
  ///   Quality control flags. They may be modified by the filter.
  /// \param obserrors
  ///   Estimates of the standard deviations of observation errors. They may be modified by the
  ///   filter.
  ObsFilter(ioda::ObsSpace &obsspace, const ObsFilterParametersBase &params,
            ObsDataPtr_<int> qcflags, ObsDataPtr_<float> obserrors);
  ObsFilter(const ObsFilter &) = delete;
  ObsFilter(ObsFilter &&) = default;
  ObsFilter& operator=(const ObsFilter &) = delete;
  ObsFilter& operator=(ObsFilter &&) = default;
  ~ObsFilter();

  /// \brief Perform any observation processing steps that do not require access to GeoVaLs or
  /// outputs produced by the observation operator.
  void preProcess();

  /// \brief Perform any observation processing steps that require access to GeoVaLs, but not to
  /// outputs produced by the observation operator.
  void priorFilter(const GeoVaLs &);

  /// \brief Perform any observation processing steps that require access to both GeoVaLs and
  /// outputs produced by the observation operator.
  ///
  /// \param gv
  ///   GeoVaLs.
  /// \param ov
  ///   Model equivalents produced by the observation operator.
  /// \param bv
  ///   Bias of departure produced by the observation operator.
  /// \param dv
  ///   Observation diagnostics produced by the observation operator.
  void postFilter(const GeoVaLs & gv,
                  const ioda::ObsVector &ov,
                  const ioda::ObsVector &bv,
                  const ObsDiagnostics &dv);

  /// \brief Check the required filter data are present prior to running this filter.
  void checkFilterData(const FilterStage filterStage);

  /// \brief Return the list of GeoVaLs required by this filter.
  oops::Variables requiredVars() const;

  /// \brief Return the list of observation diagnostics required by this filter.
  oops::ObsVariables requiredHdiagnostics() const;

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<ObsFilterBase> ofilt_;
  std::string filterName_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
