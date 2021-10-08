/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_METOFFICEBUDDYCHECK_H_
#define UFO_FILTERS_METOFFICEBUDDYCHECK_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/MetOfficeBuddyCheckParameters.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
template <typename DATATYPE> class ObsDataVector;
class ObsSpace;
}

namespace util {
class DateTime;
}

namespace ufo {

class RecursiveSplitter;
class MetOfficeBuddyPair;

/// \brief Met Office's implementation of the buddy check.
///
/// The filter cross-checks observations taken at nearby locations against each other,
/// updating their gross error probabilities (PGEs) and rejecting observations whose PGE exceeds
/// a threshold specified in the filter parameters.
///
/// Variables to be checked should be specified using the "filter variables" YAML option, supporting
/// surface (single-level) and multi-level variables. Variables can be either scalar or vector
/// (with two Cartesian components, such as the eastward and northward wind components). In the
/// latter case the two components need to specified one after the other in the "filter variables"
/// list, with the first component having the \c first_component_of_two option set to true.
/// Example:
///
/// \code{.yaml}
/// filter variables:
/// - name: air_temperature
/// - name: eastward_wind
///   options:
///     first_component_of_two: true
/// - name: northward_wind
/// \endcode
///
/// See MetOfficeBuddyCheckParameters for the documentation of the other available parameters.
///
/// This filter assumes background error estimates for each filter variable `var` can be retrieved
/// from the `var_background_error` ObsDiagnostic. These diagnostics are typically produced using
/// the BackgroundErrorVertInterp or BackgroundErrorIdentity operators. To make sure these operators
/// are applied in addition to the "main" operator calculating model equivalents of observations,
/// instruct the system to use a Composite operator with two or more components; for instance,
///
/// \code{.yaml}
/// obs operator:
///   name: Composite
///   components:
///   # operator used to evaluate H(x)
///   - name: Identity
///   # operator used to evaluate background errors
///   - name: BackgroundErrorIdentity
/// \endcode
class MetOfficeBuddyCheck : public FilterBase,
                            private util::ObjectCounter<MetOfficeBuddyCheck> {
 public:
  typedef MetOfficeBuddyCheckParameters Parameters_;
  static const std::string classname() {return "ufo::MetOfficeBuddyCheck";}

  MetOfficeBuddyCheck(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                      std::shared_ptr<ioda::ObsDataVector<int> > flags,
                      std::shared_ptr<ioda::ObsDataVector<float> > obserr);

  ~MetOfficeBuddyCheck() override;

 private:
  struct MetaData;

  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override { return QCflags::buddy; }

  /// \brief Return the name of the variable containing the background error estimate of the
  /// specified filter variable.
  Variable backgroundErrorVariable(const Variable &filterVariable) const;

  /// \brief Returns a vector of IDs of all observations that should be buddy-checked.
  std::vector<size_t> getValidObservationIds(
      const std::vector<bool> & apply,
      const boost::optional<Eigen::ArrayXXi> & profileIndex) const;

  /// \brief Collects and returns metadata of all observations.
  MetaData collectMetaData(const boost::optional<Eigen::ArrayXXi> & profileIndex) const;

  /// \brief Returns a vector of integer-valued station IDs, obtained from the source indicated by
  /// the filter parameters.
  std::vector<int> getStationIds() const;

  /// \brief Calculates and returns background error correlation scales at observation locations.
  ///
  /// Only elements with indices corresponding to those of valid observations are filled in.
  std::vector<float> calcBackgroundErrorHorizontalCorrelationScales(
      const std::vector<size_t> &validObsIds, const std::vector<float> &latitudes) const;

  /// \brief Identifies observations whose buddy checks should be logged.
  ///
  /// This function identifies observations located in one of the boxes specified in the
  /// \c tracedBoxes option, prints the metadata of these observations and sets the corresponding
  /// elements of the returned vector to true. Functions performing buddy checks are expected to
  /// log buddy checks involving these observations.
  ///
  /// \param pressures Optional -- may be null.
  std::vector<bool> flagAndPrintVerboseObservations(
      const std::vector<size_t> &validObsIds,
      const std::vector<float> &latitudes,
      const std::vector<float> &longitudes,
      const std::vector<util::DateTime> &times,
      const std::vector<float> *pressures,
      const std::vector<int> &stationIds,
      const std::vector<float> &bgErrorHorizCorrScales) const;

  /// \brief Buddy check for scalar quantities.
  ///
  /// Method: see the OPS Scientific Documentation Paper 2, sections 3.6 and 3.7
  ///
  /// \param pairs
  ///   Buddy pairs.
  /// \param flags
  ///   Observation flags associated with the variable being checked.
  /// \param verbose
  ///   Whether to log buddy checks involving particular observations.
  /// \param bgErrorHorizCorrScales
  ///   Background error horizontal correlation scales (in km).
  /// \param stationIds
  ///   Station IDs ("call signs").
  /// \param datetimes
  ///   Observation times.
  /// \param pressures
  ///   Model average pressures can be null, representing surface data, single-level (num_levels
  ///   parameter == 1) data or multi-level (num_levels > 1) data.  In all cases except
  ///   the single-level case, a vertical correlation of 1 is assumed.  For single-level data,
  ///   the estimate of the background error correlation depends upon the ratio of pressures
  ///   between each pair of observations.
  /// \param obsValues
  ///   Observed values.
  /// \param obsErrors
  ///   Estimated errors of observed values.
  /// \param bgValues
  ///   Background values.
  /// \param bgErrors
  ///   Estimated errors of background values.
  /// \param[inout] pges
  ///   Gross error probabilities. These values are updated by the buddy check.
  void checkScalarData(const std::vector<MetOfficeBuddyPair> &pairs,
                       const std::vector<int> &flags,
                       const std::vector<bool> &verbose,
                       const std::vector<float> &bgErrorHorizCorrScales,
                       const std::vector<int> &stationIds,
                       const std::vector<util::DateTime> &datetimes,
                       const Eigen::ArrayXXf *pressures,
                       const Eigen::ArrayXXf &obsValues,
                       const Eigen::ArrayXXf &obsErrors,
                       const Eigen::ArrayXXf &bgValues,
                       const Eigen::ArrayXXf &bgErrors,
                       Eigen::ArrayXXf &pges) const;

  /// \brief Buddy check for vector (two-dimensional) quantities.
  ///
  /// Method: see the OPS Scientific Documentation Paper 2, sections 3.6 and 3.7
  ///
  /// \param pairs
  ///   Buddy pairs.
  /// \param flags
  ///   Observation flags associated with the variable being checked.
  /// \param verbose
  ///   Whether to log buddy checks involving particular observations.
  /// \param bgErrorHorizCorrScales
  ///   Background error horizontal correlation scales (in km).
  /// \param stationIds
  ///   Station IDs ("call signs").
  /// \param datetimes
  ///   Observation times.
  /// \param pressures
  ///   Model average pressures can be null, represent surface data, single-level (num_levels
  ///   parameter == 1) data or multi-level (num_levels > 1) data.  In all cases except
  ///   the single-level case, a vertical correlation of 1 is assumed.  For single-level data,
  ///   the estimate of the background error correlation depends upon the ratio of pressures
  ///   between each pair of observations.
  /// \param uObsValues
  ///   Observed values of the first component, u.
  /// \param vObsValues
  ///   Observed values of the second component, v.
  /// \param obsErrors
  ///   Estimated errors of observed values (u or v).
  /// \param uBgValues
  ///   Background values of u.
  /// \param vBgValues
  ///   Background values of v.
  /// \param bgErrors
  ///   Estimated errors of background values (u or v).
  /// \param[inout] pges
  ///   Probabilities of gross error in u or v. These values are updated by the buddy check.
  void checkVectorData(const std::vector<MetOfficeBuddyPair> &pairs,
                       const std::vector<int> &flags,
                       const std::vector<bool> &verbose,
                       const std::vector<float> &bgErrorHorizCorrScales,
                       const std::vector<int> &stationIds,
                       const std::vector<util::DateTime> &datetimes,
                       const Eigen::ArrayXXf *pressures,
                       const Eigen::ArrayXXf &uObsValues,
                       const Eigen::ArrayXXf &vObsValues,
                       const Eigen::ArrayXXf &obsErrors,
                       const Eigen::ArrayXXf &uBgValues,
                       const Eigen::ArrayXXf &vBgValues,
                       const Eigen::ArrayXXf &bgErrors,
                       Eigen::ArrayXXf  &pges) const;

  /// Marks observations whose gross error probability is >= options_->rejectionThreshold
  /// as rejected by the buddy check.
  void flagRejectedObservations(
      const Variables &filtervars,
      const std::map<std::string, std::vector<float>> &grossErrProbsByVarName,
      std::vector<std::vector<bool>> &flagged) const;

 private:
  Parameters_ options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_METOFFICEBUDDYCHECK_H_
