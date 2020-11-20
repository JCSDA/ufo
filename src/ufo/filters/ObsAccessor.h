/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSACCESSOR_H_
#define UFO_FILTERS_OBSACCESSOR_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "oops/util/DateTime.h"
#include "ufo/filters/Variable.h"

namespace ioda {
class Distribution;
template <typename DATATYPE> class ObsDataVector;
class ObsSpace;
}

namespace ufo {

class RecursiveSplitter;

/// \brief This class provides access to observations that may be held on multiple MPI ranks.
///
/// It is used by filters that may be configured to
///
/// * process observations held on all MPI ranks as a single group
/// * process observations from each record (by definition, held on a single MPI rank)
///   independently from all others
/// * process observations with each distinct value of a particular variable (held on a single MPI
///   rank if this variable was used to group observations into records or on multiple MPI ranks
///   if not) independently from all others.
///
/// Depending on which of these cases applies, create an ObservationAccessor object by calling the
/// ObsAccessor::toAllObservations(),
/// ObsAccessor::toObservationsSplitIntoIndependentGroupsByRecordId() or the
/// ObsAccessor::toObservationsSplitIntoIndependentGroupsByVariable() static function. The
/// ObsAccessor will then determine whether each independent group consists of
/// observations held only on a single MPI rank. If so, methods such as getValidObservationIds() and
/// getIntVariableFromObsSpace() will return vectors constructed from data held only on the current
/// MPI rank (without any MPI communication); otherwise, these vectors will be constructed from
/// data obtained from all MPI ranks.
///
/// Call splitObservationsIntoIndependentGroups() to construct a RecursiveSplitter object whose
/// groups() method will return groups of observations that can be processed independently from
/// each other (according to the criterion specified when the ObsAccessor was constructed).
class ObsAccessor {
 public:
  ~ObsAccessor();  // defined in the .cc file
  ObsAccessor(const ObsAccessor &) = delete;
  ObsAccessor(ObsAccessor &&) = default;
  ObsAccessor & operator=(const ObsAccessor &) = delete;
  ObsAccessor & operator=(ObsAccessor &&) = default;

  /// \brief Create an accessor to observations from the observation space \p obsdb, assuming that
  /// the whole set of observations held on all MPI ranks must be processed together as a single
  /// group.
  static ObsAccessor toAllObservations(
      const ioda::ObsSpace &obsdb);

  /// \brief Create an accessor to the collection of observations held in \p obsdb, assuming that
  /// each record can be processed independently.
  static ObsAccessor toObservationsSplitIntoIndependentGroupsByRecordId(
      const ioda::ObsSpace &obsdb);

  /// \brief Create an accessor to the collection of observations held in \p obsdb, assuming that
  /// observations with different values of the variable \p variable can be processed independently.
  static ObsAccessor toObservationsSplitIntoIndependentGroupsByVariable(
      const ioda::ObsSpace &obsdb, const Variable &variable);

  /// \brief Return the IDs of observation locations that should be treated as valid by a filter.
  ///
  /// \param apply
  ///   Vector whose ith element is set to true if ith observation location held on the current
  ///   MPI rank was selected by the \c where clause in the filter's configuration.
  ///
  /// \param flags
  ///   An ObsDataVector holding the QC flags (set by any filters run previously)
  ///   of observations held on the current MPI rank.
  ///
  /// An observation location is treated as valid if (a) it has been selected by the \c where
  /// clause and (b) its QC flag for the first filtered variable is set to \c pass (see note
  /// below).
  ///
  /// If each independent group of observations is stored entirely on a single MPI rank, the
  /// returned vector contains local IDs of valid observation locations held on the current rank
  /// only. Otherwise the vector contains global IDs of valid locations held on all ranks, with IDs
  /// from 0 to nlocs(0) - 1 corresponding to locations held on rank 0, IDs from nlocs(0) to
  /// nlocs(0) + nlocs(1) - 1 corresponding to locations held on rank 1 and so on, where nlocs(i)
  /// denotes the number of locations held on ith rank.
  ///
  /// \note It is an open question what should be done if there's more than one filtered variable
  /// and their QC flags differ. Should only observation locations where none of the filtered
  /// variables have so far been rejected be treated as candidates for retaining? Or all locations
  /// where at least one of these variables has not yet been rejected? For the moment we duck the
  /// question and check only the QC flags of the first filtered variable, implicitly assuming that
  /// if multiple variables are thinned together, they all have the same QC flags. It might be best
  /// to make this user-configurable.
  std::vector<size_t> getValidObservationIds(const std::vector<bool> &apply,
                                             const ioda::ObsDataVector<int> &flags) const;

  /// \brief Return the values of the specified variable at successive observation locations.
  ///
  /// If each independent group of observations is stored entirely on a single MPI rank, the
  /// returned vector contains values observed at locations held on the current rank only.
  /// Otherwise the vector is a concatenation of vectors obtained on all ranks.
  std::vector<int> getIntVariableFromObsSpace(const std::string &group,
                                              const std::string &variable) const;
  std::vector<float> getFloatVariableFromObsSpace(const std::string &group,
                                                  const std::string &variable) const;
  std::vector<double> getDoubleVariableFromObsSpace(const std::string &group,
                                                    const std::string &variable) const;
  std::vector<std::string> getStringVariableFromObsSpace(const std::string &group,
                                                         const std::string &variable) const;
  std::vector<util::DateTime> getDateTimeVariableFromObsSpace(const std::string &group,
                                                              const std::string &variable) const;

  /// \brief Return the vector of IDs of records successive observation locations belong to.
  ///
  /// If each independent group of observations is stored entirely on a single MPI rank, the
  /// returned vector contains record IDs of observation locations held on the current rank
  /// only. Otherwise the vector is a concatenation of vectors obtained on all ranks.
  std::vector<size_t> getRecordIds() const;

  /// If each independent group of observations is stored entirely on a single MPI rank, return the
  /// number of observation locations held on the current rank. Otherwise return the total number
  /// of observation locations held on all ranks.
  size_t totalNumObservations() const;

  /// Construct a RecursiveSplitter object whose groups() method will return groups of observations
  /// that can be processed independently from each other (according to the criterion specified when
  /// the ObsAccessor was constructed).
  RecursiveSplitter splitObservationsIntoIndependentGroups(
      const std::vector<size_t> &validObsIds) const;

  /// \brief Update flags of observations held on the current MPI rank.
  ///
  /// \param isRejected
  ///   A vector of length totalNumObservations() whose ith element indicates if ith observation
  ///   should be rejected.
  ///
  /// \param[inout] flagged
  ///   A vector of vectors, each with as many elements as there are observation locations on the
  ///   current MPI rank. On output, flagged[i][j] will be set to true for each i if the element of
  ///   isRejected corresponding to jth observation location on the current rank is true.
  void flagRejectedObservations(const std::vector<bool> &isRejected,
                                std::vector<std::vector<bool> > &flagged) const;

 private:
  enum class GroupBy { NOTHING, RECORD_ID, VARIABLE };

  /// Private constructor. Construct instances of this class by calling toAllObservations(),
  /// toObservationsSplitIntoIndependentGroupsByRecordId() or
  /// toObservationsSplitIntoIndependentGroupsByVariable() instead.
  ObsAccessor(const ioda::ObsSpace &obsdb,
              GroupBy groupBy,
              boost::optional<Variable> categoryVariable);

  bool wereRecordsGroupedByCategoryVariable() const;

  void groupObservationsByRecordNumber(const std::vector<size_t> &validObsIds,
                                       RecursiveSplitter &splitter) const;

  void groupObservationsByCategoryVariable(const std::vector<size_t> &validObsIds,
                                           RecursiveSplitter &splitter) const;

 private:
  const ioda::ObsSpace *obsdb_;
  std::unique_ptr<const ioda::Distribution> inefficientDistribution_;
  const ioda::Distribution *obsDistribution_;

  GroupBy groupBy_;
  boost::optional<Variable> categoryVariable_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSACCESSOR_H_
