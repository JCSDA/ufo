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

#include "ioda/ObsDataVector.h"
#include "ufo/filters/Variable.h"

namespace ioda {
class Distribution;
template <typename DATATYPE> class ObsDataVector;
class ObsSpace;
}

namespace util {
class DateTime;
}

namespace ufo {

class Variables;
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
/// Depending on which of these cases applies, create an ObservationAccessor object by
/// calling one of the
/// ObsAccessor::toAllObservations(),
/// ObsAccessor::toObservationsSplitIntoIndependentGroupsByRecordId(),
/// ObsAccessor::toObservationsSplitIntoIndependentGroupsByVariable() or
/// ObsAccessor::toSingleObservationsSplitIntoIndependentGroupsByVariable()
/// static functions.
/// The ObsAccessor will then determine whether each independent group consists of
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
  ~ObsAccessor() = default;
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

  /// \brief Create an accessor to the collection of observations held in \p obsdb, assuming that
  /// each record is treated as a single observation.
  static ObsAccessor toSingleObservationsSplitIntoIndependentGroupsByVariable(
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
  /// \param filtervars
  ///   List of filter variables.
  ///
  /// \param retentionCandidateIfAnyFilterVariablePassedQC
  ///   Boolean switch to treat an observation as a candidate for retention if any filter variable
  ///   has not been rejected. By default this is true; if false, the observation is only treated
  ///   as a candidate for retention if all filter variables have passed QC.
  ///
  /// An observation location is treated as valid if (a) it has been selected by the \c where
  /// clause and (b) its QC flag(s) for (some/all) filtered variable(s) are set to \c pass
  /// (see below).
  ///
  /// If each independent group of observations is stored entirely on a single MPI rank, the
  /// returned vector contains local IDs of valid observation locations held on the current rank
  /// only. Otherwise the vector contains global IDs of valid locations held on all ranks, with IDs
  /// from 0 to nlocs(0) - 1 corresponding to locations held on rank 0, IDs from nlocs(0) to
  /// nlocs(0) + nlocs(1) - 1 corresponding to locations held on rank 1 and so on, where nlocs(i)
  /// denotes the number of locations held on ith rank.
  ///
  /// If there is more than one filtered variable, and their QC flags differ, there is a choice
  /// as to whether to treat observation locations as valid (i) where none of the filtered variables
  /// have so far been rejected, or (ii) where at least one of these variables has not yet been
  /// rejected. The latter choice (ii) is the default, configurable via the switch
  /// \c retentionCandidateIfAnyFilterVariablePassedQC.
  std::vector<size_t> getValidObservationIds(const std::vector<bool> &apply,
                                   const ioda::ObsDataVector<int> &flags,
                                   const Variables &filtervars,
                                   bool retentionCandidateIfAnyFilterVariablePassedQC = true) const;

  /// \brief Return the IDs of both flagged and unflagged observation locations selected by the
  /// where clause.
  ///
  /// \param apply
  ///   Vector whose ith element is set to true if ith observation location held on the current
  ///   MPI rank was selected by the \c where clause in the filter's configuration.
  ///
  /// An observation location is treated as valid if it has been selected by the \c where
  /// clause.
  ///
  /// If each independent group of observations is stored entirely on a single MPI rank, the
  /// returned vector contains local IDs of valid observation locations held on the current rank
  /// only. Otherwise the vector contains global IDs of valid locations held on all ranks, with IDs
  /// from 0 to nlocs(0) - 1 corresponding to locations held on rank 0, IDs from nlocs(0) to
  /// nlocs(0) + nlocs(1) - 1 corresponding to locations held on rank 1 and so on, where nlocs(i)
  /// denotes the number of locations held on ith rank.
  std::vector<size_t> getValidObservationIds(const std::vector<bool> &apply) const;

  /// \brief Get valid (non-missing, where-included) obs indices for a given profile.
  ///
  /// \param iProfile
  ///   Profile index, e.g. from looping through the vector returned by
  ///   \c ioda::ObsSpace::recidx_all_recnums()
  ///
  /// \param apply
  ///   Vector whose ith element is set to true if ith observation location held on the current
  ///   MPI rank was selected by the \c where clause in the filter's configuration.
  ///
  /// \param flags
  ///   An ObsDataVector holding the QC flags (set by any filters run previously)
  ///   of observations held on the current MPI rank.
  ///
  /// \param filtervars
  ///   List of filter variables.
  ///
  /// \param retentionCandidateIfAnyFilterVariablePassedQC
  ///   Boolean switch to treat an observation as a candidate for retention if any filter variable
  ///   has not been rejected. By default this is true; if false, the observation is only treated
  ///   as a candidate for retention if all filter variables have passed QC.
  ///
  /// An observation location is treated as valid if (a) it has been selected by the \c where
  /// clause and (b) its QC flag(s) for (some/all) filtered variable(s) are set to \c pass
  /// (see below). This function is only appropriate for observations that are grouped into
  /// records. Thus the returned vector contains local IDs of valid observation locations held
  /// on the current rank only.
  const std::vector<size_t> getValidObsIdsInProfile(const size_t & iProfile,
                                    const std::vector<bool> & apply,
                                    const ioda::ObsDataVector<int> &flags,
                                    const Variables & filtervars,
                                    bool retentionCandidateIfAnyFilterVariablePassedQC) const;

  /// \brief Return a boolean vector indicating whether each location was selected by the
  /// \c where clause in the filter's configuration.
  /// If each independent group of observations is stored entirely on a single MPI rank
  /// then this vector will be determined separately for each rank.
  /// Otherwise, this vector will be concatenated across all ranks.
  ///
  /// \param apply
  ///   Vector whose ith element is set to true if ith observation location held on the current
  ///   MPI rank was selected by the \c where clause in the filter's configuration.
  std::vector<bool> getGlobalApply(const std::vector<bool> &apply) const;

  /// \brief Return the values of the specified variable at successive observation locations.
  ///
  /// If each independent group of observations is stored entirely on a single MPI rank, the
  /// returned vector contains values observed at locations held on the current rank only.
  /// Otherwise the vector is a concatenation of vectors obtained on all ranks.
  std::vector<int> getIntVariableFromObsSpace(const std::string &group,
                                              const std::string &variable) const;
  std::vector<bool> getBoolVariableFromObsSpace(const std::string &group,
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
  ///
  /// \param validObsIds
  ///   Indices of valid observations.
  /// \param opsCompatibilityMode
  ///   Parameter to pass to the RecursiveSplitter's constructor.
  RecursiveSplitter splitObservationsIntoIndependentGroups(
      const std::vector<size_t> &validObsIds, bool opsCompatibilityMode = false) const;

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

  /// \brief Flags observations selected by a where clause and for which at least one filter
  /// filter variable has failed QC.
  ///
  /// \param apply
  ///   Vector whose ith element is set to true if ith observation location held on the current
  ///   MPI rank was selected by the \c where clause in the filter's configuration.
  ///
  /// \param flags
  ///   An ObsDataVector holding the QC flags (set by any filters run previously)
  ///   of observations held on the current MPI rank.
  ///
  /// \param filtervars
  ///   List of filter variables.
  ///
  /// \param[inout] flagged
  ///   A vector of vectors, each with as many elements as there are observation locations on the
  ///   current MPI rank. If any filter variable has a QC rejection for jth observation location
  ///   then, for all filter variables i, flagged[i][j] will be set to true.
  void flagObservationsForAnyFilterVariableFailingQC(const std::vector<bool> &apply,
                                                    const ioda::ObsDataVector<int> &flags,
                                                    const ufo::Variables &filtervars,
                                                    std::vector<std::vector<bool> > &flagged) const;

 private:
  /// Specifies the observation grouping by a category variable.
  /// NOTHING: no category variable used.
  /// RECORD_ID: the category variable was also used to divide the ObsSpace into records.
  /// VARIABLE: the category variable was not used to divide the ObsSpace into records.
  /// SINGLE_OBS: records are treated as single obs, in which case the category variable
  /// may or may not have been used. If it was used, the behaviour is the same regardless of whether
  /// the category variable was used to divide the ObsSpace into records.
  enum class GroupBy { NOTHING, RECORD_ID, VARIABLE, SINGLE_OBS };

  /// Private constructor. Construct instances of this class by calling toAllObservations(),
  /// toObservationsSplitIntoIndependentGroupsByRecordId(),
  /// toObservationsSplitIntoIndependentGroupsByVariable() or
  /// toSingleObservationsSplitIntoIndependentGroupsByVariable()
  /// instead.
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
  std::shared_ptr<const ioda::Distribution> obsDistribution_;

  GroupBy groupBy_;
  boost::optional<Variable> categoryVariable_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSACCESSOR_H_
