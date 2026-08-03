/*
 * (C) Copyright 2026 NOAA/OAR/GSL
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_SHAREDLISTCHECK_H_
#define UFO_FILTERS_SHAREDLISTCHECK_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/FilterParametersBase.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

/// \brief Parameters for a single entry in a compound check.
class CompoundCheckEntryParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CompoundCheckEntryParameters, Parameters)

 public:
  /// Observation variable to read for this compound match field.
  oops::RequiredParameter<std::string> obsVariable{"obs variable", this};

  /// Dictionary key in the list entries whose sublist to use for matching.
  oops::RequiredParameter<std::string> sublistToUse{"sublist to use", this};
};

/// \brief Parameters for SharedListCheck
class SharedListCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(SharedListCheckParameters, FilterParametersBase)

 public:
  /// Path to the shared list YAML file.
  oops::RequiredParameter<std::string> sharedListFile{"shared list file", this};

  /// Name of the list within the YAML file to look up.
  oops::RequiredParameter<std::string> listToUse{"list to use", this};

  /// Observation variable to match against the list entries (simple matching).
  /// Not needed when compound check is used.
  oops::OptionalParameter<std::string> variableToCheck{"variable to check", this};

  /// Compound check: match observations against list entries using multiple fields. Each entry
  ///  specifies an obs variable and the corresponding list referred to by a dictionary key
  oops::OptionalParameter<std::vector<CompoundCheckEntryParameters>> compoundCheck{
      "compound check", this};

  /// Flag mode:
  /// - "flag matched": flag observations that ARE in the list
  /// - "flag unmatched": flag observations that are NOT in the list
  oops::RequiredParameter<std::string> flagMode{"flag mode", this};
};

/// \brief Screen observations against shared lookup lists.
///
/// Supports aircraft reject lists, mesonet use lists, trusted provider lists, etc.
/// Lists are loaded from YAML files and cached in a singleton store — multiple filters
/// referencing the same file only read it once.
///
class SharedListCheck : public FilterBase {
 public:
  typedef SharedListCheckParameters Parameters_;

  SharedListCheck(ioda::ObsSpace &, const Parameters_ &,
                  ioda::ObsDataVector<int> &,
                  ioda::ObsDataVector<float> &);
  ~SharedListCheck();

 private:
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override { return QCflags::black; }
  void print(std::ostream &) const override;

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_SHAREDLISTCHECK_H_
