/*
 * (C) Copyright 2022 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_COPYFLAGSFROMEXTENDEDTOORIGINALSPACE_H_
#define UFO_FILTERS_COPYFLAGSFROMEXTENDEDTOORIGINALSPACE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Options controlling CopyFlagsFromExtendedToOriginalSpace filter
/// Parameters controlling the operation of the BackgroundCheck filter.
class CopyFlagsFromExtendedToOriginalSpaceParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(CopyFlagsFromExtendedToOriginalSpaceParameters, FilterParametersBase)

 public:
  /// Observation vertical coordinate
  oops::RequiredParameter<Variable> obsVertCoord{"observation vertical coordinate", this};
  /// Model vertical coordinate
  oops::RequiredParameter<Variable> modelVertCoord{"model vertical coordinate", this};
};

// -----------------------------------------------------------------------------

/// \brief Copies a diagnostic flag's values on model levels (in the extended space)
///  into the corresponding levels in the original (observation) space.
///
/// Example
///
///  - filter: Copy Flags From Extended To Original Space
///    where:
///    - variable:
///        name: DiagnosticFlags/flag1/variable1
///      is_false:
///    filter variables: [DiagnosticFlags/flag1/variable1]
///    model vertical coordinate: HofX/ocean_depth
///    observation vertical coordinate: DerivedObsValue/ocean_depth
///
/// For every observation level for which `DiagnosticFlags/flag1/variable1`
/// is false (unset), this will copy the value of the corresponding model level of
/// `DiagnosticFlags/flag1/variable1` into it; every observation level for which the
/// flag is set, remains unchanged; every model level remains unchanged.
///
/// Note that each model level may be linked to multiple observation levels, so one
/// model-level value may be copied into multiple observation levels. There is not yet
/// any functionality to copy values from original space to extended space; only from
/// extended to original.
///

class CopyFlagsFromExtendedToOriginalSpace : public FilterBase,
      private util::ObjectCounter<CopyFlagsFromExtendedToOriginalSpace> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef CopyFlagsFromExtendedToOriginalSpaceParameters Parameters_;

  static const std::string classname() {return "ufo::CopyFlagsFromExtendedToOriginalSpace";}

  CopyFlagsFromExtendedToOriginalSpace(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                        std::shared_ptr<ioda::ObsDataVector<int>> flags,
                        std::shared_ptr<ioda::ObsDataVector<float>> obserr);
  ~CopyFlagsFromExtendedToOriginalSpace() override;

 private:
  void applyFilter(const std::vector<bool> & apply,
                   const Variables & filtervars,
                   std::vector<std::vector<bool>> & flagged) const override;
  int qcFlag() const override {return QCflags::fguess;}
  void print(std::ostream &) const override;

 private:
  Parameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_COPYFLAGSFROMEXTENDEDTOORIGINALSPACE_H_
