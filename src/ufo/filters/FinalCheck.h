/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_FINALCHECK_H_
#define UFO_FILTERS_FINALCHECK_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/ObsProcessorBase.h"

namespace ufo {

class FinalCheckParameters : public oops::ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(FinalCheckParameters, ObsFilterParametersBase)
  // No extra parameters needed
};

/// \brief A filter run automatically at the end of the whole sequence of filters.
///
/// It does three things:
/// - verifies that all derived simulated variables have been created and if not, throws an
///   exception
/// - sets the QC flag of all observations with missing error estimates to `missing`.
/// - sets the QC flag of all observations that have been processed but are not to
///   be assimilated.

class FinalCheck : public ObsProcessorBase,
                   private util::ObjectCounter<FinalCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef FinalCheckParameters Parameters_;

  static const std::string classname() {return "ufo::FinalCheck";}

  FinalCheck(ioda::ObsSpace & obsdb, const Parameters_ & params,
             std::shared_ptr<ioda::ObsDataVector<int>> qcflags,
             std::shared_ptr<ioda::ObsDataVector<float>> obserr);
  ~FinalCheck() override;

  void doFilter() override;

 private:
  void print(std::ostream &) const override;
};

}  // namespace ufo

#endif  // UFO_FILTERS_FINALCHECK_H_
