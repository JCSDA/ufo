/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * This filter implements the super refraction check used in Naval Research Laboratory
 * Created following the advice of Ben Ruston (JCSDA/UCAR), Hui Christophersen (NRL)
 * and Neill Bowler (Met Office)
 * This method checks if the difference of the maximum and minimum values of the
 * simulated bending angles is larger than a given threshold in any 1km vertical bin
 */

#ifndef UFO_FILTERS_PROFILEMAXDIFFERENCECHECK_H_
#define UFO_FILTERS_PROFILEMAXDIFFERENCECHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// GNSSRO super refraction check
class ProfileMaxDifferenceCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ProfileMaxDifferenceCheckParameters, FilterParametersBase)

 public:
  oops::Parameter<float> threshold{"threshold",
                      "The maximum allowed difference in the checked bin. "
                      "Values greater than this value are marked as having super-refraction.",
                       0.005, this};
  oops::Parameter<Variable> variable{"variable",
                      "The variable to be checked. ",
                       Variable{"hofx/bendingAngle"},  this};
  oops::Parameter<float> binSize{"bin size",
                      "The vertical bin size (unit: meter) within which "
                      "the difference is being calculated.",
                       1000, this};
  oops::Parameter<Variable> observationHeight{"observation height",
                      "The vertical coordinate. ",
                       Variable{"MetaData/impactHeightRO"},
                       this};
  oops::Parameter<float> maxCheckHeight{"max check height",
                       "The height below what the check is applied (unit: meter). ",
                        18000, this};
};

class ProfileMaxDifferenceCheck : public FilterBase,
      private util::ObjectCounter<ProfileMaxDifferenceCheck> {
 public:
  typedef ProfileMaxDifferenceCheckParameters Parameters_;
  static const std::string classname() {return "ufo::ProfileMaxDifferenceCheck";}

  ProfileMaxDifferenceCheck(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ProfileMaxDifferenceCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::profile;}
  Parameters_ parameters_;
  int calcMaxDifference(const std::vector<float> &,
                        const std::vector<float> &)
                        const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PROFILEMAXDIFFERENCECHECK_H_
