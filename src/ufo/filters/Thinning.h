/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_THINNING_H_
#define UFO_FILTERS_THINNING_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the Thinning filter.
class ThinningParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ThinningParameters, FilterParametersBase)

 public:
  /// (Approximate) fraction of observations to be thinned.
  oops::RequiredParameter<float> amount{"amount", this};

  /// Seed used to initialize the random number generator (if it has not been initialized
  /// before). If not set, the seed is derived from the calendar time.
  oops::OptionalParameter<int> randomSeed{"random seed", this};

  /// Index of the ensemble member.
  oops::Parameter<int> member{"member", 0, this};
};

/// \brief Randomly thin a given percentage of observations.
///
/// See ThinningParameters for the documentation of the parameters controlling this filter.

class Thinning : public FilterBase,
                 private util::ObjectCounter<Thinning> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ThinningParameters Parameters_;

  static const std::string classname() {return "ufo::Thinning";}

  Thinning(ioda::ObsSpace &, const Parameters_ &,
           std::shared_ptr<ioda::ObsDataVector<int> >,
           std::shared_ptr<ioda::ObsDataVector<float> >);
  ~Thinning();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::thinned;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_THINNING_H_
