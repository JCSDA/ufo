/*
 * (C) Crown Copyright 2024 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_SUPEROB_H_
#define UFO_FILTERS_SUPEROB_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

#include "ufo/superob/SuperObBase.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

class SuperObParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SuperObParametersWrapper, oops::Parameters)
 public:
  /// Name of the superobbing algorithm.
  /// Valid names are specified using a `SuperObMaker` in subclasses of SuperObBase in the
  /// ufo/superob directory.
  oops::RequiredPolymorphicParameter<SuperObParametersBase, SuperObFactory>
    superObName{"name", this};
};

class SuperObParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(SuperObParameters, FilterParametersBase)
 public:
  /// Parameter that contains details of the algorithm to use.
  oops::RequiredParameter<SuperObParametersWrapper> algorithmParameters{"algorithm", this};
};

/// \brief Compute superobs using a specified algorithm.
///
/// This filter computes superobs from observed and simulated values.
/// The choice of algorithm can be selected using the `algorithm` parameter.
///
/// An example use of this filter is as follows:
///
///   - filter: SuperOb
///     filter variables:
///     - name: airTemperature
///     - name: windEastward
///     algorithm: mean OmB
///     action:
///       name: set
///       flag: superob
///
class SuperOb : public FilterBase,
  private util::ObjectCounter<SuperOb> {
 public:
  typedef SuperObParameters Parameters_;

  static const std::string classname() {return "ufo::SuperOb";}

  SuperOb(ioda::ObsSpace &, const Parameters_ &,
          std::shared_ptr<ioda::ObsDataVector<int> >,
          std::shared_ptr<ioda::ObsDataVector<float> >);
  ~SuperOb();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::superob;}
  Parameters_ options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_SUPEROB_H_
