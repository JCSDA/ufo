/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SUPEROB_SUPEROBBASE_H_
#define UFO_SUPEROB_SUPEROBBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/HasParameters_.h"

#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief SuperOb parameters base class.

class SuperObParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(SuperObParametersBase, Parameters)
};

class GenericSuperObParameters : public SuperObParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericSuperObParameters, SuperObParametersBase)
};


/// \brief SuperOb base class.
///
/// Subclasses of this class must implement the `computeSuperOb` method.
class SuperObBase {
 public:
  explicit SuperObBase(const SuperObParametersBase &,
                       ioda::ObsSpace &,
                       const std::vector<bool> &,
                       const Variables &,
                       const ioda::ObsDataVector<int> &,
                       std::vector<std::vector<bool>> &);
  virtual ~SuperObBase() {}

  /// Run the chosen superobbing algorithm on each record in the data set.
  void runAlgorithm() const;

  /// Compute superob for each record.
  virtual void computeSuperOb(const std::vector<std::size_t> &,
                              const std::vector<float> &,
                              const std::vector<float> &,
                              const ioda::ObsDataRow<int> &,
                              std::vector<float> &,
                              std::vector<bool> &) const = 0;

 private:
  ioda::ObsSpace & obsdb_;
  const std::vector<bool> apply_;
  const Variables & filtervars_;
  const ioda::ObsDataVector<int> & flags_;
  std::vector<std::vector<bool>> & flagged_;
};

/// SuperOb factory.
class SuperObFactory {
 public:
  static std::unique_ptr<SuperObBase> create(const std::string &,
                                             const SuperObParametersBase &,
                                             ioda::ObsSpace &,
                                             const std::vector<bool> &,
                                             const Variables &,
                                             const ioda::ObsDataVector<int> &,
                                             std::vector<std::vector<bool>> &);

  static std::unique_ptr<SuperObParametersBase> createParameters(const std::string &name);

  virtual ~SuperObFactory() = default;

 protected:
  explicit SuperObFactory(const std::string &);

 private:
  virtual std::unique_ptr<SuperObBase> make(const SuperObParametersBase &,
                                            ioda::ObsSpace &,
                                            const std::vector<bool> &,
                                            const Variables &,
                                            const ioda::ObsDataVector<int> &,
                                            std::vector<std::vector<bool>> &) = 0;

  virtual std::unique_ptr<SuperObParametersBase> makeParameters() const = 0;

  static std::map <std::string, SuperObFactory*> & getMakers() {
    static std::map <std::string, SuperObFactory*> makers_;
    return makers_;
  }
};

template<class T>
class SuperObMaker : public SuperObFactory {
  typedef oops::TParameters_IfAvailableElseFallbackType_t<T, GenericSuperObParameters>
    Parameters_;

  virtual std::unique_ptr<SuperObBase>make(const SuperObParametersBase & params,
                                           ioda::ObsSpace & obsdb,
                                           const std::vector<bool> & apply,
                                           const Variables & filtervars,
                                           const ioda::ObsDataVector<int> & flags,
                                           std::vector<std::vector<bool>> & flagged) {
    const auto & stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return std::unique_ptr<SuperObBase>
      (new T(stronglyTypedParams, obsdb, apply, filtervars, flags, flagged));
  }

  std::unique_ptr<SuperObParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit SuperObMaker(const std::string & name)
    : SuperObFactory(name) {}
};

}  // namespace ufo

#endif  // UFO_SUPEROB_SUPEROBBASE_H_
