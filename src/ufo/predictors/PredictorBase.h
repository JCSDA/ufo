/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_PREDICTORBASE_H_
#define UFO_PREDICTORS_PREDICTORBASE_H_

#include <Eigen/Core>

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>
#include <boost/noncopyable.hpp>

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;
  class ObsBias;

// -----------------------------------------------------------------------------
/// Base class for predictor parameters
class PredictorParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(PredictorParametersBase, Parameters)

 public:
  /// \brief Predictor name.
  oops::RequiredParameter<std::string> name{"name", this};
};

// -----------------------------------------------------------------------------
/// Concrete implementation of PredictorParametersBase with no new parameters. Useful for predictors
/// that don't take any options.
class EmptyPredictorParameters : public PredictorParametersBase {
  OOPS_CONCRETE_PARAMETERS(EmptyPredictorParameters, PredictorParametersBase)

  // no other parameters needed
};

// -----------------------------------------------------------------------------
/// Base class for computing predictors
///
/// Note: each concrete implementation should typedef `Parameters_` to the name of a subclass of
/// PredictorParametersBase encapsulating its configuration options. It should also provide
/// a constructor with the following signature:
///
///     PredictorBase(const Parameters_ &, const oops::ObsVariables &);
class PredictorBase : private boost::noncopyable {
 public:
  explicit PredictorBase(const PredictorParametersBase &, const oops::ObsVariables &);
  virtual ~PredictorBase() = default;

  /// compute the predictor
  virtual void compute(const ioda::ObsSpace &,
                       const GeoVaLs &,
                       const ObsDiagnostics &,
                       const ObsBias &,
                       ioda::ObsVector &) const = 0;

  /// geovars names required to compute the predictor
  const oops::Variables & requiredGeovars() const {return geovars_;}

  /// hdiags names required to compute the predictor
  const oops::ObsVariables & requiredHdiagnostics() const {return hdiags_;}

  /// predictor name
  std::string & name() {return func_name_;}
  const std::string & name() const {return func_name_;}

 protected:
  oops::ObsVariables vars_;      ///<  variables that will be bias-corrected using this predictor
  oops::Variables geovars_;      ///<  required GeoVaLs
  oops::ObsVariables hdiags_;    ///<  required ObsDiagnostics

 private:
  std::string func_name_;        ///<  predictor name
};

typedef std::vector<std::shared_ptr<PredictorBase>> Predictors;

// -----------------------------------------------------------------------------

/// Predictor Factory
class PredictorFactory {
 public:
  /// \brief Create and return a new predictor.
  ///
  /// The predictor type is determined by the \c name attribute of \p parameters.
  /// \p parameters must be an instance of the subclass of PredictorParametersBase
  /// associated with that predictor type, otherwise an exception will be thrown.
  static std::unique_ptr<PredictorBase> create(const PredictorParametersBase &parameters,
                                               const oops::ObsVariables &vars);

  /// \brief Create and return an instance of the subclass of PredictorParametersBase
  /// storing parameters of predictors of the specified type.
  static std::unique_ptr<PredictorParametersBase> createParameters(const std::string &name);

  /// \brief Return the names of all predictors that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  /// \brief Return true if a maker has been registered for a predictor of type \p name.
  static bool predictorExists(const std::string &name);

  virtual ~PredictorFactory() = default;

 protected:
  /// \brief Register a maker able to create predictors of type \p name.
  explicit PredictorFactory(const std::string &name);

 private:
  virtual std::unique_ptr<PredictorBase> make(const PredictorParametersBase &,
                                              const oops::ObsVariables &) = 0;

  virtual std::unique_ptr<PredictorParametersBase> makeParameters() const = 0;

  static std::map < std::string, PredictorFactory * > & getMakers() {
    static std::map < std::string, PredictorFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class PredictorMaker : public PredictorFactory {
  typedef typename T::Parameters_ Parameters_;

  std::unique_ptr<PredictorBase> make(const PredictorParametersBase& parameters,
                                      const oops::ObsVariables & vars) override {
    const auto &stronglyTypedParameters = dynamic_cast<const Parameters_&>(parameters);
    return boost::make_unique<T>(stronglyTypedParameters, vars);
  }

  std::unique_ptr<PredictorParametersBase> makeParameters() const override {
    return boost::make_unique<Parameters_>();
  }

 public:
  explicit PredictorMaker(const std::string & name)
    : PredictorFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_PREDICTORBASE_H_
