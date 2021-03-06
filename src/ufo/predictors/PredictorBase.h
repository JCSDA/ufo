/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_PREDICTORBASE_H_
#define UFO_PREDICTORS_PREDICTORBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// Base class for computing predictors

class PredictorBase : private boost::noncopyable {
 public:
  explicit PredictorBase(const eckit::Configuration &, const oops::Variables &);
  virtual ~PredictorBase() = default;

  /// compute the predictor
  virtual void compute(const ioda::ObsSpace &,
                       const GeoVaLs &,
                       const ObsDiagnostics &,
                       ioda::ObsVector &) const = 0;

  /// geovars names required to compute the predictor
  const oops::Variables & requiredGeovars() const {return geovars_;}

  /// hdiags names required to compute the predictor
  const oops::Variables & requiredHdiagnostics() const {return hdiags_;}

  /// predictor name
  std::string & name() {return func_name_;}
  const std::string & name() const {return func_name_;}

 protected:
  oops::Variables vars_;         ///<  variables that will be bias-corrected using this predictor
  oops::Variables geovars_;      ///<  required GeoVaLs
  oops::Variables hdiags_;       ///<  required ObsDiagnostics

 private:
  std::string func_name_;        ///<  predictor name
};

typedef std::vector<std::shared_ptr<PredictorBase>> Predictors;

// -----------------------------------------------------------------------------

/// Predictor Factory
class PredictorFactory {
 public:
  static PredictorBase * create(const eckit::Configuration &, const oops::Variables &);
  virtual ~PredictorFactory() = default;
  static bool predictorExists(const std::string &);
 protected:
  explicit PredictorFactory(const std::string &);
 private:
  virtual PredictorBase * make(const eckit::Configuration &, const oops::Variables &) = 0;
  static std::map < std::string, PredictorFactory * > & getMakers() {
    static std::map < std::string, PredictorFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class PredictorMaker : public PredictorFactory {
  virtual PredictorBase * make(const eckit::Configuration & conf, const oops::Variables & vars)
    { return new T(conf, vars); }
 public:
  explicit PredictorMaker(const std::string & name)
    : PredictorFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_PREDICTORBASE_H_
