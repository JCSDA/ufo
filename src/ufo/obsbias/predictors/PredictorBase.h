/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIAS_PREDICTORS_PREDICTORBASE_H_
#define UFO_OBSBIAS_PREDICTORS_PREDICTORBASE_H_

#include <Eigen/Dense>

#include <map>
#include <set>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"

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
  explicit PredictorBase(const eckit::Configuration &);
  virtual ~PredictorBase() {}

  /// compute the predictor
  virtual void compute(const ioda::ObsSpace &,
                       const GeoVaLs &,
                       const ObsDiagnostics &,
                       const std::vector<int> &,
                       Eigen::MatrixXd &) const = 0;

  /// geovars names required to compute the predictor
  const std::vector<std::string> & requiredGeovars() const {return geovars_;}

  /// update geovars required to compute the predictor
  void updateGeovars(const std::vector<std::string>);

  /// hdiags names required to compute the predictor
  const std::vector<std::string> & requiredHdiagnostics() const {return hdiags_;}

  /// update hdiags names required to compute the predictor
  void updateHdiagnostics(const std::vector<std::string>);

  /// predictor name
  std::string & name() {return func_name_;}
  const std::string & name() const {return func_name_;}

 private:
  std::string func_name_;         ///<  predictor name
  std::vector<std::string> geovars_;    ///<  required GeoVaLs
  std::vector<std::string> hdiags_;     ///<  required ObsDiagnostics
};

// -----------------------------------------------------------------------------

/// Predictor Factory
class PredictorFactory {
 public:
  static PredictorBase * create(const eckit::Configuration &);
  virtual ~PredictorFactory() { getMakers().clear(); }
  static bool predictorExists(const std::string &);
 protected:
  explicit PredictorFactory(const std::string &);
 private:
  virtual PredictorBase * make(const eckit::Configuration &) = 0;
  static std::map < std::string, PredictorFactory * > & getMakers() {
    static std::map < std::string, PredictorFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class PredictorMaker : public PredictorFactory {
  virtual PredictorBase * make(const eckit::Configuration & conf)
    { return new T(conf); }
 public:
  explicit PredictorMaker(const std::string & name)
    : PredictorFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_PREDICTORS_PREDICTORBASE_H_
