/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PREDICTORS_PREDICTORBASE_H_
#define UFO_PREDICTORS_PREDICTORBASE_H_

#include <map>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/obsbias_io/ObsBiasIO.h"

namespace eckit {
  class Configuration;
  class Comm;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;
  class Record;

// -----------------------------------------------------------------------------
/// Base class for computing predictors

class PredictorBase : private boost::noncopyable {
 public:
  PredictorBase(const eckit::Configuration &,
                const std::vector<int> &,
                const std::string &,
                const eckit::mpi::Comm &);
  virtual ~PredictorBase() {}

  /// Write out for saving
  virtual void write(const eckit::Configuration &,
                     ObsBiasIO< Record > &) = 0;

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
  const std::vector<int> jobs_;
  const std::string sensor_;
  oops::Variables geovars_;
  oops::Variables hdiags_;
  const eckit::mpi::Comm & comm_;

 private:
  std::string func_name_;
};

// -----------------------------------------------------------------------------

/// Predictor Factory
class PredictorFactory {
 public:
  static PredictorBase * create(const eckit::Configuration &,
                                const std::vector<int> &,
                                const std::string &,
                                const eckit::mpi::Comm &);
  virtual ~PredictorFactory() = default;
  static bool predictorExists(const std::string &);
 protected:
  explicit PredictorFactory(const std::string &);
 private:
  virtual PredictorBase * make(const eckit::Configuration &,
                               const std::vector<int> &,
                               const std::string &,
                               const eckit::mpi::Comm &) = 0;
  static std::map < std::string, PredictorFactory * > & getMakers() {
    static std::map < std::string, PredictorFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class PredictorMaker : public PredictorFactory {
  virtual PredictorBase * make(const eckit::Configuration & conf,
                               const std::vector<int> & jobs,
                               const std::string & sensor,
                               const eckit::mpi::Comm & comm)
    { return new T(conf, jobs, sensor, comm); }
 public:
  explicit PredictorMaker(const std::string & name)
    : PredictorFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_PREDICTORS_PREDICTORBASE_H_
