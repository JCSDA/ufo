/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIAS_OBSBIASBASE_H_
#define UFO_OBSBIAS_OBSBIASBASE_H_

#include <Eigen/Core>

#include <map>
#include <string>
#include <vector>

#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// Base class for observation bias operators

class ObsBiasBase : public util::Printable {
 public:
  explicit ObsBiasBase(const eckit::Configuration &);
  virtual ~ObsBiasBase() {}

  // I/O and diagnostics
  virtual void read(const std::string &) = 0;
  virtual void write(const eckit::Configuration &) const = 0;
  virtual double norm() const = 0;

  // Add increments
  virtual ObsBiasBase & operator+=(const ObsBiasIncrement &) = 0;
  virtual ObsBiasBase & operator=(const ObsBias &) = 0;

  // Bias model
  virtual void computeObsBias(ioda::ObsVector &, const Eigen::MatrixXd &) const = 0;

  // save Bias Terms for QC
  virtual void saveObsBiasTerms(ioda::ObsSpace &,
                                const std::string &,
                                const Eigen::MatrixXd &) const = 0;

  // Bias parameters interface
  virtual std::size_t size() const = 0;
  virtual double & operator[](const unsigned int) = 0;

 protected:
  std::string input_filename_;
  std::string output_filename_;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

/// Observation bias operator Factory
class ObsBiasFactory {
 public:
  static ObsBiasBase * create(const eckit::Configuration &,
                              const std::vector<std::string> &,
                              const std::vector<int> &);
  virtual ~ObsBiasFactory() { getMakers().clear(); }

 protected:
  explicit ObsBiasFactory(const std::string &);

 private:
  virtual ObsBiasBase * make(const eckit::Configuration &,
                             const std::vector<std::string> &,
                             const std::vector<int> &) = 0;
  static std::map < std::string, ObsBiasFactory * > & getMakers() {
    static std::map < std::string, ObsBiasFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class ObsBiasMaker : public ObsBiasFactory {
  virtual ObsBiasBase * make(const eckit::Configuration & conf,
                             const std::vector<std::string> & preds,
                             const std::vector<int> & jobs)
    { return new T(conf, preds, jobs); }
 public:
  explicit ObsBiasMaker(const std::string & name) : ObsBiasFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_OBSBIASBASE_H_
