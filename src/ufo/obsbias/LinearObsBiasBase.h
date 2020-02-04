/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIAS_LINEAROBSBIASBASE_H_
#define UFO_OBSBIAS_LINEAROBSBIASBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"

#include "oops/util/Printable.h"

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsVector;
  class ObsSpace;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// Base class for linear observation bias operators

class LinearObsBiasBase : public util::Printable,
                          private boost::noncopyable {
 public:
  LinearObsBiasBase() {}
  virtual ~LinearObsBiasBase() {}

/// Linear algebra operators
  virtual void diff(const ObsBias &, const ObsBias &) = 0;
  virtual void zero() = 0;
  virtual LinearObsBiasBase & operator=(const ObsBiasIncrement &) = 0;
  virtual LinearObsBiasBase & operator+=(const ObsBiasIncrement &) = 0;
  virtual LinearObsBiasBase & operator-=(const ObsBiasIncrement &) = 0;
  virtual LinearObsBiasBase & operator*=(const double) = 0;
  virtual void axpy(const double, const ObsBiasIncrement &) = 0;
  virtual double dot_product_with(const ObsBiasIncrement &) const = 0;

/// I/O and diagnostics
  virtual void read(const eckit::Configuration &) = 0;
  virtual void write(const eckit::Configuration &) const = 0;
  virtual double norm() const = 0;

/// Bias model
  virtual void computeObsBiasTL(const GeoVaLs &, const ioda::ObsDataVector<float> &,
                                ioda::ObsVector &) const = 0;

  virtual void computeObsBiasAD(GeoVaLs &, const ioda::ObsDataVector<float> &,
                                const ioda::ObsVector &) = 0;

/// Bias parameters interface
  virtual double & operator[](const unsigned int) = 0;
  virtual const double & operator[](const unsigned int) const = 0;

  virtual const ioda::ObsSpace & obspace() const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -----------------------------------------------------------------------------

/// Linear observation bias operator Factory
class LinearObsBiasFactory {
 public:
  static LinearObsBiasBase * create(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~LinearObsBiasFactory() { getMakers().clear(); }

 protected:
  explicit LinearObsBiasFactory(const std::string &);

 private:
  virtual LinearObsBiasBase * make(const ioda::ObsSpace &, const eckit::Configuration &) = 0;
  static std::map < std::string, LinearObsBiasFactory * > & getMakers() {
    static std::map < std::string, LinearObsBiasFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class LinearObsBiasMaker : public LinearObsBiasFactory {
  virtual LinearObsBiasBase * make(const ioda::ObsSpace & obs, const eckit::Configuration & conf)
    { return new T(obs, conf); }
 public:
  explicit LinearObsBiasMaker(const std::string & name) : LinearObsBiasFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_LINEAROBSBIASBASE_H_
