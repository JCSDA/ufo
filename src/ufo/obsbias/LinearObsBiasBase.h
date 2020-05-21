/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSBIAS_LINEAROBSBIASBASE_H_
#define UFO_OBSBIAS_LINEAROBSBIASBASE_H_

#include <Eigen/Core>

#include <map>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"

#include "oops/util/Printable.h"

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
  LinearObsBiasBase(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~LinearObsBiasBase() {}

  // Linear algebra operators
  virtual void diff(const ObsBias &, const ObsBias &) = 0;
  virtual void zero() = 0;
  virtual LinearObsBiasBase & operator=(const ObsBiasIncrement &) = 0;
  virtual LinearObsBiasBase & operator+=(const ObsBiasIncrement &) = 0;
  virtual LinearObsBiasBase & operator-=(const ObsBiasIncrement &) = 0;
  virtual LinearObsBiasBase & operator*=(const double) = 0;
  virtual void axpy(const double, const ObsBiasIncrement &) = 0;
  virtual double dot_product_with(const ObsBiasIncrement &) const = 0;

  // I/O and diagnostics
  virtual void read(const eckit::Configuration &) const = 0;
  virtual void write(const eckit::Configuration &) const = 0;
  virtual double norm() const = 0;

  // Bias model
  virtual void computeObsBiasTL(const GeoVaLs &,
                                const Eigen::MatrixXd &,
                                ioda::ObsVector &) const = 0;

  virtual void computeObsBiasAD(GeoVaLs &,
                                const Eigen::MatrixXd &,
                                const ioda::ObsVector &) = 0;

  // Bias parameters interface
  virtual double & operator[](const unsigned int) = 0;
  virtual const double & operator[](const unsigned int) const = 0;

  virtual const ioda::ObsSpace & obsspace() const = 0;

  // Utilities
  const eckit::mpi::Comm & mpi_comm() const {return mpi_comm_;}

 private:
  virtual void print(std::ostream &) const = 0;
  const eckit::mpi::Comm & mpi_comm_;
};

// -----------------------------------------------------------------------------

/// Linear observation bias operator Factory
class LinearObsBiasFactory {
 public:
  static LinearObsBiasBase * create(const ioda::ObsSpace &,
                                    const eckit::Configuration &,
                                    const std::vector<std::string> &,
                                    const std::vector<int> &);
  virtual ~LinearObsBiasFactory() { getMakers().clear(); }

 protected:
  explicit LinearObsBiasFactory(const std::string &);

 private:
  virtual LinearObsBiasBase * make(const ioda::ObsSpace &,
                                   const eckit::Configuration &,
                                   const std::vector<std::string> &,
                                   const std::vector<int> &) = 0;
  static std::map < std::string, LinearObsBiasFactory * > & getMakers() {
    static std::map < std::string, LinearObsBiasFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class LinearObsBiasMaker : public LinearObsBiasFactory {
  virtual LinearObsBiasBase * make(const ioda::ObsSpace & obs,
                                   const eckit::Configuration & conf,
                                   const std::vector<std::string> & preds,
                                   const std::vector<int> & jobs)
    { return new T(obs, conf, preds, jobs); }
 public:
  explicit LinearObsBiasMaker(const std::string & name) : LinearObsBiasFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_LINEAROBSBIASBASE_H_
