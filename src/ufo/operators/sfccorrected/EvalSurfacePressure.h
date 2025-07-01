/*
 * (C) Copyright 2025, UCAR Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_SFCCORRECTED_EVALSURFACEPRESSURE_H_
#define UFO_OPERATORS_SFCCORRECTED_EVALSURFACEPRESSURE_H_

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ufo/operators/sfccorrected/SurfaceOperatorBase.h"

/// Forward declarations
namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
}

namespace ufo {

class GeoVaLs;


class stationPressure_WRFDA : public SurfaceOperatorBase {
 public:
  explicit stationPressure_WRFDA(const std::string &,
                                 const Parameters_ &);
  virtual ~stationPressure_WRFDA() {}

  void simobs(const ufo::GeoVaLs &,
              const ioda::ObsSpace &,
              std::vector<float> &) const override;
  void settraj() const override;
  void TL() const override;
  void AD() const override;

 private:
  void getDataValues(const ufo::GeoVaLs &,
                     const ioda::ObsSpace &,
                     const Parameters_ &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &) const;
};

class stationPressure_UKMO : public SurfaceOperatorBase {
 public:
  explicit stationPressure_UKMO(const std::string &,
                                const Parameters_ &);
  virtual ~stationPressure_UKMO() {}

  void simobs(const ufo::GeoVaLs &,
              const ioda::ObsSpace &,
              std::vector<float> &) const override;
  void settraj() const override;
  void TL() const override;
  void AD() const override;

 private:
  void getDataValues(const ufo::GeoVaLs &,
                     const ioda::ObsSpace &,
                     const Parameters_ &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &) const;
};

class stationPressure_GSL : public SurfaceOperatorBase {
 public:
  explicit stationPressure_GSL(const std::string &,
                               const Parameters_ &);
  virtual ~stationPressure_GSL() {}

  void simobs(const ufo::GeoVaLs &,
              const ioda::ObsSpace &,
              std::vector<float> &) const override;
  void settraj() const override;
  void TL() const override;
  void AD() const override;

 private:
  void getDataValues(const ufo::GeoVaLs &,
                     const ioda::ObsSpace &,
                     const Parameters_ &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<float> &,
                     std::vector<double> &,
                     std::vector<double> &,
                     std::vector<float> &) const;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_SFCCORRECTED_EVALSURFACEPRESSURE_H_
