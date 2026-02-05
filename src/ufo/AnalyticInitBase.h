/*
 * (C) Copyright 2020-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "oops/util/AssociativeContainers.h"

#include "ufo/GeoVaLs.h"
#include "ufo/SampledLocations.h"

namespace ufo {

/// \brief Initializes GeoVaLs with analytic formula
class AnalyticInitBase {
 public:
  AnalyticInitBase() = default;
  virtual ~AnalyticInitBase() = default;

/*! \brief Fill GeoVaLs with values computed by analytic function.
 *
 * \details **AnalyticInit::fillGeoVaLs** was introduced in May, 2018 (initially
 * as a GeoVaLs constructor) for use with the interpolation test. The interpolation test
 * requires an initialization of a GeoVaLs object based on the same analytic
 * formulae used for the State initialization (see test::TestStateInterpolation()
 * for further information).  This in turn requires information about the
 * vertical profile in addition to the latitude and longitude positional
 * information in the SampledLocations object.  Currently, this information
 * about the vertical profile is obtained from an existing GeoVaLs object
 * (passed as *gvals*) that represents the output of the State::interpolate()
 * method.
 *
 * \date May, 2018: created as a constructor (M. Miesch, JCSDA)
 * \date June, 2018: moved to a method (M. Miesch, JCSDA)
 *
 * \sa test::TestStateInterpolation()
 */
  virtual void fillGeoVaLs(const SampledLocations &, GeoVaLs &) const = 0;
};

// -----------------------------------------------------------------------------

/// A factory creating instances of concrete subclasses of AnalyticInitBase.
class AnalyticInitFactory {
 public:
  static std::unique_ptr<AnalyticInitBase> create(const eckit::Configuration &);

  /// \brief Return the names of all methods that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~AnalyticInitFactory() = default;

 protected:
  explicit AnalyticInitFactory(const std::string &);

 private:
  virtual std::unique_ptr<AnalyticInitBase> make(const eckit::Configuration &) = 0;

  static std::map < std::string, AnalyticInitFactory * > & getMakers() {
    static std::map < std::string, AnalyticInitFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

/// \brief A subclass of AnalyticInitFactory able to create instances of T (a concrete subclass of
/// AnalyticInitBase).
template<typename T>
class AnalyticInitMaker : public AnalyticInitFactory {
  std::unique_ptr<AnalyticInitBase> make(const eckit::Configuration & conf) override {
    return std::make_unique<T>(conf);
  }

 public:
  explicit AnalyticInitMaker(const std::string & name) : AnalyticInitFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo
