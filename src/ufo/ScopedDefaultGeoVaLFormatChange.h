/*
 * (C) Crown copyright 2023, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SCOPEDDEFAULTGEOVALFORMATCHANGE_H_
#define UFO_SCOPEDDEFAULTGEOVALFORMATCHANGE_H_

#include "GeoVaLs.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

/// \brief Changes the default GeoVaL format within a scoped context. The original default format
/// is restored when the object goes out of scope.
class ScopedDefaultGeoVaLFormatChange {
 public:
  /// \brief Sets the default GeoVaL format for the GeoVaLs object `geovals` to `format` and
  /// stores the original default format in an internal variable.
  ScopedDefaultGeoVaLFormatChange(const GeoVaLs & geovals, GeoVaLFormat format)
    : geovals_(const_cast<GeoVaLs &>(geovals)), originalFormat_(geovals.defaultFormat()) {
    oops::Log::trace() << "ScopedDefaultGeoVaLFormatChange constructor setting default format to "
                       << format << std::endl;
    geovals_.setDefaultFormat(format);
  }

  /// \brief Restores the original default format of the GeoVaLs object passed to the constructor.
  ~ScopedDefaultGeoVaLFormatChange() {
    oops::Log::trace() << "ScopedDefaultGeoVaLFormatChange destructor resetting default format to "
                       << originalFormat_ << std::endl;
    geovals_.setDefaultFormat(originalFormat_);
  }

  ScopedDefaultGeoVaLFormatChange(const ScopedDefaultGeoVaLFormatChange &) = delete;
  ScopedDefaultGeoVaLFormatChange(ScopedDefaultGeoVaLFormatChange &&) = delete;
  ScopedDefaultGeoVaLFormatChange& operator=(const ScopedDefaultGeoVaLFormatChange &) = delete;
  ScopedDefaultGeoVaLFormatChange& operator=(ScopedDefaultGeoVaLFormatChange &&) = delete;

 private:
  GeoVaLs & geovals_;
  GeoVaLFormat originalFormat_;
};

}  // namespace ufo

#endif  // UFO_SCOPEDDEFAULTGEOVALFORMATCHANGE_H_
