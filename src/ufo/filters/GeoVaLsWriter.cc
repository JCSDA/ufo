/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/GeoVaLsWriter.h"


#include "eckit/config/LocalConfiguration.h"

#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------

GeoVaLsWriter::GeoVaLsWriter(const ioda::ObsSpace &, const Parameters_ & params,
                ObsDataPtr_<int>, ObsDataPtr_<float>)
  : config_(params.config), novars_() {}

// -----------------------------------------------------------------------------

void GeoVaLsWriter::priorFilter(const GeoVaLs & gv) {
  const double zz = sqrt(gv.dot_product_with(gv));
  oops::Log::info() << "GeoVaLsWriter norm = " << zz << std::endl;
  gv.write(config_);
}

// -----------------------------------------------------------------------------

void GeoVaLsWriter::print(std::ostream & os) const {
  os << "Filter outputting GeoVaLs" << config_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
