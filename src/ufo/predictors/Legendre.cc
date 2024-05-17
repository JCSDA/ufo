/*
 * (C) Copyright 2020-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/predictors/Legendre.h"

namespace ufo {

static PredictorMaker<Legendre> makerFuncLegendre_("legendre");

// -----------------------------------------------------------------------------

Legendre::Legendre(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    order_(parameters.order.value().value_or(1)),
    nscan_(parameters.numScanPositions) {
  if (parameters.order.value() != boost::none) {
    // override the predictor name to distinguish between Legendre predictors of different orders
    name() = name() + "_order_" + std::to_string(order_);
  }
}

// -----------------------------------------------------------------------------

void Legendre::compute(const ioda::ObsSpace & odb,
                        const GeoVaLs &,
                        const ObsDiagnostics &,
                        const ObsBias &,
                        ioda::ObsVector & out) const {
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.nlocs() == nlocs);

  // retrieve the sensor scan position and total number of scan positions (const for inst type)
  std::vector<int> scan_position(nlocs, 0);
  std::vector<double> LegPoly(order_+1, 0);

  odb.get_db("MetaData", "sensorScanPosition", scan_position);
  const std::size_t nvars = vars_.size();
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    double xscan{-1.0 + 2.0 * (scan_position[jl] - 1) / (nscan_ - 1)};
    // Transformed variable for the scan position in the range -1 to 1.
    // Calculate Legendre Polynomial for current scan position
    LegPoly[0] = 1.0;
    LegPoly[1] = xscan;
    for (std::size_t iorder=1; iorder < order_; ++iorder) {
        LegPoly[iorder+1] = ((2*iorder+1)*xscan*LegPoly[iorder]
        -iorder*LegPoly[iorder-1])/(iorder+1);
    }
    for (std::size_t iorder=0; iorder <= order_; ++iorder) {
        LegPoly[iorder] = sqrt(2*iorder+1)*LegPoly[iorder];
    }
    for (std::size_t jb = 0; jb < nvars; ++jb) {
        out[jl*nvars+jb] = LegPoly[order_];
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
