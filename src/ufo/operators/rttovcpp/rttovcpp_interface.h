/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RTTOVCPP_RTTOVCPP_INTERFACE_H_
#define UFO_OPERATORS_RTTOVCPP_RTTOVCPP_INTERFACE_H_

#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"

#include "rttov/wrapper/Profile.h"
#include "rttov/wrapper/RttovSafe.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;
}

namespace ufo {

void rttovcpp_interface(const GeoVaLs &, const ioda::ObsSpace & odb_,
                        rttov::RttovSafe & aRttov_, const std::string CoefFileName,
                        const std::vector<int> channels_, std::size_t & nlevels,
                        std::vector<bool> & skip_profile);

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RTTOVCPP_RTTOVCPP_INTERFACE_H_
