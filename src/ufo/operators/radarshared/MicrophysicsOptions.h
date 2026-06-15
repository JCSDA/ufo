/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARSHARED_MICROPHYSICSOPTIONS_H_
#define UFO_OPERATORS_RADARSHARED_MICROPHYSICSOPTIONS_H_

#include "oops/util/parameters/ParameterTraits.h"

namespace ufo {

/// Microphysics scheme options for radar operators
/// This enum is shared across all radar operators (PPRO, DirectZDA, etc.)
/// Each operator may only support a subset of these options.
enum class MicrophysicsOption {
  THOMPSON,  // Supported by: PPRO, DirectZDA
  WSM6,      // Supported by: PPRO
  NSSL,      // Supported by: PPRO, DirectZDA
  LIN,       // Supported by: DirectZDA
  GFDL,      // Supported by: DirectZDA
  TCWA2      // Supported by: PPRO
};

struct MicrophysicsOptionParameterTraitsHelper {
  typedef MicrophysicsOption EnumType;
  static constexpr char enumTypeName[] = "MicrophysicsOption";
  static constexpr util::NamedEnumerator<MicrophysicsOption> namedValues[] = {
    { MicrophysicsOption::THOMPSON, "Thompson" },
    { MicrophysicsOption::WSM6, "WSM6" },
    { MicrophysicsOption::NSSL, "NSSL" },
    { MicrophysicsOption::LIN, "Lin" },
    { MicrophysicsOption::GFDL, "GFDL" },
    { MicrophysicsOption::TCWA2, "TCWA2" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::MicrophysicsOption> :
    public EnumParameterTraits<ufo::MicrophysicsOptionParameterTraitsHelper>
{};

}  // namespace oops

#endif  // UFO_OPERATORS_RADARSHARED_MICROPHYSICSOPTIONS_H_

