/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_DIAGNOSTICFLAG_H_
#define UFO_FILTERS_DIAGNOSTICFLAG_H_

namespace ufo {

/// \brief Type representing the value of a diagnostic flag.
///
/// Currently bool, but we might switch to char in future if std::vector<bool> containers prove
/// inconvenient to use.
typedef bool DiagnosticFlag;

}  // namespace ufo

#endif  // UFO_FILTERS_DIAGNOSTICFLAG_H_
