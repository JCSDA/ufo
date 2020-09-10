/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <map>
#include <memory>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/VariableNames.h"

#include "ufo/utils/StringUtils.h"

namespace ufo {
  ProfileCheckValidator::ProfileCheckValidator(const ProfileConsistencyCheckParameters &options,
                                               ProfileDataHandler &profileDataHandler)
    : options_(options),
      profileDataHandler_(profileDataHandler)
  {
    // Set offsets due to C++ and Fortran array index starting values
    comparison_offsets_[ufo::VariableNames::StdLev] = 1;
    comparison_offsets_[ufo::VariableNames::SigBelow] = 1;
    comparison_offsets_[ufo::VariableNames::SigAbove] = 1;
    comparison_offsets_[ufo::VariableNames::IndStd] = 1;
    comparison_offsets_[ufo::VariableNames::LevErrors] = 1;
    comparison_offsets_[ufo::VariableNames::Indx] = 1;

    // List of checks performed
    std::vector <std::string> checks = options_.Checks.value();

    // Loop over each check and populate lists of integer and float values to compare
    for (const auto& check : checks) {
      if (check == "Basic") {
      } else if (check == "SamePDiffT") {
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumAnyErrors);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumSamePErrObs);
        valuesToCompare_int_.insert(ufo::VariableNames::qcflags_air_temperature);
      } else if (check == "Sign") {
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumAnyErrors);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumSignChange);
        valuesToCompare_int_.insert(ufo::VariableNames::qcflags_air_temperature);
      } else if (check == "UnstableLayer") {
        valuesToCompare_int_.insert(ufo::VariableNames::qcflags_air_temperature);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumAnyErrors);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumSuperadiabat);
        valuesToCompare_float_.insert(ufo::VariableNames::PBottom);
      } else if (check == "Interpolation") {
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumAnyErrors);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumInterpErrors);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumInterpErrObs);
        valuesToCompare_int_.insert(ufo::VariableNames::qcflags_air_temperature);
        valuesToCompare_int_.insert(ufo::VariableNames::NumStd);
        valuesToCompare_int_.insert(ufo::VariableNames::NumSig);
        valuesToCompare_int_.insert(ufo::VariableNames::StdLev);
        valuesToCompare_int_.insert(ufo::VariableNames::SigBelow);
        valuesToCompare_int_.insert(ufo::VariableNames::SigAbove);
        valuesToCompare_int_.insert(ufo::VariableNames::IndStd);
        valuesToCompare_int_.insert(ufo::VariableNames::LevErrors);
        valuesToCompare_float_.insert(ufo::VariableNames::tInterp);
        valuesToCompare_float_.insert(ufo::VariableNames::LogP);
      } else if (check == "Hydrostatic") {
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumAnyErrors);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_Num925Miss);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_Num100Miss);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumStdMiss);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumHydErrObs);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumIntHydErrors);
        valuesToCompare_int_.insert(ufo::VariableNames::qcflags_air_temperature);
        valuesToCompare_int_.insert(ufo::VariableNames::qcflags_geopotential_height);
        valuesToCompare_float_.insert(ufo::VariableNames::DC);
        valuesToCompare_float_.insert(ufo::VariableNames::ETol);
        valuesToCompare_float_.insert(ufo::VariableNames::D);
        valuesToCompare_float_.insert(ufo::VariableNames::E);
        valuesToCompare_int_.insert(ufo::VariableNames::HydError);
      } else if (check == "UInterp") {
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumSamePErrObs);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_NumInterpErrObs);
        valuesToCompare_int_.insert(ufo::VariableNames::qcflags_eastward_wind);
        valuesToCompare_float_.insert(ufo::VariableNames::uInterp);
        valuesToCompare_float_.insert(ufo::VariableNames::vInterp);
        valuesToCompare_int_.insert(ufo::VariableNames::NumStd);
        valuesToCompare_int_.insert(ufo::VariableNames::NumSig);
        valuesToCompare_int_.insert(ufo::VariableNames::StdLev);
        valuesToCompare_int_.insert(ufo::VariableNames::SigBelow);
        valuesToCompare_int_.insert(ufo::VariableNames::SigAbove);
        valuesToCompare_int_.insert(ufo::VariableNames::LevErrors);
        valuesToCompare_float_.insert(ufo::VariableNames::LogP);
      } else if (check == "RH") {
        valuesToCompare_int_.insert(ufo::VariableNames::qcflags_relative_humidity);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_TotCProfs);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_TotHProfs);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_TotCFlags);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_TotHFlags);
        valuesToCompare_int_.insert(ufo::VariableNames::counter_TotLFlags);
        valuesToCompare_float_.insert(ufo::VariableNames::Press);
        valuesToCompare_float_.insert(ufo::VariableNames::Temp);
        valuesToCompare_float_.insert(ufo::VariableNames::rh);
        valuesToCompare_float_.insert(ufo::VariableNames::td);
        valuesToCompare_float_.insert(ufo::VariableNames::tbk);
        valuesToCompare_float_.insert(ufo::VariableNames::rhbk);
        valuesToCompare_int_.insert(ufo::VariableNames::FlagH);
        valuesToCompare_int_.insert(ufo::VariableNames::Indx);
      }
    }
  }

  /// Comparison of single values
  template <typename T>
  void ProfileCheckValidator::compareOutput(const std::string &desc,
                                            const T val1,
                                            const T val2,
                                            const int offset,
                                            const float tol,
                                            int &n)
  {
    if (!differenceWithinTol(val1, val2 + offset, tol)) {
        oops::Log::debug() << "Mismatch for " << desc << " (OPS, this code): "
                           << val1 << ", " << val2 + offset << std::endl;
        n++;
    }
  }

  /// Comparison of vectors of values
  template <typename T>
  void ProfileCheckValidator::compareOutput(const std::string &desc,
                                            const std::vector <T> &vec1,
                                            const std::vector <T> &vec2,
                                            const int offset,
                                            const float tol,
                                            int &n)
  {
    // Do not compare vectors if at least one is empty
    if (oops::anyVectorEmpty(vec1, vec2))
      {
        if (vec1.empty())
          oops::Log::debug() << "Vector of " << desc << " in OPS output is empty" << std::endl;
        if (vec2.empty())
          oops::Log::debug() << "Vector of " << desc << " in this code is empty" << std::endl;
        return;
      }
    // Warn if vectors are different size but allow to continue
    if (!oops::allVectorsSameSize(vec1, vec2))
      {
        oops::Log::debug() << "Vectors to be compared for "
                             << desc << " are of different size (" << vec1.size()
                             << " and " << vec2.size() << "). "
                             << "Will compare entries until reaching the end of the shorter "
                             << "of the two." << std::endl;
      }

    // Compare vector elements up to the smaller of the two sizes.
    const size_t vecsize = std::min(vec1.size(), vec2.size());
    for (size_t jvec = 0; jvec < vecsize; ++jvec) {
      if (!differenceWithinTol(vec1[jvec], vec2[jvec] + offset, tol)) {
        oops::Log::debug() << "Mismatch for " << desc << "[" << jvec << "] "
                           << "(OPS, this code): " << vec1[jvec] << ", "
                           << vec2[jvec] + offset << std::endl;
        n++;
      }
    }
  }

  void ProfileCheckValidator::validate()
  {
    oops::Log::debug() << " Comparing values against OPS equivalents..." << std::endl;

    // Reset number of mismatches for this profile
    nMismatches_ = 0;

    float tol = options_.Comparison_Tol.value();  // Comparison tolerance

    // Compare integer values obtained in this code and OPS
    for (const auto& valueToCompare_int : valuesToCompare_int_) {
      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(valueToCompare_int, varname, groupname);
      std::string varname_OPS = "OPS_" + valueToCompare_int;
      if (groupname == "Counters") {
        /// Special case: OPS counters have one value per profile level,
        /// and are in the MetaData rather than the Counters group.
        /// This avoids the (default) treatment which assumes
        /// that variables in the Counters group have one value per profile.
        varname_OPS = "OPS_" + varname + "@MetaData";
      }

      // Obtain values for comparison
      const std::vector <int> &values_thiscode =
        profileDataHandler_.get<int>(valueToCompare_int);
      const std::vector <int> &values_OPS =
        profileDataHandler_.get<int>(varname_OPS);

      // Account for potential offset between values in this code and OPS
      int offset = 0;

      // Offsets due to C++ and Fortran array indices
      auto comparison_offsets_it = comparison_offsets_.find(valueToCompare_int);
      if (comparison_offsets_it != comparison_offsets_.end())
        offset = comparison_offsets_it->second;

      // Offsets due to particular counters being accumulated over profiles in OPS
      // (and not in this code). NumAnyErrors is not included.
      if (groupname == "Counters" && varname != "NumAnyErrors")
        offset = cumulativeCounters_[valueToCompare_int];

      // Only the first element of each counter is compared;
      // in all other cases tne entire vectors are compared.
      if (groupname == "Counters") {
        if (!oops::anyVectorEmpty(values_OPS, values_thiscode))
          compareOutput(valueToCompare_int, values_OPS[0], values_thiscode[0],
                        offset, tol, nMismatches_);
      } else {
        compareOutput(valueToCompare_int, values_OPS, values_thiscode,
                      offset, tol, nMismatches_);
      }

      // Increment cumulative counters. NumAnyErrors is not included.
      if (groupname == "Counters" && varname != "NumAnyErrors")
        cumulativeCounters_[valueToCompare_int] += values_thiscode[0];
    }

    // Compare float values obtained in this code and OPS
    for (const auto& valueToCompare_float : valuesToCompare_float_) {
      const std::vector <float> &values_thiscode =
        profileDataHandler_.get<float>(valueToCompare_float);
      const std::vector <float> &values_OPS =
        profileDataHandler_.get<float>("OPS_" + valueToCompare_float);
      compareOutput(valueToCompare_float, values_OPS, values_thiscode,
                    0, tol, nMismatches_);
    }

    oops::Log::debug() << " ... all comparisons done ("
                       << nMismatches_ << " mismatches)" << std::endl;
  }
}  // namespace ufo


