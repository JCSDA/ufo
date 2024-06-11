/* * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/BayesianBackgroundCheck.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/QCflags.h"

#include "ufo/utils/ProbabilityOfGrossError.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// BayesianBackgroundCheck: check observation closeness to background, accounting
///   for probability of gross error, i.e. that observation is bad.

BayesianBackgroundCheck::BayesianBackgroundCheck(
        ioda::ObsSpace & obsdb,
        const Parameters_ & parameters,
        std::shared_ptr<ioda::ObsDataVector<int> > flags,
        std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr, VariableNameMap(parameters.AliasFile.value())),
    parameters_(parameters)

{
  oops::Log::trace() << "BayesianBackgroundCheck constructor" << std::endl;
  allvars_ += Variables(filtervars_, "HofX");
  allvars_ += Variables(filtervars_, "QCFlags");
  for (size_t i = 0; i < filtervars_.size(); ++i) {
      allvars_ += backgrErrVariable(filtervars_[i]);
  }
}

// -----------------------------------------------------------------------------

BayesianBackgroundCheck::~BayesianBackgroundCheck() {
  oops::Log::trace() << "BayesianBackgroundCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/// Return the name of the variable containing the background error estimate of the
/// specified filter variable.

Variable BayesianBackgroundCheck::backgrErrVariable(const Variable &filterVariable) const {
  return Variable(parameters_.BkgErrGroup.value() + "/" +
                  nameMap_.convertName(filterVariable.variable()).name() +
                  parameters_.BkgErrSuffix.value());
}

// -----------------------------------------------------------------------------
/// Reduce a vector to only the elements for which j_reduced=true.

template <class T>
std::vector<T> BayesianBackgroundCheck::reduceVector(const std::vector<T> &vector_full,
                                                     const std::vector<size_t> &j_reduced) const {
  size_t reduced_size = j_reduced.size();
  std::vector<T> vector_reduced(reduced_size);
  for (size_t i_reduced=0; i_reduced < reduced_size; ++i_reduced) {
    vector_reduced[i_reduced] = vector_full[j_reduced[i_reduced]];
  }
  return vector_reduced;
}

// -----------------------------------------------------------------------------
/// Copy a reduced vector back into the correct indices of the full vector.

template <class T>
void BayesianBackgroundCheck::unreduceVector(const std::vector<T> &vector_reduced,
                                             std::vector<T> &vector_full,
                                             const std::vector<size_t> &j_reduced) const {
  size_t reduced_size = vector_reduced.size();
  for (size_t i_reduced=0; i_reduced < reduced_size; ++i_reduced) {
    vector_full[j_reduced[i_reduced]] = vector_reduced[i_reduced];
  }
}

// -----------------------------------------------------------------------------
/// Apply the Bayesian background check filter.

void BayesianBackgroundCheck::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "BayesianBackgroundCheck postFilter" << std::endl;
  const oops::ObsVariables observed = obsdb_.obsvariables();

  oops::Log::debug() << "BayesianBackgroundCheck obserr: " << *obserr_ << std::endl;

  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsObsVariables(), "ObsValue");

  Variables varhofx(filtervars_, "HofX");
  Variables varflags(filtervars_, "QCFlags");
  const float missingValueFloat = util::missingValue<float>();

  // Probability density of bad observations, PdBad:
  const std::vector<float> PdBad(obsdb_.nlocs(), parameters_.PdBad.value());

  bool previousVariableWasFirstComponentOfTwo = false;
  // Loop through all filter variables. .yaml will say if it's 2-component.
  // If so, it skips the first component and then processes both together when
  //  (previousVariableWasFirstComponentOfTwo=True).
  // Otherwise it processes the scalar variable.
  for (size_t filterVarIndex = 0; filterVarIndex < filtervars.size(); ++filterVarIndex) {
    if (filtervars[filterVarIndex].options().getBool("first_component_of_two", false)) {
      previousVariableWasFirstComponentOfTwo = true;
    } else {
      std::string varname1, varname2;
      size_t iv1, iv2;
      // H(x):
      std::vector<float> hofx1, hofx2;
      // H(x) error:
      std::vector<float> hofxerr(obsdb_.nlocs());
      if (parameters_.BkgErr.value() != boost::none) {
        // set it to a constant term
        oops::Log::debug() << "Setting a constant bg error term" << std::endl;
        std::fill(hofxerr.begin(), hofxerr.end(), *parameters_.BkgErr.value());
      } else {
        oops::Log::debug() << "Using the real background errors" << std::endl;
        data_.get(backgrErrVariable(filtervars[filterVarIndex]), hofxerr);
      }
      // PGE:
      std::vector<float> PGE1(obsdb_.nlocs());
      // Total (combined) probability distribution
      std::vector<float> TotalPd(obsdb_.nlocs(), missingValueFloat);
      // QC flags:
      std::vector<int> qcflags1(obsdb_.nlocs());
      std::vector<float> firstComponentObVal, secondComponentObVal;
      // make note of conditions on apply and flags_:
      std::vector<bool> applycondition(obsdb_.nlocs(), false);
      // index mapping between full and reduced vectors:
      std::vector<size_t> j_reduced;

      if (previousVariableWasFirstComponentOfTwo) {
        varname1 = filtervars.variable(filterVarIndex-1).variable();
        varname2 = filtervars.variable(filterVarIndex).variable();
        iv1 = observed.find(varname1);
        iv2 = observed.find(varname2);
        // H(x):
        data_.get(varhofx.variable(filterVarIndex-1), hofx1);
        data_.get(varhofx.variable(filterVarIndex), hofx2);
        // observation values:
        firstComponentObVal = obs[filterVarIndex-1];
        secondComponentObVal = obs[filterVarIndex];
        // QC flags (all zeros, or read from file if possible):
        if (obsdb_.has("QCFlags", varname1)) {
          data_.get(varflags.variable(filterVarIndex-1), qcflags1);
          oops::Log::debug() << "Got qcflags1 from file: " << qcflags1 << std::endl;
        }
        // PGE:
        obsdb_.get_db("GrossErrorProbability", varname1, PGE1);
        // Get total (combined) probability distribution if previously exists
        if (obsdb_.has("GrossErrorProbabilityTotal", varname1)) {
            obsdb_.get_db("GrossErrorProbabilityTotal", varname1, TotalPd);
          }
        for (size_t jobs=0; jobs < obsdb_.nlocs(); ++jobs) {
          if (apply[jobs]) {
            applycondition[jobs] = true;
            j_reduced.push_back(jobs);
          }
        }
      } else {
        varname1 = filtervars.variable(filterVarIndex).variable();
        iv1 = observed.find(varname1);
        // H(x):
        data_.get(varhofx.variable(filterVarIndex), hofx1);
        // observation values:
        firstComponentObVal = obs[filterVarIndex];
        // QC flags (all zeros, or read from file if possible):
        if (obsdb_.has("QCFlags", varname1)) {
          data_.get(varflags.variable(filterVarIndex), qcflags1);
          oops::Log::debug() << "Got qcflags1 from file: " << qcflags1 << std::endl;
        }
        // PGE:
        obsdb_.get_db("GrossErrorProbability", varname1, PGE1);
        // Get total (combined) probability distribution if previously exists
        if (obsdb_.has("GrossErrorProbabilityTotal", varname1)) {
          obsdb_.get_db("GrossErrorProbabilityTotal", varname1, TotalPd);
        }
        for (size_t jobs=0; jobs < obsdb_.nlocs(); ++jobs) {
          if (apply[jobs]) {
            applycondition[jobs] = true;
            j_reduced.push_back(jobs);
          }
        }
      }

      // create reduced vectors, copied from full ones, fulfilling applycondition:
      std::vector<float> firstComponentObVal_reduced = reduceVector(firstComponentObVal,
                                                                    j_reduced);
      std::vector<float> ObsErr_reduced = reduceVector((*obserr_)[varname1], j_reduced);
      std::vector<float> hofx1_reduced = reduceVector(hofx1, j_reduced);
      std::vector<float> hofxerr_reduced = reduceVector(hofxerr, j_reduced);
      std::vector<float> PdBad_reduced = reduceVector(PdBad, j_reduced);
      std::vector<int> qcflags1_reduced = reduceVector(qcflags1, j_reduced);
      std::vector<float> PGE1_reduced = reduceVector(PGE1, j_reduced);
      std::vector<float> secondComponentObVal_reduced;
      std::vector<float> hofx2_reduced;
      std::vector<float> TotalPd_reduced;
      if (previousVariableWasFirstComponentOfTwo) {
        secondComponentObVal_reduced = reduceVector(secondComponentObVal, j_reduced);
        hofx2_reduced = reduceVector(hofx2, j_reduced);
      }
      if (parameters_.SaveTotalPd) {
          TotalPd_reduced = reduceVector(TotalPd, j_reduced);
      }

      ufo::BayesianPGEUpdate(parameters_.PGEParameters,
                             firstComponentObVal_reduced,
                             ObsErr_reduced,
                             hofx1_reduced,
                             hofxerr_reduced,
                             PdBad_reduced,
                             parameters_.PerformSDiffCheck.value(),
                             qcflags1_reduced,
                             PGE1_reduced,
                             parameters_.ErrVarMax,
                             previousVariableWasFirstComponentOfTwo?
                             &secondComponentObVal_reduced : nullptr,
                             previousVariableWasFirstComponentOfTwo?
                             &hofx2_reduced : nullptr,
                             parameters_.SaveTotalPd?
                             &TotalPd_reduced : nullptr);

      // write PGE-updated values from reduced vectors back into full ones:
      unreduceVector(qcflags1_reduced, qcflags1, j_reduced);
      unreduceVector(PGE1_reduced, PGE1, j_reduced);
      unreduceVector(TotalPd_reduced, TotalPd, j_reduced);

      // Save PGE to obsdb
      obsdb_.put_db("GrossErrorProbability", varname1, PGE1);              // PGE
      if (parameters_.SaveTotalPd) {
          obsdb_.put_db("GrossErrorProbabilityTotal", varname1, TotalPd);
      }

      // Save QC flags to obsdb
      obsdb_.put_db("QCFlags", varname1, qcflags1);  // Met Office QC flags, not flagged or *flags_

      if (previousVariableWasFirstComponentOfTwo) {
        // Save PGE to obsdb
        std::vector<float> &PGE2 = PGE1;  // in old OPS, PGE same for both components of 2-vector
        obsdb_.put_db("GrossErrorProbability", varname2, PGE2);
        if (parameters_.SaveTotalPd) {
            obsdb_.put_db("GrossErrorProbabilityTotal", varname2, TotalPd);
        }
        // Save QC flags to obsdb
        std::vector<int> &qcflags2 = qcflags1;  // in old OPS, flags same for both components
        obsdb_.put_db("QCFlags", varname2, qcflags2);
        // Set flagged, for 2nd component:
        for (size_t jobs=0; jobs < obsdb_.nlocs(); ++jobs) {
          if (qcflags1[jobs] & ufo::MetOfficeQCFlags::Elem::BackRejectFlag) {
            flagged[filterVarIndex-1][jobs] = true;
          }
          if (qcflags2[jobs] & ufo::MetOfficeQCFlags::Elem::BackRejectFlag) {
            flagged[filterVarIndex][jobs] = true;
          }
          oops::Log::debug() << "flagged(1)[" << jobs << "]: "
                             << flagged[filterVarIndex-1][jobs] << std::endl;
          oops::Log::debug() << "flagged(2)[" << jobs << "]: "
                             << flagged[filterVarIndex][jobs] << std::endl;
        }
      } else {
        // Set flagged, for scalar:
        for (size_t jobs=0; jobs < obsdb_.nlocs(); ++jobs) {
          if (qcflags1[jobs] & ufo::MetOfficeQCFlags::Elem::BackRejectFlag ||
              qcflags1[jobs] & ufo::MetOfficeQCFlags::Elem::FinalRejectFlag) {
            flagged[filterVarIndex][jobs] = true;
          }
          oops::Log::debug() << "flagged[" << jobs << "]: "
                             << flagged[filterVarIndex][jobs] << std::endl;
        }
      }
      previousVariableWasFirstComponentOfTwo = false;
    }
  }
}

// -----------------------------------------------------------------------------

void BayesianBackgroundCheck::print(std::ostream & os) const {
  os << "BayesianBackgroundCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
