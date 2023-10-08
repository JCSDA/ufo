/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ProbabilityGrossErrorWholeReport.h"

#include <cmath>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"

namespace ufo {

// -----------------------------------------------------------------------------
ProbabilityGrossErrorWholeReport::ProbabilityGrossErrorWholeReport(ioda::ObsSpace & obsdb,
                                 const Parameters_ & parameters,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr),
    parameters_(parameters)
{
  oops::Log::trace() << "ProbabilityGrossErrorWholeReport contructor starting" << std::endl;
}

// -----------------------------------------------------------------------------

ProbabilityGrossErrorWholeReport::~ProbabilityGrossErrorWholeReport() {
  oops::Log::trace() << "ProbabilityGrossErrorWholeReport destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ProbabilityGrossErrorWholeReport::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ProbabilityGrossErrorWholeReport preProcess" << std::endl;
  // Missing value indicator
  const float missingValueFloat = util::missingValue<float>();
  // Dimensions
  const size_t nlocs = obsdb_.nlocs();
  const size_t nvars = filtervars.nvars();
  // Missing data indicator for stored PGEs.
  const float PGEMDI = 1.111f;
  // PGE multiplication factor used to store PGE values for later use.
  const float PGEMult = 1000.0;
  // PGE rejection limit.
  const float PGECrit = parameters_.PGEParameters.PGE_PGECrit.value();

  // Initial Probability of Gross error in whole report
  std::vector<float> PBadRep(nlocs);
  // Overall Probability density for whole report
  std::vector<float> PdReport(nlocs, 1.0f);
  // Probability density of bad observations for each variable
  std::vector<float> PdBad(nvars);
  // Probability density of bad synop observations for each variable
  std::vector<float> PdBad_synop(nvars);
  // Probability density of bad bogus observations for each variable
  std::vector<float> PdBad_bogus(nvars);
  // Probabilty density of bad observations used for each observation
  std::vector<std::vector<float>> PdBad_used(nvars, std::vector<float>(nlocs));
  // Is this variable used in the whole report error calculation?
  std::vector<bool> notUsed(nvars);
  // Vector holding filter variable names
  std::vector<std::string> varname(nvars);
  // Vector holding probability of gross error
  std::vector<std::vector<float>> varPGE(nvars, std::vector<float>(nlocs));
  // Vector holding initial probability of gross error
  std::vector<std::vector<float>> varPGEInitial(nvars, std::vector<float>(nlocs));
  // Vector holding combined probability densities
  std::vector<std::vector<float>> varPGETotal(nvars, std::vector<float>(nlocs));
  // Vector holding OPS style QCflags for each variable
  std::vector<std::vector<int>> varQCflags(nvars, std::vector<int>(nlocs));

  for (size_t ivar = 0; ivar < filtervars.nvars(); ++ivar) {
      varname[ivar] = filtervars.variable(ivar).variable();
    // Get Gross Error Probability values and Met Office QCFlags from ObsSpace
    if (obsdb_.has("GrossErrorProbability", varname[ivar]) &&
        obsdb_.has("GrossErrorProbabilityInitial", varname[ivar]) &&
        obsdb_.has("GrossErrorProbabilityTotal", varname[ivar]) &&
        obsdb_.has("QCFlags", varname[ivar])) {
    obsdb_.get_db("GrossErrorProbability", varname[ivar], varPGE[ivar]);
    obsdb_.get_db("GrossErrorProbabilityInitial", varname[ivar], varPGEInitial[ivar]);
    obsdb_.get_db("GrossErrorProbabilityTotal", varname[ivar], varPGETotal[ivar]);
    obsdb_.get_db("QCFlags", varname[ivar], varQCflags[ivar]);
    } else {
      std::stringstream errormessage;
      errormessage << "GrossErrorProbability/" + varname[ivar] + ", "
                   << "GrossErrorProbabilityInitial/" + varname[ivar] + ", "
                   << "GrossErrorProbabilityTotal/" + varname[ivar] + ", and "
                   << "QCFlags/" + varname[ivar] + " must all be present"
                   << std::endl;
      throw eckit::BadValue(errormessage.str(), Here());
    }
    // Get variable options
    PdBad[ivar] = filtervars[ivar].options().getFloat("probability_density_bad", 0.1f);
    if (filtervars[ivar].options().has("synop_probability_density_bad")) {
        PdBad_synop[ivar] = filtervars[ivar].options().getFloat("synop_probability_density_bad");
    } else {
        PdBad_synop[ivar] = PdBad[ivar];
    }
    if (filtervars[ivar].options().has("bogus_probability_density_bad")) {
        PdBad_bogus[ivar] = filtervars[ivar].options().getFloat("bogus_probability_density_bad");
    } else {
        PdBad_bogus[ivar] = PdBad[ivar];
    }
    //
    notUsed[ivar] = filtervars[ivar].options().getBool("not_used_in_whole_report", false);
  }
  // Get input probability of gross error affecting whole report
  std::vector<float> ReportPGE(nlocs);
  if (obsdb_.has("MetaData", "grossErrorProbabilityReport")) {
    obsdb_.get_db("MetaData", "grossErrorProbabilityReport", ReportPGE);
  } else {
    throw eckit::BadValue("MetaData/grossErrorProbabilityReport not present", Here());
  }
  // Get ObsType
  std::vector<int> ObsType(nlocs);
  if (obsdb_.has("MetaData", "observationTypeNum")) {
    obsdb_.get_db("MetaData", "observationTypeNum", ObsType);
  } else {
    throw eckit::BadValue("MetaData/observationTypeNum not present", Here());
  }

  // Calculate probability that whole report is affected by gross error
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    if (apply[jobs] && (ReportPGE[jobs] > 0)) {
      // Probability density given Gross error in whole report
      float PdBadRep = 1.0f;
      for (size_t ivar = 0; ivar < filtervars.nvars(); ++ivar) {
        if (!notUsed[ivar]) {
          PdReport[jobs] *= varPGETotal[ivar][jobs];
          if (ObsType[jobs] == MetOfficeObsIDs::Bogus::Bogus) {
            PdBad_used[ivar][jobs] = PdBad_bogus[ivar];
          } else if (ObsType[jobs] == MetOfficeObsIDs::Surface::SynopManual ||
                     ObsType[jobs] == MetOfficeObsIDs::Surface::SynopAuto ||
                     ObsType[jobs] == MetOfficeObsIDs::Surface::MetarManual ||
                     ObsType[jobs] == MetOfficeObsIDs::Surface::MetarAuto ||
                     ObsType[jobs] == MetOfficeObsIDs::Surface::SynopMob ||
                     ObsType[jobs] == MetOfficeObsIDs::Surface::SynopBufr ||
                     ObsType[jobs] == MetOfficeObsIDs::Surface::WOW) {
            PdBad_used[ivar][jobs] = PdBad_synop[ivar];
          } else {
            PdBad_used[ivar][jobs] = PdBad[ivar];
          }
          PdBadRep *= PdBad_used[ivar][jobs];
        }
      }
      PBadRep[jobs] = ReportPGE[jobs];
      PdBadRep *= PBadRep[jobs];
      PdReport[jobs] = PdBadRep + ((1.0f - PBadRep[jobs]) * PdReport[jobs]);
      if (PdReport[jobs] > 0.0f) {
        ReportPGE[jobs] = PdBadRep / PdReport[jobs];
      }
    }
  }

  // Calculate updated probability that individual observation OR
  // whole report is affected by gross error
  bool secondComponentOfTwo = false;
  for (size_t ivar = 0; ivar < filtervars.nvars(); ++ivar) {
    secondComponentOfTwo = filtervars[ivar].options().getBool("second_component_of_two", false);
    if (filtervars[ivar].options().getBool("no_pge_update", false))
      continue;
    if (secondComponentOfTwo) {
      varPGE[ivar] = varPGE[ivar - 1];
      varQCflags[ivar] = varQCflags[ivar - 1];
      flagged[ivar] = flagged[ivar - 1];
    } else {
      for (size_t jobs = 0; jobs < nlocs; ++jobs) {
        if (apply[jobs] && (ReportPGE[jobs] > 0) && (PdReport[jobs] > 0)) {
          if (std::abs(varPGE[ivar][jobs] - PGEMDI) > 0.001f) {
            float PdProduct = 1.0f;
            for (size_t jvar = 0; jvar < filtervars.nvars(); ++jvar) {
              if ((jvar != ivar) && !notUsed[jvar]) {
                PdProduct *= varPGETotal[jvar][jobs];
              }
            }
            // PGE in element or whole report
            float PGE1;
            if (notUsed[ivar]) {
             PGE1 = ReportPGE[jobs] + varPGE[ivar][jobs] -
                     ReportPGE[jobs] * varPGE[ivar][jobs];
            } else {
             PGE1 = ReportPGE[jobs] + (PdBad_used[ivar][jobs] * PdProduct)
                   * varPGEInitial[ivar][jobs] *
                   (1.0f - PBadRep[jobs])/PdReport[jobs];
                }
            if (PGE1 >= PGECrit) {
              varQCflags[ivar][jobs] |= ufo::MetOfficeQCFlags::Elem::BackRejectFlag;
              flagged[ivar][jobs] = true;
            }
            varPGE[ivar][jobs] = PGE1;
          }
        }
      }
    }
    // Save updated gross error probabilities and QCFlags to ObsSpace
    obsdb_.put_db("GrossErrorProbability", varname[ivar], varPGE[ivar]);
    obsdb_.put_db("QCFlags", varname[ivar], varQCflags[ivar]);
  }
  obsdb_.put_db("MetaData", "grossErrorProbabilityReport", ReportPGE);
}

// -----------------------------------------------------------------------------

void ProbabilityGrossErrorWholeReport::print(std::ostream & os) const {
  os << "ProbabilityGrossErrorWholeReport: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
