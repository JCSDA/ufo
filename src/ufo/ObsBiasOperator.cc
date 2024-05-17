/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBiasOperator.h"

#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ScopedDefaultGeoVaLFormatChange.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasOperator::ObsBiasOperator(ioda::ObsSpace & odb)
  : odb_(odb) {
  oops::Log::trace() << "ObsBiasOperator::create done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasOperator::computeObsBias(const GeoVaLs & geovals, ioda::ObsVector & ybias,
                                     const ObsBias & biascoeffs, ObsDiagnostics & ydiags,
                                     const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsBiasOperator::computeObsBias starting" << std::endl;

  const double missing = util::missingValue<double>();
  ScopedDefaultGeoVaLFormatChange change(geovals, GeoVaLFormat::REDUCED);
  const Predictors & predictors = biascoeffs.predictors();
  const std::size_t npreds = predictors.size();
  std::vector<ioda::ObsVector> predData(npreds, ioda::ObsVector(odb_));
  for (std::size_t p = 0; p < npreds; ++p) {
    predictors[p]->compute(odb_, geovals, ydiags, biascoeffs, predData[p]);
  }

  const oops::ObsVariables &simVars = biascoeffs.simVars();
  // At present we can label predictors with either the channel number or the variable name, but not
  // both. So if there are multiple channels, make sure there's only one (multi-channel) variable.
  ASSERT(simVars.channels().empty() ||
         simVars.variables().size() == simVars.channels().size());

  const std::size_t nlocs  = odb_.nlocs();
  const std::size_t nvars  = simVars.variables().size();

  const std::vector<int> & varIndexNoBC = biascoeffs.varIndexNoBC();
  for (const int jvar : varIndexNoBC) {
    for (std::size_t jp = 0; jp < npreds; ++jp) {
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        predData[jp][jl * nvars + jvar] = 0.0;
      }
    }
  }

  ybias.zero();

  /* ybias memory layout (nlocs X nvars)
   *     ch1    ch2    ch3     ch4
   * Loc --------------------------
   *  0 | 0      1      2       3
   *  1 | 4      5      6       7
   *  2 | 8      9     10      11 
   *  3 |12     13     14      15
   *  4 |16     17     18      19
   * ...|
   */

  std::vector<double> biasTerm(nlocs);
  //  For each channel: ( nlocs X 1 ) =  ( nlocs X npreds ) * (  npreds X 1 )
  for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
    if (std::find(varIndexNoBC.begin(), varIndexNoBC.end(), jvar) == varIndexNoBC.end()) {
      std::string predictorSuffix;
      if (simVars.channels().empty())
        predictorSuffix = simVars[jvar];
      else
        predictorSuffix = std::to_string(simVars.channels()[jvar]);

      for (std::size_t jp = 0; jp < npreds; ++jp) {
        // axpy
        for (std::size_t jl = 0; jl < nlocs; ++jl) {
          if (predData[jp][jl * nvars + jvar] != missing) {
            const std::size_t jrec = biascoeffs.nrecs() == 1 ? 0 : odb_.recnum()[jl];
            const double beta = biascoeffs(jrec, jvar, jp);
            biasTerm[jl] = predData[jp][jl * nvars + jvar] * beta;
            ybias[jl * nvars + jvar] += biasTerm[jl];
          }
        }
        // Save ObsBiasOperatorTerms (bias_coeff * predictor) for QC
        const std::string varname = predictors[jp]->name() + "_" + predictorSuffix;
        if (ydiags.has(varname)) {
          ydiags.allocate(1, oops::ObsVariables({varname}));
          ydiags.save(biasTerm, varname, 0);
        } else {
          oops::Log::error() << varname << " is not reserved in ydiags !" << std::endl;
          ABORT("ObsBiasOperatorTerm variable is not reserved in ydiags");
        }
      }
    }
  }

  for (std::size_t p = 0; p < npreds; ++p) {
    predData[p].save(predictors[p]->name() + "Predictor");
  }

  oops::Log::trace() << "ObsBiasOperator::computeObsBiasOperator done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasOperator::print(std::ostream & os) const {
  os << "ObsBiasOperator: linear combination of the predictors." << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
