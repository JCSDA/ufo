/*
 * (C) Copyright 2022-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/columnretrieval/ObsColumnRetrievalParameters.h"

#include "eckit/exception/Exceptions.h"

namespace ufo {

void ObsColumnRetrievalParameters::deserialize(util::CompositePath &path,
                                               const eckit::Configuration &config) {
  ObsOperatorParametersBase::deserialize(path, config);

  if (tracerVariables.value().size() != 1) {
    throw eckit::UserError(path.path() + ": ColumnRetrieval simulates a single variable, so "
                           "'tracer variables' must contain exactly one entry", Here());
  }

  if (totalNoVertice.value() && (nlayers_retrieval.value() != 1 || isAveragingKernel.value())) {
    throw eckit::UserError(path.path() + ": totalNoVertice requires nlayers_retrieval: 1 "
                           "and isAveragingKernel: false", Here());
  }
}

}  // namespace ufo
