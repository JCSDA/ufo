/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBiasPreconditioner.h"

#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasPreconditioner::ObsBiasPreconditioner(const std::vector<double> & precond) :
    precond_(precond) {}

void ObsBiasPreconditioner::multiply(const ObsBiasIncrement & dx1, ObsBiasIncrement & dx2) const {
  const Eigen::Map<const Eigen::ArrayXd> preconditioner(precond_.data(), precond_.size());
  dx2.data() = preconditioner.array()*dx1.data().array();
}



// -----------------------------------------------------------------------------

}  // namespace ufo
