/*
 * (C) Crown Copyright, 2025 Met Office UK
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/Run.h"
#include "test/base/ObsErrorCovariance.h"

#include "../ufo/DiffusionDirac.h"

#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::test::DiffusionDirac tests;
  return run.execute(tests);
}
