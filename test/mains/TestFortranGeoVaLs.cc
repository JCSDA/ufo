/*
 * (C) Crown copyright 2023, UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "../ufo/FortranGeoVaLs.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::test::FortranGeoVaLs tests;
  return run.execute(tests);
}