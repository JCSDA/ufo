#!/bin/bash
#
# (C) Crown Copyright 2024, the Met Office. All rights reserved.
#
#
# shellcheck disable=SC2317
set -euo pipefail

finally() {
    trap '' ERR
    trap '' EXIT
    if [[ -d "${WORKD:-}" ]]; then
        cd /
        rm -fr "${WORKD}"
    fi
    if [[ -d ${BASE:-} ]]; then
        cd /
        rm -fr "$BASE" 2>/dev/null || true
    fi
}

# HERE is /var/tmp/<REPONAME>/pr-<#> (cf ../.github/workflow/ci.yml)
HERE="$(cd "$(dirname "$0")" && pwd)"
THIS="$(basename "$0")"
NPROC=${NPROC:-$(nproc)}
WORKD="$(mktemp -d "${THIS}-XXXXXX" -t)"
BASE="${HERE%/*}"
TESTDIR="${BASE##*/}"
LFRICJEDI_TEST_TIER=${LFRICJEDI_TEST_TIER:-0}

trap finally ERR EXIT
cd "${WORKD}"

# -- Activate spack env if using JCSDA Docker container
if [[ -f /opt/spack-environment/activate.sh ]]; then
    # shellcheck disable=SC1091
    source /opt/spack-environment/activate.sh
fi

# -- Enable OpenMPI over subscription -----------------------------------------
if command -v ompi_info &>/dev/null; then
    echo "Check support for MPI_THREAD_MULTIPLE"
    ompi_info | grep -i 'thread support'
    ompi_vn=$(ompi_info | awk '/Ident string:/ {print $3}')
    case $ompi_vn in
        4.*) export OMPI_MCA_rmaps_base_oversubscribe=1 ;;
        5.*) export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe ;;
    esac
fi

if command -v ninja &>/dev/null; then
    GENERATOR=Ninja;
else
    GENERATOR=Unix\ Makefiles;
fi

# -- Configure
cmake -B . -S "${HERE}" -G "${GENERATOR}" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_EXE_LINKER_FLAGS='-lmpi -lmpi_cxx -lmpi_mpifh -lgfortran -lquadmath -lstdc++'

# -- Build
cmake --build . -j "${NPROC}"

# -- Test
ctest --test-dir "${TESTDIR}" -E 'coding_norms' --output-on-failure
echo "-- Run Met Office model-interface tests"
ctest -R 'orca|unifiedmodel|lfric' -E 'coding_norms|xios_server_mode|_C192$' --output-on-failure

exit
