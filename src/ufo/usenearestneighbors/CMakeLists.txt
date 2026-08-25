# (C) Crown Copyright 2025, Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set(usenearestneighbors_files
    UseNearestNeighborsAlgorithmBase.h
    UseNearestNeighborsAlgorithmBase.cc
    UseNearestNeighborsAlgorithmGatherAndMatchTimestamp.h
    UseNearestNeighborsAlgorithmGatherAndMatchTimestamp.cc
    UseNearestNeighborsAlgorithmReferencePointVariablesMean.h
    UseNearestNeighborsAlgorithmReferencePointVariablesMean.cc
    UseNearestNeighborsAlgorithmLocalPlaneFit.h
    UseNearestNeighborsAlgorithmLocalPlaneFit.cc
)

# Prefix file paths with the current folder name so they can be aggregated upstream.
PREPEND(_p_usenearestneighbors_files "usenearestneighbors" ${usenearestneighbors_files})

# Export list to parent scope; the top-level ufo CMakeLists adds this to ufo_src_files.
set(usenearestneighbors_src_files
    ${_p_usenearestneighbors_files}
    PARENT_SCOPE
)
