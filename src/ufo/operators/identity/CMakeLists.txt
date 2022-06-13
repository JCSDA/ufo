# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( identity_files
    ObsIdentity.h
    ObsIdentity.cc
    ObsIdentityParameters.h
    ObsIdentityTLAD.h
    ObsIdentityTLAD.cc
)

PREPEND( _p_identity_files       "operators/identity"       ${identity_files} )

set ( identity_src_files
      ${_p_identity_files}
      PARENT_SCOPE
)

