!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 3 (known as Metop-B), taken from processing centre 94 (DMI)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 4 (known as Metop-A), taken from processing centre 94 (DMI)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 5 (known as Metop-C), taken from processing centre 94 (DMI)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 3 (known as Metop-B), taken from processing centre 254 (EUMETSAT)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 254
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 4 (known as Metop-A), taken from processing centre 254 (EUMETSAT)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 254
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 5 (known as Metop-C), taken from processing centre 254 (EUMETSAT)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 254
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 41 (known as CHAMP), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 41 (known as CHAMP), taken from processing centre 78 (GFZ)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 42 (known as TerraSAR-X), taken from processing centre 78 (GFZ)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 42
  origc = 78
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 43 (known as TanDEM-X), taken from processing centre 78 (GFZ)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 43
  origc = 78
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 44 (known as PAZ), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 722 (known as Grace-A), taken from processing centre 78 (GFZ)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 740 (known as COSMIC-1), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 741 (known as COSMIC-2), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 742 (known as COSMIC-3), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 743 (known as COSMIC-4), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 744 (known as COSMIC-5), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 745 (known as COSMIC-6), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 786 (known as C/NOFS), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 786
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 786
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 786
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 786
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 786
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 786
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 522 (known as FY-3C), taken from processing centre 38
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 6.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.200E-01, 0.200E-01
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 6.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.200E-01, 0.200E-01
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 6.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.200E-01, 0.200E-01
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 6.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.200E-01, 0.200E-01
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 6.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.200E-01, 0.200E-01
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 6.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.200E-01, 0.200E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 523 (known as FY-3D), taken from processing centre 38
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 6.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.200E-01, 0.200E-01
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 6.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.200E-01, 0.200E-01
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 6.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.200E-01, 0.200E-01
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 6.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.200E-01, 0.200E-01
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 6.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.200E-01, 0.200E-01
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 6.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.200E-01, 0.200E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 825 (known as KOMPSAT-5), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 750 (known as COSMIC-2 E1), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 751 (known as COSMIC-2 E2), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 752 (known as COSMIC-2 E3), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 753 (known as COSMIC-2 E4), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 754 (known as COSMIC-2 E5), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 755 (known as COSMIC-2 E6), taken from processing centre 60 (UCAR)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 269 (known as Spire), taken from processing centre 178 (Spire)
!
! Each entry is for a given range of latitudes, and the nearest matrix to each
! observation will be chosen.
! The observation errors are given for a set of heights.
! If an observation is present on a height outside this range, then the value
! will be taken to be the value at the relevant end.
! If an observation is at a height between the levels indicated, then the code
! will linearly interpolate between those heights.
! The entries are given as relative errors (percentage error of the observation
! divided by the background refractivity).
!
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3.000E-06
  latitude = 75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3.000E-06
  latitude = 45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3.000E-06
  latitude = 15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3.000E-06
  latitude = -15
  clen = 1.000E+10
  heights = 0.000E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3.000E-06
  latitude = -45
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.110E+00, 0.150E-01, 0.150E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 178
  min_error = 3.000E-06
  latitude = -75
  clen = 1.000E+10
  heights = 0.000E+00, 0.025E+00, 0.100E+05, 0.700E+05
  obs_errors = 0.200E+00, 0.800E-01, 0.150E-01, 0.150E-01
/
