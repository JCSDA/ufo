!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 3 (known as Metop-B (DMI)), taken from processing centre 94
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 3
  origc = 94
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 4 (known as Metop-A (DMI)), taken from processing centre 94
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 4
  origc = 94
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 5 (known as Metop-C (DMI)), taken from processing centre 94
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 5
  origc = 94
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 740 (known as COSMIC-1 FM1 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 740
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 741 (known as COSMIC-1 FM2 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 741
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 742 (known as COSMIC-1 FM3 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 742
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 743 (known as COSMIC-1 FM4 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 743
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 744 (known as COSMIC-1 FM5 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 744
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 745 (known as COSMIC-1 FM6 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 745
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 41 (CHAMP (GFZ)), taken from processing centre 78
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 78
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 41 (CHAMP (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 41
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 722 (known as GRACE-A (GFZ)), taken from processing centre 78
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 722
  origc = 78
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 522 (known as FY-3C (CMA)), taken from processing centre 38
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 522
  origc = 38
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 523 (known as FY-3D (CMA)), taken from processing centre 38
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 523
  origc = 38
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 524 (known as FY-3E (CMA)), taken from processing centre 38
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
  satid = 524
  origc = 38
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 524
  origc = 38
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 524
  origc = 38
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 524
  origc = 38
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 524
  origc = 38
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 524
  origc = 38
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 825 (known as KOMPSAT-5 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 825
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 44 (known as PAZ (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 44
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 750 (known as COSMIC-2 E1 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 750
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 751 (known as COSMIC-2 E2 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 751
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 751 (known as COSMIC-2 E3 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 752
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 753 (known as COSMIC-2 E4 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 753
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 754 (known as COSMIC-2 E5 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 754
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 755 (known as COSMIC-2 E6 (UCAR)), taken from processing centre 60
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
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 755
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 803 (known as GRACE-C (GFZ)), taken from processing centre 78
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
  satid = 803
  origc = 78
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 803
  origc = 78
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 803
  origc = 78
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 803
  origc = 78
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 803
  origc = 78
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 803
  origc = 78
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 804 (known as GRACE-D (GFZ)), taken from processing centre 78
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
  satid = 804
  origc = 78
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 804
  origc = 78
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 804
  origc = 78
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 804
  origc = 78
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 804
  origc = 78
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 804
  origc = 78
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 66 (known as Sentinel-6A (JPL)), taken from processing centre 173
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
  satid = 66
  origc = 173
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 66
  origc = 173
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 66
  origc = 173
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 66
  origc = 173
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 66
  origc = 173
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 66
  origc = 173
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 269 (known as Spire (EUMETSAT)), taken from processing centre 254
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
  origc = 254
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 254
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 254
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 254
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 254
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 254
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 269 (known as Spire), taken from processing centre 94 (ROM SAF)
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
  origc = 94
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 94
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 94
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 94
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 94
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 94
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 269 (known as Spire), taken from processing centre 60 (UCAR)
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
  origc = 60
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 269
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 267 (known as PlanetIQ), taken from processing centre 60 (UCAR)
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
  satid = 267
  origc = 60
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 267
  origc = 60
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 267
  origc = 60
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 267
  origc = 60
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 267
  origc = 60
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 267
  origc = 60
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
!
! Observation error namelist for GPSRO refractivity observations
! This entry is for satellite 267 (known as PlanetIQ), taken from processing centre 254 (EUMETSAT)
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
  satid = 267
  origc = 254
  min_error = 0.200E-08
  latitude = 75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 267
  origc = 254
  min_error = 0.200E-08
  latitude = 45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 267
  origc = 254
  min_error = 0.200E-08
  latitude = 15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 267
  origc = 254
  min_error = 0.200E-08
  latitude = -15
  clen = 0.670E-03
  heights = 0.000E+00, 0.040E+05, 0.100E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.500E-01, 0.250E-01, 0.050E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 267
  origc = 254
  min_error = 0.200E-08
  latitude = -45
  clen = 0.670E-03
  heights = 0.000E+00, 0.070E+05, 0.250E+05, 0.400E+05, 0.500E+05  
  obs_errors = 0.300E-01, 0.070E-01, 0.070E-01, 0.400E-01, 0.500E-01
/
&GPSRO_ob_error
  satid = 267
  origc = 254
  min_error = 0.200E-08
  latitude = -75
  clen = 0.670E-03
  heights = 0.000E+00, 0.060E+05, 0.200E+05, 0.250E+05, 0.400E+05, 0.500E+05
  obs_errors = 0.100E-01, 0.060E-01, 0.060E-01, 0.100E-01, 0.400E-01, 0.500E-01
/
