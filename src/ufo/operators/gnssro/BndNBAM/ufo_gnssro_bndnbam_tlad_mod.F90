! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!> Fortran module of Gnssro NBAM (NCEP's Bending Angle Method)
!> tlad operator

module ufo_gnssro_bndnbam_tlad_mod
use fckit_configuration_module, only: fckit_configuration
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use obsspace_mod
use missing_values_mod
use lag_interp_mod
use ufo_basis_tlad_mod,only: ufo_basis_tlad
use gnssro_mod_conf
use gnssro_mod_constants
use gnssro_mod_transform, only: geop2geometric
use fckit_log_module,  only : fckit_log
use iso_c_binding, only: c_ptr, c_double
use ufo_constants_mod, only: zero, half, one, two, three, rd_over_g, rd_over_rv, rv_over_rd
use ufo_utils_mod, only: cmp_strings


implicit none
real(c_double)                             :: missing

!>  Fortran derived type for gnssro trajectory
type, extends(ufo_basis_tlad) :: ufo_gnssro_BndNBAM_tlad
  private
  integer                       :: nlev, nlev1, nlocs, iflip, nrecs
  real(kind_real), allocatable  :: jac_t(:,:), jac_prs(:,:), jac_q(:,:)
  integer, allocatable          :: nlocs_begin(:), nlocs_end(:)
  type(gnssro_conf)             :: roconf

  contains
    procedure :: setup      => ufo_gnssro_bndnbam_tlad_setup
    procedure :: delete     => ufo_gnssro_bndnbam_tlad_delete
    procedure :: settraj    => ufo_gnssro_bndnbam_tlad_settraj
    procedure :: simobs_tl  => ufo_gnssro_bndnbam_simobs_tl
    procedure :: simobs_ad  => ufo_gnssro_bndnbam_simobs_ad
end type ufo_gnssro_bndnbam_tlad

contains
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndnbam_tlad_setup(self, f_conf)
  implicit none
  class(ufo_gnssro_BndNBAM_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in)        :: f_conf

  call gnssro_conf_setup(self%roconf,f_conf)

end subroutine ufo_gnssro_bndnbam_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndnbam_tlad_settraj(self, geovals, obss)
  use gnssro_mod_transform
  use gnssro_mod_grids, only: get_coordinate_value
  use, intrinsic::  iso_c_binding
  implicit none
  class(ufo_gnssro_bndnbam_tlad), intent(inout) :: self
  type(ufo_geovals), intent(in)                 :: geovals
  type(c_ptr), value, intent(in)                :: obss
  character(len=*), parameter     :: myname  = "ufo_gnssro_bndnbam_tlad_settraj"
  character(max_string)           :: err_msg
  integer                         ::nlevAdd !num of additional levels on top of existing model levels
  integer                         :: ngrd
  type(ufo_geoval), pointer       :: t, q, gph, prs, zs
  integer                         :: iobs,k,j, klev, irec, icount
  integer                         :: nrecs
  integer                         :: nlev, nlev1, nlocs, nlevExt, nlevCheck
  real(kind_real)                 :: dw4(4), dw4_tl(4)
  real(kind_real)                 :: geomzi
  real(kind_real), allocatable    :: obsLat(:), obsImpP(:), obsLocR(:), obsGeoid(:)
  integer(c_size_t), allocatable  :: obsRecnum(:)
  real(kind_real)                 :: d_refXrad, gradRef
  real(kind_real)                 :: d_refXrad_tl
  real(kind_real), allocatable    :: grids(:)
  real(kind_real)                 :: sIndx
  integer                         :: indx
  real(kind_real)                 :: p_coef, t_coef, q_coef
  real(kind_real)                 :: fv, pw
  real(kind_real)                 :: dbetaxi, dbetan
  real(kind_real), allocatable    :: lagConst(:,:), lagConst_tl(:,:)
  real(kind_real), allocatable    :: gesT(:,:), gesQ(:,:), gesP(:,:), gesH(:,:), gesZs(:)

  real(kind_real), allocatable    :: radius(:), dzdh(:), refIndex(:)
  real(kind_real), allocatable    :: dhdp(:), dhdt(:)
  real(kind_real), allocatable    :: ref(:)
  real(kind_real), allocatable    :: refXrad(:)
  real(kind_real), allocatable    :: refXrad_s(:)
  real(kind_real), allocatable    :: refXrad_tl(:), ref_tl(:)
  real(kind_real), allocatable    :: dndp(:,:), dndt(:,:), dndq(:,:)
  real(kind_real), allocatable    :: dxidp(:,:), dxidt(:,:), dxidq(:,:)
  real(kind_real), allocatable    :: dbenddxi(:), dbenddn(:)
  integer,         allocatable    :: super_refraction_flag(:)
  integer,         allocatable    :: obsSRflag(:)
  integer                         :: hasSRflag
  integer                         :: ModelsigLevelcheck

! Make sure nothing already allocated
  call self%delete()

  missing = missing_value(missing)
  nlocs   = obsspace_get_nlocs(obss) ! number of observations
  nrecs   = obsspace_get_nrecs(obss)
  write(err_msg,*) myname, ': nlocs from gelvals and hofx, nrecs', nlocs, nrecs
  call fckit_log%debug(err_msg)

if (nlocs > 0 ) then
! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,  t)         ! air temperature
  call ufo_geovals_get_var(geovals, var_q,   q)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_sfc_geomz, zs)  ! surface geopotential height/surface altitude

  if (self%roconf%vertlayer .eq. "mass") then
    call ufo_geovals_get_var(geovals, var_prs,   prs)       ! pressure
    call ufo_geovals_get_var(geovals, var_z,     gph)       ! geopotential height
  else
    call ufo_geovals_get_var(geovals, var_prsi,  prs)       ! pressure
    call ufo_geovals_get_var(geovals, var_zi,    gph)
  end if

  nlev  = t%nval   ! number of model mass levels
  nlev1 = prs%nval
  self%nlev  = t%nval   ! number of model mass levels
  self%nlev1 = prs%nval ! number of model pressure/height levels
  self%nlocs = nlocs
  self%nrecs = nrecs

  nlevAdd = 13
  nlevCheck = int(nlev/2.0) !number of levels to check super refraction

  allocate(gesT(nlev,nlocs))
  allocate(gesQ(nlev,nlocs))
  allocate(gesP(nlev1,nlocs))
  allocate(gesH(nlev1,nlocs))
  allocate(gesZs(nlocs))

! copy geovals to local background arrays
  self%iflip = 0
  if (prs%vals(1,1) .lt. prs%vals(prs%nval,1) ) then
     self%iflip = 1
     write(err_msg,'(a)')'  ufo_gnssro_bndnbam_tlad_settraj:'//new_line('a')//                         &
                         '  Model vertical height profile is in descending order,'//new_line('a')// &
                         '  but NBAM requires it to be ascending order, need flip'
    call fckit_log%debug(err_msg)
    do k = 1, nlev
       gesT(k,:) = t%vals(nlev-k+1,:)
       gesQ(k,:) = q%vals(nlev-k+1,:)
    enddo
    do k = 1, nlev1
       gesP(k,:) = prs%vals(nlev1-k+1,:)
       gesH(k,:) = gph%vals(nlev1-k+1,:)
    enddo
  else  ! not flipping

    do k = 1, nlev
       gesT(k,:) = t%vals(k,:)
       gesQ(k,:) = q%vals(k,:)
    enddo
    do k = 1, nlev1
       gesP(k,:) = prs%vals(k,:)
       gesH(k,:) = gph%vals(k,:)
    enddo
  end if

! if all fields t and q are on mass layers,
! while p and z are on interface layers -- NBAM manner
  if ( nlev1 /= nlev ) then
     do k = nlev, 2, -1
        gesT(k,:) = half* (gesT(k,:) + gesT(k-1,:))
        gesQ(k,:) = half* (gesQ(k,:) + gesQ(k-1,:))
     enddo
  end if
        gesZs(:) = zs%vals(1,:)

! set obs space struture
  allocate(obsLat(nlocs))
  allocate(obsImpP(nlocs))
  allocate(obsLocR(nlocs))
  allocate(obsGeoid(nlocs))
  allocate(obsRecnum(nlocs))
  allocate(self%nlocs_begin(nrecs))
  allocate(self%nlocs_end(nrecs))

  call obsspace_get_db(obss, "MetaData", "latitude",         obsLat)
  call obsspace_get_db(obss, "MetaData", "impactParameterRO", obsImpP)
  call obsspace_get_db(obss, "MetaData", "earthRadiusCurvature", obsLocR)
  call obsspace_get_db(obss, "MetaData", "geoidUndulation", obsGeoid)
  call obsspace_get_recnum(obss, obsRecnum)

  if (  obsspace_has(obss,  "ObsDiag",  "superRefractionFlag") ) then
    allocate(obsSRflag(nlocs))
    call obsspace_get_db(obss,  "ObsDiag", "superRefractionFlag", obsSRflag)
    hasSRflag = 1 ! superRefractionFlag generated in hofx for real case run
  else
    hasSRflag = 0 ! superRefractionFlag does not exist in ctest
  end if

  self%nlocs_begin=1
  self%nlocs_end=1

  icount = 1
  do iobs = 1, nlocs-1
    if (obsRecnum(iobs+1) /= obsRecnum(iobs)) then
      icount = icount +1  !counting number of records
      self%nlocs_end(icount-1)= iobs
      self%nlocs_begin(icount) = iobs+1
    end if
  end do
  self%nlocs_end(nrecs)= nlocs

  if (icount /= nrecs) then
    write(err_msg,*) "record number is not consistent :", icount, nrecs
    call fckit_log%info(err_msg)
  end if
  if(cmp_strings(self%roconf%GSI_version, "GEOS")) then
     ngrd = nint(61.0/63.0 * nlev + 18)
     ModelsigLevelcheck = one
  else
     ngrd = min(80,self%roconf%ngrd) ! do not provide ngrd if EMC; default is 80km
     ModelsigLevelcheck = three
  endif

  if (trim(self%roconf%modeltopconfig) .eq. "true") then
    ! maximum = 80km for low top models refactor needed otherwise
    if ( self%roconf%modeltop > 80 ) then
       ! FAIL should either use a GSI config (NOAA or NASA)
       write(err_msg,*) myname, ' modeltopconfig not applicable for model tops over 80 km'// &
           '    ... try using default GSI GFS 16.3 setup for NOAA or NASA'
       call abor1_ftn(err_msg)
    end if
    ngrd = min(60, self%roconf%ngrd) ! maximum = 60 for low top models
    nlevAdd = (80 - self%roconf%modeltop) + 13  ! 80km is default model top basis for extrapolation
  end if

  nlevAdd = max(13, max(nlevAdd,self%roconf%nlevadd)) ! overwrite default value if present in yaml
  nlevExt = nlev + nlevAdd

  allocate(dhdp(nlev))
  allocate(dhdt(nlev))
  allocate(dzdh(nlev))
  allocate(radius(nlev))
  allocate(refIndex(nlev))
  allocate(ref(nlevExt))
  allocate(ref_tl(nlevExt))
  allocate(refXrad(0:nlevExt+1))
  allocate(refXrad_tl(0:nlevExt+1))
  allocate(refXrad_s(ngrd))
  allocate(dndp(nlev,nlev))
  allocate(dndq(nlev,nlev))
  allocate(dndt(nlev,nlev))
  allocate(dxidp(nlev,nlev))
  allocate(dxidt(nlev,nlev))
  allocate(dxidq(nlev,nlev))
  allocate(dbenddxi(nlev))
  allocate(dbenddn(nlev))
  allocate(lagConst_tl(3,nlevExt))
  allocate(lagConst(3,nlevExt))
  allocate(self%jac_t(nlev,nlocs))
  allocate(self%jac_q(nlev,nlocs))
  allocate(self%jac_prs(nlev1,nlocs))
  allocate(grids(ngrd))

! Inizialize some variables
  dxidt=zero; dxidp=zero; dxidq=zero
  dndt=zero;  dndq=zero;  dndp=zero

! tempprary manner to handle the missing hofx
  self%jac_t = missing

  do j = 1, ngrd
     grids(j) = (j-1) * ds
  end do

! calculate jacobian
  call gnssro_ref_constants(self%roconf%use_compress)

  iobs = 0
  rec_loop: do irec = 1, nrecs
    obs_loop: do icount = self%nlocs_begin(irec), self%nlocs_end(irec)

      iobs = iobs + 1

      if (hasSRflag == 1) then
         if (obsSRflag(iobs) > 0)  cycle obs_loop
      end if

      dxidt=zero; dxidp=zero; dxidq=zero
      dndt=zero;  dndq=zero;  dndp=zero

      do k = 1, nlev
!        geometric height nad dzdh jacobian
         call geop2geometric( obsLat(iobs), gesH(k,iobs)-gesZs(iobs), geomzi, dzdh(k))
!        guess radius
         radius(k) = geomzi + gesZs(iobs) + obsGeoid(iobs) + obsLocR(iobs)   ! radius r
!        guess refactivity, refactivity index,  and impact parameter
         call compute_refractivity(gesT(k,iobs), gesQ(k,iobs), gesP(k,iobs), &
                                 ref(k),self%roconf%use_compress)
         refIndex(k)= one + (r1em6*ref(k))
         refXrad(k) = refIndex(k) * radius(k)
      end do

!     data rejection based on model background !
!     (1) skip data beyond model levels
      call get_coordinate_value(obsImpP(iobs),sIndx,refXrad(1),nlev,"increasing")
      if (sIndx < ModelsigLevelcheck .or. sIndx > float(nlev))  cycle obs_loop

      do k = 1, nlev

!       jacobian for refractivity(N)
        fv    = rv_over_rd-one
        pw    = rd_over_rv+gesQ(k,iobs)*(one-rd_over_rv)
        q_coef =  n_b *gesP(k,iobs)/(gesT(k,iobs)**2*pw**2)*rd_over_rv + &
                  n_c *gesP(k,iobs)/(gesT(k,iobs)*  pw**2)*rd_over_rv
        p_coef =  n_a/gesT(k,iobs)   + &
                  n_b*gesQ(k,iobs)/(gesT(k,iobs)**2*pw) + &
                  n_c*gesQ(k,iobs)/(gesT(k,iobs)*pw)
        t_coef = -n_a*gesP(k,iobs)/gesT(k,iobs)**2 -  &
                  n_b*two*gesQ(k,iobs)*gesP(k,iobs)/(gesT(k,iobs)**3*pw) - &
                  n_c*gesQ(k,iobs)*gesP(k,iobs)/(gesT(k,iobs)**2*pw)

        dhdp=zero; dhdt=zero
        if(k > 1) then
           do j = 2, k
              dhdt(j-1)= rd_over_g*(log(gesP(j-1,iobs))-log(gesP(j,iobs)))
              dhdp(j)  = dhdp(j)-rd_over_g*(gesT(j-1,iobs)/gesP(j,iobs))
              dhdp(j-1)= dhdp(j-1)+rd_over_g*(gesT(j-1,iobs)/gesP(j-1,iobs))
           end do
        end if
        if(k == 1)then
           dndt(k,k)=dndt(k,k)+t_coef
           dndq(k,k)=dndq(k,k)+q_coef
           dndp(k,k)=dndp(k,k)+p_coef
        else
           dndt(k,k)=dndt(k,k)+half*t_coef
           dndt(k,k-1)=dndt(k,k-1)+half*t_coef
           dndq(k,k)=dndq(k,k)+half*q_coef
           dndq(k,k-1)=dndq(k,k-1)+half*q_coef
           dndp(k,k)=p_coef
        end if
        do j = 1, nlev
           dxidt(k,j)=r1em6*radius(k)*dndt(k,j) + refIndex(k)*dzdh(k)*dhdt(j)
           dxidq(k,j)=r1em6*radius(k)*dndq(k,j)
           dxidp(k,j)=r1em6*radius(k)*dndp(k,j) + refIndex(k)*dzdh(k)*dhdp(j)
        end do
      end do !nlev loop

      d_refXrad=refXrad(nlev)-refXrad(nlev-1)

      do k = 1, nlevAdd
         refXrad(nlev+k) = refXrad(nlev) + k*d_refXrad
         ref(nlev+k)     = ref(nlev+k-1)**2/ref(nlev+k-2)
      end do

      refXrad(0)=refXrad(3)
      refXrad(nlevExt+1)=refXrad(nlevExt-2)

      do klev = 1, nlev
         refXrad_tl      = zero
         refXrad_tl(klev)= one
         ref_tl          = zero
         ref_tl(klev)    = one
         lagConst        = zero
         lagConst_tl     = zero
         d_refXrad_tl    = refXrad_tl(nlev)-refXrad_tl(nlev-1)

         do k = 1, nlevAdd
            refXrad_tl(nlev+k) = refXrad_tl(nlev)+ k*d_refXrad_tl
            ref_tl(nlev+k)     = (two*ref(nlev+k-1)*ref_tl(nlev+k-1)/ref(nlev+k-2))-&
                                 (ref(nlev+k-1)**2/ref(nlev+k-2)**2)*ref_tl(nlev+k-2)
         end do
         refXrad_tl(0)=refXrad_tl(3)
         refXrad_tl(nlevExt+1)=refXrad_tl(nlevExt-2)

         do k=1,nlevExt
            call lag_interp_const_tl(lagConst(:,k),lagConst_tl(:,k),refXrad(k-1:k+1),refXrad_tl(k-1:k+1),3)
         end do

         intloop2: do j = 1, ngrd
           refXrad_s(j) = sqrt(grids(j)**2 + obsImpP(iobs)**2) !x_s^2=s^2+a^2
           call get_coordinate_value(refXrad_s(j),sIndx,refXrad(1:nlevExt),nlevExt,"increasing")
           indx=sIndx
           if (indx < nlevExt) then
             call lag_interp_smthWeights_tl(refXrad(indx-1:indx+2),refXrad_tl(indx-1:indx+2), &
                                            refXrad_s(j), lagConst(:,indx),lagConst_tl(:,indx),&
                                            lagConst(:,indx+1),lagConst_tl(:,indx+1),dw4,dw4_tl,4)
             if(indx==1) then
               dw4(4)=dw4(4)+dw4(1);dw4(1:3)=dw4(2:4);dw4(4)=zero
               dw4_tl(4)=dw4_tl(4)+dw4_tl(1);dw4_tl(1:3)=dw4_tl(2:4);dw4_tl(4)=zero
               indx=indx+1
             endif
             if(indx==nlevExt-1) then
               dw4(1)=dw4(1)+dw4(4); dw4(2:4)=dw4(1:3);dw4(1)=zero
               dw4_tl(1)=dw4_tl(1)+dw4_tl(4); dw4_tl(2:4)=dw4_tl(1:3);dw4_tl(1)=zero
               indx=indx-1
             endif

             dbetaxi=(r1em6/refXrad_s(j))*dot_product(dw4_tl,ref(indx-1:indx+2))
             dbetan =(r1em6/refXrad_s(j))*dot_product(dw4,ref_tl(indx-1:indx+2))

             if(j == 1)then
               dbenddxi(klev)=dbetaxi
               dbenddn(klev)=dbetan
             else
               dbenddxi(klev)=dbenddxi(klev)+two*dbetaxi
               dbenddn(klev)=dbenddn(klev)+two*dbetan
             end if
           else
             cycle  obs_loop
           end if ! obs inside the new "s" grids
         end do intloop2
         dbenddxi(klev)=-dbenddxi(klev)*ds*obsImpP(iobs)
         dbenddn(klev)=-dbenddn(klev)*ds*obsImpP(iobs)

       end do

       do k = 1, nlev
          self%jac_t(k,iobs)=zero
          self%jac_q(k,iobs)=zero
          self%jac_prs(k,iobs)=zero
          do j = 1, nlev
             self%jac_t(k,iobs)  = self%jac_t(k,iobs)+dbenddxi(j)*dxidt(j,k)+ &
                                                   dbenddn(j) * dndt(j,k)
             self%jac_q(k,iobs)  = self%jac_q(k,iobs)+dbenddxi(j)*dxidq(j,k)+ &
                                                     dbenddn(j) * dndq(j,k)
             self%jac_prs(k,iobs)= self%jac_prs(k,iobs)+dbenddxi(j)*dxidp(j,k)+ &
                                                    dbenddn(j) * dndp(j,k)
          end do
        end do
        if ( nlev /= nlev1)   self%jac_prs(nlev1,iobs)=  0.
    end do obs_loop
  end do rec_loop


  deallocate(obsLat)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)
  deallocate(gesT)
  deallocate(gesQ)
  deallocate(gesP)
  deallocate(gesH)
  deallocate(gesZs)
  deallocate(dhdp)
  deallocate(dhdt)
  deallocate(radius)
  deallocate(dzdh)
  deallocate(refIndex)
  deallocate(ref)
  deallocate(ref_tl)
  deallocate(refXrad)
  deallocate(refXrad_tl)
  deallocate(refXrad_s)
  deallocate(dndp)
  deallocate(dndq)
  deallocate(dndt)
  deallocate(dxidp)
  deallocate(dxidt)
  deallocate(dxidq)
  deallocate(dbenddxi)
  deallocate(dbenddn)
  deallocate(lagConst)
  deallocate(lagConst_tl)
  deallocate(obsRecnum)
  deallocate(grids)
  if (allocated(obsSRflag)) deallocate(obsSRflag)

end if
  self%ltraj = .true.

end subroutine ufo_gnssro_bndnbam_tlad_settraj
!----------------------------------------------------------------
subroutine ufo_gnssro_bndnbam_simobs_tl(self, geovals, hofx, obss)
  use gnssro_mod_transform
  implicit none
  class(ufo_gnssro_bndnbam_tlad), intent(in)    :: self
  type(ufo_geovals),             intent(in)    :: geovals
  real(kind_real),               intent(inout) :: hofx(:)
  type(c_ptr), value,            intent(in)    :: obss
  character(len=*), parameter     :: myname = "ufo_gnssro_bndnbam_simobs_tl"
  character(max_string)           :: err_msg
  integer                         :: nlev, nlev1, nlocs
  integer                         :: iobs, k, irec, icount
  type(ufo_geoval), pointer       :: t_tl, prs_tl, q_tl
  real(kind_real), allocatable    :: gesT_tl(:,:), gesP_tl(:,:), gesQ_tl(:,:)
  real(kind_real)                 :: sumIntgl

! check if trajectory was set
  if (.not. self%ltraj) then
      write(err_msg,*) myname, ' trajectory was not set!'
      call abor1_ftn(err_msg)
  endif


if (geovals%nlocs > 0 ) then
  hofx = missing_value(missing)

! check if nlocs is consistent in geovals & hofx
  if (self%nlocs /= size(hofx)) then
      write(err_msg,*) myname, ' error: nlocs inconsistent!'
      call abor1_ftn(err_msg)
  endif

  call ufo_geovals_get_var(geovals, var_ts,  t_tl)         ! air temperature
  call ufo_geovals_get_var(geovals, var_q,   q_tl)         ! specific humidity
  if (self%roconf%vertlayer .eq. "mass") then
    call ufo_geovals_get_var(geovals, var_prs,   prs_tl)       ! pressure
  else
    call ufo_geovals_get_var(geovals, var_prsi,  prs_tl)       ! pressure
  end if

  nlocs = self%nlocs
  nlev  = self%nlev
  nlev1 = self%nlev1
  allocate(gesT_tl(nlev, nlocs))
  allocate(gesQ_tl(nlev, nlocs))
  allocate(gesP_tl(nlev1,nlocs))

  if( self%iflip == 1 ) then
    do k = 1, nlev
       gesT_tl(k,:) = t_tl%vals(nlev-k+1,:)
       gesQ_tl(k,:) = q_tl%vals(nlev-k+1,:)
    enddo
    do k = 1, nlev1
       gesP_tl(k,:) = prs_tl%vals(nlev1-k+1,:)
    enddo
  else  ! not flipping
    do k =1, nlev
       gesT_tl(k,:) = t_tl%vals(k,:)
       gesQ_tl(k,:) = q_tl%vals(k,:)
    enddo
    do k =1, nlev1
       gesP_tl(k,:) = prs_tl%vals(k,:)
    enddo
  end if

! if all fields t and q are on mass layers,
! while p and z are on interface layers -- NBAM manner
  if ( nlev1 /= nlev ) then
     do k = nlev, 2, -1
        gesT_tl(k,:) = half* (gesT_tl(k,:) + gesT_tl(k-1,:))
        gesQ_tl(k,:) = half* (gesQ_tl(k,:) + gesQ_tl(k-1,:))
     enddo
  end if

  iobs = 0
  rec_loop: do irec = 1, self%nrecs
     obs_loop: do icount = self%nlocs_begin(irec), self%nlocs_end(irec)
        iobs = iobs + 1
        if (self%jac_t(1,iobs) /= missing ) then
        sumIntgl = 0.0
        do k = 1, nlev
            sumIntgl = sumIntgl + self%jac_t(k,iobs)*gesT_tl(k,iobs) &
                     + self%jac_q(k,iobs)*gesQ_tl(k,iobs) &
                     + self%jac_prs(k,iobs)*gesP_tl(k,iobs)

        end do
        hofx(iobs)=sumIntgl
        end if
     end do obs_loop
  end do rec_loop

  deallocate(gesT_tl)
  deallocate(gesP_tl)
  deallocate(gesQ_tl)

end if

end subroutine ufo_gnssro_bndnbam_simobs_tl

!----------------------------------------------------------------
subroutine ufo_gnssro_bndnbam_simobs_ad(self, geovals, hofx, obss)
  implicit none
  class(ufo_gnssro_bndnbam_tlad), intent(in)    :: self
  type(ufo_geovals),             intent(inout) :: geovals
  real(kind_real),               intent(in)    :: hofx(:)
  type(c_ptr), value,            intent(in)    :: obss
  type(ufo_geoval), pointer       :: t_ad, q_ad, prs_ad
  real(kind_real), allocatable    :: gesT_ad(:,:), gesP_ad(:,:), gesQ_ad(:,:)
  character(len=*), parameter   :: myname = "ufo_gnssro_bndnbam_simobs_ad"
  character(max_string)         :: err_msg
  integer                       :: nlocs, iobs, k, nlev,nlev1, icount, irec

! check if trajectory was set
  if (.not. self%ltraj) then
     write(err_msg,*) myname, ' trajectory wasnt set!'
     call abor1_ftn(err_msg)
  endif

if (self%nlocs > 0 ) then

! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
      write(err_msg,*) myname, ' error: nlocs inconsistent!'
      call abor1_ftn(err_msg)
  endif

  call ufo_geovals_get_var(geovals, var_ts,  t_ad)         ! air temperature
  call ufo_geovals_get_var(geovals, var_q,   q_ad)         ! specific humidity

  if (self%roconf%vertlayer .eq. "mass") then
    call ufo_geovals_get_var(geovals, var_prs,   prs_ad)       ! pressure
  else
    call ufo_geovals_get_var(geovals, var_prsi,  prs_ad)       ! pressure
  end if

  nlocs = self%nlocs
  nlev  = self%nlev
  nlev1 = self%nlev1
  allocate(gesT_ad(nlev,nlocs))
  allocate(gesQ_ad(nlev,nlocs))
  allocate(gesP_ad(nlev1,nlocs))

  gesT_ad = 0.0_kind_real
  gesQ_ad = 0.0_kind_real
  gesP_ad = 0.0_kind_real

  iobs = 0
  rec_loop: do irec = 1, self%nrecs
    obs_loop: do icount = self%nlocs_begin(irec), self%nlocs_end(irec)
      iobs = iobs + 1
      if (self%jac_t(1,iobs) /= missing .and. hofx(iobs) /= missing) then

          do k = 1,nlev1
             if (k == nlev + 1) then
                gesP_ad(k,iobs) = gesP_ad(k,iobs) + hofx(iobs)*self%jac_prs(k,iobs)
             else
                gesT_ad(k,iobs) = gesT_ad(k,iobs) + hofx(iobs)*self%jac_t(k,iobs)
                gesQ_ad(k,iobs) = gesQ_ad(k,iobs) + hofx(iobs)*self%jac_q(k,iobs)
                gesP_ad(k,iobs) = gesP_ad(k,iobs) + hofx(iobs)*self%jac_prs(k,iobs)
             end if
          end do

          if ( nlev1 /= nlev ) then
            do k = 2, nlev
               gesT_ad(k-1,iobs) = half*gesT_ad(k,iobs) +  gesT_ad(k-1,iobs)
               gesT_ad(k,iobs)   = half*gesT_ad(k,iobs)
               gesQ_ad(k-1,iobs) = half*gesQ_ad(k,iobs) +  gesQ_ad(k-1,iobs)
               gesQ_ad(k,iobs)   = half*gesQ_ad(k,iobs)
            enddo
          end if

      end if
    end do  obs_loop
  end do   rec_loop

  if( self%iflip == 1 ) then
    do k = 1, nlev
       t_ad%vals(nlev-k+1,:) = gesT_ad(k,:)
       q_ad%vals(nlev-k+1,:) = gesQ_ad(k,:)
    enddo
    do k = 1, nlev1
       prs_ad%vals(nlev1-k+1,:) = gesP_ad(k,:)
    enddo
  else  ! not flipping
    do k =1, nlev
       t_ad%vals(k,:) = gesT_ad(k,:)
       q_ad%vals(k,:) = gesQ_ad(k,:)
    enddo
    do k =1, nlev1
       prs_ad%vals(k,:) = gesP_ad(k,:)
    enddo
  end if

  deallocate(gesT_ad)
  deallocate(gesP_ad)
  deallocate(gesQ_ad)
end if

end subroutine ufo_gnssro_bndnbam_simobs_ad
! ------------------------------------------------------------------------------

subroutine ufo_gnssro_bndnbam_tlad_delete(self)
  implicit none
  class(ufo_gnssro_bndnbam_tlad), intent(inout) :: self
  character(len=*), parameter :: myname="ufo_gnssro_bndnbam_tlad_delete"

  self%nlocs = 0
  self%nrecs = 0
  self%nlev  = 0
  self%nlev1 = 0
  if (allocated(self%nlocs_begin)) deallocate(self%nlocs_begin)
  if (allocated(self%nlocs_end)) deallocate(self%nlocs_end)
  if (allocated(self%jac_t)) deallocate(self%jac_t)
  if (allocated(self%jac_prs)) deallocate(self%jac_prs)
  if (allocated(self%jac_q)) deallocate(self%jac_q)

  self%ltraj = .false.

end subroutine ufo_gnssro_bndnbam_tlad_delete

! ------------------------------------------------------------------------------
end module ufo_gnssro_bndnbam_tlad_mod
