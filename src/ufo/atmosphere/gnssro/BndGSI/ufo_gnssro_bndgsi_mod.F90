!  
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro bending angle observations following 
!> the NCEP/GSI (2018 Aug) implementation

module ufo_gnssro_bndgsi_mod
  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use vert_interp_mod
  use ufo_basis_mod,      only: ufo_basis
  use lag_interp_mod,     only: lag_interp_const, lag_interp_smthWeights
  use obsspace_mod

  implicit none
  public             :: ufo_gnssro_BndGSI
  private

  !> Fortran derived type for gnssro trajectory
  type, extends(ufo_basis) :: ufo_gnssro_BndGSI
  contains
   procedure :: simobs    => ufo_gnssro_bndgsi_simobs
  end type ufo_gnssro_BndGSI

  contains
! ------------------------------------------------------------------------------
  subroutine ufo_gnssro_bndgsi_simobs(self, geovals, hofx, obss)
      use gnssro_mod_constants
      use gnssro_mod_transform
      use gnssro_mod_grids, only: get_coordinate_value
      use obsspace_mod, only: obsspace_get_db
      implicit none
      class(ufo_gnssro_bndGSI), intent(in)    :: self
      type(ufo_geovals),        intent(in)    :: geovals
      real(kind_real),          intent(inout) :: hofx(:)
      type(c_ptr), value,       intent(in)    :: obss
    
      character(len=*), parameter     :: myname_ ="ufo_gnssro_bndgsi_simobs"
      logical, parameter              :: use_compress = .true.
      real(kind_real), parameter      :: r1em6 = 1.0e-6_kind_real
      real(kind_real), parameter      :: r1em3 = 1.0e-3_kind_real
      real(kind_real), parameter      :: six   = 6.0_kind_real
      real(kind_real), parameter      :: ds    = 10000.0_kind_real
      integer, parameter              :: nlevAdd = 13 !num of additional levels on top of exsiting model levels
      integer, parameter              :: newAdd  = 20 !num of additional levels on top of extended levels
      integer, parameter              :: ngrd    = 80 !num of new veritcal grids for bending angle computation
      integer, parameter              :: max_string    = 800
      real(kind_real), parameter      :: miss_values   = -99999.00
      real(kind_real), parameter      :: crit_gradRefr = 157.0_kind_real !criteria for the refractivity gradient

      character(max_string)           :: err_msg
      integer                         :: iobs, ilev, igrd, klev
      integer                         :: nlev, nlev1, nobs, nlevExt, nlevCheck
      real(kind_real)                 :: rnlevExt
      real(kind_real)                 :: w4(4), dw4(4)
      type(ufo_geoval), pointer       :: t, gphi, prsi, q
      real(kind_real), allocatable    :: gesT(:,:), gesZi(:,:), gesPi(:,:), gesQ(:,:)
      real(kind_real), allocatable    :: refr(:,:), radius(:,:)
      real(kind_real), allocatable    :: refrIndex(:), refrXrad(:), geomzi(:), refrXrad_new(:)
      real(kind_real), allocatable    :: lagConst(:,:)
      real(kind_real), allocatable    :: obsLat(:), obsImpP(:), obsLocR(:), obsGeoid(:)
      real(kind_real)                 :: specHmean, tmean
      real(kind_real)                 :: d_refrXrad
      real(kind_real)                 :: derivRefr_s(ngrd),grids(ngrd),refrXrad_s(ngrd)
      real(kind_real)                 :: sIndx, sIndxExt
      integer                         :: indx
      real(kind_real)                 :: bndIntgd, bendingAngle, obsImpH, gradRefr
      logical                         :: obs_check, qc_layer_SR
      integer                         :: nobs_outIntgl
      integer                         :: count_SR, top_layer_SP, top_layer_SR, bot_layer_SR !for super refraction
      integer                         :: count_rejection
       real(kind_real)                 :: jacob    
      nobs = 0
      nlev = 0
      nlev1 = 0
      nlevExt = 0
      sIndxExt = one
      hofx(:) = miss_values

      ! check if nobs is consistent in geovals & hofx
      if (geovals%nobs /= size(hofx)) then
        write(err_msg,*) myname_, ' error: nobs inconsistent!'
        call abor1_ftn(err_msg)
      endif
      ! check if prsi (pressure at model interface levels) variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_prsi, prsi)
      ! check if t (temperature) variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_t, t)
      ! check if q(specific humidity) variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_q, q)
      ! check if gphi (geopotential height at model interface levels) variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_zi, gphi)

 
      nlev  = t%nval ! number of model levels
      nlev1 = prsi%nval ! number of model interface levels 
      if (nlev1 /= nlev+1) then
         write(err_msg,*) myname_, ' Numbers of vertical profiles, nlev1 and nlev, dont match'
         call abor1_ftn(err_msg)
      endif
      nlevExt = nlev + nlevAdd
      nlevCheck = min(23, nlev) !number of levels to check super refraction
      nobs  = geovals%nobs ! number of observations

      allocate(gesPi(nlev1,nobs)) 
      allocate(gesZi(nlev1,nobs)) 
      allocate(gesT(nlev,nobs)) 
      allocate(gesQ(nlev,nobs)) 

      !FV3 background is from top to bottom. Reserse the vertical order 
      do ilev=1, nlev
         gesT(ilev,:) = t%vals(nlev-ilev+1,:)
         gesQ(ilev,:) = q%vals(nlev-ilev+1,:)
      enddo

      do ilev=1, nlev1
         gesPi(ilev,:) = prsi%vals(nlev1-ilev+1,:)
         gesZi(ilev,:) = gphi%vals(nlev1-ilev+1,:)
      enddo

      do igrd = 0, ngrd-1
          grids(igrd+1) = igrd * ds
      end do 

      allocate(geomzi(nlev))      !geometric height at interface model levels
      allocate(radius(nlev,nobs)) !distance between earth center to model interface level

      allocate(refr(nlevExt,nobs))    !refractivity N 
      allocate(refrIndex(nlev))       !refractivity index n
      allocate(refrXrad(0:nlevExt+1)) !x=nr, r: radius
      allocate(lagConst(3,nlevExt))   !x=nr, r: radius
      allocate(refrXrad_new(nlevExt+newAdd))

      allocate(obsLat(nobs))
      allocate(obsImpP(nobs))
      allocate(obsLocR(nobs))
      allocate(obsGeoid(nobs))

      call obsspace_get_db(obss, "Metadata", "Latitude", obsLat)
      call obsspace_get_db(obss, "Metadata", "IMPP", obsImpP)   !observed impact parameter; meter
      call obsspace_get_db(obss, "Metadata", "ELRC", obsLocR)   !local radius of earth; meter
      call obsspace_get_db(obss, "Metadata", "GEODU", obsGeoid) !Geoid; meter

      nobs_outIntgl = 0 !initialize count of observations out of integral grids  
      count_rejection = 0

      obs_loop: do iobs = 1, nobs 
         do ilev = 1,nlev 

         ! compute guess geometric height from geopotential height at model interface levels
            call geop2geometric( obsLat(iobs), gesZi(ilev,iobs), geomzi(ilev), jacob)

         ! compute guess radius 
            radius(ilev,iobs) = geomzi(ilev) + obsGeoid(iobs) + obsLocR(iobs)   ! radius r

         ! compute guess refractivity and refractivity index at model interface levels
            if(ilev > 1) then
               specHmean = (gesQ(ilev,iobs) + gesQ(ilev-1,iobs))/two
               tmean = (gesT(ilev,iobs) + gesT(ilev-1,iobs) )/two
            else
               specHmean = gesQ(1,iobs)
               tmean = gesT(1,iobs)
            endif
            call compute_refractivity(tmean, specHmean, gesPi(ilev,iobs), refr(ilev,iobs),use_compress) !refr N
            refrIndex(ilev) = one + (r1em6*refr(ilev,iobs))         ! refr index n
            refrXrad(ilev)  = refrIndex(ilev) * radius(ilev,iobs)   ! x=nr

         end do 

         ! data rejection based on model background !
         ! (1) skip data below the model levels
         call get_coordinate_value(obsImpP(iobs), sIndx,refrXrad(1),nlev,"increasing")
         if (sIndx < one .or. sIndx > float(nlev)) then 
             count_rejection = count_rejection + 1
             cycle
         end if

         ! (2) super-refraction
         qc_layer_SR=.false.
         count_SR=0
         top_layer_SR=0
         bot_layer_SR=0
         obsImpH = (obsImpP(iobs) - obsLocR(iobs)) * r1em3 !impact heigt: a-r_earth
         if (obsImpH <= six) then
            do klev=nlevCheck,1,-1

               ! check for model SR layer 
               gradRefr = 1000.0_kind_real * (refr(klev+1,iobs)-refr(klev,iobs)) / (radius(klev+1,iobs)-radius(klev,iobs))

               if (abs(gradRefr)>= half*crit_gradRefr) then  !Super refractivity - likely, to be used in obs SR qc
                  qc_layer_SR=.true.   !SR-likely layer detected
               endif

               !if (((ref_rad(klev+1)-ref_rad(klev))/(radius(klev+1,iobs)-radius(klev,iobs))) < zero) then
               if (abs(gradRefr) >= 0.75_kind_real*crit_gradRefr) then  !relax to close-to-SR conditions
                  count_SR=count_SR+1 ! layers of SR

                  if (count_SR > 1 ) then
                     !if(abs(bot_layer_SR-klev) > 1) write(6,*) 'WARNING GPSRO: non-consecutive SR layers'
                     bot_layer_SR=klev
                  else
                     top_layer_SR=klev
                     bot_layer_SR=top_layer_SR
                  endif

               endif
            end do
            if (top_layer_SR >= 1 .and. obsImpP(iobs) <= refrXrad(top_layer_SR+2)) then !obs inside model SR layer
               count_rejection = count_rejection + 1
               cycle
            end if
          endif


         ! Extend atmosphere above interface level nlev
         d_refrXrad = refrXrad(nlev) - refrXrad(nlev-1)
         do ilev=1,nlevAdd
            refrXrad(nlev+ilev)=refrXrad(nlev)+ ilev*d_refrXrad ! extended x_i
            refr(nlev+ilev,iobs)=refr(nlev+ilev-1,iobs)**2/refr(nlev+ilev-2,iobs) ! exended N_i
         end do

         refrXrad(0)=refrXrad(3)
         refrXrad(nlevExt+1)=refrXrad(nlevExt-2)

         do ilev = 1,nlevExt
            call lag_interp_const(lagConst(:,ilev),refrXrad(ilev-1:ilev+1),3)
         enddo

         ! Set up a new equally-spaced vertical grid for integral 

         derivRefr_s = zero
         grids_loop: do igrd =1,ngrd
         !use the new grids (s) for bending angle computation
           refrXrad_s(igrd)=sqrt(grids(igrd)**2 + obsImpP(iobs)**2) !x_s^2=s^2+a^2

           call get_coordinate_value(refrXrad_s(igrd), sIndx,refrXrad(1:nlevExt),nlevExt,"increasing")

           rnlevExt = float(nlevExt)

           if (sIndx < rnlevExt) then  !obs inside the new grid
           !HS if (sIndx > zero .and. sIndx < rnlevExt) then  !obs inside the new grid
              indx=sIndx

!             Compute derivative at new grids (dN/ds) using Lagrange interpolators
              call lag_interp_smthWeights(refrXrad(indx-1:indx+2),refrXrad_s(igrd),&
                   lagConst(:,indx),lagConst(:,indx+1),&
                   w4,dw4,4)
              if(indx==1) then
                 w4(4)=w4(4)+w4(1); w4(1:3)=w4(2:4);w4(4)=zero
                 dw4(4)=dw4(4)+dw4(1);dw4(1:3)=dw4(2:4);dw4(4)=zero
                 indx=indx+1
              endif
              if(indx==nlevExt-1) then
                 w4(1)=w4(1)+w4(4); w4(2:4)=w4(1:3);w4(1)=zero
                 dw4(1)=dw4(1)+dw4(4); dw4(2:4)=dw4(1:3);dw4(1)=zero
                 indx=indx-1
              endif
              !need make sure dot_product is available or add the code
              derivRefr_s(igrd)=dot_product(dw4,refr(indx-1:indx+2,iobs)) !derivative dN/dx_s
              derivRefr_s(igrd)=max(zero,abs(derivRefr_s(igrd)))

           else
              obs_check=.true.
              nobs_outIntgl=nobs_outIntgl+1
              d_refrXrad=refrXrad(nlev)-refrXrad(nlev-1)
              do klev=1,newAdd
                 refrXrad_new(nlevExt+klev)=refrXrad(nlevExt)+ klev*d_refrXrad ! extended x_i
              end do
              do klev=1,nlevExt
                 refrXrad_new(klev)=refrXrad(klev)
              enddo
              call get_coordinate_value(refrXrad_s(igrd), sIndx,refrXrad_new(1:nlevExt+newAdd),nlevExt+newAdd,"increasing")
              sIndxExt=max(sIndx,sIndxExt)
           endif !obs in new grid
         end do grids_loop

!        bending angle (radians)
         bendingAngle = ds*derivRefr_s(1)/refrXrad_s(1)
         do igrd = 2,ngrd
            bndIntgd     = ds*derivRefr_s(igrd)/refrXrad_s(igrd)
            bendingAngle = bendingAngle + two*bndIntgd
         end do
         bendingAngle=r1em6 * obsImpP(iobs) * bendingAngle
         hofx(iobs) = bendingAngle

      end do obs_loop

      write(6,*) 'bndGSI: hofx ', &
                 'min = ', minval(hofx, mask=hofx > miss_values), 'min index = ', minloc(hofx), &
                 'max = ', maxval(hofx, mask=hofx > miss_values), 'max index = ', maxloc(hofx)
      write(6,*) 'bndGSI: ', count_rejection, ' out of ', nobs, ' rejected due to model vertical range and super refraction'

      !for tuning the nlevExt. New grids (s) should be in range  with nlevExt. If not, adjust the hardwired 
      if (nobs_outIntgl>=1) then
         write(6,*)'bndGSI: Warning',nobs_outIntgl,'obs outside integration grid. Increase nlevExt to',&
         int(sIndxExt)
      endif

      deallocate(obsLat)
      deallocate(obsImpP)
      deallocate(obsLocR)
      deallocate(obsGeoid)

      deallocate(gesPi) 
      deallocate(gesZi) 
      deallocate(gesT) 
      deallocate(gesQ) 
      deallocate(refr) 
      deallocate(refrIndex) 
      deallocate(refrXrad) 
      deallocate(geomzi) 
      deallocate(radius) 
      deallocate(lagConst) 
      deallocate(refrXrad_new) 

   end subroutine ufo_gnssro_bndgsi_simobs
! ------------------------------------------------------------------------------
end module ufo_gnssro_bndgsi_mod
