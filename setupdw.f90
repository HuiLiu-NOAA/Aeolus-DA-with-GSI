module dw_setup
  implicit none
  private
  public:: setup
        interface setup; module procedure setupdw; end interface

contains
!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  setupdw --- Compute rhs of oi for doppler lidar winds
!
! !INTERFACE:
!
subroutine setupdw(obsLL,odiagLL,lunin,mype,bwork,awork,nele,nobs,is,conv_diagsave)

! !USES:

  use mpeu_util, only: die,perr,getindex
  use kinds, only: r_kind,r_single,r_double,i_kind


  use qcmod, only: dfact,dfact1,npres_print,ptop,pbot

  use gridmod, only: nsig,get_ijk

  use guess_grids, only: hrdifsig,geop_hgtl,ges_lnprsl,&
       nfldsig,sfcmod_gfs,sfcmod_mm5,comp_fact10

  use constants, only: grav_ratio,flattening,grav,zero,rad2deg,deg2rad, &
       grav_equator,one,two,somigliana,semi_major_axis,eccentricity,r1000,&
       wgtlim, r10
  use constants, only: tiny_r_kind,half,cg_term,huge_single

  use obsmod, only: rmiss_single,lobsdiag_forenkf
  use obsmod, only: netcdf_diag, binary_diag, dirname, ianldate
  use nc_diag_write_mod, only: nc_diag_init, nc_diag_header, nc_diag_metadata, &
       nc_diag_write, nc_diag_data2d
  use nc_diag_read_mod, only: nc_diag_read_init, nc_diag_read_get_dim, nc_diag_read_close
  use m_obsdiagNode, only: obs_diag
  use m_obsdiagNode, only: obs_diags
  use m_obsdiagNode, only: obsdiagLList_nextNode
  use m_obsdiagNode, only: obsdiagNode_set
  use m_obsdiagNode, only: obsdiagNode_get
  use m_obsdiagNode, only: obsdiagNode_assert

  use obsmod, only: lobsdiagsave,nobskeep,lobsdiag_allocated,time_offset
  use m_obsNode, only: obsNode
  use m_dwNode, only: dwNode
  use m_dwNode, only: dwNode_appendto
  use m_obsLList, only: obsLList
  use obsmod, only: luse_obsdiag
  use gsi_4dvar, only: nobs_bins,hr_obsbin
  use state_vectors, only: svars3d, levels, nsdim

  use jfunc, only: last, jiter, miter, jiterstart
  use convinfo, only: nconvtype,cermin,cermax,cgross,cvar_b,cvar_pg,ictype
  use convinfo, only: icsubtype

  use m_dtime, only: dtime_setup, dtime_check

  use gsi_bundlemod, only : gsi_bundlegetpointer
  use gsi_metguess_mod, only : gsi_metguess_get,gsi_metguess_bundle
  use sparsearr, only: sparr2, new, size, writearray, fullarray

  implicit none

! !INPUT PARAMETERS:

  type(obsLList ),target,dimension(:),intent(in):: obsLL
  type(obs_diags),target,dimension(:),intent(in):: odiagLL

  integer(i_kind)                                  ,intent(in   ) :: lunin   ! unit from which to read observations
  integer(i_kind)                                  ,intent(in   ) :: mype    ! mpi task id
  integer(i_kind)                                  ,intent(in   ) :: nele    ! number of data elements per observation
  integer(i_kind)                                  ,intent(in   ) :: nobs    ! number of observations
  integer(i_kind)                                  ,intent(in   ) :: is      ! ndat index
  logical                                          ,intent(in   ) :: conv_diagsave ! logical to save innovation dignostics

! !INPUT/OUTPUT PARAMETERS:
                                                  ! array containing information about ...
  real(r_kind),dimension(100+7*nsig)               ,intent(inout) :: awork !  data counts and gross checks
  real(r_kind),dimension(npres_print,nconvtype,5,3),intent(inout) :: bwork !  obs-ges stats 

! !DESCRIPTION:  For doppler lidar wind observations, this routine
!  \begin{enumerate}
!         \item reads obs assigned to given mpi task (geographic region),
!         \item simulates obs from guess,
!         \item apply some quality control to obs,
!         \item load weight and innovation arrays used in minimization
!         \item collects statistics for runtime diagnostic output
!         \item writes additional diagnostic information to output file
!  \end{enumerate}
!
! !REVISION HISTORY:
!   1998-05-15  yang, weiyu
!   1999-08-24  derber, j., treadon, r., yang, w., first frozen mpp version
!   2004-06-17  treadon - update documentation
!   2004-07-15  todling - protex-compliant prologue; added intent/only's
!   2004-10-06  parrish - increase size of dwork array for nonlinear qc
!   2004-11-22  derber - remove weight, add logical for boundary point
!   2004-12-22  treadon - move logical conv_diagsave from obsmod to argument list
!   2005-03-02  dee - remove garbage from diagnostic file
!   2005-03-09  parrish - nonlinear qc change to account for inflated obs error
!   2005-05-27  derber - level output change
!   2005-07-27  derber  - add print of monitoring and reject data
!   2005-09-28  derber  - combine with prep,spr,remove tran and clean up
!   2005-10-14  derber  - input grid location and fix regional lat/lon
!   2005-11-03  treadon - correct error in ilone,ilate data array indices
!   2005-11-29  derber - remove psfcg and use ges_lnps instead
!   2006-01-31  todling/treadon - store wgt/wgtlim in rdiagbuf(6,ii)
!   2006-02-02  treadon - rename lnprsl as ges_lnprsl
!   2006-02-24  derber  - modify to take advantage of convinfo module
!   2006-05-30  derber,treadon - modify diagnostic output
!   2006-06-06  su - move to wgtlim to constants module
!   2006-07-28  derber  - modify to use new inner loop obs data structure
!                       - unify NL qc
!   2006-07-31  kleist - use ges_ps
!   2006-08-28      su - fix a bug in variational qc
!   2007-03-19  tremolet - binning of observations
!   2007-06-05  tremolet - add observation diagnostics structure
!   2007-08-28      su - modify the gross check error 
!   2008-05-23  safford - rm unused vars
!   2008-12-03  todling - changed handle of tail%time
!   2009-03-19  mccarty/brin - set initial obs error to that from bufr
!   2009-08-19  guo     - changed for multi-pass setup with dtime_check().
!   2010-08-01  woollen  - add azmth and elevation angle check in duplication (denoted as jsw)
!   2010-09-01  masutani - remove repe_dw and get representativeness error from coninfo  (msq)
!   2010-11-20  woollen -  dpress is adjusted by zsges  (denoted as jsw)
!   2010-12-03  woollen -  fix low level adjust ment to factw (denoted as jsw)
!   2010-12-06  masutani - pass subtype kx to identify KNMI product  (msq)
!   2011-04-18  mccarty - updated kx determination for ADM, modified presw calculation
!   2011-05-05  mccarty - re-removed repe_dw, added +1 conditional for reproducibility on ADM
!   2011-05-26  mccarty - moved MSQ error logic from read_lidar
!   2013-01-26  parrish - change from grdcrd to grdcrd1, tintrp2a to tintrp2a1, tintrp2a11,
!                           tintrp3 to tintrp31 (to allow successful debug compile on WCOSS)
!   2013-10-19  todling - metguess now holds background
!   2014-01-28  todling - write sensitivity slot indicator (ioff) to header of diagfile
!   2014-12-30  derber - Modify for possibility of not using obsdiag
!   2015-10-01  guo   - full res obvsr: index to allow redistribution of obsdiags
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!   2016-06-24  guo     - fixed the default value of obsdiags(:,:)%tail%luse to luse(i)
!                       . removed (%dlat,%dlon) debris.
!   2016-11-29  shlyaeva - save linearized H(x) for EnKF
!   2017-02-06  todling - add netcdf_diag capability; hidden as contained code
!   2017-02-09  guo     - Remove m_alloc, n_alloc.
!                       . Remove my_node with corrected typecast().
!   2019-07-26  hliu  - add Bias correction, QCs, and errors of Aeolus L2B HLOS wind component
!   2019-11-16  hliu  - updates of blacklist of FM-B data
!   2020-04-20  hliu  - updates for FM-B data (thinning of Mie winds, error specification, bias correction)
!   2020-06-22  hliu  - Add diag output of Aeolus winds in netcdf format
!   2020-08-11  hliu  - Implement a speed-dependent bias correction of Aeolus - GFS
!   2020-12-06  hliu  - implement TLS regression bias correction of Aeolus - GFS

!
! !REMARKS:
!   language: f90
!   machine:  ibm RS/6000 SP; SGI Origin 2000; Compaq/HP
!
! !AUTHOR: 
!   yang             org: np20                date: 1998-05-15
!
!EOP
!-------------------------------------------------------------------------

! Declare external calls for code analysis
  external:: tintrp3
  external:: grdcrd1
  external:: stop2

! Declare local parameters
  real(r_kind),parameter:: r0_001 = 0.001_r_kind
  real(r_kind),parameter:: r8 = 8.0_r_kind
  real(r_kind),parameter:: ten = 10.0_r_kind
  character(len=*),parameter:: myname="setupdw"
  real(r_kind),parameter:: dmiss = 9.0e+10_r_kind !missing value for msq error adj - wm

  integer:: kikx, likx, nnn
  real(r_kind) :: lat1, lon1, hhh1, lat2, lon2, hhh2, dist_km, dist_hh
!hliu

! Declare local variables
  
  real(r_double) rstation_id
  real(r_kind) sinazm,cosazm,scale
  real(r_kind) ratio_errors,dlat,dlon,dtime,error,dpres,zsges    !jsw
  real(r_kind) dlnp,pobl,rhgh,rsig,rlow
! hliu  real(r_kind) zob,termrg,dz,termr,sin2,termg
  real(r_kind) zob,termrg,dz,termr,sin2,termg, zobt, zobb,zobt0, zobb0, zobt2, zobb2 
  real(r_kind) ugesindwt,vgesindwt,ugesindwb,vgesindwb,wshear
  real(r_kind) ugesindwt2, vgesindwt2, ugesindwb2, vgesindwb2
! hliu

  real(r_kind) sfcchk,slat,psges,dwwind
  real(r_kind) ugesindw,vgesindw,factw,presw, dwges, dwobs
  real(r_kind) residual,obserrlm,obserror,ratio,val2
  real(r_kind) ress,ressw
  real(r_kind) val,valqc,ddiff,rwgt,sfcr,skint
  real(r_kind) cg_dw,wgross,wnotgross,wgt,arg,term,exp_arg,rat_err2
  real(r_kind) errinv_input,errinv_adjst,errinv_final
  real(r_kind) err_input,err_adjst,err_final,tfact
  real(r_kind),dimension(nele,nobs):: data
  real(r_kind),dimension(nobs):: dup
  real(r_kind),dimension(nsig):: hges,zges,prsltmp
  real(r_single),allocatable,dimension(:,:)::rdiagbuf

  integer(i_kind) mm1,ikxx,nn,isli,ibin,ioff,ioff0
  integer(i_kind) jsig
  integer(i_kind) i,nchar,nreal,k,j,k1,jj,l,ii,k2
  integer(i_kind) ier,ilon,ilat,ihgt,ilob,id,itime,ikx,iatd,inls,incls
  integer(i_kind) iazm,ielva,iuse,ilate,ilone, idsat
  integer(i_kind) idomsfc,isfcr,iff10,iskint

  real(r_kind) :: delz
  type(sparr2) :: dhx_dx
  real(r_single), dimension(nsdim) :: dhx_dx_array
  integer(i_kind) :: iz, u_ind, v_ind, nind, nnz

  character(8) station_id
  character(8),allocatable,dimension(:):: cdiagbuf

  logical,dimension(nobs):: luse,muse
  integer(i_kind),dimension(nobs):: ioid ! initial (pre-distribution) obs ID
  logical proceed

  logical:: in_curbin,in_anybin, save_jacobian
  type(dwNode),pointer:: my_head
  type(obs_diag),pointer:: my_diag
  type(obs_diags),pointer:: my_diagLL
  
  equivalence(rstation_id,station_id)

  real(r_kind),allocatable,dimension(:,:,:  ) :: ges_ps
  real(r_kind),allocatable,dimension(:,:,:  ) :: ges_z
  real(r_kind),allocatable,dimension(:,:,:,:) :: ges_u
  real(r_kind),allocatable,dimension(:,:,:,:) :: ges_v

  type(obsLList),pointer,dimension(:):: dwhead

!hliu ---------- for Aeolus TLS bias correction ----------
  integer(i_kind), parameter:: nlat=19, nray= 22, nmie=21
  real(r_kind), Parameter :: pi = 3.1415926, d2r= pi/180.0

real(r_single) rayasc(2,nray,nlat), raydes(2,nray,nlat)
real(r_single) mieasc(2,nmie,nlat), miedes(2,nmie,nlat)
real(r_single) coefray(2,nlat), coefmie(2,nlat)
real(r_single) hght_ray(nray), hght_mie(nmie), latint(nlat)

real(r_single) dwwindt, dwwindb, lat0, obsdif, latdif, bc6, dwwindt2, dwwindb2
integer(i_kind) n, lnd, lnd2, ireg, kk
real(r_single)  x1, x2, f1, f2, x0

real(r_single)  rayerr(3,nray), mieerr(3,nmie)
real(r_single)  raydel(3,nray), miedel(3,nmie),delta   !ratio of errob/errbg HL
real(r_single)  beta0, beta1, gesopt, obsopt, dinv, del

 data hght_ray /24.5, 22.5, 20.5, 18.5, 17.5, 16.5, 15.5, 14.5, 13.5, &
                12.5, 11.5, 10.5,  9.5,  8.5,  7.5,  6.5,  5.5,  4.5, &
                 3.5,  2.5,  1.5,  0.5/

 data hght_mie /22.5, 20.5, 18.5, 17.5, 16.5, 15.5, 14.5, 13.5, 12.5, 11.5, &
                10.5,  9.5,  8.5,  7.5,  6.5,  5.5,  4.5,  3.5,  2.5,  1.5, 0.5/

! hliu-----------------------------------
!  print*, 'jiter0,nobs== ', jiter, nobs
   if( jiter == 1) then
      hght_ray = hght_ray *1000.0   !(m)
      hght_mie = hght_mie *1000.0   !(m)
   endif

      do j=1, nlat
       latint(j) = -90.0 + 10*(j-1)
      enddo

!hliu -----------------------------------
      call read_bias_correctionTLSonline
!hliu -----------------------------------

  dwhead => obsLL(:)

  save_jacobian = conv_diagsave .and. jiter==jiterstart .and. lobsdiag_forenkf

! Check to see if required guess fields are available
  call check_vars_(proceed)
  if(.not.proceed) return  ! not all vars available, simply return

! If require guess vars available, extract from bundle ...
  call init_vars_

!*********************************************************************************
! Read and reformat observations in work arrays.  
  read(lunin)data,luse,ioid

!    index information for data array (see reading routine)
  ikxx=1      ! index of ob type
  ilon=2      ! index of grid relative obs location (x)
  ilat=3      ! index of grid relative obs location (y)
  itime=4     ! index of observation time in data array
  ihgt=5      ! index of obs vertical coordinate in data array(height-m)
  ielva=6     ! index of elevation angle(radians)
  iazm=7      ! index of azimuth angle(radians) in data array
  inls=8      ! index of accumulation length (m)   !hliu
  incls=9     ! index of number of cloud laser shots  == 0 for now
  iatd=10     ! index of layer depth     
  ilob=11     ! index of lidar observation
  ier=12      ! index of obs error
  id=13       ! index of station id
!hliu --------------------------------------------
  idsat=  14  ! satellite id   hliu
  iuse=   15  ! index of use parameter
  idomsfc=16  ! index of dominate surface type
  iskint= 17  ! index of skin temperature
  iff10 = 18  ! index of 10 m wind factor
  isfcr = 19  ! index of surface roughness
  ilone=  20  ! index of longitude (degrees)
  ilate=  21  ! index of latitude (degrees)

! iuse=14     ! index of use parameter
! idomsfc=15  ! index of dominate surface type
! iskint=16   ! index of skin temperature
! iff10 = 17  ! index of 10 m wind factor
! isfcr = 18  ! index of surface roughness
! ilone=19    ! index of longitude (degrees)
! ilate=20    ! index of latitude (degrees)
!hliu ---------------------

  do i=1,nobs
!hliu for monitoring data (iuse = -1 in convinfo), usage set to 100 in read_lidar.f90
! and  muse = .false. but error> 0.0. So, the thining is not done for monitored OBS
     data(incls, i) = zero       ! use it for setting monitoring of aeolus wind for now

     muse(i)=nint(data(iuse,i)) <= jiter
  end do

!------------------------------------------------
!hliu  thin Mie winds to 90km interval (for every layer)
!------------------------------------------------
!  if( jiter == 1) then      ! jiter= 1 ->3
    nnn = 0
    do k=1,nobs
     kikx=nint(data(ikxx, k))
     if( icsubtype(kikx)==11 .and. muse(k) ) nnn = nnn + 1
    enddo
    print*, 'total mie winds= ', nnn
 
   do k=1,nobs
      do l=k+1,nobs
         kikx=nint(data(ikxx, k))
         likx=nint(data(ikxx, l))
 
      if( icsubtype(kikx)==11 .and. icsubtype(likx)==11) then   ! Mie cloudy-sky
         lat1 = data(ilate, k) *d2r
         lon1 = data(ilone, k) *d2r
         hhh1 = data(ihgt,  k)      ! m

         lat2 = data(ilate, l) *d2r
         lon2 = data(ilone, l) *d2r
         hhh2 = data(ihgt,  l)      ! m
!! -------------  distance between the two locations --------------------------
!!   Distance,d = 6377.0 * acos((sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2)
!!                * cos(lon2-lon1))
!!-----------------------------------------------------------------------------
         dist_km = 6377.0 * acos((sin(lat1) * sin(lat2)) + &
                   cos(lat1) * cos(lat2) * cos((lon2-lon1)) )
         dist_hh = abs(hhh1-hhh2)     ! m
 
!! -----------   the minimum layer thickness is 250m ---------
   if( dist_km < 90.0  .and. dist_hh < 500.0  .and. muse(k) .and. muse(l) &
            .and. data(ier,k) < r1000 )then
                   muse(l) =.false.            ! still monitored
             data( iuse,l) = 206
             data(incls,l) = 101       ! use it for identifier for thinning for now
         end if
       end if
      end do
   end do

   nnn = 0
   do k=1,nobs
     kikx=nint(data(ikxx,k))
     if( icsubtype(kikx)==11 .and. muse(k) ) nnn = nnn + 1
   enddo
    print*, 'total mie winds after thinning= ', nnn
!  endif
!!hliu -------- thinning of Mie winds ------


  dup=one
  do k=1,nobs
     do l=k+1,nobs
        if(data(ilat,k) == data(ilat,l) .and.  &
           data(ilon,k) == data(ilon,l) .and.  &
           data(ihgt,k) == data(ihgt,l) .and. &
           data(iazm,k) == data(iazm,l) .and. &     ! jsw check azmth angle
           data(ielva,k) == data(ielva,l) .and. &   ! jsw check eleveaiton angle
           data(ier,k) < r1000 .and. data(ier,l) < r1000 .and. &
           muse(k) .and. muse(l))then
           tfact=min(one, abs(data(itime,k)-data(itime,l))/dfact1)
           dup(k)=dup(k)+one-tfact*tfact*(one-dfact)
           dup(l)=dup(l)+one-tfact*tfact*(one-dfact)
        end if
     end do
  end do


! If requested, save select data for output to diagnostic file
  if(conv_diagsave)then
     ii=0
     nchar=1
     ioff0=27            !hliu rdiagbuf length definition
     nreal=ioff0
     if (lobsdiagsave) nreal=nreal+4*miter+1
     if (save_jacobian) then
       nnz   = 2                   ! number of non-zero elements in dH(x)/dx profile
       nind   = 1
       call new(dhx_dx, nnz, nind)
       nreal = nreal + size(dhx_dx)
     endif
     allocate(cdiagbuf(nobs),rdiagbuf(nreal,nobs))
     if(netcdf_diag) call init_netcdf_diag_
  end if

  scale=one
  rsig=float(nsig)
  mm1=mype+1

  call dtime_setup()
!---------------------------
!---------------------------
  do i=1,nobs
!---------------------------
!---------------------------
! Convert obs lats and lons to grid coordinates
     dtime=data(itime,i)
     call dtime_check(dtime, in_curbin, in_anybin)
     if(.not.in_anybin) cycle

     if(in_curbin) then
        dlat=data(ilat,i)
        dlon=data(ilon,i)
        dpres=data(ihgt,i)
 
        ikx=nint(data(ikxx,i))
     endif

!    Link observation to appropriate observation bin
     if (nobs_bins>1) then
        ibin = NINT( dtime/hr_obsbin ) + 1
     else
        ibin = 1
     endif
     IF (ibin<1.OR.ibin>nobs_bins) write(6,*)mype,'Error nobs_bins,ibin= ',nobs_bins,ibin

     if (luse_obsdiag) my_diagLL => odiagLL(ibin)

!    Link obs to diagnostics structure
     if (luse_obsdiag) then
        my_diag => obsdiagLList_nextNode(my_diagLL      ,&
                create = .not.lobsdiag_allocated        ,&
                   idv = is             ,&
                   iob = ioid(i)        ,&
                   ich = 1              ,&
                  elat = data(ilate,i)  ,&
                  elon = data(ilone,i)  ,&
                  luse = luse(i)        ,&
                 miter = miter          )

        if(.not.associated(my_diag)) call die(myname, &
                'obsdiagLList_nextNode(), create =', .not.lobsdiag_allocated)
     endif

     if(.not.in_curbin) cycle

! Save observation latitude.  This is needed when converting 
! geopotential to geometric height (hges --> zges below)
     slat=data(ilate,i)*deg2rad

! Interpolate log(surface pressure), model terrain, 
! and log(pres) at mid-layers to observation location.
     factw=data(iff10,i)
     if(sfcmod_gfs .or. sfcmod_mm5) then
        sfcr = data(isfcr,i)
        skint = data(iskint,i)
        isli = data(idomsfc,i)
        call comp_fact10(dlat,dlon,dtime,skint,sfcr,isli,mype,factw)
     end if

     call tintrp2a11(ges_ps,psges,dlat,dlon,dtime,hrdifsig,&
          mype,nfldsig)
     call tintrp2a1(ges_lnprsl,prsltmp,dlat,dlon,dtime,hrdifsig,&
          nsig,mype,nfldsig)
     call tintrp2a1(geop_hgtl,hges,dlat,dlon,dtime,hrdifsig,&
          nsig,mype,nfldsig)

     call tintrp2a11(ges_z,zsges,dlat,dlon,dtime,hrdifsig,&      ! jsw
          mype,nfldsig)                                    ! jsw
          dpres=dpres-zsges              !jsw need to adjust dpres by zsges


! Convert geopotential height at layer midpoints to geometric height using
! equations (17, 20, 23) in MJ Mahoney's note "A discussion of various
! measures of altitude" (2001).  Available on the web at
! http://mtp.jpl.nasa.gov/notes/altitude/altitude.html
!
! termg  = equation 17
! termr  = equation 21
! termrg = first term in the denominator of equation 23
! zges   = equation 23

     sin2  = sin(slat)*sin(slat)
     termg = grav_equator * &
           ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
     termr = semi_major_axis /(one + flattening + grav_ratio -  &
                 two*flattening*sin2)
     termrg = (termg/grav)*termr
     do k=1,nsig
        zges(k) = (termr*hges(k)) / (termrg-hges(k))  ! eq (23)
     end do

! Given observation height, (1) adjust 10 meter wind factor if
! necessary, (2) convert height to grid relative units, (3) compute
! compute observation pressure (for diagnostic purposes only), and
! (4) compute location of midpoint of first model layer above surface
! in grid relative units

! Adjust 10m wind factor if necessary.  Rarely do we have a
! lidar obs within 10 meters of the surface.  Almost always,
! the code below resets the 10m wind factor to 1.0 (i.e., no
! reduction in wind speed due to surface friction).

!    adjust wind near surface        jsw
  if (dpres<zges(1)) then            !jsw
      if(zges(1)>10)then
         term = (zges(1)-dpres)/(zges(1)-ten)
         term = min(max(term,zero),one)
         if(zges(1)<10) term=1
         factw = one-term+factw*term
      endif
    else
      factw=one
    endif


! Convert observation height (in dpres) from meters to grid relative
! units.  Save the observation height in zob for later use.
     zob = dpres
     call grdcrd1(dpres,zges,nsig,1)


!! hliu -------------------------------------------------
!! get top and bottom heights of L2B layers in grid relative
!! need to read out these heights directly from the data later
!
!     zobt= zob + 0.25*data(iatd,i)        ! top of L2B layers
!     zobb= zob - 0.25*data(iatd,i)        ! bottom of L2B layers

!     zobt2= zob + 0.5*data(iatd,i)        ! top of L2B layers
!     zobb2= zob - 0.5*data(iatd,i)        ! bottom of L2B layers

!      zobt0 = zobt2   ! save zobt in (m) for calculate shear
!      zobb0 = zobb2   ! save zobt in (m)

!     call grdcrd1(zobt,zges,nsig,1)
!     call grdcrd1(zobb,zges,nsig,1)

!     call grdcrd1(zobt2,zges,nsig,1)
!     call grdcrd1(zobb2,zges,nsig,1)
!!hliu --------------------------------------


! Set indices of model levels below (k1) and above (k2) observation.
! wm - updated so {k1,k2} are at min {1,2} and at max {nsig-1,nsig}
     k=dpres
     k1=min(max(1,k),nsig-1)
     k2=min(k1+1,nsig)
!    k1=max(1,k)         - old method
!    k2=min(k+1,nsig)    - old method

! Compute observation pressure (only used for diagnostics)
     dz       = zges(k2)-zges(k1)
     dlnp     = prsltmp(k2)-prsltmp(k1)
     pobl     = prsltmp(k1) + (dlnp/dz)*(zob-zges(k1))
     presw   = ten*exp(pobl)

! Determine location in terms of grid units for midpoint of
! first layer above surface
     sfcchk=log(psges)
     call grdcrd1(sfcchk,prsltmp,nsig,-1)

! Check to see if observation is below midpoint of first
! above surface layer.  If so, set rlow to that difference

     rlow=max(sfcchk-dpres,zero)

! Check to see if observation is above midpoint of layer
! at the top of the model.  If so, set rhgh to that difference.
     rhgh=max(dpres-r0_001-nsig,zero)

! Increment obs counter along with low and high obs counters
     if(luse(i))then
        awork(1)=awork(1)+one
        if(rhgh/=zero) awork(2)=awork(2)+one
        if(rlow/=zero) awork(3)=awork(3)+one
     end if

! Set initial obs error to that supplied in BUFR stream.
     error = data(ier,i)


!hliu-------------------------------------------------------------
! Use Hollingthworth error estimate of Aeolus winds 

 if (ictype(ikx)==48 ) then           ! Aeolus winds

    lat0 = data(ilate,i)              !deg
    if( lat0  >   25.0) then
       ireg = 1   !NH
    else if( lat0  <= -25.0) then
       ireg = 2   !SH
    else
       ireg = 3   !TR
    end if

      if( icsubtype(ikx)==20) then      ! Ray winds
           kk = minloc( abs(hght_ray-zob),1)
        error = rayerr(ireg, kk)
        delta = raydel(ireg, kk)        ! from HL estimate, errob/errbg

       if ( delta > 0.0) then
        delta = delta **2
        else
        delta = 4.0                     ! ray default
       endif

!  add 2021.01.29
       if( zob > 15000.0 .and. abs(lat0)>50.0 ) error = error * 2
!new   if( zob > 15000.0 .and. abs(lat0)>50.0 ) error = error * 2

      else                              ! mie winds
           kk = minloc( abs(hght_mie-zob),1)
        error = mieerr(ireg, kk)
        delta = miedel(ireg, kk)

       if ( delta > 0.0) then
        delta = delta **2
        else
        delta = 2.0                     ! mie default
       endif

!  inflate Mie winds with arge error in upper troposphere of SH 
!new   if( zob > 11000.0 ) error = error * 2.5
       if( zob > 12000.0 .and. lat0 < 0.0 ) error = error * 2

      endif

! inflate Aeolus error in polar areas for high density of Ray + Mie
       if( abs(lat0) > 80.0 ) error = error * 2

       if(error <= -99.0)  muse(i)=.false.     ! for rayleigh winds <= 2.5km
 endif
! hliu-----------------------


! Removed repe_dw, but retained the "+ one" for reproducibility
!hliu!  for ikx=100 or 101 - wm
!hliu     if (ictype(ikx)==100 .or. ictype(ikx)==101)error = error + one
!hliu! msq error change moved from read_lidar, wrapped to avoid changing 
!hliu!  ADM values
!hliu     if (ictype(ikx)==200 .or. ictype(ikx)==201) then 
!hliu        if (data(ier,i) > dmiss) then                  
!hliu           error = 3.0_r_kind                                
!hliu        else
!hliu           error = data(ier,i) / cos(data(ielva,i))
!hliu        endif
!hliu     endif    

     ratio_errors = error/abs(error + 1.0e6_r_kind*rhgh + r8*rlow)

!hliu   error, or data(ier,i) can not be zero ------------
     error = one/error

     if(dpres < zero .or. dpres > rsig)ratio_errors = zero
 
! Simulate dw wind from guess (forward model)
! First, interpolate u,v guess to observation location

     call tintrp31(ges_u,ugesindw,dlat,dlon,dpres,dtime,&
        hrdifsig,mype,nfldsig)
     call tintrp31(ges_v,vgesindw,dlat,dlon,dpres,dtime,&
        hrdifsig,mype,nfldsig)


!! hliu-----------------------------
!! hliu  simulation at top of L2B layers
!     call tintrp31(ges_u,ugesindwt,dlat,dlon,zobt,dtime,&
!        hrdifsig,mype,nfldsig)
!     call tintrp31(ges_v,vgesindwt,dlat,dlon,zobt,dtime,&
!        hrdifsig,mype,nfldsig)
!
!! hliu  simulation at bottom of L2B layers
!     call tintrp31(ges_u,ugesindwb,dlat,dlon,zobb,dtime,&
!        hrdifsig,mype,nfldsig)
!     call tintrp31(ges_v,vgesindwb,dlat,dlon,zobb,dtime,&
!        hrdifsig,mype,nfldsig)

!! hliu-----------------------------
!     call tintrp31(ges_u,ugesindwt2,dlat,dlon,zobt2,dtime,&
!        hrdifsig,mype,nfldsig)
!     call tintrp31(ges_v,vgesindwt2,dlat,dlon,zobt2,dtime,&
!        hrdifsig,mype,nfldsig)
 
!! hliu  simulation at bottom of L2B layers
!     call tintrp31(ges_u,ugesindwb2,dlat,dlon,zobb2,dtime,&
!        hrdifsig,mype,nfldsig)
!     call tintrp31(ges_v,vgesindwb2,dlat,dlon,zobb2,dtime,&
!        hrdifsig,mype,nfldsig)
!! hliu-----------------------------



! Next, convert wind components to line of sight value
! wm     if (nint(data(isubtype,i))==100.or.nint(data(isubtype,i))==101) then
!hliu --------------------------------------------------------------------
!    if (ictype(ikx)==100 .or. ictype(ikx)==101) then
!     KNMI  product  msq
     if (ictype(ikx)==100 .or. ictype(ikx)==101 .or. ictype(ikx)==48) then
!hliu --------------------------------------------------------------------

        cosazm  = -cos(data(iazm,i))  ! cos(azimuth)  ! mccarty msq 
        sinazm  = -sin(data(iazm,i))  ! sin(azimuth)  ! mccarty msq
     else
        cosazm  = cos(data(iazm,i))  ! cos(azimuth)
        sinazm  = sin(data(iazm,i))  ! sin(azimuth)
     endif

     dwwind=(ugesindw*sinazm+vgesindw*cosazm)*factw

!!hliu -----------------------------------------------------
!     dwwindt =(ugesindwt*sinazm + vgesindwt*cosazm)*factw    
!     dwwindb =(ugesindwb*sinazm + vgesindwb*cosazm)*factw    

!     dwwindt2=(ugesindwt2*sinazm + vgesindwt2*cosazm)*factw    
!     dwwindb2=(ugesindwb2*sinazm + vgesindwb2*cosazm)*factw    
!!hliu -----------------------------------------------------


     iz = max(1, min( int(dpres), nsig))
     delz = max(zero, min(dpres - float(iz), one))

     if (save_jacobian) then
        u_ind = getindex(svars3d, 'u')
        if (u_ind < 0) then
           print *, 'Error: no variable u in state vector. Exiting.'
           call stop2(1300)
        endif
        v_ind = getindex(svars3d, 'v')
        if (v_ind < 0) then
           print *, 'Error: no variable v in state vector. Exiting.'
           call stop2(1300)
        endif

        dhx_dx%st_ind(1)  = iz               + sum(levels(1:u_ind-1))
        dhx_dx%end_ind(1) = min(iz + 1,nsig) + sum(levels(1:u_ind-1))

        dhx_dx%val(1) = (one - delz) * sinazm * factw
        dhx_dx%val(2) = delz * sinazm * factw

        dhx_dx%st_ind(2)  = iz               + sum(levels(1:v_ind-1))
        dhx_dx%end_ind(2) = min(iz + 1,nsig) + sum(levels(1:v_ind-1))

        dhx_dx%val(3) = (one - delz) * cosazm * factw
        dhx_dx%val(4) = delz * cosazm * factw
     endif


!hliu ------------------------------------------------
! QCs for the simulations/observations pair associated
! with large vertical wind shear of GFS 

!   wind profile roughness (w_i+1 + w_i-1 - 2*w_i) /depth

!    wshear = ( 0.5*(dwwindt2 + dwwindb2) - dwwind) /( (zobt0-zobb0)**2 )    ! m/s/m/m

!    if( abs(wshear) > 5.0e-3 ) then
!      muse(i) = .false.
!      ratio_errors = zero
!      data(iuse,i) = 206
!    end if


     ddiff = data(ilob,i) - dwwind


!hliu--------------------------------------------------------
! Apply the TLS Bias corrections in bins of lat and levels
!------------------------------------------------------------
     dwobs = data(ilob,i)
     dwges = dwwind
     bc6 = 0.0

 if(ictype(ikx)==48) then                 ! for Aeolus
           lat0 = data(ilate,i)           !deg
            lnd = minloc( abs(latint - lat0 ), 1 )
         latdif = lat0 - latint(lnd)

          if(latdif > 0.0) then
           lnd2 = min(lnd+1, nlat)
           else
           lnd2 = max(lnd-1, 2)
          endif

           x1 = latint(lnd)
           x2 = latint(lnd2)
           x0 = lat0

!---------------------------------------------------------
! --- interpolate BC coeffients to Aeolus latitude 
!---------------------------------------------------------
  if( icsubtype(ikx)==20) then             ! Rayleigh

       kk= minloc( abs(hght_ray - zob), 1)
!-----------------------------------------------------
        if( data(iazm,i) > pi ) then       ! ascending >180deg or pi 
          coefray =  rayasc(:, kk, :)
         else                              ! descending
          coefray =  raydes(:, kk, :)
        endif

       if( lnd .ne. lnd2) then
           f1 = coefray(1, lnd)
           f2 = coefray(1, lnd2)
           beta0 = f1 + (f2-f1)* (x0-x1)/(x2-x1)

           f1 = coefray(2, lnd)
           f2 = coefray(2, lnd2)
           beta1 = f1 + (f2-f1)* (x0-x1)/(x2-x1)
        else
           beta0 = coefray(1, lnd)
           beta1 = coefray(2, lnd)
       endif

!-----------------------------------------
   else                       ! Mie cloudy
!-----------------------------------------
       kk= minloc( abs(hght_mie - zob), 1)
!-----------------------------------------
        if( data(iazm,i) > pi ) then       ! ascending >180deg or pi 
          coefmie =  mieasc(:, kk, :)
         else                              ! descending
          coefmie =  miedes(:, kk, :)
        endif

       if( lnd .ne. lnd2) then
           f1 = coefmie(1, lnd)
           f2 = coefmie(1, lnd2)
           beta0 = f1 + (f2-f1)* (x0-x1)/(x2-x1)

           f1 = coefmie(2, lnd)
           f2 = coefmie(2, lnd2)
           beta1 = f1 + (f2-f1)* (x0-x1)/(x2-x1)
        else
           beta0 = coefmie(1, lnd)
           beta1 = coefmie(2, lnd)
       endif
   endif            !if( icsubtype(ikx)==20)

!-------------------------------------------------
!  print*, 'delta= ', delta, beta0, beta1
!-----------------------------------------------------------
!  TLS regression of O over B,  O(y) = beta0 + beta1 * B(x)
!-----------------------------------------------------------
!   here y = OBS, x = ges GFS
!  yopt and xopt: mean estimates for y and x
!
!    dinv = y-beta0-beta1*x      ! innovation
!    xopt = x + beta1*dinv/(beta1^2+del)
!    yopt = beta0 + beta1*xopt
!         = beta0 + beta1*x +beta1^2 *di/(beta1^2 + del)
!         = y - di + beta1^2 * di/(beta1^2 +del)
!         = y - del*di/(beta1^2 + del)


          dinv = dwobs - beta0 - beta1 * dwges
        gesopt = dwges + beta1 * dinv/(beta1*beta1 + delta)

!  Derive bias correction --------------------
!   mean value estimate for mean of O-B vs B:
!        bc6 = beta0 + (beta1-1) * gesopt
! 
!   mean value estimate for mean of O-B vs O:
!       obsopt = dwobs - delta * dinv/(beta1*beta1 + delta)    or
        obsopt = beta0 + beta1*gesopt

     if( abs(beta1) > 0.0 ) then
         bc6 = beta0/beta1 + (beta1-1)/beta1 * obsopt
      else
         bc6 = 0.0
     endif

!------------------------------------------------------------
        ddiff = ddiff - bc6       ! subtract bias correction
!------------------------------------------------------------
 endif          !if(ictype(ikx)==48)
!hliu -------------------------------

!    Gross check using innovation normalized by error
     obserror = one/max(ratio_errors*error,  tiny_r_kind)
     obserrlm = max(cermin(ikx),min(cermax(ikx), obserror))
     residual = abs(ddiff)
     ratio    = residual/obserrlm

     if (ratio > cgross(ikx) .or. ratio_errors < tiny_r_kind) then
        if(luse(i))awork(4) = awork(4) + one
        error = zero
        ratio_errors=zero
     else
        ratio_errors=ratio_errors/sqrt(dup(i))
     endif

     if (ratio_errors*error <= tiny_r_kind) muse(i) = .false.
     if (nobskeep>0 .and. luse_obsdiag) call obsdiagNode_get(my_diag, jiter=nobskeep, muse=muse(i))
 
!----hliu-----------------------------------------------------
! Discard L2B Rayleigh winds below 2km (with large errors)
!---------------------------------------------------------
!  OBS rejected:  muse = false and error == 0
!  OBS monitored: muse = false and error \= 0
!---------------------------------------------------------

! QC Rayleigh clear:
   if(icsubtype(ikx)==20 ) then
       if( data(ihgt,i) <= 2000.0 .or. data(ihgt,i) > 19000.0) then
        muse(i)=.false.
        ratio_errors = zero
        data(iuse,i) = 200
        bc6 = 0.0
       endif
   endif

! QC Mie winds: skip few Mie winds with large error
   if(icsubtype(ikx)==11 ) then
       if( zob > 16000.0 ) then
        muse(i)=.false.
        ratio_errors = zero
        data(iuse,i) = 200
        bc6 = 0.0
       endif
   endif

! skip negative height of L2B data
     if( data(ihgt,i) < 0.0) then
       muse(i)=.false.
       ratio_errors = zero
       data(iuse,i) = 200
     endif
!hliu-------------------------------------------------------

!hliu-------------------------------------------
! Apply ESA recommended QCs
! Discard L2B winds with large L2B uncertainty
!-----------------------------------------------
! Rayleigh clear:
     if(icsubtype(ikx)==20 .and. data(ier,i)>12.0) then
       muse(i)=.false.
       ratio_errors = zero
       data(iuse,i) = 201
        bc6 = 0.0
     endif
! Mie cloudy:
     if(icsubtype(ikx)==11 .and. data(ier,i)>5.0) then
       muse(i)=.false.
       ratio_errors = zero
       data(iuse,i) = 201
        bc6 = 0.0
     endif

!hliu -------------------------------------------------------------------------
 ! MR> Winds within 20 hPa of the model orography were discarded because it was
 ! noted that ground return winds contaminating the L2B dataset (this is
 ! improved in the next version of the L2B processor (v3.10)).

     if ( abs(psges*r10 - presw) < 20.0_r_kind) then
       muse(i)=.false.
       ratio_errors = zero
       data(iuse,i) = 203
     endif

!hliu ------------------------------------------------------
! Discard winds with short horizontal accumulation lengths:
! Rayleigh winds < 60 km and Mie winds < 5 km.
!  if( mod(i, 100)==0) print*, 'erraccu= ', data(inls,i)  ! (m)

     if ( icsubtype(ikx) == 20 .and. data(inls,i) < 60000.0_r_kind) then
       muse(i)=.false.
       ratio_errors = zero
       data(iuse,i) = 208
     endif
     if ( icsubtype(ikx) == 11 .and. data(inls,i) < 5000.0_r_kind) then
       muse(i)=.false.
       ratio_errors = zero
       data(iuse,i) = 208
     endif

  if( muse(i) == .false. ) bc6 = 0.0     ! hliu for output


!    Compute penalty terms
     val   = error*ddiff
     if(luse(i))then
        val2     = val*val
        exp_arg  = -half*val2
        rat_err2 = ratio_errors**2
        if (cvar_pg(ikx) > tiny_r_kind .and. error > tiny_r_kind) then
           arg  = exp(exp_arg)
           wnotgross= one-cvar_pg(ikx)
           cg_dw=cvar_b(ikx)
           wgross = cg_term*cvar_pg(ikx)/(cg_dw*wnotgross)
           term = log((arg+wgross)/(one+wgross))
           wgt  = one-wgross/(arg+wgross)
           rwgt = wgt/wgtlim
        else
           term = exp_arg
           wgt  = wgtlim
           rwgt = wgt/wgtlim
        endif
        valqc = -two*rat_err2*term


!       Accumulate statistics for obs belonging to this task
        if(muse(i))then
           if(rwgt < one) awork(21) = awork(21)+one
           jsig = dpres
           jsig=max(1,min(jsig,nsig))
           awork(jsig+6*nsig+100)=awork(jsig+6*nsig+100)+val2*rat_err2
           awork(jsig+5*nsig+100)=awork(jsig+5*nsig+100)+one
           awork(jsig+3*nsig+100)=awork(jsig+3*nsig+100)+valqc   
        endif


! Loop over pressure level groupings and obs to accumulate statistics
! as a function of observation type.

        do k = 1,npres_print
           if(presw > ptop(k) .and. presw <= pbot(k)) then   
              ress =scale*ddiff
              ressw=ress*ress
              val2 =val*val
              rat_err2 = ratio_errors**2
              nn=1
              if (.not. muse(i)) then
                 nn=2                                            ! hliu rejected OBS
!hliu            if(ratio_errors*error >=tiny_r_kind) nn=3       ! monitored OBS
                 if(ratio_errors*error >=tiny_r_kind .and. data(incls,i) ==zero) nn=3       ! monitored OBS
              end if

              bwork(k,ikx,1,nn) = bwork(k,ikx,1,nn)+one             ! count
              bwork(k,ikx,2,nn) = bwork(k,ikx,2,nn)+ddiff           ! bias    
              bwork(k,ikx,3,nn) = bwork(k,ikx,3,nn)+ressw           ! (o-g)**2
              bwork(k,ikx,4,nn) = bwork(k,ikx,4,nn)+val2*rat_err2   ! penalty
              bwork(k,ikx,5,nn) = bwork(k,ikx,5,nn)+valqc           ! nonlin qc penalty
           end if
  
        end do
     end if

     if (luse_obsdiag) then
        call obsdiagNode_set(my_diag,wgtjo=(error*ratio_errors)**2, &
                jiter=jiter,muse=muse(i),nldepart=ddiff)
     endif

!    If obs is "acceptable", load array with obs info for use
!    in inner loop minimization (int* and stp* routines)
     if (.not. last .and. muse(i)) then
 
        allocate(my_head)
        call dwNode_appendto(my_head,dwhead(ibin))

        my_head%idv = is
        my_head%iob = ioid(i)
        my_head%elat= data(ilate,i)
        my_head%elon= data(ilone,i)

!       Set (i,j,k) indices of guess gridpoint that bound obs location
        my_head%dlev = dpres
        my_head%factw= factw
        call get_ijk(mm1,dlat,dlon,dpres,my_head%ij,my_head%wij)

        do j=1,8
           my_head%wij(j)=factw*my_head%wij(j)  
        end do                                                 

        my_head%res    = ddiff
        my_head%err2   = error**2
        my_head%raterr2=ratio_errors**2    
        my_head%time   = dtime
        my_head%b      = cvar_b(ikx)
        my_head%pg     = cvar_pg(ikx)
        my_head%cosazm = cosazm                  ! v factor
        my_head%sinazm = sinazm                  ! u factor
        my_head%luse   = luse(i)

        if(luse_obsdiag) then
           call obsdiagNode_assert(my_diag, my_head%idv,my_head%iob,1,myname,'my_diag:my_head')
           my_head%diags => my_diag
        endif
        my_head => null()
     endif

! Save select output for diagnostic file  
     if(conv_diagsave)then
        ii=ii+1
        rstation_id = data(id,i)
        err_input   = data(ier,i)
        err_adjst   = data(ier,i)
        if (ratio_errors*error>tiny_r_kind) then
           err_final = one/(ratio_errors*error)
        else
           err_final = huge_single
        endif

        errinv_input = huge_single
        errinv_adjst = huge_single
        errinv_final = huge_single
        if (err_input>tiny_r_kind) errinv_input=one/err_input
        if (err_adjst>tiny_r_kind) errinv_adjst=one/err_adjst
        if (err_final>tiny_r_kind) errinv_final=one/err_final

        if (binary_diag) call contents_binary_diag_(my_diag)
        if (netcdf_diag) call contents_netcdf_diag_(my_diag)

     end if

  end do

! Release memory of local guess arrays
  call final_vars_

! Write information to diagnostic file
  if(conv_diagsave) then
    if(netcdf_diag) call nc_diag_write
    if(binary_diag .and. ii>0)then
       write(7)' dw',nchar,nreal,ii,mype,ioff0
       write(7)cdiagbuf(1:ii),rdiagbuf(:,1:ii)
       deallocate(cdiagbuf,rdiagbuf)
     end if
  end if
  
! End of routine

  return
  contains

  subroutine check_vars_ (proceed)
  logical,intent(inout) :: proceed
  integer(i_kind) ivar, istatus
! Check to see if required guess fields are available
  call gsi_metguess_get ('var::ps', ivar, istatus )
  proceed=ivar>0
  call gsi_metguess_get ('var::z' , ivar, istatus )
  proceed=proceed.and.ivar>0
  call gsi_metguess_get ('var::u', ivar, istatus )
  proceed=proceed.and.ivar>0
  call gsi_metguess_get ('var::v', ivar, istatus )
  proceed=proceed.and.ivar>0
  end subroutine check_vars_ 

  subroutine init_vars_

  real(r_kind),dimension(:,:  ),pointer:: rank2=>NULL()
  real(r_kind),dimension(:,:,:),pointer:: rank3=>NULL()
  character(len=5) :: varname
  integer(i_kind) ifld, istatus

! If require guess vars available, extract from bundle ...
  if(size(gsi_metguess_bundle)==nfldsig) then
!    get ps ...
     varname='ps'
     call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank2,istatus)
     if (istatus==0) then
         if(allocated(ges_ps))then
            write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
            call stop2(999)
         endif
         allocate(ges_ps(size(rank2,1),size(rank2,2),nfldsig))
         ges_ps(:,:,1)=rank2
         do ifld=2,nfldsig
            call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank2,istatus)
            ges_ps(:,:,ifld)=rank2
         enddo
     else
         write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
         call stop2(999)
     endif
!    get z ...
     varname='z'
     call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank2,istatus)
     if (istatus==0) then
         if(allocated(ges_z))then
            write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
            call stop2(999)
         endif
         allocate(ges_z(size(rank2,1),size(rank2,2),nfldsig))
         ges_z(:,:,1)=rank2
         do ifld=2,nfldsig
            call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank2,istatus)
            ges_z(:,:,ifld)=rank2
         enddo
     else
         write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
         call stop2(999)
     endif
!    get u ...
     varname='u'
     call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank3,istatus)
     if (istatus==0) then
         if(allocated(ges_u))then
            write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
            call stop2(999)
         endif
         allocate(ges_u(size(rank3,1),size(rank3,2),size(rank3,3),nfldsig))
         ges_u(:,:,:,1)=rank3
         do ifld=2,nfldsig
            call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank3,istatus)
            ges_u(:,:,:,ifld)=rank3
         enddo
     else
         write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
         call stop2(999)
     endif
!    get v ...
     varname='v'
     call gsi_bundlegetpointer(gsi_metguess_bundle(1),trim(varname),rank3,istatus)
     if (istatus==0) then
         if(allocated(ges_v))then
            write(6,*) trim(myname), ': ', trim(varname), ' already incorrectly alloc '
            call stop2(999)
         endif
         allocate(ges_v(size(rank3,1),size(rank3,2),size(rank3,3),nfldsig))
         ges_v(:,:,:,1)=rank3
         do ifld=2,nfldsig
            call gsi_bundlegetpointer(gsi_metguess_bundle(ifld),trim(varname),rank3,istatus)
            ges_v(:,:,:,ifld)=rank3
         enddo
     else
         write(6,*) trim(myname),': ', trim(varname), ' not found in met bundle, ier= ',istatus
         call stop2(999)
     endif
  else
     write(6,*) trim(myname), ': inconsistent vector sizes (nfldsig,size(metguess_bundle) ',&
                 nfldsig,size(gsi_metguess_bundle)
     call stop2(999)
  endif
  end subroutine init_vars_

  subroutine init_netcdf_diag_
  character(len=80) string
  character(len=128) diag_conv_file
  integer(i_kind) ncd_fileid,ncd_nobs
  logical append_diag
  logical,parameter::verbose=.false. 
     write(string,900) jiter
900  format('conv_dw_',i2.2,'.nc4')
     diag_conv_file=trim(dirname) // trim(string)

     inquire(file=diag_conv_file, exist=append_diag)

     if (append_diag) then
        call nc_diag_read_init(diag_conv_file,ncd_fileid)
        ncd_nobs = nc_diag_read_get_dim(ncd_fileid,'nobs')
        call nc_diag_read_close(diag_conv_file)

        if (ncd_nobs > 0) then
           if(verbose) print *,'file ' // trim(diag_conv_file) // ' exists.  Appending.  nobs,mype=',ncd_nobs,mype
        else
           if(verbose) print *,'file ' // trim(diag_conv_file) // ' exists but contains no obs.  Not appending. nobs,mype=',ncd_nobs,mype
           append_diag = .false. ! if there are no obs in existing file, then do not try to append
        endif
     end if

     call nc_diag_init(diag_conv_file, append=append_diag)

     if (.not. append_diag) then ! don't write headers on append - the module will break?
        call nc_diag_header("date_time",ianldate )
        call nc_diag_header("Number_of_state_vars", nsdim          )
     endif
  end subroutine init_netcdf_diag_

!=====================================================
  subroutine contents_binary_diag_(odiag)
  type(obs_diag),pointer,intent(in):: odiag
        cdiagbuf(ii)    = station_id         ! station id

        rdiagbuf(1,ii)  = ictype(ikx)        ! observation type
        rdiagbuf(2,ii)  = icsubtype(ikx)     ! observation subtype

        rdiagbuf(3,ii)  = data(ilate,i)      ! observation latitude (degrees)
        rdiagbuf(4,ii)  = data(ilone,i)      ! observation longitude (degrees)
        rdiagbuf(5,ii)  = rmiss_single       ! station elevation (meters)
        rdiagbuf(6,ii)  = presw              ! observation pressure (hPa)
        rdiagbuf(7,ii)  = data(ihgt,i)       ! observation height (meters)
        rdiagbuf(8,ii)  = dtime-time_offset  ! obs time (hours relative to analysis time)

        rdiagbuf(9,ii)  = rmiss_single       ! input prepbufr qc or event mark
        rdiagbuf(10,ii) = rmiss_single       ! setup qc or event mark
        rdiagbuf(11,ii) = data(iuse,i)       ! read_prepbufr data usage flag
        if(muse(i)) then
           rdiagbuf(12,ii) = one             ! analysis usage flag (1=use, -1=not used)
        else
           rdiagbuf(12,ii) = -one
        endif

        rdiagbuf(13,ii) = rwgt                 ! nonlinear qc relative weight
        rdiagbuf(14,ii) = errinv_input         ! prepbufr inverse obs error
        rdiagbuf(15,ii) = errinv_adjst         ! read_prepbufr inverse obs error
        rdiagbuf(16,ii) = errinv_final         ! final inverse observation error

        rdiagbuf(17,ii) = data(ilob,i)         ! observation
        rdiagbuf(18,ii) = ddiff                ! obs-ges used in analysis (with bias correction)
        rdiagbuf(19,ii) = dwwind               ! ges 
!       rdiagbuf(19,ii) = data(ilob,i)-dwwind  ! obs-ges w/o bias correction (future slot)
 
        rdiagbuf(20,ii) = factw                ! 10m wind reduction factor
        rdiagbuf(21,ii) = data(ielva,i)*rad2deg! elevation angle (degrees)
        rdiagbuf(22,ii) = data(iazm,i)*rad2deg ! bearing or azimuth (degrees)
        rdiagbuf(23,ii) = data(inls,i)         ! number of laser shots
        rdiagbuf(24,ii) = data(incls,i)        ! number of cloud laser shots
        rdiagbuf(25,ii) = data(iatd,i)         ! atmospheric depth
! hliu 2020.04.02
!       rdiagbuf(26,ii) = data(ilob,i)         ! line of sight component of wind orig.
        rdiagbuf(26,ii) = data(ier,i)          ! L2B uncertainty
! hliu 2020.04.02

        rdiagbuf(27,ii) = 1.e+10_r_single      ! ges ensemble spread (filled in by EnKF)

        ioff=ioff0
        if (lobsdiagsave) then
           do jj=1,miter 
              ioff=ioff+1 
              if (odiag%muse(jj)) then
                 rdiagbuf(ioff,ii) = one
              else
                 rdiagbuf(ioff,ii) = -one
              endif
           enddo
           do jj=1,miter+1
              ioff=ioff+1
              rdiagbuf(ioff,ii) = odiag%nldepart(jj)
           enddo
           do jj=1,miter
              ioff=ioff+1
              rdiagbuf(ioff,ii) = odiag%tldepart(jj)
           enddo
           do jj=1,miter
              ioff=ioff+1
              rdiagbuf(ioff,ii) = odiag%obssen(jj)
           enddo
        endif
        if (save_jacobian) then
           call writearray(dhx_dx, rdiagbuf(ioff+1:nreal,ii))
           ioff = ioff + size(dhx_dx)
        endif

  end subroutine contents_binary_diag_


!=====================================================
  subroutine contents_netcdf_diag_(odiag)
!=====================================================
  type(obs_diag),pointer,intent(in):: odiag
! Observation class
  character(7),parameter     :: obsclass = '     dw'
  real(r_single),parameter::     missing = -9.99e9_r_single
  real(r_kind),dimension(miter) :: obsdiag_iuse
           call nc_diag_metadata("Station_ID",              station_id             )
           call nc_diag_metadata("Observation_Class",       obsclass               )
           call nc_diag_metadata("Observation_Type",        ictype(ikx)            )
           call nc_diag_metadata("Observation_Subtype",     icsubtype(ikx)         )
           call nc_diag_metadata("Latitude",                sngl(data(ilate,i))    )
           call nc_diag_metadata("Longitude",               sngl(data(ilone,i))    )
           call nc_diag_metadata("Station_Elevation",       missing                )
           call nc_diag_metadata("Pressure",                sngl(presw)            )
           call nc_diag_metadata("Height",                  sngl(data(ihgt,i))     )
           call nc_diag_metadata("Time",                    sngl(dtime-time_offset))
           call nc_diag_metadata("Prep_QC_Mark",            missing                )
           call nc_diag_metadata("Prep_Use_Flag",           sngl(data(iuse,i))     )
!          call nc_diag_metadata("Nonlinear_QC_Var_Jb",     var_jb                 )
           call nc_diag_metadata("Nonlinear_QC_Rel_Wgt",    sngl(rwgt)             )                 
           if(muse(i)) then
              call nc_diag_metadata("Analysis_Use_Flag",    sngl(one)              )
           else
              call nc_diag_metadata("Analysis_Use_Flag",    sngl(-one)             )              
           endif

           call nc_diag_metadata("Errinv_Input",            sngl(errinv_input)     )
           call nc_diag_metadata("Errinv_Adjust",           sngl(errinv_adjst)     )
           call nc_diag_metadata("Errinv_Final",            sngl(errinv_final)     )

           call nc_diag_metadata("Observation",                   sngl(data(ilob,i)))
           call nc_diag_metadata("Obs_Minus_Forecast_adjusted",   sngl(ddiff)      )
           call nc_diag_metadata("Obs_Minus_Forecast_unadjusted", sngl(data(ilob,i)-dwwind))

!hliu 2020.06.21
!          call nc_diag_metadata("Bias_Correction",               sngl(bc6))
           call nc_diag_metadata("Forecast", sngl(dwwind))

!          call nc_diag_metadata("Layer_Vert_Shear", sngl(wshear*1.0e6)) !m/s/km/km
!          call nc_diag_metadata("Layer_top", sngl(dwwindt))
!          call nc_diag_metadata("Layer_bom", sngl(dwwindb))
!          call nc_diag_metadata("Layer_top2", sngl(dwwindt2))
!          call nc_diag_metadata("Layer_bom2", sngl(dwwindb2))

!          call nc_diag_metadata("Layer_depth", sngl(data(iatd,i)*0.001))

           call nc_diag_metadata("Azim_Angle_deg", sngl(data(iazm,i)*rad2deg))
           call nc_diag_metadata("Hori_Accumulation_km", sngl(0.001*data(inls,i)))
           call nc_diag_metadata("L2B_uncertainty", sngl(data(ier,i)))

!_RT_NC4_TODO
!_RT    rdiagbuf(20,ii) = factw                ! 10m wind reduction factor
!_RT    rdiagbuf(21,ii) = data(ielva,i)*rad2deg! elevation angle (degrees)
!_RT    rdiagbuf(22,ii) = data(iazm,i)*rad2deg ! bearing or azimuth (degrees)
!_RT    rdiagbuf(23,ii) = data(inls,i)         ! number of laser shots
!_RT    rdiagbuf(24,ii) = data(incls,i)        ! number of cloud laser shots
!_RT    rdiagbuf(25,ii) = data(iatd,i)         ! atmospheric depth
!_RT    rdiagbuf(26,ii) = data(ilob,i)         ! line of sight component of wind orig.
 
           if (lobsdiagsave) then
              do jj=1,miter
                 if (odiag%muse(jj)) then
                       obsdiag_iuse(jj) =  one
                 else
                       obsdiag_iuse(jj) = -one
                 endif
              enddo
   
              call nc_diag_data2d("ObsDiagSave_iuse",     obsdiag_iuse                             )
              call nc_diag_data2d("ObsDiagSave_nldepart", odiag%nldepart )
              call nc_diag_data2d("ObsDiagSave_tldepart", odiag%tldepart )
              call nc_diag_data2d("ObsDiagSave_obssen",   odiag%obssen   )             
           endif

           if (save_jacobian) then
              call fullarray(dhx_dx, dhx_dx_array)
              call nc_diag_data2d("Observation_Operator_Jacobian", dhx_dx_array)
           endif
   
  end subroutine contents_netcdf_diag_

  subroutine final_vars_
    if(allocated(ges_v )) deallocate(ges_v )
    if(allocated(ges_u )) deallocate(ges_u )
    if(allocated(ges_z )) deallocate(ges_z )
    if(allocated(ges_ps)) deallocate(ges_ps)
  end subroutine final_vars_


!hliu -------------------------------------
subroutine read_bias_correctionTLSonline

   rayasc = 0.0; raydes = 0.0;  mieasc = 0.0; miedes = 0.0

   open(961, file='BC_TLS_Ray.asc',form='formatted')
   open(962, file='BC_TLS_Ray.des',form='formatted')
   open(963, file='BC_TLS_Mie.asc',form='formatted')
   open(964, file='BC_TLS_Mie.des',form='formatted')

 do i=1, 2
   do k=1, nray
    read(961, '(19f8.3)') (rayasc(i, k, n), n=1, nlat)
    read(962, '(19f8.3)') (raydes(i, k, n), n=1, nlat)
   enddo

   do k=1, nmie
    read(963, '(19f8.3)') (mieasc(i, k, n), n=1, nlat)
    read(964, '(19f8.3)') (miedes(i, k, n), n=1, nlat)
   enddo
 enddo

    close(961);close(962);close(963);close(964)

!   if( ianldate <= 2019082406 ) then   ! for Aug. 2019

!----------- rad in HL error estimates -------------
   open(971, file='Ray_errNH.txt',form='formatted')
   open(972, file='Ray_errSH.txt',form='formatted')
   open(973, file='Ray_errTR.txt',form='formatted')

   open(981, file='Mie_errNH.txt',form='formatted')
   open(982, file='Mie_errSH.txt',form='formatted')
   open(983, file='Mie_errTR.txt',form='formatted')
   
    do k=1, nray
     read(971, '(3f8.2)')  rayerr(1,k), raydel(1,k)    ! NH
     read(972, '(3f8.2)')  rayerr(2,k), raydel(2,k)    ! SH
     read(973, '(3f8.2)')  rayerr(3,k), raydel(3,k)    ! TR
    end do

    do k=1, nmie
     read(981, '(3f8.2)')  mieerr(1,k), miedel(1,k)    ! NH
     read(982, '(3f8.2)')  mieerr(2,k), miedel(2,k)    ! SH
     read(983, '(3f8.2)')  mieerr(3,k), miedel(3,k)    ! TR
    end do

    close(971);close(972);close(973);close(981);close(982);close(983)
!
! print*, 'rayasc= ', rayasc
! print*, 'raydes= ', raydes
! print*, 'mieasc= ', mieasc
! print*, 'miedes= ', miedes
 
!print*, 'rayerr= ', rayerr
!print*, 'mieerr= ', mieerr
!
!print*, 'raydel= ', raydel
!print*, 'miedel= ', miedel

 end subroutine read_bias_correctionTLSonline


end subroutine setupdw
end module dw_setup

