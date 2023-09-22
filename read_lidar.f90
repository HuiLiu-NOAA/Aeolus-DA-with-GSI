subroutine read_lidar(nread,ndata,nodata,infile,obstype,lunout,twind,sis,nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_lidar                   read doppler lidar winds
!   prgmmr: yang             org: np20                date: 1998-05-15
!
! abstract:  This routine reads doppler lidar wind files.  
!
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! program history log:
!   1998-05-15  yang, weiyu
!   1999-08-24  derber, j., treadon, r., yang, w., first frozen mpp version
!   2004-06-16  treadon - update documentation
!   2004-07-29  treadon - add only to module use, add intent in/out
!   2005-08-02  derber - modify to use convinfo file
!   2005-09-08  derber - modify to use input group time window
!   2005-10-11  treadon - change convinfo read to free format
!   2005-10-17  treadon - add grid and earth relative obs location to output file
!   2005-10-18  treadon - remove array obs_load and call to sumload
!   2005-10-26  treadon - add routine tag to convinfo printout
!   2006-02-03  derber  - add new obs control
!   2006-02-08  derber  - modify to use new convinfo module
!   2006-02-24  derber  - modify to take advantage of convinfo module
!   2006-07-27  msq/terry - removed cosine factor for line of sight winds and obs err
!   2007-03-01  tremolet - measure time from beginning of assimilation window
!   2008-04-18  safford - rm unused vars
!   2010-08-01  woollen -  change bufr table (denoted as jsw)
!   2010-08-01  woollen -  change kx to ikx (bug) (denoted as jsw)
!   2010-09-01  masutani -  remove statements related to old cos(lat) correction !msq
!   2010-10-06  masutani -- use ikx, ikx=999 for missing type Bufrtable was updated
!   2010-11-05  mccarty/woollen -  add level to dwld
!   2010-11-30  masutani - add kx to cdata_all(21), change maxdat to 21  (denoted msq)
!   2011-04-15  mccarty - change maxdat back to 20, kx in setupdw taken from ictype
!   2011-05-05  mccarty - cleaned up unnecessary print statement
!   2011-05-26  mccarty - remove dwlerror logic (moved to setupdw) 
!   2011-08-01  lueken  - added module use deter_sfc_mod
!   2013-01-26  parrish - change from grdcrd to grdcrd1 (to allow successful debug compile on WCOSS)
!   2015-02-23  Rancic/Thomas - add l4densvar to time window logical
!   2015-10-01  guo     - consolidate use of ob location (in deg
!   2019-02-22  mccarty - A number of Aeolus-centric updates, including handling of subtypes for 
!                         different wind retrieval methodology
!   2019-07-12  huiliu  - various updates of Aeolus L2B HLOS wind
!   2020-01-28  huiliu  - one update of Aeolus wind on the 250m altitude adjustment
!
!   input argument list:
!     infile   - unit from which to read BUFR data
!     obstype  - observation type to process
!     lunout   - unit to which to write data for further processing
!     twind    - input group time window (hours)
!
!   output argument list:
!     nread    - number of doppler lidar wind observations read
!     ndata    - number of doppler lidar wind profiles retained for further processing
!     nodata   - number of doppler lidar wind observations retained for further processing
!     sis      - satellite/instrument/sensor indicator
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,r_double,i_kind
  use gridmod, only: nlat,nlon,regional,tll2xy,rlats,rlons
  use convinfo, only: nconvtype,ctwind, &  !added mccarty
      ncmiter,ncgroup,ncnumgrp,icuse,ictype,ioctype,icsubtype  !mccarty
  use constants, only: deg2rad,zero,r60inv ! check the usage   msq
  use obsmod, only: iadate,offtime_data
  use gsi_4dvar, only: l4dvar,l4densvar,time_4dvar,winlen,iwinbgn
  use deter_sfc_mod, only: deter_sfc2
  use mpimod, only: npe
  implicit none

! Declare passed variables
  character(len=*),intent(in   ) :: obstype,infile
  character(len=20),intent(in  ) :: sis
  integer(i_kind) ,intent(in   ) :: lunout
  integer(i_kind) ,intent(inout) :: nread,ndata,nodata
  integer(i_kind),dimension(npe),intent(inout) :: nobs
  real(r_kind)    ,intent(in   ) :: twind

! Declare local parameters
  integer(i_kind) :: maxobs
  integer(i_kind),parameter:: maxdat=27     

  real(r_double),parameter:: r360 = 360.0_r_double
  real(r_kind),parameter::   r90  =  90.0_r_kind

! Declare local variables
  logical dwl,outside

  character(80) hdstr,dwstr 
  character(10) date
  character(8) subset, station_id, horiz_seq_str, vert_seq_str
!hliu
  character(len=16),allocatable,dimension(:):: dw_ctype

  integer(i_kind) lunin,i,kx,ilat,ikx,idomsfc
  integer(i_kind) jdate,ihh,idd,idate,iret,im,iy,k,levs
  integer(i_kind) nmrecs,ilon,nreal,nchanl,nmsgmax
  integer(i_kind) nhdr, ndwl
!hliu
  integer(i_kind):: ndws
  integer(i_kind),allocatable,dimension(:):: dw_itype,dw_stype, dw_ikx

  real(r_kind) time,usage,dlat,dlon,dlat_earth,dlon_earth
  real(r_kind) dlat_earth_deg,dlon_earth_deg
  real(r_kind) hloswind,sfcr,tsavg,ff10,toff,t4dv,layer_depth ! msq changed to hloswind
  real(r_kind),allocatable,dimension(:,:):: cdata_all

  real(r_double) rstation_id
  real(r_double) rkx                        !msq
  real(r_double),dimension(:),allocatable  :: hdr
  real(r_double),dimension(:,:),allocatable :: dwld 
  real(r_double),dimension(:),allocatable  :: aeolusd
  integer(i_kind),parameter                 :: n_horiz_seq = 4
  real(r_double),dimension(n_horiz_seq)     :: horiz_seq
  integer(i_kind),parameter                 :: n_vert_seq  = 5
  real(r_double),dimension(n_vert_seq)      :: vert_seq
  integer(i_kind)                           :: subtype


  integer(i_kind) idate5(5),minobs,minan,nmind
  real(r_kind) time_correction

  integer(i_kind):: ilev        ! mccarty
  equivalence(rstation_id,station_id)

  data lunin / 10 /

!**************************************************************************
! Initialize variables
  nmrecs=0
  nreal=maxdat
  nchanl=0
  ilon=2
  ilat=3

!hliu------------------------------------------------
! Check convinfo file to see if requesting to process DW data
  ikx = 0
  do i=1,nconvtype
      if ( trim(sis)==trim(ioctype(i))) ikx=ikx+1
  end do

! If no dw data requested to be process, exit routine
  if(ikx==0) then
   write(6,*)'READ LIDAR: CONVINFO DOES NOT INCLUDE ANY ',trim(sis),' DATA'
   return
  end if

! Allocate and load arrays to contain DW types only.
  ndws=ikx
  allocate(dw_ctype(ndws), dw_itype(ndws), dw_stype(ndws), dw_ikx(ndws))
  ikx=0
  do i=1,nconvtype
      if (trim(sis)==trim(ioctype(i))) then
        ikx=ikx+1
        dw_ctype(ikx)=ioctype(i)
        dw_itype(ikx)=ictype(i)
        dw_stype(ikx)=icsubtype(i)
        dw_ikx(ikx)  =i
     endif
  end do
!hliu -------------------

  call getcount_bufr(trim(infile),nmsgmax,maxobs)
  allocate(cdata_all(maxdat,maxobs))
!print*, 'maxobs= ', maxobs

! Open, then read date from bufr data
  open(lunin,file=trim(infile),form='unformatted')
  call openbf(lunin,'IN',lunin)
  call datelen(10)
  call readmg(lunin,subset,idate,iret)
  if(iret/=0) then
      print*,' failed to dw read data from ',lunin    ! msq
      call closbf(lunin)
      return
  endif

! Time offset
  call time_4dvar(idate,toff)

! If date in lidar file does not agree with analysis date, 
! print message and stop program execution.
  write(date,'( i10)') idate
  read (date,'(i4,3i2)') iy,im,idd,ihh

  if(offtime_data) then
!       in time correction for observations to account for analysis
!                    time being different from obs file time.
!hliu  this time_correction is to move all times of OBS to analysis time, thish
!hliu  is only good for 3DVar, not needed for 4dvar ------

     write(date,'( i10)') idate
     read (date,'(i4,3i2)') iy,im,idd,ihh
     idate5(1)=iy
     idate5(2)=im
     idate5(3)=idd
     idate5(4)=ihh
     idate5(5)=0
     call w3fs21(idate5,minobs)  ! obs ref time in minutes relative to historic date
     idate5(1)=iadate(1)
     idate5(2)=iadate(2)
     idate5(3)=iadate(3)
     idate5(4)=iadate(4)
     idate5(5)=0
     call w3fs21(idate5,minan)   ! analysis ref time in minutes relative to historic date

!    add obs reference time, then subtract analysis time to get obs time relative to analysis

     time_correction=float(minobs-minan)*r60inv

  else
     time_correction=zero
  end if

  write(6,*)'READ_LIDAR dw or aeolus: time offset is ',toff,' hours.'

! Big loop over bufr file	

  obsloop: do
     call readsb(lunin,iret) 
     if(iret/=0) then
        call readmg(lunin,subset,jdate,iret)
        if(iret/=0) exit obsloop
        cycle obsloop
     end if
     nmrecs=nmrecs+1

     if (subset=='DWLDAT') then
         call read_dwldat_
         nread = nread + 1
     else if (subset=='FN023000') then 
         call read_aeolus_
         nread = nread + 1
     else 
         cycle obsloop
     endif
! End of bufr read loop
  end do obsloop

! Write observations to scratch file
  call count_obs(ndata,maxdat,ilat,ilon,cdata_all,nobs)
  write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
  write(lunout) ((cdata_all(k,i),k=1,maxdat),i=1,ndata)

!!print*, 'ndata in read_lidar= ', ndata, maxobs

! Close unit to bufr file
  deallocate(cdata_all)

  if (allocated(hdr))     deallocate(hdr)
  if (allocated(dwld))    deallocate(dwld)
  if (allocated(aeolusd)) deallocate(aeolusd)

  call closbf(lunin)

! End of routine
  return

  contains




!-------------------------------------------------------------
!-------------------------------------------------------------
  subroutine read_dwldat_
!-------------------------------------------------------------
!-------------------------------------------------------------

     if (.not. allocated(hdr))  allocate(hdr(5))
     if (.not. allocated(dwld)) allocate(dwld(8,24))

!    Extract type, date, and location information
! 
     call ufbint(lunin,rkx,1,1,iret,'TYP')           !msq
     kx=nint(rkx)                                    !msq
!  data dwstr  /'HEIT ELEV BEARAZ NOLS NOLC ADPL LOSC LOSCU'/  !msq  jsw
!  data dwstr2  /'ADWL ELEV BORA NOLS NOLC ADPL LOSC SDLE'/ !msq  used for KNMI data prepared by GMAO

     if (kx==100.or.kx==101) then
!        ADM data 
         hdstr = 'SID CLON CLAT DHR TYP'
         dwstr = 'ADWL ELEV BORA NOLS NOLC ADPL LOSC SDLE'
     else if (kx==201.or.kx==202) then
!        GWOS data
         hdstr = 'SID XOB YOB DHR TYP'
         dwstr = 'HEIT ELEV BEARAZ NOLS NOLC ADPL LOSC LOSCU'
     else
!        undefined dwl data
         hdstr = 'SID XOB YOB DHR TYP'
         dwstr = 'HEIT ELEV BEARAZ NOLS NOLC ADPL LOSC LOSCU'
         kx=999
     endif
     call ufbint(lunin,hdr,5,1,iret,hdstr)

     
     ikx=0
     do i=1,nconvtype
        if(trim(obstype) == trim(ioctype(i)) .and. kx == ictype(i))ikx = i
     end do
!    Determine if this is doppler wind lidar report
     dwl= (ikx /= 0) .and. (subset=='DWLDAT')  ! jsw chenge kx to ikx (bug)
     if(.not. dwl) then
          return
     endif


     t4dv = toff + hdr(4)
     if (l4dvar.or.l4densvar) then
        if (t4dv<zero .OR. t4dv>winlen) return
     else   
        time=hdr(4) + time_correction
        if (abs(time) > ctwind(ikx) .or. abs(time) > twind) return
     endif

     hdr(2)=mod(hdr(2),r360)  ! msq
     if (hdr(2) < zero)  hdr(2)=hdr(2)+r360

     dlat_earth_deg = hdr(3)
     dlon_earth_deg = hdr(2)
  write(6,*)'READ_LIDAR dw1d/nasa: ', dlat_earth_deg, dlon_earth_deg, dwld(1,ilev)

     dlat_earth = hdr(3) * deg2rad
     dlon_earth = hdr(2) * deg2rad

     if(regional)then
        call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
        if (outside)  return
     else
        dlat = dlat_earth
        dlon = dlon_earth
        call grdcrd1(dlat,rlats,nlat,1)
        call grdcrd1(dlon,rlons,nlon,1)
     endif

     call ufbint(lunin,dwld,8,24,levs,dwstr) !mccarty, msq

     do ilev=1,levs !mccarty, jsw

!       If wind data, extract observation.
        nodata=min(nodata+1,maxobs)
        ndata=min(ndata+1,maxobs)
        usage = zero

! hliu for dw1 -------- for minitoring --------------
        if(icuse(ikx) < 0)  usage=100._r_kind
! hliu -------- for minitoring --------------

        if(ncnumgrp(ikx) > 0 )then                     ! cross validation on
           if(mod(ndata,ncnumgrp(ikx))== ncgroup(ikx)-1)usage=ncmiter(ikx)
        end if

        station_id=subset

        hloswind=dwld(7,ilev)/(cos(dwld(2,ilev)*deg2rad))    ! obs wind (line of sight component)
        call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr)

        cdata_all(1,ndata)=ikx                    ! obs type
        cdata_all(2,ndata)=dlon                   ! grid relative longitude
        cdata_all(3,ndata)=dlat                   ! grid relative latitude
        cdata_all(4,ndata)=t4dv                   ! obs time (analyis relative hour)
        cdata_all(5,ndata)=dwld(1,ilev)           ! obs height (altitude) (m), NASA DW1D
        cdata_all(6,ndata)=dwld(2,ilev)*deg2rad   ! elevation angle (radians)
        cdata_all(7,ndata)=dwld(3,ilev)*deg2rad   ! bearing or azimuth (radians)
        cdata_all(8,ndata)=dwld(4,ilev)           ! number of laser shots
        cdata_all(9,ndata)=dwld(5,ilev)           ! number of cloud laser shots
        cdata_all(10,ndata)=dwld(6,ilev)          ! obs layer depth
        cdata_all(11,ndata)=hloswind               ! obs wind (line of sight component) msq
        cdata_all(12,ndata)=dwld(8,ilev)          ! standard deviation (obs error) msq
        cdata_all(13,ndata)=rstation_id           ! station id
        cdata_all(14,ndata)=hdr(1)                ! satellite id
        cdata_all(15,ndata)=usage                 ! usage parameter
        cdata_all(16,ndata)=idomsfc+0.001_r_kind  ! dominate surface type
        cdata_all(17,ndata)=tsavg                 ! skin temperature      
        cdata_all(18,ndata)=ff10                  ! 10 meter wind factor  
        cdata_all(19,ndata)=sfcr                  ! surface roughness     
! -------------------   NASA dw1d dataset only --------------
        cdata_all(20,ndata)=dlon_earth_deg        ! earth relative longitude (degrees)
        cdata_all(21,ndata)=dlat_earth_deg        ! earth relative latitude (degrees)
        cdata_all(22,ndata)=zero                  ! reference Pressure
        cdata_all(23,ndata)=zero                  ! retrieval derivative of wind w.r.t. Pressure
        cdata_all(24,ndata)=zero                  ! reference Temperature
        cdata_all(25,ndata)=zero                  ! retrieval derivative of wind w.r.t. Temperature
        cdata_all(26,ndata)=zero                  ! retrieval Backscatter
        cdata_all(27,ndata)=zero                  ! retrieval derivative of wind w.r.t. Backscatter
     enddo   ! ilev

     return
!-------------------------------------------------------------
!-------------------------------------------------------------
  end subroutine read_dwldat_
!-------------------------------------------------------------



!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine read_aeolus_
!----------------------------------------------------------------------
!----------------------------------------------------------------------
     hdstr = 'SAID SIID YEAR MNTH DAYS HOUR MINU SECW RCVCH LL2BCT'
     nhdr  = 10
! hliu
!    dwstr = 'HLSW HLSWEE CONFLG PRES TMDBST BKSTR DWPRS DWTMP DWBR'
!    ndwl  = 9
     dwstr = 'HLSW HLSWEE CONFLG PRES TMDBST BKSTR DWPRS DWTMP DWBR HOIL'
     ndwl  = 10

     if (.not. allocated(hdr))     allocate(hdr(nhdr))     ! allocate according to # of variables in hdstr
     if (.not. allocated(aeolusd)) allocate(aeolusd(ndwl)) ! aalocate according to # of variables in dwstr
     !read header

     call ufbint(lunin,hdr,nhdr,1,iret,hdstr)

!          Check obs time
     idate5(1) = nint(hdr(3)) ! year
     idate5(2) = nint(hdr(4)) ! month
     idate5(3) = nint(hdr(5)) ! day
     idate5(4) = nint(hdr(6)) ! hour
     idate5(5) = nint(hdr(7)) ! minute

     if( idate5(1) < 1900 .or. idate5(1) > 3000 .or. &
         idate5(2) < 1    .or. idate5(2) >   12 .or. &
         idate5(3) < 1    .or. idate5(3) >   31 .or. &
         idate5(4) <0     .or. idate5(4) >   24 .or. &
         idate5(5) <0     .or. idate5(5) >   60 )then

        write(6,*)'READ_LIDAR aeolus:  ### ERROR IN READING AEOLUS BUFR DATA:', &
        ' STRANGE OBS TIME (YMDHM):', idate5(1:5)
        return
     endif

!    Retrieve obs time
     call w3fs21(idate5,nmind)
!  add in seconds      ! determined outside of this routine

     if (l4dvar.or.l4densvar) then
!       t4dv = (real(nmind-iwinbgn,r_kind) + real(hdr(8),r_kind)*r60inv)*r60inv + time_correction
! hliu  (a bug for aeolus, hdr(6)= hour)
        t4dv = (real(nmind-iwinbgn,r_kind) + real(hdr(8),r_kind)*r60inv)*r60inv 
!       t4dv is defined next for 3dvar option
     endif
!hliu

!    determine slot in convinfo, including setting type to satid and subtype
!    For type, WMO Satellite ID (SAID)

     kx = nint(hdr(1))

!    For subtype:
!        in bufr:
!           Channel (RCVCH):  0==Mie, 1==Rayleigh
!           Classification Type (LL2BCT):  0==Clear, 1==Cloudy
!        for subtype:
!           (RCVCH+1)*10 + (LL2BCT)
!          ...or...
!           subtype:  10==Mie,clear, 11==Mie,cloudy, 20==Rayleigh,Clear, 21==Rayleigh,Cloudy

     subtype = (nint(hdr(9)) + 1) * 10 + nint(hdr(10))

     ikx=0
!hliu -------------------------------------
! do i=1,nconvtype
!  if(trim(obstype) == trim(ioctype(i)) .and. kx == ictype(i) .and. subtype == icsubtype(i)) ikx = i
! end do

! print*,'satid= ',obstype,kx,subtype,dw_ctype(i),dw_itype(i),dw_stype(i)

  floop: do i=1,ndws
   if(trim(obstype) == trim(dw_ctype(i)) .and. kx == dw_itype(i) .and. &
            subtype == dw_stype(i)) then
     ikx=dw_ikx(i)
     exit floop
   endif
  end do floop

!-------- skip the other subtypes not interested -----
    if (ikx==0)  then
!!    print*, ' subtype skipped', subtype
     return
    endif
! hliu ----------------------------------

! check time window in subset
     if (l4dvar.or.l4densvar) then
        if (t4dv<zero .OR. t4dv>winlen) return
     else
!       time=hdr(4) + time_correction     !(a bug, hdr(4)=month for Aeolus, hliu)
        time = float(nmind - minan)*r60inv
        t4dv = time             ! for save in cdata_all for later 3dvar use
        if (abs(time) > ctwind(ikx) .or. abs(time) > twind) return
     endif

!    Read geolocation sequences (WM)
!     - the bufr table that I developed can access the horizontal and vertical information
!        as sequences.  This was the most reasonable appraoch that I came up with to the 
!        repeated sequences for the horizontal beginning, end, and centroid of the ob (HBEG,
!        HEND, and HCENT, respectively) and the vertical top, bottom, and centroid of the ob
!        (VTOP, VBOT, and VCENT, respectively)
!        Note - for first implementation, I am treating the ob as a point valid at the centroid
!
!        Order of variables in horizontal sequence: CRDSIG CLATH CLONH TISE
!        Order of variables in vertical sequence:   CRDSIG HEITH BEARAZ ELEV SATRG

     call ufbseq(lunin,horiz_seq,n_horiz_seq,1,iret,horiz_seq_str)

!    determine layer depth
!    First - read top of layer sequence, save height
     vert_seq_str = 'VTOP'
     call ufbseq(lunin,vert_seq,n_vert_seq,1,iret,vert_seq_str)
     layer_depth = vert_seq(2)

!    Second - read bottom of layer sequence, subtract bottom height from top height
     vert_seq_str = 'VBOT'
     call ufbseq(lunin,vert_seq,n_vert_seq,1,iret,vert_seq_str)
     layer_depth = layer_depth - vert_seq(2)

!    read vertical weighted centroid information
     vert_seq_str = 'VCENT'
     call ufbseq(lunin,vert_seq,n_vert_seq,1,iret,vert_seq_str)

!    read horizontal weighted centroid information
!    read vertical weighted centroid information
     horiz_seq_str = 'HCENT'
     call ufbseq(lunin,horiz_seq,n_horiz_seq,1,iret,horiz_seq_str)

     call ufbint(lunin,aeolusd,ndwl,1,iret,dwstr)

     if (abs(aeolusd(1)) > 1000.0_r_kind) return ! check if observation is realistic 
                                                 !  - 1000 ms-1 chosen somewhat arbitrarily, it appears 
                                                 !    unrealistic obs are reported as 1.0e11.  If the
                                                 !    magnitude is larger than threshhold; cycle loop on 
                                                 !    unrealistic observation

!    Do lat/lon handling
     dlat_earth = horiz_seq(2)
     dlon_earth = horiz_seq(3)

!  hliu ------------
!  write(6,*)'lll= ', dlat_earth, dlon_earth, vert_seq(2)

     if( abs(dlat_earth) > r90  .or. abs(dlon_earth) > r360 .or. &
        (abs(dlat_earth) == r90 .and. dlon_earth /= zero) )then
        write(6,*)'READ_LIDAR aeolus:  ### ERROR IN READING AEOLUS BUFR DATA:', &
           ' STRANGE OBS POINT (LAT,LON):', dlat, dlon
        return
     endif

     dlon_earth = mod(dlon_earth,r360)  ! msq
     if (dlon_earth < zero)  dlon_earth=dlon_earth+r360

     dlat_earth_deg = dlat_earth
     dlon_earth_deg = dlon_earth
     dlat_earth = dlat_earth * deg2rad
     dlon_earth = dlon_earth * deg2rad

     if(regional)then
        call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
        if (outside)  return
     else
        dlat = dlat_earth
        dlon = dlon_earth
        call grdcrd1(dlat,rlats,nlat,1)
        call grdcrd1(dlon,rlons,nlon,1)
     endif

     station_id=subset

!    determine surface parameters from lat/lon
     call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr)

!    set usage flag:
!       in bufr, conflg == 0: valid
!                       == 1: invalid
!       here, to setupdw:
!                usage == 0: use
!                usage == 100: don't use - off in convinfo
!                usage == 101: don't use - invalid in bufr
     usage = zero

     if(icuse(ikx) < 0)usage=100._r_kind
     if(aeolusd(3) > 0)usage=101._r_kind        ! confidence flag=1, invalid
! print*, 'uuu= ', ikx,  usage , icuse(ikx), aeolusd(3) 

!     hdstr = 'SAID SIID YEAR MNTH DAYS HOUR MINU SECW RCVCH LL2BCT'
!     dwstr = 'HLSW HLSWEE CONFLG PRES TMDBST BKSTR DWPRS DWTMP DWBR'
!     vert_seq = 'CRDSIG HEITH BEARAZ ELEV SATRG'
!     horiz_seq = 'CRDSIG CLATH CLONH TISE'

     nodata=min(nodata+1,maxobs)
     ndata=min(ndata+1,maxobs)

! hliu ikx is the order of DW data type in convinfo table (180, 181 so far)
     cdata_all(1,ndata)=ikx                    ! obs type
     cdata_all(2,ndata)=dlon                   ! grid relative longitude
     cdata_all(3,ndata)=dlat                   ! grid relative latitude
     cdata_all(4,ndata)=t4dv                   ! obs time (analyis relative hour)
     cdata_all(5,ndata)=vert_seq(2)            ! obs height (altitude) (m)

     cdata_all(6,ndata)=vert_seq(4)*deg2rad    ! elevation angle (radians)
     cdata_all(7,ndata)=vert_seq(3)*deg2rad    ! bearing or azimuth (radians)

! hliu  
     cdata_all(8,ndata)=aeolusd(10)            ! accumulation length (m), hliu
!    cdata_all(8,ndata)=zero                   ! number of laser shots - ZERO FOR NOW

     cdata_all(9,ndata)=zero                   ! number of cloud laser shots - ZERO FOR NOW
     cdata_all(10,ndata)=layer_depth           ! obs layer depth
     cdata_all(11,ndata)=aeolusd(1)            ! obs wind (line of sight component) msq
     cdata_all(12,ndata)=aeolusd(2)            ! standard deviation (obs error) msq
     cdata_all(13,ndata)=rstation_id           ! station id
     cdata_all(14,ndata)=hdr(1)                ! satellite id

     cdata_all(15,ndata)=usage                 ! usage parameter         
     cdata_all(16,ndata)=idomsfc+0.001_r_kind  ! dominate surface type
     cdata_all(17,ndata)=tsavg                 ! skin temperature      
     cdata_all(18,ndata)=ff10                  ! 10 meter wind factor  
     cdata_all(19,ndata)=sfcr                  ! surface roughness     
     cdata_all(20,ndata)=dlon_earth_deg        ! earth relative longitude (degrees)
     cdata_all(21,ndata)=dlat_earth_deg        ! earth relative latitude (degrees)
     cdata_all(22,ndata)=aeolusd(4)            ! reference Pressure (ECMWF ges)
     cdata_all(23,ndata)=aeolusd(7)            ! reference derivative of wind w.r.t. Pressure
     cdata_all(24,ndata)=aeolusd(5)            ! reference Temperature (ECMWF ges)
     cdata_all(25,ndata)=aeolusd(8)            ! reference derivative of wind w.r.t. Temperature
     cdata_all(26,ndata)=aeolusd(6)            ! retrieval Backscatter
     cdata_all(27,ndata)=aeolusd(9)            ! retrieval derivative of wind w.r.t. Backscatter

     return
     
  end subroutine read_aeolus_

end subroutine read_lidar
