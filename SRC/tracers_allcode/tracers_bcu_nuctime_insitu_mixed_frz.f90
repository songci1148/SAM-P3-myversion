module tracers


! This module serves as a template for adding tracer transport in the model. The tracers can be 
! chemical tracers, or bin microphysics drop/ice categories, etc. 
! The number of tracers is set by the parameter ntracers which is set in domain.f90.
! Also, the logical flag dotracers should be set to .true. in namelist (default is .false.).
! The model will transport the tracers around automatically (advection and SGS diffusion).
! The user must supply the initialization in the subroutine tracers_init() in this module.
! By default, the surface flux of all tracers is zero. Nonzero values can be set in tracers_flux().
! The local sinks/sources of tracers should be supplied in tracers_physics().



 use grid
 use vars, only: w, tabs, qv, qcl, qci, qpl, qpi !BG for tracers
 use params, only: epsv,dotracers
 use microphysics, only: micro_proc_rates
 implicit none

 real tracer  (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, 0:ntracers) 
 real fluxbtr (nx, ny, 0:ntracers) ! surface flux of tracers
 real fluxttr (nx, ny, 0:ntracers) ! top boundary flux of tracers
 real trwle(nz,0:ntracers)  ! resolved vertical flux 
 real trwsb(nz,0:ntracers)  ! SGS vertical flux 
 real tradv(nz,0:ntracers)  ! tendency due to vertical advection 
 real trdiff(nz,0:ntracers)  ! tendency due to vertical diffusion 
 real trphys(nz,0:ntracers)  ! tendency due to physics 
 real trlsadv(nz,0:ntracers)  ! tendency due to large-scale vertical advection 
 character *4 tracername(0:ntracers)
 character *40 tracerlongname(0:ntracers)
 character *10 tracerunits(0:ntracers)

CONTAINS

 subroutine tracers_init()

  integer k,ntr
  character *2 ntrchar
  character ntrchar_single
  integer, external :: lenstr

 tracer = 0.
 fluxbtr = 0.
 fluxttr = 0.

! Add your initialization code here. Default is to set to 0 in setdata.f90.

 if(nrestart.eq.0) then

   ! Add your initialization code here.
   tracer = 0.

 end if

! Specify te tracers' default names:

   ! Default names are TRACER01, TRACER02, etc:

 do ntr = 1,ntracers
   if(ntracers.lt.10) then
     write(ntrchar_single,'(i1)') ntr
     tracername(ntr) = 'TR'//ntrchar_single
   else
     write(ntrchar,'(i2)') ntr
     do k=1,3-lenstr(ntrchar)-1
       ntrchar(k:k)='0'
     end do
     tracername(ntr) = 'TR'//ntrchar(1:2)
   end if
   tracerlongname(ntr) = TRIM(tracername(ntr))
   tracerunits(ntr) = '[TR]'
 end do

 ntr =0
 ntr = ntr+1
 tracerlongname(ntr) = 'Mixed-phase ice nucleation time tracer (set to 1 where mixed-phase ice nucleation occurs, elsewhere decays w/tau=2 hours)'
 ntr = ntr+1
 tracerlongname(ntr) = 'In-situ ice nucleation time tracer (set to 1 where in-situ ice nucleation occurs, elsewhere decays w/tau=2 hours)'
 ntr = ntr+1
 tracerlongname(ntr) = 'Security ice nucleation time tracer (set to 1 where security ice nucleation occurs, elsewhere decays w/tau=2 hours)'
 ntr = ntr+1
 !tracerlongname(ntr) = 'Total ice nucleation time tracer (set to 1 where all ice nucleation occurs, elsewhere decays w/tau=2 hours)'
 !ntr = ntr+1
 tracerlongname(ntr) = 'Buoyant cloudy updraft (BCU) tracer (set to 1 in BCU, elsewhere decays w/tau=2 hours)'

 end subroutine tracers_init



 subroutine tracers_flux()

! Set surface and top fluxes of tracers. Default is 0 set in setdata.f90

 end subroutine tracers_flux


 
 subroutine tracers_physics()
 ! add here a call to a subroutine that does something to tracers besides advection and diffusion.
 ! The transport is done automatically. 
  integer :: i, j, k, ntr
  real :: fac_decay
  real :: tmpz(nzm), tv0(nzm)

  trphys = 0. ! Default tendency due to physics. You code should compute this to output statistics.

  do ntr = 1,ntracers
    do k = 1,nzm
      trphys(k,ntr) = SUM(tracer(1:nx,1:ny,k,ntr))
    end do

    select case (ntr)
    case (1)
      ! Tracer 1: Mixphase freezing tracer (==1 if freezing in mix phase, decays as 1/(2hours) elsewhere)

      fac_decay = 1. - dtn/7200.

      do k = 1,nzm
        ! apply decay everywhere (relax to zero with a decay timescale of one day)
        tracer(:,:,k,ntr) = fac_decay*tracer(:,:,k,ntr)
      !aaaa
        ! overwrite with a value of one where active (w>1, qc+qi>1e-5)
        do j = 1,ny
          do i = 1,nx
            if (( micro_proc_rates(i,j,k,1) + micro_proc_rates(i,j,k,4) + micro_proc_rates(i,j,k,5)  ) .gt.1.e-8)  then !mix phase frz, no secondary ice
             !sort of deposition freezing + immersion freezing of cloud droplets and rain
              tracer(i,j,k,ntr) = 1.
            end if
          end do
        end do
      end do

     case (2)
    ! Tracer 2: In situ nucleation tracer

      fac_decay = 1. - dtn/7200.

      do k = 1,nzm
        ! apply decay everywhere (relax to zero with a decay timescale of one day)
        tracer(:,:,k,ntr) = fac_decay*tracer(:,:,k,ntr)

        ! overwrite with a value of one where active
        do j = 1,ny
          do i = 1,nx
            !Mohler freezign (if enabled) + Liu&Penner freezing
            if (( micro_proc_rates(i,j,k,2) + micro_proc_rates(i,j,k,3)) .gt.1.e-8)  then
              tracer(i,j,k,ntr) = 1.
            end if
          end do
        end do
      end do

    case (3)
      ! Tracer 3: Security statement nucleation tracer

      fac_decay = 1. - dtn/7200.

      do k = 1,nzm
        ! apply decay everywhere (relax to zero with a decay timescale of one day)
        tracer(:,:,k,ntr) = fac_decay*tracer(:,:,k,ntr)

        ! overwrite with a value of one where active
        do j = 1,ny
          do i = 1,nx
            if ( (micro_proc_rates(i,j,k,10)) .gt.1.e-8)  then
            !security freezing: if droplets present at T=-37, freeze them
              tracer(i,j,k,ntr) = 1.
            end if
          end do
        end do
      end do


    case (4)
      ! Tracer 5: Active tracer (==1 for w>1, qc+qi>1e-5 kg/kg and positively buoyant; decays as 1/(2 hours) elsewhere)
      !   Similar to Romps (2010, JAS, 10.1175/2010JAS3371.1), though his tracer was
      !   set to zero elsewhere and did not account for positive buoyancy.
      fac_decay = 1. - dtn/7200.

      do k = 1,nzm
        tv0(k) = 0.
        do j = 1,ny
          do i = 1,nx
            tv0(k) = tv0(k) + tabs(i,j,k)*( 1. + epsv*qv(i,j,k) - qcl(i,j,k) - qci(i,j,k) - qpl(i,j,k) - qpi(i,j,k) )
          end do
        end do
      end do
      tv0(:) = tv0(:) / float(nx_gl*ny_gl) ! normalize over whole domain
      if(dompi) then
        ! sum across processors, if necessary
        call task_sum_real(tv0,tmpz,nzm)
        tv0(1:nzm) = tmpz(1:nzm)
      end if

      do k = 1,nzm
        ! apply decay everywhere (relax to zero with a decay timescale of one day)
        tracer(:,:,k,ntr) = fac_decay*tracer(:,:,k,ntr)

        ! overwrite with a value of one where active (w>1, qc+qi>1e-5)
        do j = 1,ny
          do i = 1,nx
            if( ( (w(i,j,k) + w(i,j,k+1)).gt.2. ) .AND. (qcl(i,j,k)+qci(i,j,k).gt.1.e-5) ) then
              ! if w>1 and qcl+qci>1e-5, check for positive buoyancy
              if( tabs(i,j,k)*( 1. + epsv*qv(i,j,k) - qcl(i,j,k) - qci(i,j,k) - qpl(i,j,k) - qpi(i,j,k) ) .gt. tv0(k) ) then
                tracer(i,j,k,ntr) = 1.
              end if
            end if
          end do
        end do
      end do

    case default

      write(*,*) 'No physics tendency specified for tracer number ', ntr
      call task_abort()

    end select

    do k = 1,nzm
      ! compute physics tendency for this level/tracer
      trphys(k,ntr) = SUM(tracer(1:nx,1:ny,k,ntr)) - trphys(k,ntr)
    end do

  end do ! do ntr = 1,ntracers

 end subroutine tracers_physics



 subroutine tracers_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

! Initialize the list of tracers statistics variables written in statistics.f90

   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount
   integer ntr


   do ntr=1,ntracers

     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))
     deflist(count) = trim(tracerlongname(ntr))
     !deflist(count) = trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'FLX'
     deflist(count) = 'Total flux of '//trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))//' kg/m2/s'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'FLXS'
     deflist(count) = 'SGS flux of '//trim(tracername(ntr))
     unitlist(count) = trim(tracerunits(ntr))//' kg/m2/s'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'ADV'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to vertical advection')
     unitlist(count) = trim(tracerunits(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'DIFF'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to vertical SGS transport')
     unitlist(count) = trim(tracername(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
     count = count + 1
     trcount = trcount + 1
     namelist(count) = trim(tracername(ntr))//'PHYS'
     deflist(count) = 'Tendency of '//trim(tracername(ntr)//'due to physics')
     unitlist(count) = trim(tracername(ntr))//'/day'
     status(count) = 1
     average_type(count) = 0
   end do

 end subroutine tracers_hbuf_init

end module tracers
