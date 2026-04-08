subroutine setperturb

!  Random noise

use vars
use params
use microphysics, only: micro_field, index_water_vapor, iqit, inci, ReffIce_P3, IceMassMixingRatio_P3, reffi, nCat_ice_P3
use sgs, only: setperturb_sgs

implicit none

integer i,j,k,ptype,it,jt
real rrr,ranf_
real xxx,yyy,zzz
! parameters
real :: iwc_target_kgm3
real :: reff_top_um, reff_base_um
real :: nL_base, nL_top
real :: zbase, ztop, radius
integer :: ii

! --- variables used by case(30) and case(31) (must be declared up here)
real xc, yc, rxy, noise
real sum_noise, mean_noise, fracz, reff_um, nL, invrho, qice_mmr, n_mix, Tlocal
integer npts

call ranset_(3*rank)

ptype = perturb_type

call setperturb_sgs(ptype)  ! set sgs fields

select case (ptype)

  case(0)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(k.le.5) then
            t(i,j,k)=t(i,j,k)+0.02*rrr*(6-k)
         endif
       end do
      end do
     end do

  case(1)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.6.e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
         endif
       end do
      end do
     end do

  case(2) ! warm bubble

     if(masterproc) then
       print*, 'initialize with warm bubble:'
       print*, 'bubble_x0=',bubble_x0
       print*, 'bubble_y0=',bubble_y0
       print*, 'bubble_z0=',bubble_z0
       print*, 'bubble_radius_hor=',bubble_radius_hor
       print*, 'bubble_radius_ver=',bubble_radius_ver
       print*, 'bubble_dtemp=',bubble_dtemp
       print*, 'bubble_dq=',bubble_dq
     end if

     call task_rank_to_index(rank,it,jt)
     do k=1,nzm
       zzz = z(k)
       do j=1,ny
         yyy = dy*(j+jt)
         do i=1,nx
          xxx = dx*(i+it)
           if((xxx-bubble_x0)**2+YES3D*(yyy-bubble_y0)**2.lt.bubble_radius_hor**2 &
            .and.(zzz-bubble_z0)**2.lt.bubble_radius_ver**2) then
              rrr = cos(pi/2.*(xxx-bubble_x0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(yyy-bubble_y0)/bubble_radius_hor)**2 &
               *cos(pi/2.*(zzz-bubble_z0)/bubble_radius_ver)**2
              t(i,j,k) = t(i,j,k) + bubble_dtemp*rrr
              micro_field(i,j,k,index_water_vapor) = &
                  micro_field(i,j,k,index_water_vapor) + bubble_dq*rrr
           end if
         end do
       end do
     end do

  case(3)   ! gcss wg1 smoke-cloud case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.0.5e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
         endif
       end do
      end do
     end do

  case(4)  ! gcss wg1 arm case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(z(k).le.200.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr*(1.-z(k)/200.)
         endif
       end do
      end do
     end do

  case(5)  ! gcss wg1 BOMEX case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(z(k).le.1600.) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+0.025e-3*rrr
         endif
       end do
      end do
     end do

  case(6)  ! GCSS Lagragngian ASTEX

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.6.e-3) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)+2.5e-5*rrr
         endif
       end do
      end do
     end do

  case(22) !bloss: Try to make a general perturbation for boundary layer cloud simulations.
           ! Add noise everywhere that the water mass mixing ratio is more than half
           !   the value at the surface.  The noise has amplitude of 0.1K and 2% of the initial q0(k).

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         rrr=1.-2.*ranf_()
         if(q0(k).gt.0.5*q0(1)) then
            t(i,j,k)=t(i,j,k)+0.1*rrr
            micro_field(i,j,k,index_water_vapor)= &
                 (1. + 0.02*rrr)*micro_field(i,j,k,index_water_vapor)
         endif
       end do
      end do
     end do

  case(30)
    !============================================================
    ! Isolated anvil cloud perturbation (SAM-consistent)
    !
    ! Adds weak random temperature noise (±0.01 K) ONLY in a
    ! prescribed anvil region:
    !   base = 8 km, top = 13 km, radius = 30 km,
    ! centered in the 256-km square domain.
    !
    ! NOTE: This does NOT initialize ice condensate. In your SAM build,
    ! ice mass/number indices are not exposed here, so the cloud must
    ! already exist (typically via restart).
    !============================================================

    ! Domain center (global coordinates)
    xc = 0.5 * nx * dx
    yc = 0.5 * ny * dy

    call task_rank_to_index(rank,it,jt)

    do k = 1, nzm
      zzz = z(k)

      if (zzz .ge. 8000.0 .and. zzz .le. 13000.0) then

        do j = 1, ny
          yyy = dy * (j + jt)

          do i = 1, nx
            xxx = dx * (i + it)

            rxy = sqrt( (xxx - xc)**2 + YES3D*(yyy - yc)**2 )
            if (rxy .le. 30000.0) then
            ! Weak temperature perturbation (±0.01 K)
            noise = 0.01 * (2.0*ranf_() - 1.0)
            t(i,j,k) = t(i,j,k) + noise

            ! VERY weak water vapor perturbation (safe, SAM-consistent)
            micro_field(i,j,k,index_water_vapor) = &
                  micro_field(i,j,k,index_water_vapor) * (1.0 + 1.0e-6*noise)
            endif
          end do
        end do
      endif
    end do

  case(31)
    !============================================================
    ! Isolated anvil cloud initialization (Case 31)
    !
    ! Initializes an optically thick anvil cloud with:
    !  - Base = 8 km, Top = 13 km, Radius = 30 km (centered domain)
    !  - Ice water content (IWC) = 0.3 g/m3 (uniform target)
    !  - Effective radius: 20 um at top -> 40 um at base (linear)
    !  - Ice number concentration: 900 -> 3600 L^-1 (base->top),
    !      enforce >3000 L^-1 above 10.5 km
    !  - Random potential temperature perturbations ~ ±0.01 K (zero mean)
    !  - Attempts a simple buoyant compensation so cloud is near-neutral
    !============================================================

    if(masterproc) then
      print*, 'Initialize isolated anvil cloud (case 31): IWC=0.3 g/m3, reff 20->40um, N 900->3600 L^-1'
    end if

    xc = 0.5 * nx_gl * dx
    yc = 0.5 * ny_gl * dy

    iwc_target_kgm3 = 0.3e-3 ! 0.3 g/m3 -> kg/m3
    reff_top_um = 20.0
    reff_base_um = 40.0
    nL_base = 900.0
    nL_top = 3600.0
    zbase = 8000.0
    ztop = 13000.0
    radius = 30000.0
    ii = 1

    ! First pass: set microphysics fields and accumulate statistics for mean noise
    sum_noise = 0.0
    npts = 0

    call task_rank_to_index(rank,it,jt)

    do k = 1, nzm
      zzz = z(k)
      if (zzz .ge. zbase .and. zzz .le. ztop) then
        ! fractional height from base (0) to top (1)
        fracz = (zzz - zbase) / (ztop - zbase)

        ! determine target reff and number concentration at this level
        reff_um = reff_base_um + (reff_top_um - reff_base_um) * fracz
        nL = nL_base + (nL_top - nL_base) * fracz
        if (zzz .ge. 10500.0) nL = max(nL, 3000.0)

        do j = 1, ny
          yyy = dy * (j + jt)
          do i = 1, nx
            xxx = dx * (i + it)
            rxy = sqrt( (xxx - xc)**2 + YES3D*(yyy - yc)**2 )
            if (rxy .le. radius) then

              ! convert to mixing ratios (#/kg or kg/kg)
              invrho = 1.0 / rho(k)

              qice_mmr = iwc_target_kgm3 * invrho  ! kg/kg
              n_mix = (nL / 1000.0) * invrho       ! convert L^-1 -> m^-3 then /rho -> #/kg

              ! set total ice mass and number in micro_field (category ii)
              micro_field(i,j,k, iqit(ii)) = qice_mmr
              micro_field(i,j,k, inci(ii)) = n_mix

              ! set P3 arrays (per-category)
              IceMassMixingRatio_P3(i,j,k,ii) = qice_mmr
              ReffIce_P3(i,j,k,ii) = reff_um
              reffi(i,j,k) = reff_um

              ! random small temperature perturbation (±0.01 K)
              noise = 0.01 * (2.0*ranf_() - 1.0)
              t(i,j,k) = t(i,j,k) + noise
              sum_noise = sum_noise + noise
              npts = npts + 1

            end if
          end do
        end do

      end if
    end do

    ! remove mean of noise in the cloud so the perturbations have zero mean
    if (npts .gt. 0) then
      mean_noise = sum_noise / dble(npts)

      do k = 1, nzm
        zzz = z(k)
        if (zzz .ge. zbase .and. zzz .le. ztop) then
          do j = 1, ny
            do i = 1, nx
              xxx = dx * (i + it)
              yyy = dy * (j + jt)
              rxy = sqrt( (xxx - xc)**2 + YES3D*(yyy - yc)**2 )
              if (rxy .le. radius) then
                t(i,j,k) = t(i,j,k) - mean_noise
              end if
            end do
          end do
        end if
      end do
    end if

    ! Buoyant compensation: simple heuristic to offset virtual temperature decrease
    ! due to added condensate mass: add ~ T * qice (small correction)
    do k = 1, nzm
      zzz = z(k)
      if (zzz .ge. zbase .and. zzz .le. ztop) then
        do j = 1, ny
          do i = 1, nx
            xxx = dx * (i + it)
            yyy = dy * (j + jt)
            rxy = sqrt( (xxx - xc)**2 + YES3D*(yyy - yc)**2 )
            if (rxy .le. radius) then
              qice_mmr = micro_field(i,j,k, iqit(ii))
              Tlocal = tabs(i,j,k)
              ! add a small warming to offset the virtual temperature reduction
              t(i,j,k) = t(i,j,k) + Tlocal * qice_mmr
            end if
          end do
        end do
      end if
    end do


  case(-1)
    !bloss: no perturbation

  case default

       if(masterproc) print*,'perturb_type is not defined in setperturb(). Exitting...'
       call task_abort()

end select

end subroutine setperturb
