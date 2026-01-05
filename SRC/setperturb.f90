subroutine setperturb

!  Random noise

use vars
use params
use microphysics, only: micro_field, index_water_vapor
use sgs, only: setperturb_sgs

implicit none

integer i,j,k,ptype,it,jt
real rrr,ranf_
real xxx,yyy,zzz

! --- variables used by case(30) (must be declared up here)
real xc, yc, rxy, noise

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

  case(-1)
    !bloss: no perturbation

  case default

       if(masterproc) print*,'perturb_type is not defined in setperturb(). Exitting...'
       call task_abort()

end select

end subroutine setperturb
