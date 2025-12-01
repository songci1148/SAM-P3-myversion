! Wrapper call to avoid circular dependency between /microphysics and SGS module BGBGBG

subroutine micro_sgs_init(tke)
  use sgs, only sgs_micro_init
  implicit none
  real :: tke(nx,ny,nz)
  call sgs_micro_init(tke)

end
