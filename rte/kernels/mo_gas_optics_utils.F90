!> Physical constants, planetary/atmospheric parameters, and utility functions for
!>    low-level gas optics calculations including Planck source functions.
!>
!> layer mass for each species
!> layer number density for each species (TK)
!> The latter two don't have C bindings
!
! Copyright 2026-, Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
module mo_gas_optics_utils
  use mo_rte_kind,       only : wp, wl
  use mo_gas_optics_constants, &
                         only: boltzmann_k, planck_h, lightspeed, &
                               m_h2o, m_dry, avogad, R_univ_gconst, grav
  use mo_rte_util_array, only : zero_array

  implicit none

  interface compute_Planck_source
    module procedure compute_Planck_source_2D, compute_Planck_source_1D
  end interface

  private :: B_nu
  public  :: compute_Planck_source, get_layer_number, interp_tlev_from_tlay

contains
  ! -------------------------------------------------------------------------------------------------
  ! Planck source functions
  ! -------------------------------------------------------------------------------------------------
  !
  ! Planck function (gets wrapped by 1D, 2D codes)
  !
  elemental function B_nu(T, nu)
    real(wp), intent(in) :: T, nu
    real(wp)             :: B_nu
    B_nu = 100._wp*2._wp*planck_h*((nu*100._wp)**3)*(lightspeed**2) / &
         (exp((planck_h * lightspeed * nu * 100._wp) / (boltzmann_k * T)) - 1._wp)
  end function
  ! -----------------------------------------
  subroutine compute_Planck_source_2D(&
      ncol, nlay, nnu, &
      nus, dnus, T, &
      source) bind(C, name="rte_compute_Planck_source_2D")
    integer,  &
      intent(in ) :: ncol, nlay, nnu
    real(wp), dimension(nnu), &
      intent(in ) :: nus, dnus
    real(wp), dimension(ncol, nlay), &
      intent(in ) :: T
    real(wp), dimension(ncol, nlay, nnu), &
      intent(out) :: source

     ! Local variables
    integer :: icol, ilay, inu

   !$acc                         parallel loop    collapse(3) copyin(nus,dnus,T) copyout(source)
   !$omp target teams distribute parallel do simd collapse(3)
   do inu = 1, nnu
      do ilay = 1, nlay
        do icol = 1, ncol
          source(icol, ilay, inu) = B_nu(T(icol, ilay), nus(inu)) * dnus(inu)
        end do
      end do
    end do

  end subroutine compute_Planck_source_2D
  ! -----------------------------------------
  subroutine compute_Planck_source_1D(&
      ncol, nnu, &
      nus, dnus, T, &
      source) bind(C, name="rte_compute_Planck_source_1D")
    integer,  &
      intent(in ) :: ncol, nnu
    real(wp), dimension(nnu), &
      intent(in ) :: nus, dnus
    real(wp), dimension(ncol), &
      intent(in ) :: T
    real(wp), dimension(ncol, nnu), &
      intent(out) :: source

     ! Local variables
     integer :: icol, ilay, inu

    !$acc                         parallel loop    collapse(2) copyin(nus,dnus,T)  copyout(source)
    !$omp target teams distribute parallel do simd collapse(2) map(to:nus,dnus,T) map(from:source)
    do inu = 1, nnu
      do icol = 1, ncol
        source(icol, inu) = B_nu(T(icol), nus(inu)) * dnus(inu)
      end do
    end do

  end subroutine compute_Planck_source_1D
  ! -------------------------------------------------------------------------------------------------
  ! layer number density
  ! -------------------------------------------------------------------------------------------------
  function get_layer_number(ncol, nlay, vmr_h2o, plev) result(col_dry)
    !
    !> Number density (#/m^-2) of dry air molecules
    !>    "col_dry" in RRTMGP
    ! input
    integer, intent(in)                           :: ncol, nlay
    real(wp), dimension(ncol, nlay  ), intent(in) :: vmr_h2o  ! volume mixing ratio of water vapor to dry air
    real(wp), dimension(ncol, nlay+1), intent(in) :: plev     ! Layer boundary pressures [Pa]
    ! output
    real(wp), dimension(ncol, nlay) :: col_dry ! Column dry amount
    ! ------------------------------------------------
    real(wp):: delta_plev, m_air, fact
    integer :: icol, ilev
    ! ------------------------------------------------
    !$acc                parallel loop gang vector collapse(2) copyin(plev,vmr_h2o)  copyout(col_dry)
    !$omp target teams distribute parallel do simd collapse(2) map(to:plev,vmr_h2o) map(from:col_dry)
    do ilev = 1, nlay
      do icol = 1, ncol
        delta_plev = abs(plev(icol,ilev) - plev(icol,ilev+1))
        ! Get average mass of moist air per mole of moist air
        fact = 1._wp / (1.+vmr_h2o(icol,ilev))
        m_air = (m_dry + m_h2o * vmr_h2o(icol,ilev)) * fact
        col_dry(icol,ilev) = delta_plev/grav * avogad * fact/m_air
      end do
    end do
  end function get_layer_number
  ! -------------------------------------------------------------------------------------------------
  ! interpolate temperature at levels from value at layer centers
  ! -------------------------------------------------------------------------------------------------
  function interp_tlev_from_tlay(ncol, nlay, tlay, play, plev) result(tlev)
    integer,  intent(in) :: ncol, nlay
    real(wp), intent(in) :: tlay(ncol, nlay), play(ncol, nlay), plev(ncol, nlay+1)
    real(wp) :: tlev(ncol, nlay+1)

    integer :: icol, ilay
    !$acc                parallel loop gang vector
    !$omp target teams distribute parallel do simd
    do icol = 1, ncol
      tlev(icol,1)      = tlay(icol,1) &
                         + (plev(icol,1)-play(icol,1))*(tlay(icol,2)-tlay(icol,1))  &
                                                        / (play(icol,2)-play(icol,1))
      tlev(icol,nlay+1) = tlay(icol,nlay)                                                             &
                        + (plev(icol,nlay+1)-play(icol,nlay))*(tlay(icol,nlay)-tlay(icol,nlay-1))  &
                                                  / (play(icol,nlay)-play(icol,nlay-1))
      end do
     !$acc                parallel loop gang vector collapse(2)
     !$omp target teams distribute parallel do simd collapse(2)
     do ilay = 2, nlay
        do icol = 1, ncol
           tlev(icol,ilay) = (play(icol,ilay-1)*tlay(icol,ilay-1)*(plev(icol,ilay  )-play(icol,ilay)) &
                           +  play(icol,ilay  )*tlay(icol,ilay  )*(play(icol,ilay-1)-plev(icol,ilay))) /  &
                             (plev(icol,ilay)*(play(icol,ilay-1) - play(icol,ilay)))
        end do
      end do
    end function interp_tlev_from_tlay
  ! -------------------------------------------------------------------------------------------------
end module mo_gas_optics_utils
