! Copyright 2026-
!    Trustees of Columbia University in the City of New York
! All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! ----------------------------------------------------------------------------
!!
! Gas, cloud, and aerosol optics classes need to be initialized with data; ddq data is distributed as a netCDF file.
!    The gas optics classes themselves don't include methods for reading the data so we don't conflict with users'
!    local environment. This module provides a straight-forward serial implementation of reading the data
!    and calling gas_optics%load().
!
!
module mo_optics_ddq_utils
  use mo_rte_kind,           only: wp, wl
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_gas_optics_ddq,     only: ty_gas_optics_ddq
  use mo_testing_utils,      only: stop_on_err
  use mo_gas_optics_ddq_kernels, &
                             only: fax_norder, fax_nterms, xsec_nterms
  ! --------------------------------------------------
  use mo_simple_netcdf, only: read_field, read_char_vec, var_exists, get_dim_size
  use netcdf
  implicit none

  private
  public :: load_gas_optics
  integer, private :: gas_name_len = 8
contains
  !--------------------------------------------------------------------------------------------------------------------
  ! read optical coefficients from NetCDF file
  subroutine load_gas_optics(gas_optics, filename)
    class(ty_gas_optics_ddq), intent(inout) :: gas_optics
    character(len=*),         intent(in   ) :: filename
    ! --------------------------------------------------
    ! -------------------------------------
    ! Spectral discretization
    real(wp), allocatable :: nus(:), weights(:)  ! (nnu) Wavenumbers and weights for the DDQ scheme
    ! -------------------------------------
    ! Functional approximations to cross-sections (fax) for line absorption
    character(len=gas_name_len), allocatable :: fax_species_names(:)
    real(wp), allocatable :: fax_a(:,:,:), fax_b(:,:,:) ! (0:2, nspecies, nnu)
    real(wp), allocatable :: fax_T0(:)       ! (nspecies)
    real(wp), allocatable :: fax_c(:,:,:)    ! (0:3, nspecies, nnu)
    real(wp), allocatable :: fax_p0(:)       ! (nspecies)
    real(wp), allocatable :: fax_sigma0(:,:) ! (     nspecies, nnu), reference absorption coefficient at p_0, T_0
    real(wp), allocatable :: fax_S(:)        ! (     nspecies), self-broadening coefficients
    ! -------------------------------------
    ! cross-section fits (xsec)
    character(len=gas_name_len), allocatable :: xsec_species_names(:)
    real(wp),         allocatable :: xsec_p(:,:,:) ! (0:3, nspecies, nnu), (const, T, T^2, p scaling of cross-section)
    ! -------------------------------------
    ! MT_CKD continuum absorption (mtckd)
    character(len=gas_name_len), allocatable :: mtckd_species_names(:)
    real(wp),         allocatable :: mtckd_cself(:, :), mtckd_cfrgn(:, :), mtckd_n(:, :) ! self- and foreign continuua
    real(wp),         allocatable :: mtckd_T0, mtckd_p0
    ! -------------------------------------
    ! Splar source function
    real(wp), allocatable :: rayleigh_xsec(:) ! (nnu)
    ! Splar source function
    real(wp), allocatable :: solar_source(:) ! (nnu)

    ! -----------------
    !
    ! Book-keeping variables
    !
    integer :: ncid
    integer :: nnu, &
               fax_nspecies,  &
               xsec_nspecies, &
               mtckd_nspecies
    ! --------------------------------------------------
    !
    ! How big are the various arrays?
    !
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("load_gas_optics(): can't open file " // trim(fileName))
    nnu            = get_dim_size(ncid,'nu')
    fax_nspecies   = get_dim_size(ncid,'fax_nspecies')
    xsec_nspecies  = get_dim_size(ncid,'xsec_nspecies')
    mtckd_nspecies = get_dim_size(ncid,'mtckd_nspecies')

    ! -----------------
    !
    ! Read the many arrays
    !
    nus     = read_field(ncid, 'nu',      nnu)
    weights = read_field(ncid, 'weights', nnu)
    fax_species_names  = read_char_vec(ncid, 'fax_species_names', fax_nspecies)
    fax_a  = read_field(ncid, 'fax_a', fax_norder+1, fax_nspecies, nnu)
    fax_b  = read_field(ncid, 'fax_b', fax_norder+1, fax_nspecies, nnu)
    fax_c  = read_field(ncid, 'fax_c', fax_nterms+1, fax_nspecies, nnu)
    fax_T0 = read_field(ncid, 'fax_T0', fax_nspecies)
    fax_p0 = read_field(ncid, 'fax_p0', fax_nspecies)
    fax_S  = read_field(ncid, 'fax_S',  fax_nspecies)
    fax_sigma0 = read_field(ncid, 'fax_sigma0', fax_nspecies, nnu)

    xsec_species_names  = read_char_vec(ncid, 'xsec_species_names',  xsec_nspecies)
    xsec_p = read_field(ncid, 'xsec_p', xsec_nterms+1, xsec_nspecies, nnu)

    mtckd_species_names = read_char_vec(ncid, 'mtckd_species_names', mtckd_nspecies)
    mtckd_cself = read_field(ncid, 'mtckd_cself', mtckd_nspecies, nnu)
    mtckd_cfrgn = read_field(ncid, 'mtckd_cfrgn', mtckd_nspecies, nnu)
    mtckd_n     = read_field(ncid, 'mtckd_n',     mtckd_nspecies, nnu)
    mtckd_T0    = as_scalar(read_field(ncid, 'mtckd_T0', 1))
    mtckd_p0    = as_scalar(read_field(ncid, 'mtckd_p0', 1))
    ! --------------------------------------------------
    !
    ! Initialize the gas optics class with data. The calls look slightly different depending
    !   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
    ! gas_optics%load() returns a string; a non-empty string indicates an error.
    !
    if(var_exists(ncid, 'solar_spectral_radiance')) then
      !
      ! Cheating a bit, expecting rayleigh_xsec and solar_spectral_radiance to be present together
      !
      rayleigh_xsec = read_field(ncid, 'rayleigh_xsec', nnu)
      solar_source = read_field(ncid, 'solar_source', nnu)
      call stop_on_err(gas_optics%load( &
                      nus, weights,     &
                      fax_species_names, fax_a, fax_b, fax_T0, fax_c, fax_p0, fax_sigma0, fax_S, &
                      xsec_species_names, xsec_p, &
                      mtckd_species_names, mtckd_cself, mtckd_cfrgn, mtckd_n, mtckd_T0, mtckd_p0, &
                      rayleigh_xsec, solar_source))
    else
      call stop_on_err(gas_optics%load( &
                      nus, weights,     &
                      fax_species_names, fax_a, fax_b, fax_T0, fax_c, fax_p0, fax_sigma0, fax_S, &
                      xsec_species_names, xsec_p, &
                      mtckd_species_names, mtckd_cself, mtckd_cfrgn, mtckd_n, mtckd_T0, mtckd_p0))
    end if
    ! --------------------------------------------------
    ncid = nf90_close(ncid)
  end subroutine load_gas_optics
    ! --------------------------------------------------
  function as_scalar(x)
    real(wp), dimension(1), intent(in):: x
    real(wp) :: as_scalar

    as_scalar = x(1)
  end function as_scalar
    ! --------------------------------------------------
end module mo_optics_ddq_utils
