module mo_gas_optics_ddq
! This code implements data-driven quadrature
!
! Copyright 2026-, Trustees of Columbia University.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
!> ## Class implementing data-driven quadrature
!>
!> Implements a class for computing spectrally-resolved gas optical properties and source functions
!>   given atmopsheric physical properties (profiles of temperature, pressure, and gas concentrations)
!>   The class must be initialized with data (provided as a netCDF file) before being used.
!>
!> Two variants apply to internal Planck sources (longwave radiation in the Earth's atmosphere) and to
!>   external stellar radiation (shortwave radiation in the Earth's atmosphere).
!>   The variant is chosen based on what information is supplied during initialization.
!   (It might make more sense to define two sub-classes)
!
! -------------------------------------------------------------------------------------------------module gas_optics_ddq
  use mo_rte_kind,           only: wp, wl
  use mo_rte_config,         only: check_extents, check_values
  use mo_rte_util_array,     only: zero_array
  use mo_rte_util_array_validation, &
                             only: any_vals_less_than, any_vals_outside, extents_are
  use mo_optical_props,      only: ty_optical_props
  use mo_source_functions,   only: ty_source_func_lw
  use mo_gas_optics_constants,   only: boltzmann_k, lightspeed, planck_h
  use mo_gas_optics_util_string, only: lower_case, string_in_array, string_loc_in_array
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props_arry, &
                                   ty_optical_props_1scl, &
                                   ty_optical_props_2str, &
                                   ty_optical_props_nstr
  use mo_gas_optics,         only: ty_gas_optics
  use mo_gas_optics_utils,   only: compute_Planck_source, get_layer_number
  use mo_gas_optics_ddq_kernels, &
                             only: compute_tau_absorption
  implicit none
  private
  public :: ty_gas_optics_ddq

  ! -------------------------------------------------------------------------------------------------
  type, extends(ty_gas_optics), public :: ty_gas_optics_ddq
    private
    ! Spectral discretization
    ! -------------------------------------
    real(wp), allocatable :: nus(:), weights(:)  ! (nnu) Wavenumbers and weights for the DDQ scheme
    !
    real(wp), allocatable :: solar_source(:)     ! (nnu)
    !
    ! -------------------------------------
    ! Functional approximations to cross-sections (fax) for line absorption
    ! -------------------------------------
    ! Which gases?
    character(len=8), allocatable :: fax_species_names(:)
    !
    ! Temperature dependence of absorption coefficients
    !
    real(wp), allocatable :: fax_a(:,:,:), fax_b(:,:,:) ! (0:2, nspecies, nnu)
    real(wp)              :: fax_T0
    !
    ! Pressure dependence
    !
    real(wp), allocatable :: fax_c(:,:,:)    ! (0:3, nspecies, nnu)
    real(wp)              :: fax_p0
    !
    ! Reference cross-section, self-broadening factor
    !
    real(wp), allocatable :: fax_sigma0(:,:) ! (     nspecies, nnu), reference absorption coefficient at p_0, T_0
    real(wp), allocatable :: fax_S(:)        ! (     nspecies), self-broadening coefficients
    ! -------------------------------------
    ! cross-section fits (xsec)
    ! -------------------------------------
    character(len=8), allocatable :: xsec_species_names(:)
    real(wp),         allocatable :: xsec_p(:,:,:) ! (0:3, nspecies, nnu), (const, T, T^2, p scaling of cross-section)
    ! -------------------------------------
    ! MT_CKD continuum absorption (mtckd)
    character(len=8), allocatable :: mtckd_species_names(:)
    real(wp),         allocatable :: mtckd_cself(:, :), mtckd_cfrgn(:, :), mtckd_n(:, :) ! self- and foreign continuua
    real(wp)                      :: mtckd_T0, mtckd_p0
    ! -------------------------------------
    ! Unique list of all gases known to the class
    character(len=8), allocatable :: species_names(:)
  contains
    procedure, public :: source_is_internal
    procedure, public :: source_is_external
    procedure, public :: is_loaded
    procedure, public :: get_press_min
    procedure, public :: get_press_max
    procedure, public :: get_temp_min
    procedure, public :: get_temp_max
    ! procedure, public :: finalize

    procedure, public :: get_ngas
    procedure, public :: get_gases

    procedure, public :: load

    procedure, public  :: gas_optics_int
    procedure, public  :: gas_optics_ext

  end type
contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Compute gas optical depth given temperature, pressure, and composition
  !
  function gas_optics_ext(this,                         &
                          play, plev, tlay, gas_desc,   & ! mandatory inputs
                          optical_props, toa_src,       & ! mandatory outputs
                          col_dry) result(error_msg)
    class(ty_gas_optics_ddq), intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   !! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  !! Gas volume mixing ratios
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props !!
    real(wp), dimension(:,:), intent(  out) :: toa_src   !! Incoming solar irradiance(ncol,ngpt)
    character(len=128)                      :: error_msg !! Empty if successful
    ! Optional inputs
    real(wp), dimension(:,:), intent(in   ), &
                           optional, target :: col_dry !! Column dry amount (molecules/cm^2); dim(ncol,nlay)
  end function gas_optics_ext
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Compute gas optical depth and Planck source functions,
  !  given temperature, pressure, and composition
  !
  function gas_optics_int(this,                             &
                          play, plev, tlay, tsfc, gas_desc, &
                          optical_props, sources,           &
                          col_dry, tlev) result(error_msg)
    class(ty_gas_optics_ddq), intent(in   ) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   !! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    real(wp), dimension(:),   intent(in   ) :: tsfc      !! surface skin temperatures [K]; (ncol)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  !! Gas volume mixing ratios
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props !! Optical properties
    class(ty_source_func_lw    ),  &
                              intent(inout) :: sources    !! Planck sources
    character(len=128)                      :: error_msg  !! Empty if successful
    real(wp), dimension(:,:), intent(in   ), &
                          optional, target :: col_dry, &  !! Column dry amount (molecules/cm^2); dim(ncol,nlay)
                                                 tlev     !! level temperatures [K]l (ncol,nlay+1)
    ! --------------
    integer :: ncol, nlay, nnu
    ! --------------
    error_msg = ""
    ncol = size(play,1)
    nlay = size(play,2)
    nnu  = this%get_ngpt() ! How does this work?

    ! Compute layer number here
    !   get_col_dry returns number density per cm^2
    ! What to do if a gas the DDQ knows about isn't provided?
    !    Use zeros? Choose a subset of the coefficients?

    ! Absoption optical depth

    ! fill in the optical properties
    call optical_props%set_top_at_1(play(1,1) < play(1, nlay))

    select type(optical_props)
      type is (ty_optical_props_2str)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(ncol, nlay, nnu, optical_props%g)
      type is (ty_optical_props_nstr)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(optical_props%get_nmom(), &
                      ncol, nlay, nnu, optical_props%p)
    end select

    !
    ! Planck function sources
    !
    call compute_Planck_source(ncol, nlay, nnu,              &
                               this%nus, this%weights, tlay, &
                               sources%lay_source)
    ! This will fail if Tlev isn't provided
    !   There's interpolation code in RRTMGP gas optics -
    !   should we make this generic and package it with the gas optics type?
    if(.not. present(tlev)) then
      error_msg = "tlev required for DDQ! (someone should fix this)"
      return
    end if
    call compute_Planck_source(ncol, nlay+1, nnu,            &
                               this%nus, this%weights, tlev, &
                               sources%lev_source)
    call compute_Planck_source(ncol,         nnu,            &
                               this%nus, this%weights, tsfc, &
                               sources%sfc_source)

    call zero_array(ncol, nnu, sources%sfc_source_Jac)

  end function gas_optics_int
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Initialiation
  !
  function load(this,                       &
                nus, weights, solar_source, &
                fax_species_names, fax_a, fax_b, fax_T0, fax_c, fax_p0, fax_sigma0, fax_S, &
                xsec_species_names, xsec_p, &
                mtckd_species_names, mtckd_cself, mtckd_cfrgn, mtckd_n, mtckd_T0, mtckd_p0) &
                result(error_msg)
    class(ty_gas_optics_ddq), intent(inout) :: this
    ! -------------------------------------
    ! Spectral discretization
    real(wp), intent(in) :: nus(:), weights(:)  ! (nnu) Wavenumbers and weights for the DDQ scheme
    real(wp), intent(in) :: solar_source(:)     ! (nnu)
    ! -------------------------------------
    ! Functional approximations to cross-sections (fax) for line absorption
    character(len=*), intent(in) :: fax_species_names(:)
    real(wp), intent(in) :: fax_a(:,:,:), fax_b(:,:,:) ! (0:2, nspecies, nnu)
    real(wp), intent(in) :: fax_T0
    real(wp), intent(in) :: fax_c(:,:,:)    ! (0:3, nspecies, nnu)
    real(wp), intent(in) :: fax_p0
    real(wp), intent(in) :: fax_sigma0(:,:) ! (     nspecies, nnu), reference absorption coefficient at p_0, T_0
    real(wp), intent(in) :: fax_S(:)        ! (     nspecies), self-broadening coefficients
    ! -------------------------------------
    ! cross-section fits (xsec)
    character(len=*), intent(in) :: xsec_species_names(:)
    real(wp),         intent(in) :: xsec_p ! (0:3, nspecies, nnu), (const, T, T^2, p scaling of cross-section)
    ! -------------------------------------
    ! MT_CKD continuum absorption (mtckd)
    character(len=*), intent(in) :: mtckd_species_names(:)
    real(wp),         intent(in) :: mtckd_cself(:, :), mtckd_cfrgn(:, :), mtckd_n(:, :) ! self- and foreign continuua
    real(wp),         intent(in) :: mtckd_T0, mtckd_p0
    character(len = 128) :: error_msg

    ! -------------------------------------
    integer :: nnu, fax_nspecies, xsec_nspecies, mtckd_nspecies, i
    integer, parameter :: fax_norder = 2, fax_nterms = 3, xsec_nterms = 3
    ! -------------------------------------
    error_msg = this%ty_optical_props%init(band_lims_wvn = &
                                             spread(nus, dim=1, ncopies=2))
    if(len_trim(error_msg) /= 0) return

    nnu = size(nus)
    fax_nspecies   = size(  fax_species_names)
    xsec_nspecies  = size( xsec_species_names)
    mtckd_nspecies = size(mtckd_species_names)

    ! Check all dimensions?
    ! Make a consolidated list of all gases?

    ! Species names!
    allocate(this%fax_a(0:fax_norder, fax_nspecies, nnu), &
             this%fax_b(0:fax_norder, fax_nspecies, nnu), &
             this%fax_c(0:fax_nterms, fax_nspecies, nnu), &
             this%fax_sigma0(         fax_nspecies, nnu), &
             this%fax_S (nnu))
    this%fax_a = fax_a
    this%fax_b = fax_b
    this%fax_c = fax_c
    this%fax_sigma0 = fax_sigma0
    this%fax_S  = fax_S
    this%fax_T0 = fax_T0
    this%fax_p0 = fax_p0

    ! Species names!
    allocate(this%xsec_p(0:xsec_nterms, xsec_nspecies, nnu))
    this%xsec_p = xsec_p

    ! species names!
    allocate(this%mtckd_cself(mtckd_nspecies, nnu), &
             this%mtckd_cfrgn(mtckd_nspecies, nnu), &
             this%mtckd_n    (mtckd_nspecies, nnu))
    this%mtckd_cself = mtckd_cself
    this%mtckd_cfrgn = mtckd_cfrgn
    this%mtckd_n     = mtckd_n

  end function load
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if class has been initialized
  !
  pure function is_loaded(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    logical                              :: is_loaded
    is_loaded = allocated(this%nus)
  end function is_loaded
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if initialized for internal sources/longwave, false otherwise
  !
  pure function source_is_internal(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    logical                          :: source_is_internal
    source_is_internal = .not. allocated(this%solar_source)
  end function source_is_internal
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if initialized for external sources/shortwave, false otherwise
  !
  pure function source_is_external(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    logical                          :: source_is_external
    source_is_external = allocated(this%solar_source)
  end function source_is_external

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> How man gases are known to the optics?
  !
  pure function get_ngas(this)
    class(ty_gas_optics_ddq), intent(in)  :: this
    integer :: get_ngas !! names of the gases for which absorption coefficents

    get_ngas = size(this%species_names)
  end function get_ngas
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> For which species are absorption coefficients available?
  !
  pure function get_gases(this)
    class(ty_gas_optics_ddq), intent(in)  :: this
    character(32), dimension(get_ngas(this)) :: get_gases !! names of the gases for which absorption coefficents

    get_gases = this%species_names
  end function get_gases
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the minimum pressure for which absorption coefficient fits are valid
  !
  pure function get_press_min(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    real(wp)                             :: get_press_min !! minimum pressure for which fits are valid

    get_press_min = 1._wp
  end function get_press_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the maximum pressure for which absorption coefficient fits are valid
  !
  pure function get_press_max(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    real(wp)                             :: get_press_max !! maximum pressure for which fits are valid

    get_press_max = 1100.e2_wp
  end function get_press_max

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the minimum temparature for which absorption coefficient fits are valid
  !
  pure function get_temp_min(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    real(wp)                             :: get_temp_min !! minimum temperature for which fits are valid

    get_temp_min = 180._wp
  end function get_temp_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the maximum temparature for which absorption coefficient fits are valid
  !
  pure function get_temp_max(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    real(wp)                             :: get_temp_max !! maximum temperature for which fits are valid

    get_temp_max = 320._wp
  end function get_temp_max
  !--------------------------------------------------------------------------------------------------------------------

end module mo_gas_optics_ddq
