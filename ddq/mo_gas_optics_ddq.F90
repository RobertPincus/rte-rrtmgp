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
  use mo_gas_optics_utils,   only: compute_Planck_source, get_layer_number, interp_tlev_from_tlay
  use mo_gas_optics_ddq_kernels, &
                             only: tau_absorption_from_fits
  implicit none
  private
  public :: ty_gas_optics_ddq

  integer, parameter :: gas_name_len = 8
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
    character(len=gas_name_len), &
             allocatable :: fax_species_names(:)
    !
    ! Temperature dependence of absorption coefficients
    !
    real(wp), allocatable :: fax_a(:,:,:), fax_b(:,:,:) ! (0:2, nspecies, nnu)
    real(wp), allocatable :: fax_T0(:) ! (nspecies)
    !
    ! Pressure dependence
    !
    real(wp), allocatable :: fax_c(:,:,:)    ! (0:3, nspecies, nnu)
    real(wp), allocatable :: fax_p0(:)       ! (nspecies)
    !
    ! Reference cross-section, self-broadening factor
    !
    real(wp), allocatable :: fax_sigma0(:,:) ! (     nspecies, nnu), reference absorption coefficient at p_0, T_0
    real(wp), allocatable :: fax_S(:)        ! (     nspecies), self-broadening coefficients
    ! -------------------------------------
    ! cross-section fits (xsec)
    ! -------------------------------------
    character(len=gas_name_len), &
              allocatable :: xsec_species_names(:)
    real(wp), allocatable :: xsec_p(:,:,:) ! (0:3, nspecies, nnu), (const, T, T^2, p scaling of cross-section)
    ! -------------------------------------
    ! MT_CKD continuum absorption (mtckd)
    character(len=gas_name_len), &
              allocatable :: mtckd_species_names(:)
    real(wp), allocatable :: mtckd_cself(:, :), mtckd_cfrgn(:, :), mtckd_n(:, :) ! self- and foreign continuua
    real(wp)              :: mtckd_T0, mtckd_p0
    ! -------------------------------------
    ! Unique list of all gases known to the class
    character(len=gas_name_len), &
              allocatable :: species_names(:)
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
    ! --------------
    integer :: ncol, nlay, nnu
    integer :: icol, inu
    ! ----------------------------------------------------------
    error_msg = ""
    !
    ! Source function needs temperature at interfaces/levels and at layer centers
    !   Allocate small local array for tlev unconditionally
    !
    ncol = size(play,1)
    nlay = size(play,2)
    nnu  = this%get_ngpt()

    ! --------------
    call optical_props%set_top_at_1(play(1,1) < play(1, nlay))

    ! Absoption optical depth
    error_msg =  compute_tau_absorption(this,   &
                    play, plev, tlay, gas_desc, &
                    optical_props%tau, col_dry)

    ! Revisit after Rayleigh scattering is complete
    select type(optical_props)
      type is (ty_optical_props_2str)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(ncol, nlay, nnu, optical_props%g)
      type is (ty_optical_props_nstr)
        call zero_array(ncol, nlay, nnu, optical_props%ssa)
        call zero_array(optical_props%get_nmom(), &
                      ncol, nlay, nnu, optical_props%p)
    end select

    ! ----------------------------------------------------------
    !
    ! External source function is constant
    !
    !$acc enter data create(toa_src)
    !$omp target enter data map(alloc:toa_src)
    if(check_extents) then
      if(.not. extents_are(toa_src, ncol, nnu)) &
        error_msg = "gas_optics(): array toa_src has wrong size"
    end if
    if(error_msg  /= '') return

    !$acc parallel loop collapse(2)
    !$omp target teams distribute parallel do simd collapse(2)
    do inu = 1,nnu
       do icol = 1,ncol
          toa_src(icol,inu) = this%solar_source(inu)
       end do
    end do
    !$acc exit data copyout(toa_src)
    !$omp target exit data map(from:toa_src)

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
    real(wp), dimension(size(plev,1),size(plev,2)), target  :: tlev_arr
    real(wp), dimension(:,:),                       pointer :: tlev_wk
    ! ----------------------------------------------------------
    error_msg = ""
    !
    ! Source function needs temperature at interfaces/levels and at layer centers
    !   Allocate small local array for tlev unconditionally
    !
    ncol = size(play,1)
    nlay = size(play,2)
    nnu  = this%get_ngpt()

    ! --------------
    call optical_props%set_top_at_1(play(1,1) < play(1, nlay))

    ! Absoption optical depth
    error_msg =  compute_tau_absorption(this,   &
                    play, plev, tlay, gas_desc, &
                    optical_props%tau, col_dry)

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
    ! Planck function sources, interpolating tlev from tlay if necessary
    !
    if (present(tlev)) then
      tlev_wk => tlev
    else
      tlev_arr = interp_tlev_from_tlay(ncol, nlay, tlay, play, plev)
      tlev_wk => tlev_arr
    end if
    call compute_Planck_source(ncol, nlay, nnu,              &
                               this%nus, this%weights, tlay, &
                               sources%lay_source)
    call compute_Planck_source(ncol, nlay+1, nnu,            &
                               this%nus, this%weights, tlev_wk, &
                               sources%lev_source)
    call compute_Planck_source(ncol,         nnu,            &
                               this%nus, this%weights, tsfc, &
                               sources%sfc_source)

    call zero_array(ncol, nnu, sources%sfc_source_Jac)

  end function gas_optics_int
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Absorption optical depth
  !
  function compute_tau_absorption(this,             &
                          play, plev, tlay, gas_desc, &
                          tau_abs, col_dry) result(error_msg)
    class(ty_gas_optics_ddq),   intent(in ) :: this
    real(wp), dimension(:,:),   intent(in ) :: play, &   !! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),         intent(in ) :: gas_desc  !! Gas volume mixing ratios
    real(wp), dimension(:,:,:), intent(out) :: tau_abs   !! Cabsorption optical depth; (ncol,nlay,nnu)
    real(wp), dimension(:,:),   intent(in   ), &
                           optional, target :: col_dry     !! Column dry amount (molecules/cm^2); dim(ncol,nlay)
    character(len=128)                      :: error_msg
    ! -----------------------
    real(wp), dimension(size(play, 1),size(play, 2)), &
                                    target  :: dry_num_arr
    real(wp), dimension(:,:),       pointer :: dry_num
    real(wp), dimension(:,:,:), allocatable :: vmrs
    integer,  dimension(:),     allocatable :: temp
    integer,  dimension(size(this%fax_species_names))   ::   fax_num_index
    integer,  dimension(size(this%xsec_species_names))  ::  xsec_num_index
    integer,  dimension(size(this%mtckd_species_names)) :: mtckd_num_index
    character(len=32), &
              dimension(:), allocatable     :: temp_gas_names
    character(len=gas_name_len), &             ! Gases provided by users; union of these with known gases
              dimension(:), allocatable     :: provided_gases, gases_to_use
    integer :: ncol, nlay, ngas, nnu
    integer :: igas, idx_h2o
    ! -----------------------------------------------------------
    ncol = size(play, 1)
    nlay = size(play, 2)
    nnu = this%get_ngpt()
    ! -----------------------------------------------------------
    ! Data validation
    !
    ! initialization
    if (.not. this%is_loaded()) then
      error_msg = 'ERROR: spectral configuration not loaded'
      return
    end if
    !
    ! Check input data sizes and values
    !
    !$acc        data copyin(play,plev,tlay) create(   vmr,col_gas)
    !$omp target data map(to:play,plev,tlay) map(alloc:vmr,col_gas)
    if(check_extents) then
      if(.not. extents_are(play, ncol, nlay  )) &
        error_msg = "gas_optics(): array play has wrong size"
      if(.not. extents_are(tlay, ncol, nlay  )) &
        error_msg = "gas_optics(): array tlay has wrong size"
      if(.not. extents_are(plev, ncol, nlay+1)) &
        error_msg = "gas_optics(): array plev has wrong size"
      if(.not. extents_are(tau_abs, ncol, nlay, nnu)) &
        error_msg = "gas_optics(): optical depth have the wrong extents"
      if(present(col_dry)) then
        if(.not. extents_are(col_dry, ncol, nlay)) &
          error_msg = "gas_optics(): array col_dry has wrong size"
      end if
    end if

    if(error_msg == '') then
      if(check_values) then
        if(any_vals_outside(play, this%get_press_min(),this%get_press_max())) &
          error_msg = "gas_optics(): array play has values outside range"
        if(any_vals_less_than(plev, 0._wp)) &
          error_msg = "gas_optics(): array plev has values outside range"
        if(any_vals_outside(tlay, this%get_temp_min(),  this%get_temp_max())) &
          error_msg = "gas_optics(): array tlay has values outside range"
        if(present(col_dry)) then
          if(any_vals_less_than(col_dry, 0._wp)) &
            error_msg = "gas_optics(): array col_dry has values outside range"
        end if
      end if
    end if

    if(error_msg /= '') return
    ! -----------------------------------------------------------

    ! allocate on assignment
    temp_gas_names(:) = gas_desc%get_gas_names()
    ! Which gases does the user provide?
    provided_gases(:) = [(trim(temp_gas_names(igas)), &
                          igas = 1, size(temp_gas_names))]
    ! Which gases does the user provide that the scheme knows about?
    gases_to_use = pack(provided_gases, &
                        mask = [(string_in_array(provided_gases(igas), &
                                                 this%species_names),  &
                                 igas = 1, size(provided_gases))])
    ! vmr array, index 0 is the
    allocate(vmrs( 0:size(gases_to_use), ncol, nlay))
    call zero_array(ncol, nlay, vmrs(0,:,:))
    do igas = 1, ngas
      error_msg = gas_desc%get_vmr(gases_to_use(igas), vmrs(igas,:,:))
      if (error_msg /= "") return
    end do

    fax_num_index = [(max(string_loc_in_array(this%fax_species_names(igas), &
                                              gases_to_use),                &
                          0), &
                     igas = 1, size(this%fax_species_names) &
                    )]
    xsec_num_index = [(max(string_loc_in_array(this%xsec_species_names(igas), &
                                               gases_to_use),                &
                          0), &
                       igas = 1, size(this%xsec_species_names) &
                     )]
    mtckd_num_index = [(max(string_loc_in_array(this%mtckd_species_names(igas), &
                                                gases_to_use),                &
                            0), &
                        igas = 1, size(this%mtckd_species_names) &
                      )]

    if (present(col_dry)) then
      !$acc        enter data copyin(col_dry)
      !$omp target enter data map(to:col_dry)
      dry_num => col_dry
    else
      !$acc        enter data create(   col_dry_arr)
      !$omp target enter data map(alloc:col_dry_arr)
      dry_num => dry_num_arr
      idx_h2o = string_loc_in_array("h2o", gases_to_use)
      dry_num_arr = get_layer_number(ncol, nlay,       &
                                     vmrs(idx_h2o,:,:), &
                                     plev) ! dry air column amounts computation
    end if

    call tau_absorption_from_fits(ncol, nlay, this%get_ngpt(), ngas, &
                  this%nus, &
                  play, tlay, dry_num, vmrs, &
                  size(fax_num_index), fax_num_index,   &
                  this%fax_a, this%fax_b, this%fax_T0, this%fax_c, this%fax_p0, this%fax_sigma0, this%fax_S, &
                  size(xsec_num_index), xsec_num_index,  &
                  this%xsec_p,                           &
                  size(mtckd_num_index), mtckd_num_index, &
                  this%mtckd_cself, this%mtckd_cfrgn, this%mtckd_n, this%mtckd_T0, this%mtckd_p0, &
                  tau_abs)
  end function compute_tau_absorption
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
    real(wp), intent(in) :: fax_T0(:)       ! (nspecies)
    real(wp), intent(in) :: fax_c(:,:,:)    ! (0:3, nspecies, nnu)
    real(wp), intent(in) :: fax_p0(:)       ! (nspecies)
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

    character(len=gas_name_len), dimension(:), &
      allocatable :: all_names
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

    allocate(this%fax_species_names(  fax_nspecies),      &
             this%fax_a(0:fax_norder, fax_nspecies, nnu), &
             this%fax_b(0:fax_norder, fax_nspecies, nnu), &
             this%fax_c(0:fax_nterms, fax_nspecies, nnu), &
             this%fax_sigma0(         fax_nspecies, nnu), &
             this%fax_S (nnu))
    this%fax_species_names = fax_species_names
    this%fax_a = fax_a
    this%fax_b = fax_b
    this%fax_c = fax_c
    this%fax_sigma0 = fax_sigma0
    this%fax_S  = fax_S
    this%fax_T0 = fax_T0
    this%fax_p0 = fax_p0

    allocate(this%xsec_species_names(   xsec_nspecies),    &
             this%xsec_p(0:xsec_nterms, xsec_nspecies, nnu))
    this%xsec_species_names = xsec_species_names
    this%xsec_p = xsec_p

    allocate(this%mtckd_species_names(mtckd_nspecies),    &
             this%mtckd_cself(mtckd_nspecies, nnu), &
             this%mtckd_cfrgn(mtckd_nspecies, nnu), &
             this%mtckd_n    (mtckd_nspecies, nnu))
    this%mtckd_species_names = mtckd_species_names
    this%mtckd_cself = mtckd_cself
    this%mtckd_cfrgn = mtckd_cfrgn
    this%mtckd_n     = mtckd_n

    ! Make a consolidated list of all gases?
    allocate(all_names(fax_nspecies + xsec_nspecies + mtckd_nspecies))
    all_names(:fax_nspecies              ) = &
      fax_species_names(:)
    all_names(fax_nspecies+1:fax_nspecies+xsec_nspecies) = &
      xsec_species_names(:)
    all_names(fax_nspecies+xsec_nspecies+1:) = &
      mtckd_species_names(:)
    do i = 2, size(all_names)
      if (string_in_array(all_names(i), all_names(:i-1))) all_names(i) = ""
    end do

    this%species_names(:) = pack(all_names, mask = all_names /= "")
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

    get_press_min = 0.01_wp
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

    get_temp_min = 150._wp
  end function get_temp_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the maximum temparature for which absorption coefficient fits are valid
  !
  pure function get_temp_max(this)
    class(ty_gas_optics_ddq), intent(in) :: this
    real(wp)                             :: get_temp_max !! maximum temperature for which fits are valid

    get_temp_max = 350._wp
  end function get_temp_max
  !--------------------------------------------------------------------------------------------------------------------

end module mo_gas_optics_ddq
