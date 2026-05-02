! This code is part of rte-fjx
!
! Contacts: Robert Pincus and George Milly
! email:  rrtmgp@aer.com
!
! Copyright 2022-  Trustees of Columbia University. All rights reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
!> ## Class implementing the RRTMGP correlated-_k_ distribution 
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
! -------------------------------------------------------------------------------------------------
module mo_gas_optics_fjx
  use mo_rte_kind,           only: wp, wl
  use mo_rte_config,         only: check_extents, check_values
  use mo_rte_util_array,     only: zero_array
  use mo_rte_util_array_validation,     only: any_vals_less_than, any_vals_outside, extents_are !currently lt & extents_are
  use mo_optical_props,      only: ty_optical_props
  use mo_source_functions,   only: ty_source_func_lw
  use mo_gas_optics_util_string, only: lower_case, string_in_array, string_loc_in_array
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_gas_optics,         only: ty_gas_optics
  use mo_rte_fjx_interpolation_tables
  use mo_gas_optics_constants,   only: avogad, m_dry, grav
  implicit none
  private
  real(wp), parameter :: pi = acos(-1._wp)

  ! -------------------------------------------------------------------------------------------------
  type, extends(ty_gas_optics), public :: ty_gas_optics_fjx
    private
    !---- Variables in file 'FJX_spec.dat' (RD_XXX)
    !> will we want any of this wavelength information?
    ! WL: Centres of wavelength bins - 'effective wavelength'  (nm)
!   real(wp)  WL(S_)
    ! WBIN: Boundaries of wavelength bins                  (microns)
!   real(wp)  WBIN(S_+1)
    ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
    ! FW: Solar flux in W/m2
    ! FP: PAR quantum action spectrum
    ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
    real(wp), dimension(:), allocatable :: qrayl, fl, fw

   !type(ty_rte_fjx_interp_table), dimension(num_cross_sections_gas_optics) :: QQQQ
    type(ty_rte_fjx_interp_table), dimension(:), allocatable :: QQQQ !is dim necessary?

    integer :: num_cross_sections_gas_optics
   !character(len=*), dimension(num_cross_sections_gas_optics), parameter :: cross_section_names=['O2','O3']
   !character(len=*), dimension(:), allocatable :: cross_section_names !is dim necessary?
    character(len=3), dimension(:), allocatable :: cross_section_names !(len=*) not allowed, other solution?

  contains
    ! Type-bound procedures
    ! Public procedures
    ! public interface
    procedure, public :: source_is_internal
    procedure, public :: source_is_external
    procedure, public :: get_press_min
    procedure, public :: get_press_max
    procedure, public :: get_temp_min
    procedure, public :: get_temp_max
    procedure, public :: load
    ! Internal procedures
    procedure, public  :: gas_optics_int
    procedure, public  :: gas_optics_ext
  end type ty_gas_optics_fjx
  ! -------------------------------------------------------------------------------------------------

contains
  ! --------------------------------------------------------------------------------------
  !
  ! Public procedures
  !
  ! --------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> Compute gas optical depth and Planck source functions,
  !>  given temperature, pressure, and composition
  !
  function gas_optics_int(this,                             &
                          play, plev, tlay, tsfc, gas_desc, &
                          optical_props, sources,           &
                          col_dry, tlev) result(error_msg)
    ! inputs
    class(ty_gas_optics_fjx), intent(in) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   !! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    real(wp), dimension(:),   intent(in   ) :: tsfc      !! surface skin temperatures [K]; (ncol)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  !! Gas volume mixing ratios
    ! output
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props !! Optical properties
    class(ty_source_func_lw    ),  &
                              intent(inout) :: sources       !! Planck sources
    character(len=128)                      :: error_msg     !! Empty if succssful 
    ! Optional inputs
    real(wp), dimension(:,:),   intent(in   ), &
                           optional, target :: col_dry, &  !! Column dry amount; dim(ncol,nlay)
                                               tlev        !! level temperatures [K]; (ncol,nlay+1)
    ! ----------------------------------------------------------

  error_msg="gas optics int should not have been called"

  end function gas_optics_int
  !------------------------------------------------------------------------------------------
  !
  !> Compute gas optical depth given temperature, pressure, and composition
  !>    Top-of-atmosphere stellar insolation is also reported 
  !
  function gas_optics_ext(this,                         &
                          play, plev, tlay, gas_desc,   & ! mandatory inputs
                          optical_props, toa_src,       & ! mandatory outputs
                          col_dry) result(error_msg)      ! optional input

    class(ty_gas_optics_fjx), intent(in) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   !! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &   !! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay      !! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),       intent(in   ) :: gas_desc  !! Gas volume mixing ratios
    ! output
    class(ty_optical_props_arry),  &
                              intent(inout) :: optical_props 
    real(wp), dimension(:,:), intent(  out) :: toa_src     !! Incoming solar irradiance(ncol,ngpt)
    character(len=128)                      :: error_msg   !! Empty if successful

    ! Optional inputs
    real(wp), dimension(:,:), intent(in   ), &
                           optional, target :: col_dry ! Column dry amount; dim(ncol,nlay)
    ! ----------------------------------------------------------
    ! Local variables
    integer :: L,C,ncol,nlay,n,idx
    real(wp) :: TTTX,massfact
    real(wp), dimension(size(play,dim=1),size(play,dim=2)) :: dry_col
    real(wp), dimension(size(play,dim=1),size(play,dim=2),optical_props%get_ngpt()) :: OD,SSA
    real(wp), dimension(optical_props%get_ngpt()) :: ODRAY,ODABS,xsec
    real(wp), dimension(size(play,dim=1),size(play,dim=2),this%num_cross_sections_gas_optics) :: gas_conc

    ncol  = size(play,dim=1)
    nlay  = size(play,dim=2)
    massfact = 100._wp*avogad/(m_dry*1000._wp*grav*10._wp) !remove factor of 100 for Pascals? m_dry*1000 -> grams
    error_msg=''

    !input sanitizing here
    if (check_values) then
      if (any_vals_less_than(tlay, 0._wp)) & 
        error_msg="all temperature values must be non-negative"
      if (any_vals_less_than(plev, 0._wp)) &   ! we don't use play right?
        error_msg="all pressure values must be non-negative"
      if(error_msg /= "") return
    endif
    if (check_extents) then
      if (.not. extents_are(plev,ncol,nlay+1)) &
        error_msg="plev incorrectly sized"
      if (.not. extents_are(tlay,ncol,nlay)) &
        error_msg="tlay incorrectly sized"
      if(error_msg /= "") return
    endif

    !output sanitizing
    select type (optical_props)
      type is (ty_optical_props_1scl)
        error_msg='rte-fjx: must provide scattering optical properties'
        return
    end select
        
    !check that data is loaded...
    if (.not. allocated(this%qrayl)) then
      error_msg='load function must be called first'
      return
    endif

    !fill in gas vmr from gas_desc type
    do n = 1, this%num_cross_sections_gas_optics
      error_msg=gas_desc%get_vmr(this%cross_section_names(n),gas_conc(:,:,n)) !no assumption on order of names?
      if (error_msg /= "") return
    end do

    ! probably check if col_dry provided, else calc
    ! should be subroutine outside of function?
    do L = 1,nlay
      do C = 1,ncol
        dry_col(C,L)=(plev(C,L)-plev(C,L+1))*massfact
        ODRAY  = dry_col(C,L)*this%qrayl(:)
        ! we're gonna assume that all species get interpolated in temperature. 
        TTTX = tlay(C,L)
        odabs(:) = 0._wp ! absorption optical depth vs g-point 
        do n = 1, this%num_cross_sections_gas_optics
          idx = string_loc_in_array(this%cross_section_names(n),this%QQQQ%title) !redundant?
          xsec(:) = interp_x(TTTX, this%QQQQ(idx))
          odabs(:) = odabs(:) + dry_col(C,L) * gas_conc(C,L,n) * xsec(:)
        end do 
        OD(C,L,:)  = ODRAY(:) + ODABS(:)       !total optical depth
        SSA(C,L,:) = ODRAY/OD(C,L,:)
      enddo
    enddo

    optical_props%tau(:,:,:) = OD(:,:,:)
    select type(optical_props)
      type is (ty_optical_props_2str)
        optical_props%ssa(:,:,:) = SSA(:,:,:)
        call zero_array(ncol,nlay,optical_props%get_ngpt(),optical_props%g)
      type is (ty_optical_props_nstr)
        optical_props%ssa(:,:,:) = SSA(:,:,:)
        call zero_array(optical_props%get_nmom(),ncol,nlay,optical_props%get_ngpt(),optical_props%p)
        optical_props%p(2,:,:,:) = 0.1_wp ! rayleigh scattering
    end select

    toa_src(1,:)=this%fl(:)
    !SPhot !solar #/cm2/s   ! SUSIM average 11Nov94(low) + 29Mar92(med-high)

   !toa_src(1,:)=this%fw(:)
    !SWatt |solarheat W/m2  | v75a+ (17:18 @ 485nm) scaled to 1360.8 W/m2

    print *,'gas optics ext has been called'
    ! ----------------------------------------------------------
  end function gas_optics_ext
  !------------------------------------------------------------------------------------------
  !
  !
  ! Initialization
  !
  !--------------------------------------------------------------------------------------------------------------------

  function load (this, in_table, wavebounds, gpt_lims, qrayl_in, fl_in, fw_in) result(error_msg)
    class(ty_gas_optics_fjx), intent(inout) :: this
    type(ty_rte_fjx_interp_table), dimension(:), intent(in) :: in_table
    real(wp), dimension(:,:), intent(in) :: wavebounds
    integer, dimension(:,:), intent(in) :: gpt_lims
    real(wp), dimension(:), intent(in) :: qrayl_in, fl_in, fw_in   ! make target?
    character(len=128)           :: error_msg     !! Empty if successful

    integer :: gas_loc, n, ngpt !get rid of gas_loc now?
    error_msg = ""

    ngpt=gpt_lims(2,1) !previously set from opt_props%get_ngpt() where opt_props was passed in
    ! Sanitize everything! Ensure are same size, non-negative  
    ! Ensure that you have the same number of g-points for qrayl_in, fl_in,
    if(check_extents) then
      if (.not. extents_are(qrayl_in,ngpt)) error_msg="qrayl incorrectly sized"
      if (.not. extents_are(fl_in,ngpt)) error_msg="flux incorrectly sized"
      if (.not. extents_are(fw_in,ngpt)) error_msg="flux incorrectly sized"
      if(error_msg /= "") return
    end if
    if (check_values) then
      if (any_vals_less_than(qrayl_in, 0._wp)) & 
        error_msg="all rayleigh values must be non-negative"
      if (any_vals_less_than(fl_in, 0._wp)) &
        error_msg="all flux values must be non-negative"
      if (any_vals_less_than(fw_in, 0._wp)) &
        error_msg="all flux values must be non-negative"
      if(error_msg /= "") return
    endif

    this%num_cross_sections_gas_optics=size(in_table)
    allocate(this%QQQQ(this%num_cross_sections_gas_optics))
    allocate(this%cross_section_names(this%num_cross_sections_gas_optics)) !redundant but maintain?

    do n = 1, this%num_cross_sections_gas_optics
     !print *, n
      ! do we not want to check against a hard-coded list of expected absorbers anymore?
      ! for now, we let the input data determine the process...
      this%cross_section_names(n)=in_table(n)%title
    ! if (.not. string_in_array(cross_section_names(n),in_table%title)) then
    !   error_msg= 'cross-section data for ' // cross_section_names(n) //' not found'
    !   return
    ! else
    !   gas_loc = string_loc_in_array(cross_section_names(n),in_table%title)
        ! Ensure that the table has the right number of g-points
        if (check_extents) then
         !print *, shape(in_table(n)%cross_section)
          if (.not. extents_are(in_table(n)%cross_section,ngpt,size(in_table(n)%coordinate))) then !ok? h2o 1 coord
            error_msg="cross section incorrectly sized"
            return
          endif
        end if
        if (check_values) then
          if (any_vals_less_than(in_table(n)%cross_section,0._wp)) then
            error_msg="cross section values must be non-negative"
            return
          endif
        end if
        ! check that interpolation is in temperature
        if (in_table(n)%if_p) then
          error_msg="gas optics requires species interpolated in temperature"
          return
        end if 
        this%QQQQ(n) = in_table(n)
    ! end if
    end do 

    allocate(this%qrayl(ngpt))
    this%qrayl=qrayl_in ! allocate on assignment?
    allocate(this%fl(ngpt))
    this%fl=fl_in
    allocate(this%fw(ngpt))
    this%fw=fw_in

    error_msg=this%init(wavebounds,gpt_lims) !?
    if(len_trim(error_msg) /= 0) return

  end function load 

  !--------------------------------------------------------------------------------------------------------------------

  !
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if initialized for internal sources/longwave, false otherwise
  !
  pure function source_is_internal(this)
    class(ty_gas_optics_fjx), intent(in) :: this
    logical                          :: source_is_internal
    source_is_internal = .FALSE.
  end function source_is_internal
  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return true if initialized for external sources/shortwave, false otherwise
  !
  pure function source_is_external(this)
    class(ty_gas_optics_fjx), intent(in) :: this
    logical                          :: source_is_external
    source_is_external = .TRUE.
  end function source_is_external

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the minimum pressure on the interpolation grids
  !
  pure function get_press_min(this)
    class(ty_gas_optics_fjx), intent(in) :: this
    real(wp)                                :: get_press_min !! minimum pressure for which the k-dsitribution is valid

    get_press_min = 0.
  end function get_press_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the maximum pressure on the interpolation grids
  !
  pure function get_press_max(this)
    class(ty_gas_optics_fjx), intent(in) :: this
    real(wp)                                :: get_press_max !! maximum pressure for which the k-dsitribution is valid

    get_press_max = 1.0E6
  end function get_press_max

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the minimum temparature on the interpolation grids
  !
  pure function get_temp_min(this)
    class(ty_gas_optics_fjx), intent(in) :: this
    real(wp)                                :: get_temp_min !! minimum temperature for which the k-dsitribution is valid 

    get_temp_min = 0.
  end function get_temp_min

  !--------------------------------------------------------------------------------------------------------------------
  !
  !> return the maximum temparature on the interpolation grids
  !
  pure function get_temp_max(this)
    class(ty_gas_optics_fjx), intent(in) :: this
    real(wp)                                :: get_temp_max !! maximum temperature for which the k-dsitribution is valid

    get_temp_max = 1.0E6
  end function get_temp_max
  !--------------------------------------------------------------------------------------------------------------------

end module mo_gas_optics_fjx
