! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Program to read fixed-format netCDF files and compute broadband fluxes
!
! Program is invoked as rte_examples num_cols scheme_name input_file output_file
! -------------------------------------------------------------------------------------------------
program rte_examples
  ! --------------------------------------------------
  ! Modules for working with RTE
  !
  ! Working precision for real variables
  use mo_rte_kind,           only: wp
  !
  ! Optical properties of the atmosphere as array of values
  !   In the longwave we include only absorption optical depth (_1scl)
  !   Shortwave calculations would use optical depth, single-scattering albedo, asymmetry parameter (_2str)
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, &
                                   ty_optical_props_1scl, &
                                   ty_optical_props_2str
  !
  ! Gas optics uses a derived type to represent gas concentrations compactly...
  use mo_gas_concentrations, only: ty_gas_concs
  !
  ! ... and another type to encapsulate the longwave source functions.
  use mo_source_functions,   only: ty_source_func_lw
  !
  ! RTE longwave and shortwave solvers
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  !
  ! RTE driver uses a derived type to reduce spectral fluxes to whatever the user wants
  !   Here we're just reporting broadband fluxes
  use mo_fluxes,             only: ty_fluxes_broadband
  ! --------------------------------------------------
  ! Gas optics: maps physical state of the atmosphere to optical properties
  !    The optics that gets used is chosen at run time an argument
  use mo_gas_optics,         only: ty_gas_optics
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp ! Correlated-k
  use mo_optics_ssm,         only: ty_optics_ssm        ! Simple spectral models
  ! --------------------------------------------------
  !
  ! modules for reading and writing files
  !
  use mo_testing_utils,      only: stop_on_err
  !
  ! Some gas optics classes need to be initialized with data read from files
  !
  use mo_optics_utils_rrtmgp,only: load_gas_optics_rrtmgp => load_gas_optics

  use mo_rte_examples_io,    only: inquire_rte_example, read_rte_example, write_rte_example
  ! --------------------------------------------------
  implicit none
  !
  ! Command line arguments
  !
  integer            :: nUserArgs
  character(len=16)  :: char_input, scheme_name
  character(len=256) :: scheme_file, problem_file, solution_file
  integer            :: block_size
  integer            :: ncol, nlay
  !
  ! Output variables
  !
  real(wp), dimension(:,:), target, &
                            allocatable :: flux_up, flux_dn, flux_dir
  !
  ! Derived types from the RTE and RRTMGP libraries
  !
  type(ty_gas_concs)           :: gas_concs
  class(ty_optical_props_arry), &
                 allocatable   :: optical_props
  class(ty_gas_optics), &
                 allocatable   :: gas_optics
  type(ty_source_func_lw)      :: source
  type(ty_fluxes_broadband)    :: fluxes

  real(wp), dimension(:,:),  allocatable &
	   :: pres_layer, pres_level, temp_layer, temp_level
  real(wp), dimension(:),    allocatable &
	   :: surface_albedo, solar_zenith_angle, total_solar_irradiance, &
	               surface_emissivity, surface_temperature
  real(wp), dimension(:, :), allocatable  &
	   :: sfc_props_spectral, toa_flux


  logical :: do_rrtmgp = .false., &
             do_ssm    = .false., &
             do_ddq    = .false.
  logical :: do_lw
  ! --------------------------------------------------
  nUserArgs = command_argument_count()
  if (nUserArgs <  5) call stop_on_err("Usage: rte_examples block_size scheme_name scheme_data problem_file solution_file")

  call get_command_argument(1, char_input)
  read(char_input, '(i8)') block_size
  if(block_size <= 0) call stop_on_err("Specify positive block_size.")
  call get_command_argument(2, scheme_name)
  call get_command_argument(3, scheme_file)
  call get_command_argument(4, problem_file)
  call get_command_argument(5, solution_file)

  do_rrtmgp = trim(scheme_name) == "rrtmgp"
  do_ssm    = trim(scheme_name) == "ssm"
  do_ddq    = trim(scheme_name) == "ddq"
  if(.not. any([do_rrtmgp, do_ssm, do_ddq])) &
    call stop_on_err("Unknown optics scheme")

  !
  ! Read problem file (from col, lay/lev, variant to cols, lay/lev)
  !   determine which gas concetentrations are available (needed by RRTMGP)
  !
  call inquire_rte_example(problem_file, ncol, nlay, gas_concs)
  !
  ! Initialize optics schemes
  !
  if(do_rrtmgp) then
    allocate(ty_gas_optics_rrtmgp::gas_optics)
  else if (do_ssm) then
    allocate(ty_optics_ssm::gas_optics)
  end if

  select type (gas_optics)
    type is (ty_gas_optics_rrtmgp)
      call load_gas_optics_rrtmgp(gas_optics, trim(scheme_file), gas_concs)
    type is (ty_optics_ssm)
      ! SSM requires no files, instead use "lw" or "sw" as the file name to pick
      call stop_on_err(gas_optics%configure(do_sw = trim(scheme_file) == "sw"))
  end select
  !
  ! Read longwave/shortwave example, expand boundary conditions to spectral dimension
  !
  allocate(sfc_props_spectral(gas_optics%get_nband(), ncol))
  if (gas_optics%source_is_internal()) then
    call read_rte_example(problem_file, do_lw, &
  	               pres_layer, pres_level, temp_layer, temp_level, gas_concs, &
  	               surface_emissivity = surface_emissivity, &
                   surface_temperature = surface_temperature)
    allocate(ty_optical_props_1scl::optical_props)
    call stop_on_err(source%alloc(ncol, nlay, gas_optics))
    !
    ! From broadband to spectral surface emissivity
    !
    sfc_props_spectral(:,:) = spread(surface_emissivity, &
    	                             dim=1, ncopies = gas_optics%get_nband())

  else
    call read_rte_example(problem_file, do_lw, &
  	               pres_layer, pres_level, temp_layer, temp_level, gas_concs, &
  	               surface_albedo     = surface_albedo,     &
                   solar_zenith_angle = solar_zenith_angle, &
                   total_solar_irradiance = total_solar_irradiance)
    allocate(ty_optical_props_2str::optical_props)
    allocate(toa_flux(ncol, gas_optics%get_ngpt()))
    sfc_props_spectral(:,:) = spread(surface_albedo, &
    	                               dim=1, ncopies = gas_optics%get_nband())
  end if
  select type (optical_props)
    type is (ty_optical_props_1scl)
      call stop_on_err(optical_props%alloc_1scl(ncol, nlay, gas_optics))
    type is (ty_optical_props_2str)
      call stop_on_err(optical_props%alloc_2str(ncol, nlay, gas_optics))
  end select
  ! ------------------------------------------------------------------------
  !
  ! Cycle over blocks, compute fluxes
  ! We'll add the cycling over blocks later
  !
  allocate(fluxes%flux_up(ncol, nlay+1), &
           fluxes%flux_dn(ncol, nlay+1))
  if (gas_optics%source_is_internal()) then
    !
    ! Longwave calculations
    !
  	call stop_on_err( &
	    gas_optics%gas_optics(pres_layer,  &
                            pres_level,  &
                            temp_layer,  &
                            surface_temperature, &
                            gas_concs,     &
                            optical_props, &
                            source,        &
                            tlev = temp_level) &
	  )
	  call stop_on_err(         &
	    rte_lw(optical_props,   &
	           source,          &
	           sfc_props_spectral,   &
	           fluxes) &
	  )
  else
    !
    ! Shortwave calculations
    !
    allocate(fluxes%flux_dn_dir(ncol, nlay+1))
    call stop_on_err( &
      gas_optics%gas_optics(pres_layer,  &
                            pres_level,  &
                            temp_layer,  &
                            gas_concs,     &
                            optical_props, &
                            toa_flux)      &
      )
    call stop_on_err( &
      rte_sw(optical_props,   &
             cos(solar_zenith_angle * acos(-1._wp)/180._wp), &
             toa_flux,        &
             sfc_props_spectral, &
             sfc_props_spectral, &
             fluxes)          &
    )
  end if
  ! ------------------------------------------------------------------------
  call write_rte_example(ncol, nlay, fluxes, solution_file)
end program rte_examples
