! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Example program to demonstrate the calculation of longwave radiative fluxes in clear, aerosol-free skies.
!   The example files come from the Radiative Forcing MIP (https://www.earthsystemcog.org/projects/rfmip/)
!   The large problem (1800 profiles) is divided into blocks
!
! Program is invoked as rrtmgp_rfmip_lw [block_size input_file  coefficient_file upflux_file downflux_file]
!   All arguments are optional but need to be specified in order.
! Or as ssm_rfmip_lw [block size] to test the simple spectral model
! -------------------------------------------------------------------------------------------------
program rrtmgp_rfmip_lw
  ! --------------------------------------------------
  !
  ! Modules for working with rte and rrtmgp
  !
  ! Working precision for real variables
  !
  use mo_rte_kind,           only: wp
  !
  ! Optical properties of the atmosphere as array of values
  !   In the longwave we include only absorption optical depth (_1scl)
  !   Shortwave calculations would use optical depth, single-scattering albedo, asymmetry parameter (_2str)
  !
  use mo_optical_props,      only: ty_optical_props_1scl
  !
  !
  ! Gas optics uses a derived type to represent gas concentrations compactly...
  !
  use mo_gas_concentrations, only: ty_gas_concs
  !
  ! ... and another type to encapsulate the longwave source functions.
  !
  use mo_source_functions,   only: ty_source_func_lw
  !
  ! RTE longwave driver
  !
  use mo_rte_lw,             only: rte_lw
  !
  ! RTE driver uses a derived type to reduce spectral fluxes to whatever the user wants
  !   Here we're just reporting broadband fluxes
  !
  use mo_fluxes,             only: ty_fluxes_broadband
  ! --------------------------------------------------
  ! Gas optics: maps physical state of the atmosphere to optical properties
  !    This example can use either a k-distribution from RRTMGP or a simple spectral model
  !    The optics that gets used is chosen at run time from the program name
  !
  use mo_gas_optics,         only: ty_gas_optics
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_optics_ssm,         only: ty_optics_ssm
  !
  ! modules for reading and writing files
  !
  use mo_rfmip_io,           only: read_size, read_and_block_pt, read_and_block_gases_ty, unblock_and_write, &
                                   read_and_block_lw_bc, determine_gas_names
  use mo_testing_utils,      only: stop_on_err
  !
  ! RRTMGP's gas optics class needs to be initialized with data read from a netCDF files
  !
  use mo_optics_utils_rrtmgp,only: load_gas_optics

  implicit none
  ! --------------------------------------------------
  !
  ! Local variables
  !
  character(len=512) :: invoked_name
  character(len=512) :: rfmip_file = 'multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc'
  character(len=512) :: kdist_file = 'coefficients_lw.nc'
  character(len=132) :: flxdn_file, flxup_file
  integer            :: nargs, ncol, nlay, nbnd, nexp, nblocks, block_size, forcing_index, physics_index, n_quad_angles = 1
  integer            :: b, icol, ibnd
  character(len=4)   :: block_size_char, forcing_index_char = '1', physics_index_char = '1'
  logical            :: do_rrtmgp, do_ssm

  character(len=32 ), &
            dimension(:),             allocatable :: kdist_gas_names, rfmip_gas_games
  real(wp), dimension(:,:,:),         allocatable :: p_lay, p_lev, t_lay, t_lev ! block_size, nlay, nblocks
  real(wp), dimension(:,:,:), target, allocatable :: flux_up, flux_dn
  real(wp), dimension(:,:  ),         allocatable :: sfc_emis, sfc_t  ! block_size, nblocks (emissivity is spectrally constant)
  real(wp), dimension(:,:  ),         allocatable :: sfc_emis_spec    ! nbands, block_size (spectrally-resolved emissivity)

  !
  ! Classes used by rte+rrtmgp
  !
  type(ty_source_func_lw)     :: source
  type(ty_optical_props_1scl) :: optical_props
  type(ty_fluxes_broadband)   :: fluxes

  !
  ! Optics is determined at run time
  !
  class(ty_gas_optics), allocatable :: gas_optics
  !
  ! ty_gas_concentration holds multiple columns; we make an array of these objects to
  !   leverage what we know about the input file
  !
  type(ty_gas_concs), dimension(:), allocatable  :: gas_conc_array

  ! -------------------------------------------------------------------------------------------------
  !
  ! Code starts
  !
  ! Determine which gas optics to use based on the name by which the program is evoked
  !   (possibly fragile)
  ! Based on the possibilities: rrtmgp_rfmip_lw, ssm_rfmip_lw (with or without .exe on Windows)
  call get_command_argument(0, invoked_name)
  do_rrtmgp = (invoked_name(len_trim(invoked_name)-14:len_trim(invoked_name)-8 ) == "rrtmgp_" .or. &
               invoked_name(len_trim(invoked_name)-18:len_trim(invoked_name)-12) == "rrtmgp_")
  do_ssm    = (invoked_name(len_trim(invoked_name)-11:len_trim(invoked_name)-8 ) == "ssm_"    .or. &
               invoked_name(len_trim(invoked_name)-15:len_trim(invoked_name)-12) == "ssm_")
  if (.not. (do_rrtmgp .or. do_ssm)) call stop_on_err("Don't recogize which optics to use")

  nargs = command_argument_count()
  if(nargs >= 2) call get_command_argument(2, rfmip_file)
  call read_size(rfmip_file, ncol, nlay, nexp)
  if(nargs >= 1) then
    call get_command_argument(1, block_size_char)
    read(block_size_char, '(i4)') block_size
  else
    block_size = ncol
  end if

  if(do_rrtmgp) then
    allocate(ty_gas_optics_rrtmgp::gas_optics)

    print *, "Usage: rrtmgp_rfmip_lw [block_size] [rfmip_file] [k-distribution_file] [forcing_index (1,2,3)] [physics_index (1,2)]"
    if(nargs >= 3) call get_command_argument(3, kdist_file)
    if(nargs >= 4) call get_command_argument(4, forcing_index_char)
    if(nargs >= 5) call get_command_argument(5, physics_index_char)

    read(forcing_index_char, '(i4)') forcing_index
    if(forcing_index < 1 .or. forcing_index > 3) &
      stop "Forcing index is invalid (must be 1,2 or 3)"

    read(physics_index_char, '(i4)') physics_index
    if(physics_index < 1 .or. physics_index > 2) &
      stop "Physics index is invalid (must be 1 or 2)"
    if(physics_index == 2) n_quad_angles = 3

    flxdn_file = 'rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p' //  &
                 trim(physics_index_char) // 'f' // trim(forcing_index_char) // '_gn.nc'
    flxup_file = 'rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p' // &
                 trim(physics_index_char) // 'f' // trim(forcing_index_char) // '_gn.nc'
    !
    ! Identify the set of gases used in the calculation based on the forcing index
    !   A gas might have a different name in the k-distribution than in the files
    !   provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')
    !
    call determine_gas_names(rfmip_file, kdist_file, forcing_index, kdist_gas_names, rfmip_gas_games)
    print *, "Calculation uses RFMIP gases: ", (trim(rfmip_gas_games(b)) // " ", b = 1, size(rfmip_gas_games))
  else if (do_ssm) then
    allocate(ty_optics_ssm::gas_optics)
    print *, "Usage: ssm_rfmip_lw [block_size] [rfmip_file]"
    flxdn_file = 'rld_ssm_rfmip-rad-irf.nc'
    flxup_file = 'rlu_ssm_rfmip-rad-irf.nc'
    !
    ! These variables are needed for the fragile RFMIP IO
    !
    kdist_gas_names = ["co2"]
    rfmip_gas_games = ["carbon_dioxide"]
  end if

  !
  ! How big is the problem? Does it fit into blocks of the size we've specified?
  !
  if(mod(ncol*nexp, block_size) /= 0 ) call stop_on_err("rrtmgp_rfmip_lw: number of columns doesn't fit evenly into blocks.")
  nblocks = (ncol*nexp)/block_size
  print *, "Doing ",  nblocks, "blocks of size ", block_size

  ! --------------------------------------------------
  !
  ! Prepare data for use in rte+rrtmgp
  !
  !
  ! Allocation on assignment within reading routines
  !
  call read_and_block_pt(rfmip_file, block_size, p_lay, p_lev, t_lay, t_lev)
  !
  ! Read the gas concentrations and surface properties
  !
  call read_and_block_gases_ty(rfmip_file, block_size, kdist_gas_names, rfmip_gas_games, gas_conc_array)
  call read_and_block_lw_bc(rfmip_file, block_size, sfc_emis, sfc_t)

  select type (gas_optics)
    type is (ty_gas_optics_rrtmgp)
      !
      ! Read k-distribution information. load_gas_optics() reads data from netCDF and calls
      !   gas_optics%init(); users might want to use their own reading methods
      !
      call load_gas_optics(gas_optics, trim(kdist_file), gas_conc_array(1))
      if(.not. gas_optics%source_is_internal()) &
        stop "rrtmgp_rfmip_lw: k-distribution file isn't LW"
      !
      ! RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
      !   is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
      !   This introduces an error but shows input sanitizing.
      !
      ! Are the arrays ordered in the vertical with 1 at the top or the bottom of the domain?
      if(p_lay(1, 1, 1) < p_lay(1, nlay, 1)) then
        p_lev(:,1,:) = gas_optics%get_press_min() + epsilon(gas_optics%get_press_min())
      else
        p_lev(:,nlay+1,:) &
                     = gas_optics%get_press_min() + epsilon(gas_optics%get_press_min())
      end if
    type is (ty_optics_ssm)
      call stop_on_err(gas_optics%configure())
  end select
  nbnd = gas_optics%get_nband()
  print *, "number of bands is", nbnd

  !
  ! Allocate space for output fluxes (accessed via pointers in ty_fluxes_broadband),
  !   gas optical properties, and source functions. The %alloc() routines carry along
  !   the spectral discretization from the k-distribution.
  !
  allocate(flux_up(    block_size, nlay+1, nblocks), &
           flux_dn(    block_size, nlay+1, nblocks))
  allocate(sfc_emis_spec(nbnd, block_size))
  call stop_on_err(source%alloc            (block_size, nlay, gas_optics))
  call stop_on_err(optical_props%alloc_1scl(block_size, nlay, gas_optics))
  !
  ! OpenACC directives put data on the GPU where it can be reused with communication
  ! NOTE: these are causing problems right now, most likely due to a compiler
  ! bug related to the use of Fortran classes on the GPU.
  !
  !$acc enter data create(sfc_emis_spec)
  !$omp target enter data map(alloc:sfc_emis_spec)
  !$acc enter data create(optical_props, optical_props%tau)
  !$omp target enter data map(alloc:optical_props%tau)
  !$acc enter        data create(source, source%lay_source, source%lev_source, source%sfc_source)
  !$omp target enter data map(alloc:source%lay_source, source%lev_source, source%sfc_source)
  ! --------------------------------------------------
  !
  ! Loop over blocks
  !
  do b = 1, nblocks
    fluxes%flux_up => flux_up(:,:,b)
    fluxes%flux_dn => flux_dn(:,:,b)
    !
    ! Expand the spectrally-constant surface emissivity to a per-band emissivity for each column
    !   (This is partly to show how to keep work on GPUs using OpenACC)
    !
    !$acc parallel loop collapse(2) copyin(sfc_emis)
    !$omp target teams distribute parallel do simd collapse(2) map(to:sfc_emis)
    do icol = 1, block_size
      do ibnd = 1, nbnd
        sfc_emis_spec(ibnd,icol) = sfc_emis(icol,b)
      end do
    end do
    !
    ! Compute the optical properties of the atmosphere and the Planck source functions
    !    from pressures, temperatures, and gas concentrations...
    !
    call stop_on_err(gas_optics%gas_optics(p_lay(:,:,b), &
                                           p_lev(:,:,b),       &
                                           t_lay(:,:,b),       &
                                           sfc_t(:  ,b),       &
                                           gas_conc_array(b),  &
                                           optical_props,      &
                                           source,             &
                                           tlev = t_lev(:,:,b)))
    !
    ! ... and compute the spectrally-resolved fluxes, providing reduced values
    !    via ty_fluxes_broadband
    !
    call stop_on_err(rte_lw(optical_props,   &
                            source,          &
                            sfc_emis_spec,   &
                            fluxes, n_gauss_angles = n_quad_angles))
  end do
  !$acc exit data delete(sfc_emis_spec)
  !$omp target exit data map(release:sfc_emis_spec)
  !$acc exit data delete(optical_props%tau, optical_props)
  !$omp target exit data map(release:optical_props%tau)
  !$acc exit data delete(source%lay_source, source%lev_source, source%sfc_source)
  !$omp target exit data map(release:source%lay_source, source%lev_source, source%sfc_source)
  !$acc exit data delete(source)
  ! --------------------------------------------------m
  call unblock_and_write(trim(flxup_file), 'rlu', flux_up)
  call unblock_and_write(trim(flxdn_file), 'rld', flux_dn)
end program rrtmgp_rfmip_lw
