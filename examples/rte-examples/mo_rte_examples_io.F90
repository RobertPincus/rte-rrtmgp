module mo_rte_examples_io
  use mo_rte_kind,           only: wp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_gas_optics_util_string, &
                            only: lower_case, string_in_array
  use mo_testing_utils,     only: stop_on_err
  use mo_simple_netcdf,     only: get_dim_size, read_field, &
                                  create_dim, create_var, write_field
  use netcdf
  implicit none

  private
  public :: inquire_rte_example, read_rte_example, write_rte_example
  ! h2o, co2, ch4, n2o, co, n2, o2, cfc11, cfc12, ! DDQ
  ! h2o co2 o3 n2o co ch4 o2 n2 ccl4 cfc11 cfc12 cfc22 hfc143a hfc125 hfc23 hfc32 hfc134a cf4 no2 ! all gases in RRTMGP
  character(len=24), parameter, &
    dimension(9) :: state_variables = [ &
    "pres_layer            ", &
    "pres_level            ", &
    "temp_layer            ", &
    "temp_level            ", &
    "surface_albedo        ", &
    "solar_zenith_angle    ", &
    "total_solar_irradiance", &
    "surface_emissivity    ", &
    "surface_temperature   "  ]

  ! Module variable - set when data is read, used when writing out resuls
  integer :: nvar = 0
contains
  ! --------------------------------------------------
  !
  ! Get problem sizes and the list of available gas concentrations (which RRTMGP needs)
  ! Anything that isn't a state variable is assumed to be a gas concentration
  ! State variables are
  !    pres_layer, pres_leve, temp_layer, temp_level, ps?
  !    surface_albedo, solar_zenith_angle, total_solar_irradiance
  !    surface_emissivity, surface_temperature
  !
  subroutine inquire_rte_example(filename, &
  	               ncol, nlay, available_gases)
    character(len=*),   intent(in)  :: filename
    integer,            intent(out) :: ncol, nlay
    type(ty_gas_concs), intent(out) :: available_gases

    integer :: ncid, ivar, nvar, ngases, var_dim, ndims, xtype
    integer, parameter :: max_vars = 64
    character(len=32) :: varnames(max_vars) = ""
    logical           :: is_gas(max_vars)
    integer           :: dimids(3)
    ! -----
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("inquire_rte_example: can't find file " // trim(fileName))

    ! The following should never happen - get_dim_size will fail first
    if(nf90_inq_dimid(ncid, "variant", var_dim) /= NF90_NOERR) return

    do ivar = 1, max_vars
      if (nf90_inquire_variable(ncid, ivar, name = varnames(ivar), xtype = xtype, &
                                ndims = ndims, dimids = dimids) /= NF90_NOERR) then
        exit
      else
        !
        ! A field might be a gas if it depends only on variant and is real valued
        !
        is_gas(ivar) = (string_in_array(varnames(ivar), ["h2o", "o3 "])) .or. &
                       ((ndims == 1) .and. (dimids(1) == var_dim) &
                        .and. (xtype == NF90_FLOAT .or. xtype == NF90_DOUBLE))
      end if
    end do
    ! write out varnames and is_gas for debugging
    ngases = 0
    do ivar = 1, max_vars
      if (trim(varnames(ivar)) == "") exit
      ngases = ngases + 1
    end do

    is_gas(ivar:max_vars) = .false.
    call stop_on_err(available_gases%init(gas_names=pack(varnames(1:ngases), mask=is_gas(1:ngases))))

    !
    ! The problem sizes as the calling routines will see them
    !
    ncol = get_dim_size(ncid, 'col') * get_dim_size(ncid, 'variant')
    nlay = get_dim_size(ncid, 'layer')

    ncid = nf90_close(ncid)

  end subroutine inquire_rte_example
  ! --------------------------------------------------
  subroutine read_rte_example(filename, is_lw, &
  	               pres_layer, pres_level, temp_layer, temp_level, gas_concs, &
  	               surface_albedo, solar_zenith_angle, total_solar_irradiance, &
  	               surface_emissivity, surface_temperature)
    character(len=*), intent(in) :: filename
    logical,          intent(in) :: is_lw
    real(wp), intent(out), dimension(:,:), allocatable &
       :: pres_layer, pres_level, temp_layer, temp_level
    type(ty_gas_concs), intent(inout) :: gas_concs
    real(wp), intent(out), dimension(:), allocatable , &
      optional :: surface_albedo, solar_zenith_angle, total_solar_irradiance, &
                   surface_emissivity, surface_temperature

    integer :: ncid
    integer :: ncol, nlay, ngas, igas
    character(len=32) :: gas_names(gas_concs%get_num_gases())
    real(wp), dimension(:), allocatable :: temp_conc

    real(wp), dimension(:,:,:), allocatable :: temp
    ! -----
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_rte_example: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'layer')
    nvar = get_dim_size(ncid, 'variant')

    !
    ! State variables
    !
    pres_layer = to_cols(spread(read_field(ncid, "pres_layer", ncol, nlay  ), &
                                dim = 3, ncopies = nvar))
    pres_level = to_cols(spread(read_field(ncid, "pres_level", ncol, nlay+1), &
                                dim = 3, ncopies = nvar))
    temp_layer = to_cols(read_field(ncid, "temp_layer", &
    	                              ncol, nlay,   nvar))
    temp_level = to_cols(read_field(ncid, "temp_level", &
    	                              ncol, nlay+1, nvar))
    if(present(surface_emissivity)) &
	    surface_emissivity = reshape(read_field(ncid, "surface_emissivity", &
	    	                                      ncol, nvar), &
	    	                           shape = [ncol * nvar])
    if(present(surface_temperature)) &
	    surface_temperature = reshape(read_field(ncid, "surface_temperature", &
	    	                                       ncol, nvar), &
	    	                            shape = [ncol * nvar])
    if(present(surface_albedo)) &
	    surface_albedo = reshape(read_field(ncid, "surface_albedo", &
                                              ncol, nvar), &
                                   shape = [ncol * nvar])
    if(present(solar_zenith_angle)) &
	    solar_zenith_angle = reshape(read_field(ncid, "solar_zenith_angle", &
                                              ncol, nvar), &
                                   shape = [ncol * nvar])
    if(present(total_solar_irradiance)) &
	    total_solar_irradiance = reshape(read_field(ncid, "total_solar_irradiance", &
                                              ncol, nvar), &
                                   shape = [ncol * nvar])

    gas_names = gas_concs%get_gas_names()
    do igas = 1, gas_concs%get_num_gases()
      if (string_in_array(gas_names(igas), ["h2o", "o3 "])) then
          call stop_on_err(                    &
          	gas_concs%set_vmr(gas_names(igas), &
          	                  to_cols(read_field(ncid, gas_names(igas), &
    	                                            ncol, nlay, nvar)) )  &
          )
      else
          temp_conc = read_field(ncid, gas_names(igas), nvar)
          call stop_on_err( &
            gas_concs%set_vmr(gas_names(igas), &
                              to_cols(spread(spread(temp_conc, dim=1, ncopies=nlay), &
                                             dim=1, ncopies=ncol)))                  &
          )
      end if
    end do

    ncid = nf90_close(ncid)

  end subroutine read_rte_example
  ! --------------------------------------------------
  subroutine write_rte_example(ncol, nlay, fluxes, solution_file)
    integer,                   intent(in) :: ncol, nlay
    type(ty_fluxes_broadband), intent(in) :: fluxes
    character(len=*),          intent(in) :: solution_file

    integer :: ncid
    real(wp), dimension(ncol, nlay+1, nvar) :: temp

    if(nf90_create(trim(solution_file), NF90_NETCDF4, ncid) /= NF90_NOERR) &
      call stop_on_err("Can't create file " // trim(solution_file))

    call create_dim(ncid, "col",    ncol/nVar)
    call create_dim(ncid, "layer",  nlay)
    call create_dim(ncid, "level",  nlay+1)
    call create_dim(ncid, "variant", nvar)
    call create_var(ncid, "flux_up",                         &
                          ["col    ", "level  ", "variant"], &
                          [ncol/nvar, nlay+1, nvar])
    call create_var(ncid, "flux_dn",                         &
                          ["col    ", "level  ", "variant"], &
                          [ncol/nvar, nlay+1, nvar])
    call stop_on_err( &
      write_field(ncid, "flux_up", to_variant(fluxes%flux_up)) &
    )
    call stop_on_err( &
      write_field(ncid, "flux_dn", to_variant(fluxes%flux_dn)) &
    )

    if(associated(fluxes%flux_dn_dir)) then
      call create_var(ncid, "flux_dir",                      &
                        ["col    ", "level  ", "variant"],   &
                          [ncol/nvar, nlay+1, nvar])
      call stop_on_err( &
        write_field(ncid, "flux_dir", to_variant(fluxes%flux_dn_dir)) &
      )
    end if

  end subroutine write_rte_example
  ! --------------------------------------------------
  !
  ! Convert from ncol, nlay, nvariant to ncol*nvariant, nlay
  !
  function to_cols(array_in) result(array_out)
    real(wp), dimension(:,:,:), intent(in) :: array_in
    real(wp), dimension(size(array_in, 1) * size(array_in, 3), &
    	                  size(array_in, 2)) :: array_out

    integer :: nx, ny, nv, v

    nx = size(array_in, 1)
    ny = size(array_in, 2)
    nv = size(array_in, 3)
    do v = 1, nv
      array_out(nx*(v-1)+1:nx*v, :) = array_in(:, :, v)
    end do
  end function to_cols
  ! --------------------------------------------------
  !
  ! Convert from ncol, nlay, nvariant to ncol*nvariant, nlay
  !
  function to_variant(array_in) result(array_out)
    real(wp), dimension(:,:), intent(in) :: array_in
    real(wp), dimension(size(array_in, 1)/nvar, &
                        size(array_in, 2),      &
                        nvar) :: array_out

    integer :: nx, v
    nx = size(array_in, 1)/nvar
    do v = 1, nvar
      array_out(:,:,v) = array_in(nx*(v-1)+1:nx*v, :)
    end do
  end function to_variant
  ! --------------------------------------------------

end module mo_rte_examples_io
