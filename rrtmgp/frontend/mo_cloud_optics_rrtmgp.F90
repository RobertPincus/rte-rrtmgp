! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
! Provides cloud optical properties as a function of effective radius for the RRTMGP bands
!   or by g-point. Based on Mie calculations for liquid and results from Yang et al. (2013)
!     (doi:10.1175/JAS-D-12-039.1) for ice with variable surface roughness.
!   Can use either look-up tables by spectral band or by g-point.
!
! The class can be used as-is but is also intended as an example of how to extend the RTE framework
! -------------------------------------------------------------------------------------------------

module mo_cloud_optics_rrtmgp
  use mo_rte_kind,      only: wp, wl
  use mo_rte_config,    only: check_values, check_extents
  use mo_rte_util_array_validation,&
                        only: any_vals_less_than, any_vals_outside, extents_are
  use mo_optical_props, only: ty_optical_props,      &
                              ty_optical_props_arry, &
                              ty_optical_props_1scl, &
                              ty_optical_props_2str, &
                              ty_optical_props_nstr
  use mo_cloud_optics_rrtmgp_kernels, &
                        only: compute_cld_from_table
  implicit none
  private
  ! -----------------------------------------------------------------------------------
  type, extends(ty_optical_props), public :: ty_cloud_optics_rrtmgp
    private
    !
    ! Ice surface roughness category - needed for Yang (2013) ice optics parameterization
    !
    integer            :: icergh = 0  ! (1 = none, 2 = medium, 3 = high)
    !
    ! Lookup table information
    !
    ! Upper and lower limits of the tables
    real(wp) :: radliq_lwr = 0._wp, radliq_upr = 0._wp
    real(wp) :: diamice_lwr = 0._wp, diamice_upr = 0._wp
    ! How many steps in the table? (for convenience)
    integer  :: liq_nsteps = 0,        ice_nsteps = 0
    ! How big is each step in the table?
    real(wp) :: liq_step_size = 0._wp, ice_step_size = 0._wp
    !
    ! Cloud optics lookup tables  - by g-point or by band (with ngpt=nbnd)
    !
    real(wp), dimension(:,:  ), allocatable :: extliq, ssaliq, asyliq ! (nsize_liq, ngpt)
    real(wp), dimension(:,:,:), allocatable :: extice, ssaice, asyice ! (nsize_ice, ngpt, nrghice)
    !
    ! -----
  contains
    procedure, public :: load
    procedure, public :: finalize
    procedure, public :: cloud_optics
    procedure, public :: get_min_radius_liq
    procedure, public :: get_min_diameter_ice
    procedure, public :: get_max_radius_liq
    procedure, public :: get_max_diameter_ice
    procedure, public :: get_num_ice_roughness_types
    procedure, public :: set_ice_roughness
  end type ty_cloud_optics_rrtmgp

contains
  ! ------------------------------------------------------------------------------
  !
  ! Routines to load lookup table data needed for cloud optics calculations either
  !    by spectral band or by g-point.
  !
  ! ------------------------------------------------------------------------------
  function load(this, band_lims_wvn, &
                radliq_lwr, radliq_upr, &
                diamice_lwr, diamice_upr, &
                extliq, ssaliq, asyliq, &
                extice, ssaice, asyice, &
                band_lims_gpt) result(error_msg)

    class(ty_cloud_optics_rrtmgp), intent(inout) :: this
    real(wp), dimension(:,:),   intent(in   ) :: band_lims_wvn ! beginning and ending wavenumbers for each band
    ! Lookup table interpolation constants
    ! Lower and upper bounds of the tables; also the constant for calculating interpolation indices for liquid
    real(wp),                   intent(in   ) :: radliq_lwr, radliq_upr
    real(wp),                   intent(in   ) :: diamice_lwr, diamice_upr
    ! LUT coefficients
    ! Extinction, single-scattering albedo, and asymmetry parameter for liquid and ice respectively
    real(wp), dimension(:,:),   intent(in   ) :: extliq, ssaliq, asyliq
    real(wp), dimension(:,:,:), intent(in   ) :: extice, ssaice, asyice
    integer,  dimension(:,:), optional,&
                                intent(in   ) :: band_lims_gpt ! beginning and ending g-points for each band
    character(len=128)    :: error_msg
    ! -------
    !
    ! Local variables
    !
    integer :: nrghice, nsize_liq, nsize_ice
    integer :: ngpt,i,j,k

    error_msg = this%init(band_lims_wvn, band_lims_gpt, name="RRTMGP cloud optics")
    !
    ! LUT coefficient dimensions
    !
    nsize_liq = size(extliq,dim=1)
    nsize_ice = size(extice,dim=1)
    nrghice   = size(extice,dim=3)
    ngpt     = this%get_ngpt() ! Same as the number of bands if defined by-band
    !
    ! Error checking
    !   Can we check for consistency between table bounds
    !
    if(.not. extents_are(extliq, nsize_liq, ngpt)) &
      error_msg = "cloud_optics%init(): array extliq isn't consistently sized"
    if(.not. extents_are(ssaliq, nsize_liq, ngpt)) &
      error_msg = "cloud_optics%init(): array ssaliq isn't consistently sized"
    if(.not. extents_are(asyliq, nsize_liq, ngpt)) &
      error_msg = "cloud_optics%init(): array asyliq isn't consistently sized"
    if(.not. extents_are(extice, nsize_ice, ngpt, nrghice)) &
      error_msg = "cloud_optics%init(): array extice isn't consistently sized"
    if(.not. extents_are(ssaice, nsize_ice, ngpt, nrghice)) &
      error_msg = "cloud_optics%init(): array ssaice isn't consistently sized"
    if(.not. extents_are(asyice, nsize_ice, ngpt, nrghice)) &
      error_msg = "cloud_optics%init(): array asyice isn't consistently sized"
    if(error_msg /= "") return

    this%liq_nsteps = nsize_liq
    this%ice_nsteps = nsize_ice
    this%liq_step_size = (radliq_upr - radliq_lwr)/real(nsize_liq-1,wp)
    this%ice_step_size = (diamice_upr - diamice_lwr)/real(nsize_ice-1,wp)
    ! Allocate LUT coefficients
    allocate(this%extliq(nsize_liq, ngpt), &
             this%ssaliq(nsize_liq, ngpt), &
             this%asyliq(nsize_liq, ngpt), &
             this%extice(nsize_ice, ngpt, nrghice), &
             this%ssaice(nsize_ice, ngpt, nrghice), &
             this%asyice(nsize_ice, ngpt, nrghice))

    !$acc enter data create(this)                                               &
    !$acc            create(this%extliq, this%ssaliq, this%asyliq)  &
    !$acc            create(this%extice, this%ssaice, this%asyice)
    !$omp target enter data &
    !$omp map(alloc:this%extliq, this%ssaliq, this%asyliq) &
    !$omp map(alloc:this%extice, this%ssaice, this%asyice)
    ! Load band LUT constants
    this%radliq_lwr = radliq_lwr
    this%radliq_upr = radliq_upr
    this%diamice_lwr = diamice_lwr
    this%diamice_upr = diamice_upr

    ! Load LUT coefficients
#ifndef AMDFLANG_WORKAROUND
    !$acc kernels
    !$omp target
    this%extliq = extliq
    this%ssaliq = ssaliq
    this%asyliq = asyliq
    this%extice = extice
    this%ssaice = ssaice
    this%asyice = asyice
    !$acc end kernels
    !$omp end target
#else
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, size(extliq, 2)
       do i = 1, size(extliq, 1)
          this%extliq(i,j) = extliq(i,j)
       end do
    end do
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, size(ssaliq, 2)
       do i = 1, size(ssaliq, 1)
          this%ssaliq(i,j) = ssaliq(i,j)
       end do
    end do
    !$omp target teams distribute parallel do collapse(2)
    do j = 1, size(asyliq, 2)
       do i = 1, size(asyliq, 1)
          this%asyliq(i,j) = asyliq(i,j)
       end do
    end do
    !$omp target teams distribute parallel do collapse(3)
    do k = 1, size(extice, 3)
       do j = 1, size(extice, 2)
          do i = 1, size(extice, 1)
             this%extice(i,j,k) = extice(i,j,k)
          end do
       end do
    end do
    !$omp target teams distribute parallel do collapse(3)
    do k = 1, size(ssaice, 3)
       do j = 1, size(ssaice, 2)
          do i = 1, size(ssaice, 1)
             this%ssaice(i,j,k) = ssaice(i,j,k)
          end do
       end do
    end do
    !$omp target teams distribute parallel do collapse(3)
    do k = 1, size(asyice, 3)
       do j = 1, size(asyice, 2)
          do i = 1, size(asyice, 1)
             this%asyice(i,j,k) = asyice(i,j,k)
          end do
       end do
    end do
#endif
    !
    ! Set default ice roughness - min values
    !
    error_msg = this%set_ice_roughness(1)
  end function load
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Finalize
  !
  !--------------------------------------------------------------------------------------------------------------------
  subroutine finalize(this)
    class(ty_cloud_optics_rrtmgp), intent(inout) :: this

    this%radliq_lwr = 0._wp
    this%radliq_upr = 0._wp
    this%diamice_lwr = 0._wp
    this%diamice_upr = 0._wp

    ! Lookup table cloud optics coefficients
    if(allocated(this%extliq)) then

      !$acc exit data delete(this%extliq, this%ssaliq, this%asyliq)  &
      !$acc           delete(this%extice, this%ssaice, this%asyice)  &
      !$acc           delete(this)
      !$omp target exit data map(release:this%extliq, this%ssaliq, this%asyliq) &
      !$omp map(release:this%extice, this%ssaice, this%asyice)


      deallocate(this%extliq, this%ssaliq, this%asyliq, &
                 this%extice, this%ssaice, this%asyice)
      this%liq_nsteps = 0
      this%ice_nsteps = 0
      this%liq_step_size = 0._wp
      this%ice_step_size = 0._wp
    end if

  end subroutine finalize
  ! ------------------------------------------------------------------------------
  !
  ! Derive cloud optical properties from provided cloud physical properties for
  ! either band or g-point discretization
  !
  ! ------------------------------------------------------------------------------
  !
  ! Compute single-scattering properties
  !
  function cloud_optics(this,                     &
                        clwp, ciwp, reliq, dgice, &
                        optical_props) result(error_msg)
    class(ty_cloud_optics_rrtmgp), &
              intent(in   ) :: this
    real(wp), intent(in   ) :: clwp  (:,:), & ! cloud liquid water path (g/m2)
                               ciwp  (:,:), & ! cloud ice water path    (g/m2)
                               reliq (:,:), & ! cloud liquid particle effective radius (microns)
                               dgice (:,:)    ! cloud ice particle effective diameter (microns)
    class(ty_optical_props_arry), &
              intent(inout) :: optical_props
                                               ! Dimensions: (ncol,nlay,ngpt/nbnd)

    character(len=128)      :: error_msg
    ! ------- Local -------
    logical(wl), dimension(size(clwp,1), size(clwp,2)) :: liqmsk, icemsk
    real(wp),    dimension(size(clwp,1), size(clwp,2), this%get_ngpt()) :: &
                ltau, ltaussa, ltaussag, itau, itaussa, itaussag
                ! Optical properties: tau, tau*ssa, tau*ssa*g
                ! liquid and ice separately
    integer  :: ncol, nlay, ngpt
    integer  :: nsizereg
    integer  :: icol, ilay, igpt
    ! scalars for total tau, tau*ssa
    real(wp) :: tau, taussa

    ! ----------------------------------------
    !
    ! Error checking
    !
    ! ----------------------------------------

    error_msg = ''
    if(.not.(allocated(this%extliq))) then
      error_msg = 'cloud optics: no data has been initialized'
      return
    end if

    ncol = size(clwp,1)
    nlay = size(clwp,2)
    ngpt = this%get_ngpt()
    !
    ! Array sizes
    !
    if (check_extents) then
      if(size(liqmsk,1) /= ncol .or. size(liqmsk,2) /= nlay) &
        error_msg = "cloud optics: liqmask has wrong extents"
      if(size(icemsk,1) /= ncol .or. size(icemsk,2) /= nlay) &
        error_msg = "cloud optics: icemsk has wrong extents"
      if(size(ciwp,  1) /= ncol .or. size(ciwp,  2) /= nlay) &
        error_msg = "cloud optics: ciwp has wrong extents"
      if(size(reliq, 1) /= ncol .or. size(reliq, 2) /= nlay) &
        error_msg = "cloud optics: reliq has wrong extents"
      if(size(dgice, 1) /= ncol .or. size(dgice, 2) /= nlay) &
        error_msg = "cloud optics: dgice has wrong extents"
      if(optical_props%get_ncol() /= ncol .or. optical_props%get_nlay() /= nlay) &
        error_msg = "cloud optics: optical_props have wrong extents"
      if(error_msg /= "") return
    end if

    !
    ! Spectral consistency
    !
    if(check_values) then
      if(.not. this%bands_are_equal(optical_props)) &
        error_msg = "cloud optics: optical properties don't have the same band structure"
      if(error_msg /= "") return
    end if

    !$acc data copyin(clwp, ciwp, reliq, dgice)                         &
    !$acc      create(ltau, ltaussa, ltaussag, itau, itaussa, itaussag) &
    !$acc      create(liqmsk,icemsk)
    !$omp target data map(to:clwp, ciwp, reliq, dgice) &
    !$omp map(alloc:ltau, ltaussa, ltaussag, itau, itaussa, itaussag) &
    !$omp map(alloc:liqmsk, icemsk)
    !
    ! Cloud masks; don't need value re values if there's no cloud
    !
    !$acc parallel loop gang vector default(present) collapse(2)
    !$omp target teams distribute parallel do simd collapse(2)
    do ilay = 1, nlay
      do icol = 1, ncol
        liqmsk(icol,ilay) = clwp(icol,ilay) > 0._wp
        icemsk(icol,ilay) = ciwp(icol,ilay) > 0._wp
      end do
    end do

    !
    ! Particle size, liquid/ice water paths
    !
    if(check_values) then
      if(any_vals_outside(reliq, liqmsk, this%radliq_lwr, this%radliq_upr)) &
        error_msg = 'cloud optics: liquid effective radius is out of bounds'
      if(any_vals_outside(dgice, icemsk, this%diamice_lwr, this%diamice_upr)) &
        error_msg = 'cloud optics: ice effective diameter is out of bounds'
      if(any_vals_less_than(clwp, liqmsk, 0._wp) .or. any_vals_less_than(ciwp, icemsk, 0._wp)) &
        error_msg = 'cloud optics: negative clwp or ciwp where clouds are supposed to be'
    end if
    if(error_msg == "") then
      !
      !
      ! ----------------------------------------
      !
      ! The band or g-point tables determining extinction coeffient, single-scattering
      !   albedo, and asymmetry parameter (g) as a function of effective radius.
      ! We compute the optical depth tau (=exintinction coeff * condensed water path)
      !   and the products tau*ssa and tau*ssa*g for liquid and ice cloud separately.
      ! These are used to determine the optical properties of ice and water cloud together.
      ! We could compute the properties for liquid and ice separately and
      !    use ty_optical_props_arry%increment but this involves substantially more division.
      !
      if (allocated(this%extliq)) then
        !
        ! Cloud optical properties by spectral band or g-point
        !
        ! Liquid
        !
        call compute_cld_from_table(ncol, nlay, ngpt, liqmsk, clwp, reliq,              &
                                    this%liq_nsteps,this%liq_step_size,this%radliq_lwr, &
                                    this%extliq, this%ssaliq, this%asyliq,              &
                                    ltau, ltaussa, ltaussag)
        !
        ! Ice
        !
        call compute_cld_from_table(ncol, nlay, ngpt, icemsk, ciwp, dgice,              &
                                    this%ice_nsteps,this%ice_step_size,this%diamice_lwr,&
                                    this%extice(:,:,this%icergh),                       &
                                    this%ssaice(:,:,this%icergh),                       &
                                    this%asyice(:,:,this%icergh),                       &
                                    itau, itaussa, itaussag)
      endif

      !
      ! Combine liquid and ice contributions into total cloud optical properties
      !   See also the increment routines in mo_optical_props_kernels
      !
      select type(optical_props)
      type is (ty_optical_props_1scl)
        !$acc parallel loop gang vector default(present) collapse(3) &
        !$acc               copyin(optical_props) copyout(optical_props%tau)
        !$omp target teams distribute parallel do simd collapse(3) &
        !$omp map(from:optical_props%tau)

        do igpt = 1, ngpt
          do ilay = 1, nlay
            do icol = 1,ncol
              ! Absorption optical depth  = (1-ssa) * tau = tau - taussa
              optical_props%tau(icol,ilay,igpt) = (ltau(icol,ilay,igpt) - ltaussa(icol,ilay,igpt)) + &
                                                   (itau(icol,ilay,igpt) - itaussa(icol,ilay,igpt))
            end do
          end do
        end do
      type is (ty_optical_props_2str)
        !$acc parallel loop gang vector default(present) collapse(3) &
        !$acc               copyin(optical_props) copyout(optical_props%tau, optical_props%ssa, optical_props%g)
        !$omp target teams distribute parallel do simd collapse(3) &
        !$omp map(from:optical_props%tau, optical_props%ssa, optical_props%g)
        do igpt = 1, ngpt
          do ilay = 1, nlay
            do icol = 1,ncol
              tau    = ltau    (icol,ilay,igpt) + itau   (icol,ilay,igpt)
              taussa = ltaussa (icol,ilay,igpt) + itaussa(icol,ilay,igpt)
              optical_props%g  (icol,ilay,igpt) = (ltaussag(icol,ilay,igpt) + itaussag(icol,ilay,igpt)) / &
                                                        max(epsilon(tau), taussa)
              optical_props%ssa(icol,ilay,igpt) = taussa/max(epsilon(tau), tau)
              optical_props%tau(icol,ilay,igpt) = tau
            end do
          end do
        end do
      type is (ty_optical_props_nstr)
        error_msg = "cloud optics: n-stream calculations not yet supported"
      end select
    end if
    !$acc end data
    !$omp end target data
  end function cloud_optics
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Inquiry functions
  !
  !--------------------------------------------------------------------------------------------------------------------
  function set_ice_roughness(this, icergh) result(error_msg)
    class(ty_cloud_optics_rrtmgp), intent(inout) :: this
    integer,                intent(in   ) :: icergh
    character(len=128)                    :: error_msg

    error_msg = ""
    if(.not. allocated(this%extice)) &
      error_msg = "cloud_optics%set_ice_roughness(): can't set before initialization"
    if (icergh < 1 .or. icergh > this%get_num_ice_roughness_types()) &
       error_msg = 'cloud optics: cloud ice surface roughness flag is out of bounds'
    if(error_msg /= "") return

    this%icergh = icergh
  end function set_ice_roughness
  !-----------------------------------------------
  function get_num_ice_roughness_types(this) result(i)
    class(ty_cloud_optics_rrtmgp), intent(in   ) :: this
    integer                               :: i

    i = 0
    if(allocated(this%extice)) i = size(this%extice, dim=3)
  end function get_num_ice_roughness_types
  !-----------------------------------------------
  function get_min_radius_liq(this) result(r)
    class(ty_cloud_optics_rrtmgp), intent(in   ) :: this
    real(wp)                              :: r

    r = this%radliq_lwr
  end function get_min_radius_liq
  !-----------------------------------------------
  function get_max_radius_liq(this) result(r)
    class(ty_cloud_optics_rrtmgp), intent(in   ) :: this
    real(wp)                              :: r

    r = this%radliq_upr
  end function get_max_radius_liq
  !-----------------------------------------------
  function get_min_diameter_ice(this) result(d)
    class(ty_cloud_optics_rrtmgp), intent(in   ) :: this
    real(wp)                              :: d

    d = this%diamice_lwr
  end function get_min_diameter_ice
  !-----------------------------------------------
  function get_max_diameter_ice(this) result(d)
    class(ty_cloud_optics_rrtmgp), intent(in   ) :: this
    real(wp)                              :: d

    d = this%diamice_upr
  end function get_max_diameter_ice
end module mo_cloud_optics_rrtmgp
