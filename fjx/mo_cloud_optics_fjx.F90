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
!> ## Compute cloud optical properties from cloud physical properties
!>
!>
!
! -------------------------------------------------------------------------------------------------
module mo_cloud_optics_fjx
  use mo_rte_kind,           only: wp, wl
  use mo_rte_config,         only: check_extents
 !use mo_rte_util_array,     only: zero_array
  use mo_rte_util_array_validation,     only: extents_are
 !use mo_gas_optics_util_string, only: lower_case, string_in_array, string_loc_in_array
  use mo_optical_props,      only: ty_optical_props, & 
                                   ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  implicit none
  private

  ! -------------------------------------------------------------------------------------------------
  type, extends(ty_optical_props), public :: ty_cloud_optics_fjx
    private

    real(wp), allocatable, dimension(:,:,:)   :: QCC, SCC  !(wavelength, size-based classification, particle type)
    real(wp), allocatable, dimension(:,:,:,:) :: PCC  !(Legendre moment, wavelength, size-based classification, particle type)
    !need DCC too for density, and RCC
    real(wp), allocatable, dimension(:)       :: DCC
    real(wp), allocatable, dimension(:,:)     :: RCC

  contains
    ! Type-bound procedures
    ! Public procedures
    procedure, public :: load
    procedure, public :: cloud_optics
    ! Internal procedures
    procedure, private  :: optici
    procedure, private  :: opticl
  end type ty_cloud_optics_fjx
  ! -------------------------------------------------------------------------------------------------

contains
  ! --------------------------------------------------------------------------------------
  function cloud_optics(this,                             &
                        lwpx, iwpx, refflx, reffix,       &
                        tlay, optical_props) result(error_msg)
    ! inputs
    class(ty_cloud_optics_fjx), intent(in) :: this
    real(wp), dimension(:,:), intent(in)   :: lwpx, & ! liquid water path; (ncol,nlay)
                                              iwpx, & ! ice water path; (ncol, nlay)
                                              refflx, reffix, &
                                              tlay    ! temperature, for particle classification
    ! output
    class(ty_optical_props_arry),  &
                             intent(inout) :: optical_props !! Optical properties
    character(len=128)                     :: error_msg
   ! ----------------------------------------------------------
    integer :: L,C,nlay,ncol,P,ngpt,nmom
    real(wp)               :: dens_l,dens_i !should be same actually, just get from this%DCC, not liq/ice subroutine?
    real(wp), allocatable, dimension(:)       :: ext_l,ssa_l,ext_i,ssa_i
    real(wp), allocatable, dimension(:,:)     :: leg_l,leg_i
    real(wp), allocatable, dimension(:,:,:)   :: ext,ssa
    real(wp), allocatable, dimension(:,:,:,:) :: leg
    logical :: cld_present

    error_msg = ""

    ncol=size(lwpx,dim=1)
    nlay=size(lwpx,dim=2)
    ngpt=size(this%QCC,dim=1)!check against optical_props? use to allocate ext, ssa, leg, arrays?
    nmom=size(this%PCC,dim=1)

    allocate(ext(ncol,nlay,ngpt))
    allocate(ssa(ncol,nlay,ngpt))
    allocate(leg(nmom,ncol,nlay,ngpt))
    !the below vars allocate correctly on assignment upon return from optici/opticl
    !but if optici/opticl isn't called they are not allocated, and ext_i can not be added to ext_l for instance
    allocate(ext_l(ngpt))
    allocate(ext_i(ngpt))
    allocate(ssa_l(ngpt))
    allocate(ssa_i(ngpt))
    allocate(leg_l(nmom,ngpt))
    allocate(leg_i(nmom,ngpt))
    !these values likely need to be initialized as zero too, as Prather is summing them together and expecting this
    leg_l=0._wp
    leg_i=0._wp
    !that alone doesn't eliminate NaNs as leg also includes ext_l/i and ssa_l/i terms
   !ext_l=0._wp  !must be done inside layer&column loop below to reset for each iteration
   !ext_i=0._wp
    ssa_l=0._wp
    ssa_i=0._wp
    !i don't think it should be necessary to initialize (ice+liq sum) ext
   !ext=0._wp
    !i think it should be necessary to initialize (ice+liq sum) ssa & leg now, though it seems to work without
    ssa=0._wp
    leg=0._wp
    !do we need to worry about possibility of a layer with clouds, but not all gpts have nonzero ext?

    !sanitize inputs

    !sanitize outputs
    select type (optical_props)
      type is (ty_optical_props_1scl)
        error_msg='rte-fjx: must provide scattering optical properties'
        return
    end select

    !check that data is loaded
    if (.not. allocated(this%DCC)) then
      error_msg='load function must be called first'
      return
    endif

    do L=1,nlay
      do C=1,ncol
        cld_present = .false. !must be done inside loop to reset for each iteration
        ext_l=0._wp  !Prather used two different variables for the array returned from opticl/optici
        ext_i=0._wp  !and the one multiplied by the additional factors
        if (lwpx(C,L) > 1.e-5_wp .and. refflx(C,L) > 0.1_wp) then
          cld_present=.true.
          error_msg = this%opticl(refflx(C,L),          dens_l,ext_l,ssa_l,leg_l)
          ext_l = lwpx(C,L) * 0.75_wp * ext_l / (refflx(C,L) * dens_l)
        endif
        if (iwpx(C,L) > 1.e-5_wp .and. reffix(C,L) > 0.1_wp) then
          cld_present=.true.
          error_msg = this%optici(reffix(C,L),tlay(C,L),dens_i,ext_i,ssa_i,leg_i)
          ext_i = iwpx(C,L) * 0.75_wp * ext_i / (reffix(C,L) * dens_i)
        endif
        ext(C,L,:) = ext_l + ext_i
        if (cld_present) then
          ssa(C,L,:) = ((ssa_l * ext_l) + (ssa_i * ext_i)) / (ext_l + ext_i)
          do P=1,nmom
            leg(P,C,L,:) = ((ext_l * ssa_l * leg_l(P,:)) + (ext_i * ssa_i * leg_i(P,:))) / ((ssa_l * ext_l) + (ssa_i * ext_i))
          enddo
        endif
      enddo
    enddo  

    !fill optical_props
    optical_props%tau(:,:,:) = ext(:,:,:)
   !print *, "ext", ext
    select type(optical_props)
      type is (ty_optical_props_2str)
        optical_props%ssa(:,:,:) = ssa(:,:,:)
   !    print *, "ssa", ssa
   !    print *, "leg", leg(2,:,:,:)/3._wp
        optical_props%g(:,:,:) = leg(2,:,:,:)/3._wp
      type is (ty_optical_props_nstr)
        optical_props%ssa(:,:,:) = ssa(:,:,:)
        do P=1,nmom
          optical_props%p(P,:,:,:) = leg(P,:,:,:)/(2._wp*(real(P,wp)-1._wp)+1._wp) !right? convert P to float?
        enddo
    end select

  end function cloud_optics
  ! --------------------------------------------------------------------------------------
  function optici(this,                            &
                  reff,teff,ddens,qqext,ssalb,ssleg) result(error_msg) !need to pass temperature into module above
   !implicit none (if defined as subroutine?)
    class(ty_cloud_optics_fjx), intent(in) :: this
    real(wp), intent(in)                   :: reff, teff
    real(wp), intent(out)                  :: ddens
    real(wp), dimension(:), allocatable, intent(out)    :: qqext
    real(wp), dimension(:), allocatable, intent(out)    :: ssalb
    real(wp), dimension(:,:), allocatable, intent(out)  :: ssleg
    character(len=128)                     :: error_msg

    integer  :: I,J,K,L,NR,nsize,ngpt,nmom
    real(wp) :: FNR
    error_msg = ""

    nsize=size(this%QCC,dim=2)
    ngpt=size(this%QCC,dim=1)
    nmom=size(this%PCC,dim=1)

    allocate(qqext(ngpt))
    allocate(ssalb(ngpt))
    allocate(ssleg(nmom,ngpt))

    if (teff .ge. 233.15_wp) then
      K = 2  ! ice irreg (warm)
    else
      K = 3  ! ice hexag (cold)
    endif
    ddens = this%DCC(K)
    I = 1      !must have at least 2 Reff bins, interpolate in Reff
    do NR = 2,nsize-1
      if (reff .gt. this%RCC(NR,K)) then
        I = NR
      endif
    enddo
    FNR = (reff - this%RCC(I,K)) / (this%RCC(I+1,K) - this%RCC(I,K))
    FNR = min(1._wp, max(0._wp, FNR))

    ! each wavelength S-bins J has its own indexed optical properties
    do J=1,ngpt
      qqext(J) = this%QCC(J,I,K) + FNR*(this%QCC(J,I+1,K)-this%QCC(J,I,K))
      ssalb(J) = this%SCC(J,I,K) + FNR*(this%SCC(J,I+1,K)-this%SCC(J,I,K))
      do L=1,8
        ssleg(L,J) = this%PCC(L,J,I,K) + FNR*(this%PCC(L,J,I+1,K)-this%PCC(L,J,I,K))
      enddo
    enddo

  end function optici
  ! --------------------------------------------------------------------------------------
  function opticl(this,                             &
                  reff,ddens,qqext,ssalb,ssleg) result(error_msg) !temperature isn't actually used for liquid cloud
   !implicit none (if defined as subroutine?)
    class(ty_cloud_optics_fjx), intent(in) :: this
    real(wp), intent(in)                   :: reff
    real(wp), intent(out)                  :: ddens
   !real(wp), dimension(size(this%QCC,dim=1)), intent(out)    :: qqext !this doesn't work, not sure why
    real(wp), dimension(:), allocatable, intent(out)    :: qqext
   !real(wp), dimension(size(this%QCC,dim=1)), intent(out)    :: ssalb
    real(wp), dimension(:), allocatable, intent(out)    :: ssalb
   !real(wp), dimension(size(this%PCC,dim=1),size(this%QCC,dim=1)), intent(out)  :: ssleg
    real(wp), dimension(:,:), allocatable, intent(out)  :: ssleg
    character(len=128)                     :: error_msg

    integer  :: I,J,K,L,NR,nsize,ngpt,nmom
    real(wp) :: FNR
    error_msg = ""

    nsize=size(this%QCC,dim=2)
    ngpt=size(this%QCC,dim=1)
    nmom=size(this%PCC,dim=1)

    allocate(qqext(ngpt))
    allocate(ssalb(ngpt))
    allocate(ssleg(nmom,ngpt))

    K = 1   ! liquid water Mie clouds
    ddens = this%DCC(K)
    I = 1      !must have at least 2 Reff bins, interpolate in Reff
    do NR = 2,nsize-1
      if (reff .gt. this%RCC(NR,K)) then
        I = NR
      endif
    enddo
    FNR = (reff - this%RCC(I,K)) / (this%RCC(I+1,K) - this%RCC(I,K))
    FNR = min(1._wp, max(0._wp, FNR))

    !  each wavelength S-bins J has its own indexed optical properties
    do J=1,ngpt
      qqext(J) = this%QCC(J,I,K) + FNR*(this%QCC(J,I+1,K)-this%QCC(J,I,K))
      ssalb(J) = this%SCC(J,I,K) + FNR*(this%SCC(J,I+1,K)-this%SCC(J,I,K))
      do L=1,8
        ssleg(L,J) = this%PCC(L,J,I,K) + FNR*(this%PCC(L,J,I+1,K)-this%PCC(L,J,I,K))
      enddo
    enddo

  end function opticl
  ! --------------------------------------------------------------------------------------
  !
  !
  ! Initialization
  !
  !--------------------------------------------------------------------------------------------------------------------
  function load(this,ntype,nsize,ngpt,nmom,dcc,rcc,qcc,scc,pcc,wavebounds,gpt_lims) result(error_msg)
    class(ty_cloud_optics_fjx),   intent(inout) :: this
    integer,                      intent(in)    ::  nmom, ngpt, nsize, ntype
    real(wp), dimension(:) ,      intent(in) :: DCC
    real(wp), dimension(:,:),     intent(in) :: RCC
    real(wp), dimension(:,:,:),   intent(in) :: QCC, SCC  !(wavelength, size-based classification, particle type)
    real(wp), dimension(:,:,:,:), intent(in) :: PCC       !(Legendre moment, wavelength, size-based classification, particle type)
    real(wp), dimension(:,:), intent(in) :: wavebounds
    integer, dimension(:,:), intent(in) :: gpt_lims
    character(len=128)                       :: error_msg     !! Empty if successful

    error_msg = ""
    ! input sanitizing
    ! maybe also check in cloud_optics that ngpt matches that of optical_props passed in?
    if (check_extents) then
      if (.not. extents_are(dcc,ntype)) &
        error_msg="dcc incorrectly sized"
      if (.not. extents_are(rcc,nsize,ntype)) &
        error_msg="rcc incorrectly sized"
      if (.not. extents_are(qcc,ngpt,nsize,ntype)) &
        error_msg="qcc incorrectly sized"
      if (.not. extents_are(scc,ngpt,nsize,ntype)) &
        error_msg="scc incorrectly sized"
      if (.not. extents_are(pcc,nmom,ngpt,nsize,ntype)) &
        error_msg="pcc incorrectly sized"
      if(error_msg /= "") return
    endif

    allocate(this%DCC(ntype)) !didn't think this was necessary?
    this%DCC(:) = DCC
    allocate(this%RCC(nsize,ntype)) 
    this%RCC(:,:) = RCC
    allocate(this%QCC(ngpt,nsize,ntype))
    this%QCC(:,:,:) = QCC
    allocate(this%SCC(ngpt,nsize,ntype))
    this%SCC(:,:,:) = SCC
    allocate(this%PCC(nmom,ngpt,nsize,ntype))
    this%PCC(:,:,:,:) = PCC

    error_msg=this%init(wavebounds,gpt_lims) !?
    if(len_trim(error_msg) /= 0) return

  end function load

  !--------------------------------------------------------------------------------------------------------------------

end module mo_cloud_optics_fjx
