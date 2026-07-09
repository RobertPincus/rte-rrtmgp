module mo_gas_optics_ddq_kernels
  use mo_rte_kind,             only: wp, wl
  use mo_gas_optics_constants, only: boltzmann_k, lightspeed, planck_h

  implicit none
  private
  public :: tau_absorption_from_fits, add_tau_rayleigh
  integer, parameter, public :: fax_norder = 2, fax_nterms = 3, xsec_nterms = 3

contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Compute absorption optical depth from second-order polynomial approximations
  !    for absorption cross-section
  !
  subroutine tau_absorption_from_fits(ncol, nlay, nnu, ngas, &
                  nus, &
                  play, tlay, dry_num, vmrs, &
                  fax_ngas, fax_num_index, fax_a, fax_b, fax_T0, fax_c, fax_p0, fax_sigma0, fax_S, &
                  xsec_ngas, xsec_num_index, xsec_p, &
                  mtckd_ngas, mtckd_num_index, mtckd_cself, mtckd_cfrgn, mtckd_n, mtckd_T0, mtckd_p0, &
                  tau) bind(C, name="ddq_compute_tau_absorption")
    integer,  intent(in)  :: ncol, nlay, nnu, ngas
    real(wp), intent(in)  :: nus(nnu)
    real(wp), dimension(     ncol, nlay), &
              intent(in)  :: play, tlay, dry_num
! VMRs might need to start at ngas = 0, with element 0 having value 0, to stand in for gases users haven't supplied
!   otherwise we need to filter the
    real(wp), intent(in)  :: vmrs(0:ngas, ncol, nlay)

    ! Functional approximations to cross-sections
    integer,  intent(in) :: fax_ngas
    integer,  intent(in) :: fax_num_index(fax_ngas)
    real(wp), dimension(0:2, fax_ngas, nnu), &
              intent(in)  :: fax_a, fax_b
    real(wp), intent(in)  :: fax_c(0:3, fax_ngas, nnu)
    real(wp), intent(in)  :: fax_sigma0(fax_ngas, nnu)
    real(wp), dimension(fax_ngas) &
                          :: fax_S, fax_T0, fax_p0
    ! Cross-section fits
    integer,  intent(in) :: xsec_ngas
    integer,  intent(in) :: xsec_num_index(xsec_ngas)
    real(wp), intent(in) :: xsec_p(0:3, xsec_ngas, nnu)

    integer,  intent(in) :: mtckd_ngas
    integer,  intent(in) :: mtckd_num_index(mtckd_ngas)
    real(wp), dimension(mtckd_ngas, nnu), &
              intent(in) :: mtckd_cself, mtckd_cfrgn, mtckd_n
    real(wp), dimension(mtckd_ngas), &
              intent(in) :: mtckd_T0, mtckd_p0

    real(wp), intent(out) :: tau(ncol, nlay, nnu)
    ! -----------------
    integer  :: igas, icol, ilay, inu
    real(wp) :: vmr, num_density
    real(wp) :: t, x, P_scale, T_scale, pres, delta_T
    real(wp) :: cself, cfrgn, R ! MT_CKD

    do inu = 1, nnu
      do ilay = 1, nlay
        do icol = 1, ncol
          num_density = dry_num(icol, ilay)
          t = 0
          !
          ! Functional approximation to cross-sections
          !
          do igas = 1, fax_ngas
            vmr = vmrs(fax_num_index(igas), icol, ilay)
            ! Increase pressure to account for self-broadening
            pres = play(icol, ilay) * (1 + vmr * fax_S(igas))
            x = log(pres/fax_p0(igas))
            ! fax_c(3,:,:) is the hinge point x_h
            P_scale =                  &
                   fax_c(0, igas, inu) &
                +  fax_c(1, igas, inu) * x &
                + (fax_c(2, igas, inu) - fax_c(1, igas, inu)) &
                  * max(x - fax_c(3, igas, inu), 0._wp)
            delta_T = tlay(icol, ilay) - fax_T0(igas)
            T_scale =                            &
                (fax_a(0, igas, inu)             &
                +fax_a(1, igas, inu)*delta_T     &
                +fax_a(2, igas, inu)*delta_T**2) &
              / (fax_b(0, igas, inu)             &
                +fax_b(1, igas, inu)*delta_T     &
                +fax_b(2, igas, inu)*delta_T**2)
            t = t &
              + (fax_sigma0(igas, inu) * exp(P_scale + T_scale)) & ! cross-section [m**2/mol]
              * (vmr * num_density)                    ! Integrated number density [mol/m**2]
          end do
          !
          ! Cross-sections pressure and temperature dependence following doi:10.1029/2022MS003239
          !
          do igas = 1, xsec_ngas
            vmr = vmrs(xsec_num_index(igas), icol, ilay)
            t = t &
              + (xsec_p(0, igas, inu) &
                + xsec_p(1, igas, inu) * tlay(icol, ilay)    &
                + xsec_p(2, igas, inu) * tlay(icol, ilay)**2 &
                + xsec_p(3, igas, inu) * play(icol, ilay))    &
              * (vmr * num_density)
          end do
          !
          ! MT_CKD continuum
          !
          do igas = 1, mtckd_ngas
            vmr = vmrs(mtckd_num_index(igas), icol, ilay)
            cself = (mtckd_T0(igas)/tlay(icol, ilay))**(mtckd_n(igas, inu)) &
                  * (play(icol, ilay)/mtckd_p0(igas)) * vmr                 &
                  * mtckd_cself(igas, inu)
            cfrgn = (mtckd_T0(igas)/tlay(icol, ilay))                       &
                  * (play(icol, ilay)/mtckd_p0(igas)) * (1._wp - vmr)       &
                  * mtckd_cfrgn(igas, inu)
            ! nu supplied in kaysers (cm^-1); convert to MKS
            R = 100._wp * nus(inu) &
              * tanh((planck_h * lightspeed * 100._wp * nus(inu)) &
                     / (2._wp * boltzmann_k * tlay(icol, ilay)))
           t = t                   &
             + R * (cself + cfrgn) &
             * (vmr * num_density)
          end do
          tau(icol, ilay, inu) = t
        end do
      end do
    end do
  end subroutine tau_absorption_from_fits
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Compute absorption optical depth from second-order polynomial approximations
  !    for absorption cross-section
  !
  subroutine add_tau_rayleigh(ncol, nlay, nnu,  &
                              dry_num,          &
                              rayleigh_xsec,    &
                              tau, ssa) bind(C, name="ddq_add_tau_rayleigh")
    integer,  intent(in)    :: ncol, nlay, nnu
    real(wp), intent(in)    :: dry_num(ncol, nlay)
    real(wp), intent(in)    :: rayleigh_xsec(nnu)
    real(wp), intent(inout) :: tau(ncol, nlay, nnu), ssa(ncol, nlay, nnu)
    ! -----------------
    integer  :: icol, ilay, inu
    real(wp) :: t, t_r

    do inu = 1, nnu
      do ilay = 1, nlay
        do icol = 1, ncol
          t_r = dry_num(ncol, nlay) * rayleigh_xsec(inu)
          t = tau(icol, ilay, inu)
          tau(icol, ilay, inu) = t + t_r
          ssa(icol, ilay, inu) = t_r/(t + t_r)
        end do
      end do
    end do
  end subroutine add_tau_rayleigh
  !--------------------------------------------------------------------------------------------------------------------
end module mo_gas_optics_ddq_kernels
