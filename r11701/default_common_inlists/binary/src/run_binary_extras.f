! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
      module run_binary_extras

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use crlibm_lib
      use utils_lib

      implicit none

      ! Nova regime limits: set here once, used by in_nova_regime, the
      ! eta/ignition-mass fits and the nova mass-transfer gates below.
      real(dp), parameter :: nova_m_wd_min      = 0.25d0  ! Msun
      real(dp), parameter :: nova_m_wd_max      = 1.44d0  ! Msun
      real(dp), parameter :: nova_log10mdot_min = -18d0   ! log10(Msun/yr)

      contains

      ! True if WD mass [Msun] and log10(Mdot [Msun/yr]) fall inside the nova regime
      logical function in_nova_regime(m_wd_msun, log10mdot)
         real(dp), intent(in) :: m_wd_msun, log10mdot
         in_nova_regime = (m_wd_msun >= nova_m_wd_min .and. m_wd_msun <= nova_m_wd_max &
                           .and. log10mdot >= nova_log10mdot_min)
      end function in_nova_regime

      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
         b% warn_binary_extra =.false.

         b% other_sync_spin_to_orbit => my_sync_spin_to_orbit
         b% other_tsync => my_tsync
         b% other_mdot_edd => my_mdot_edd !pointer 
         b% other_rlo_mdot => my_rlo_mdot
         
         b% other_jdot_mb => mb_torque_selector
         b% other_jdot_ls => jdot_ls_with_mb

      end subroutine extras_binary_controls

      ! This is the MESA default routine except that 
      ! it also accounts for additional AML due to e.g., magnetic braking
      subroutine jdot_ls_with_mb(binary_id, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: delta_J, MOI
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         b% jdot_ls = 0
         ! ignore in first step, or if not doing rotation
         if (b% doing_first_model_of_run) &
            return

         ! bulk change in spin angular momentum takes tides into account
         delta_J = b% s_donor% total_angular_momentum_old - &
                   b% s_donor% total_angular_momentum

         ! ignore angular momentum lost through winds
         if (b% s_donor% mstar_dot < 0) &
            delta_J = delta_J - b% s_donor% angular_momentum_removed * &
            abs(b% mdot_system_wind(b% d_i) / b% s_donor% mstar_dot)

         ! Ignore angular momentum lost through magnetic braking
         ! (this assumes run_star_extras places MB AML in extra_omegadot)
         if (b% s_donor% extra_omegadot(1) < 0) then
            MOI = dot_product(b% s_donor% dm_bar(1:b% s_donor% nz), &
                              b% s_donor% i_rot(1:b% s_donor% nz))
            delta_J = delta_J + (b% s_donor% extra_omegadot(1) * MOI) &
                                 * b% s_donor% dt
         end if

         b% jdot_ls = b% jdot_ls + delta_J

         ! Repeat for accretor
         if (b% point_mass_i == 0) then
            delta_J = b% s_accretor% total_angular_momentum_old - &
                      b% s_accretor% total_angular_momentum

            if (b% s_accretor% mstar_dot < 0) then
               ! all AM lost via wind from the accretor is lost from the system
               delta_J = delta_J - b% s_accretor% angular_momentum_removed
            end if
            ! Ignore angular momentum lost through magnetic braking
            if (b% s_accretor% extra_omegadot(1) < 0) then
               MOI = dot_product(b% s_accretor% dm_bar(1:b% s_accretor% nz), &
                                 b% s_accretor% i_rot(1:b% s_accretor% nz))
               delta_J = delta_J + (b% s_accretor% extra_omegadot(1) * MOI) &
                                    * b% s_accretor% dt
            end if

            b% jdot_ls = b% jdot_ls + delta_J

         else if (b% model_twins_flag) then
            b% jdot_ls = b% jdot_ls + b% jdot_ls
         end if

         b% jdot_ls = b% jdot_ls / b% s_donor% dt

      end subroutine jdot_ls_with_mb

      subroutine my_tsync(id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
         integer, intent(in) :: id
         character (len=strlen), intent(in) :: sync_type !synchronization timescale
         real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ).
         real(dp), intent(in) :: qratio !mass_other_star/mass_this_star
         real(dp), intent(in) :: m
         real(dp), intent(in) :: r_phot
         real(dp), intent(in) :: osep ! orbital separation (cm)
         real(dp), intent(out) :: t_sync
         integer, intent(out) :: ierr
         real(dp) :: rGyr_squared , moment_of_inertia
         real(dp) :: one_div_t_sync_conv, one_div_t_sync_rad, one_div_t_sync ! t_sync_rad, t_sync_conv
         type (binary_info), pointer :: b
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
           write(*,*) 'failed in star_ptr'
           return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s%nz))
         rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))

         ! Implemented the option for both equilibrium and dynamical tides
         if (sync_type == "Hut_conv") then
            !sync_type .eq. "Hut_conv"!Convective envelope + Radiative core
            ! eq. (11) of Hut, P. 1981, A&A, 99, 126
            t_sync = 3.0d0*k_div_T(b, s,.true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
            ! invert it.
            !write(*,*) 'star id', s% id
            if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
               b% s1% xtra2 = 1d0/t_sync
               !write(*,*) 'two timescales ', b% s1% xtra1, b% s1% xtra2
            else if (b% point_mass_i /= 2 .and. b% s2% id == s% id) then
               b% s2% xtra2 = 1d0/t_sync
               !write(*,*) 'two timescales ', b% s2% xtra1, b% s2% xtra2
            else
               write(*,*) 'something is not going well with the stars IDs '
            end if

            t_sync = 1d0/t_sync
            !write(*,*) 'Hut_conv ', t_sync
         else if (sync_type == "Hut_rad") then
            !sync_type .eq. "Hut_rad"! Radiative envelope + convective core
            ! eq. (11) of Hut, P. 1981, A&A, 99, 126
            t_sync = 3.0*k_div_T(b, s,.false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
            ! invert it.
            !write(*,*) 'star id', s% id
            if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
               b% s1% xtra1 = 1d0/t_sync
               !write(*,*) 'two timescales ', b% s1% xtra1, b% s1% xtra2
            else if (b% point_mass_i /= 2 .and. b% s2% id == s% id) then
               b% s2% xtra1 = 1d0/t_sync
               !write(*,*) 'two timescales ', b% s2% xtra1, b% s2% xtra2
            else
               write(*,*) 'something is not going well with the stars IDs '
            end if

            t_sync = 1d0/t_sync
            !write(*,*) 'Hut_rad ', t_sync
         else if (sync_type == "structure_dependent") then !  Calculates both timescales from "Hut_rad" and "Hut_conv" and picks the shortest
            one_div_t_sync_conv = 3.0d0*k_div_T_posydon(b, s, .true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
            one_div_t_sync_rad = 3.0d0*k_div_T_posydon(b, s, .false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
            !write(*,*) 'star id', s% id
            if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
               b% s1% xtra1 = 1d0/one_div_t_sync_rad
               b% s1% xtra2 = 1d0/one_div_t_sync_conv
               !write(*,*) 'two timescales ', b% s1% xtra1, b% s1% xtra2
            else if (b% point_mass_i /= 2 .and. b% s2% id == s% id) then
               b% s2% xtra1 = 1d0/one_div_t_sync_rad
               b% s2% xtra2 = 1d0/one_div_t_sync_conv
               !write(*,*) 'two timescales ', b% s2% xtra1, b% s2% xtra2
            else
               write(*,*) 'something is not going well with the stars IDs '
            end if
            !write(*,*) 'two 1/timescales ', one_div_t_sync_conv , one_div_t_sync_rad
            !write(*,*) 'two timescales ', b% s1% ixtra1, b% s1% ixtra2
            one_div_t_sync = MAX(one_div_t_sync_conv,one_div_t_sync_rad)
            !one_div_t_sync = one_div_t_sync_conv1 + one_div_t_sync_conv2 + one_div_t_sync_rad ! if we want to combine them
            t_sync = 1d0/one_div_t_sync
            !write(*,*) 't_tides in years', t_sync / secyer
         else if (sync_type == "Orb_period") then ! sync on timescale of orbital period
            t_sync = b% period ! synchronize on timescale of orbital period
         else
            ierr = -1
            write(*,*) 'unrecognized sync_type', sync_type
            return
        end if
        t_sync = t_sync / Ftid
      end subroutine my_tsync

      subroutine get_tsync(id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
         integer, intent(in) :: id
         character (len=strlen), intent(in) :: sync_type ! synchronization timescale
         real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ).
         real(dp), intent(in) :: qratio ! mass_other_star/mass_this_star
         real(dp), intent(in) :: m
         real(dp), intent(in) :: r_phot
         real(dp), intent(in) :: osep ! orbital separation (cm)
         real(dp), intent(out) :: t_sync
         integer, intent(out) :: ierr
         real(dp) :: rGyr_squared, moment_of_inertia
         type (binary_info), pointer :: b
         type (star_info), pointer :: s

         include 'formats'

         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         ! calculate the gyration radius squared
         moment_of_inertia = dot_product(s% i_rot(:s% nz), s% dm_bar(:s% nz))
         rGyr_squared = (moment_of_inertia/(m*r_phot*r_phot))
         if (sync_type == "Hut_conv") then
            ! eq. (11) of Hut, P. 1981, A&A, 99, 126
            t_sync = 3.0*k_div_T(b, s, .true.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
            ! invert it.
            t_sync = 1d0/t_sync
         else if (sync_type == "Hut_rad") then
            ! eq. (11) of Hut, P. 1981, A&A, 99, 126
            t_sync = 3.0*k_div_T(b, s,.false.)*(qratio*qratio/rGyr_squared)*pow6(r_phot/osep)
            ! invert it.
            t_sync = 1d0/t_sync
         else if (sync_type == "Orb_period") then ! sync on timescale of orbital period
            t_sync = b% period ! synchronize on timescale of orbital period
         else
            ierr = -1
            write(*,*) 'unrecognized sync_type', sync_type
            return
         end if
         t_sync = t_sync / Ftid
      end subroutine get_tsync

      subroutine my_sync_spin_to_orbit(id, nz, osep, qratio, rl, dt_next, Ftid,sync_type, sync_mode, ierr)
         use const_def, only: dp, strlen
         integer, intent(in) :: id
         integer, intent(in) :: nz
         real(dp), intent(in) :: osep ! orbital separation (cm)
         real(dp), intent(in) :: qratio ! mass_other_star/mass_this_star
         real(dp), intent(in) :: rl ! roche lobe radius (cm)
         real(dp), intent(in) :: dt_next ! next timestep
         real(dp), intent(in) :: Ftid ! efficiency of tidal synchronization. (time scale / Ftid ).

         character (len=strlen), intent(in) :: sync_type ! synchronization timescale
         character (len=strlen), intent(in) :: sync_mode ! where to put/take angular momentum
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         type (binary_info), pointer :: b

         integer :: k
         real(dp), dimension(nz) :: j_sync, delta_j
         real(dp) :: t_sync, m, r_phot, omega_orb
         real(dp) :: a1,a2
         real(dp) :: omegadot_mb, dj_mb_k

         include 'formats'
         ierr = 0

         t_sync = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_ptr'
            return
         end if

         call binary_ptr(s% binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if

         omega_orb = 2d0*pi/b% period
         do k=1,nz
            j_sync(k) = omega_orb*s% i_rot(k)
         end do

         if (.not. b% use_other_tsync) then !Default tidal synchronization timescale calculation
            call get_tsync(s% id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
            if (ierr/=0) return
         else
            call b% other_tsync(s% id, sync_type, Ftid, qratio, m, r_phot, osep, t_sync, ierr)
            if (ierr/=0) return
         end if
         a1 = f2(b% eccentricity)
         a2 = pow_cr(1-pow2(b% eccentricity), 1.5d0)*f5(b% eccentricity)

         ! Option for tides to apply only to the envelope. (Qin et al. 2018 implementation)
         !if (.not. b% have_radiative_core(id)) then ! convective core
         !    !write(*,*) 'applying tides only in radiative envelope'
         !    do k=1,nz
         !       if (s% mixing_type(k) /= convective_mixing) then
         !           delta_j(k) = (1d0 - exp_cr(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
         !       else
         !           delta_j(k) = 0.0
         !       end if
         !    end do
         !else
         !    !write(*,*) 'applying tides only in convective regions'
         !    do k=1,nz
         !       if (s% mixing_type(k) == convective_mixing) then
         !           delta_j(k) = (1d0 - exp_cr(-a2*dt_next/t_sync))*(s% j_rot(k) - a1/a2*j_sync(k))
         !       else
         !           delta_j(k) = 0.0
         !       end if
         !    end do
         !end if

         ! Tides apply in all layers
         ! write(*,*) 'applying tides in all layers'
         do k=1,nz
            ! correct for spin down due to magnetic braking:
            omegadot_mb = s% extra_omegadot(k)
            dj_mb_k = s% extra_omegadot(k) * s% i_rot(k) * dt_next
            delta_j(k) = (1d0 - exp_cr(-a2*dt_next/t_sync))*&
                         ((s% j_rot(k) + dj_mb_k) - a1/a2*j_sync(k)) 
         end do


         if (b% point_mass_i /= 1 .and. b% s1% id == s% id) then
            b% t_sync_1 = t_sync
         else
            b% t_sync_2 = t_sync
         end if

         if (.not. b% doing_first_model_of_run) then
            do k=1,nz
               s% extra_jdot(k) = s% extra_jdot(k) - delta_j(k)/dt_next
            end do
         end if
      end subroutine my_sync_spin_to_orbit

      real(dp) function f2(e)
         real(dp), intent(in) :: e

         f2 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f2 after eq. 11
         if (e > 0d0) then
            f2 = 1d0 + 15d0/2d0 * pow2(e) + 45d0/8d0 * pow4(e) + 5d0/16d0 * pow6(e)
         end if

      end function f2

      real(dp) function f3(e)
         real(dp), intent(in) :: e

         f3 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f3 after eq. 11
         if (e > 0d0) then
            f3 = 1d0 + 15d0/4d0*pow2(e) + 15d0/8d0 * pow4(e) + 5d0/64d0 * pow6(e)
         end if

      end function f3


      real(dp) function f4(e)
         real(dp), intent(in) :: e

         f4 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f4 after eq. 11
         if (e > 0d0) then
            f4 = 1d0 + 3d0/2d0 * pow2(e) + 1d0/8d0 * pow4(e)
         end if

      end function f4


      real(dp) function f5(e)
         real(dp), intent(in) :: e

         f5 = 1d0

         ! Hut 1981, A&A, 99, 126, definition of f5 after eq. 11
         if (e > 0d0) then
            f5 = 1d0 + 3d0*pow2(e) + 3d0/8d0 * pow4(e)
         end if

      end function f5

      real(dp) function mass_conv_core(s)
         type (star_info), pointer :: s
         integer :: j, nz, k
         real(dp) :: dm_limit
         include 'formats'
         mass_conv_core = 0.0d0
         dm_limit = s% conv_core_gap_dq_limit*s% xmstar
         nz = s% nz
         do j = 1, s% n_conv_regions
            ! ignore possible small gap at center
            if (s% cz_bot_mass(j) <= s% m(nz) + dm_limit) then
               mass_conv_core = s% cz_top_mass(j)/Msun
               ! jump over small gaps
               do k = j+1, s% n_conv_regions
                  if (s% cz_bot_mass(k) - s% cz_top_mass(k-1) >= dm_limit) exit
                  mass_conv_core = s% cz_top_mass(k)/Msun
               end do
               exit
            end if
         end do
      end function mass_conv_core


      real(dp) function k_div_T(b, s, has_convective_envelope)
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         logical, intent(in) :: has_convective_envelope

         integer :: k,i, h1
         real(dp) osep, qratio, m, r_phot,porb, m_env, r_env, tau_conv, P_tid, f_conv,E2, Xs

         ! k/T computed as in Hurley, J., Tout, C., Pols, O. 2002, MNRAS, 329, 897
         ! Kudos to Francesca Valsecchi for help implementing and testing this

         k_div_T = 0d0

         osep = b% separation
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            qratio = 1.0d0/qratio
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if
         porb = b% period

         if (has_convective_envelope) then
            m_env = 0d0
            r_env = 0d0
            do k=1, s% nz
               if (s% mixing_type(k) /= convective_mixing .and. &
                    s% rho(k) > 1d5*s% rho(1)) then
                  r_env = (r_phot - s% r(k))/Rsun
                  m_env = (s% m(1) - s% m(k))/Msun
                  exit
               end if
            end do
            tau_conv = 0.431d0*pow_cr(m_env*r_env* &
                (r_phot/Rsun-r_env/2d0)/3d0/s% L_phot,one_third) * secyer
            P_tid = 1d0/abs(1d0/porb-s% omega_avg_surf/(2d0*pi))
            f_conv = min(1.0d0, pow_cr(P_tid/(2d0*tau_conv),b% tidal_reduction))

            k_div_T = 2d0/21d0*f_conv/tau_conv*m_env/(m/Msun)
         else ! radiative envelope
            ! New fitting E2 (Qin et al. 2018)
            do i = s% nz, 1, -1
               if (s% brunt_N2(i) >= 0) exit
            end do
            !write(*,*) i
            h1 = s% net_iso(ih1)
            Xs = s% xa(h1,1)
            ! E2 is different for H-rich and He stars (Qin et al. 2018)
            if (Xs < 0.4d0) then ! HeStar
               E2 = exp10_cr(-0.93_dp)*pow_cr(s% r(i)/r_phot,6.7_dp)! HeStars
            else
               E2 = exp10_cr(-0.42_dp)*pow_cr(s% r(i)/r_phot,7.5_dp)! H-rich stars
            !write(*,*) E2, s% r(i)
            end if
            if (isnan(E2)) then  !maybe this won't be used.
               k_div_T = 1d-20
            else
               k_div_T = sqrt(standard_cgrav*m*r_phot*r_phot/pow5(osep)/(Msun/pow3(Rsun)))
               k_div_T = k_div_T*pow_cr(1d0+qratio,5d0/6d0)
               k_div_T = k_div_T * E2
            end if
         end if

         ! protect against NaN, this can happen near sync
         if (isnan(k_div_T)) k_div_T = 1d-99

      end function k_div_T

      subroutine loop_conv_layers(s, n_conv_regions_posydon, n_zones_of_region, bot_bdy, top_bdy, &
                                     cz_bot_mass_posydon, cz_bot_radius_posydon, &
                                     cz_top_mass_posydon, cz_top_radius_posydon)
         type (star_info), pointer :: s
         ! integer, intent(out) :: ierr

         logical :: in_convective_region
         integer :: k, j, nz
         logical, parameter :: dbg = .false.
         integer, intent(out) :: n_conv_regions_posydon
         !integer :: max_num_mixing_regions
         !max_num_mixing_regions = 100
         !integer, intent(out), dimension (:), allocatable :: n_zones_of_region, bot_bdy, top_bdy
         !real(dp),intent(out), dimension (:), allocatable :: cz_bot_mass_posydon, cz_bot_radius_posydon
         !real(dp),intent(out), dimension (:), allocatable :: cz_top_mass_posydon, cz_top_radius_posydon
         integer :: min_zones_for_convective_tides
         integer ::  pot_n_zones_of_region, pot_bot_bdy, pot_top_bdy
         real(dp) :: pot_cz_bot_mass_posydon, pot_cz_bot_radius_posydon
         integer, intent(out), dimension (max_num_mixing_regions) :: n_zones_of_region, bot_bdy, top_bdy
         real(dp),intent(out), dimension (max_num_mixing_regions) :: cz_bot_mass_posydon
         real(dp),intent(out) :: cz_bot_radius_posydon(max_num_mixing_regions)
         real(dp),intent(out), dimension (max_num_mixing_regions) :: cz_top_mass_posydon, cz_top_radius_posydon

         include 'formats'
         !ierr = 0
         min_zones_for_convective_tides = 10
         nz = s% nz
         n_zones_of_region(:)=0
         bot_bdy(:)=0
         top_bdy(:)=0
         cz_bot_mass_posydon(:)=0.0d0
         cz_bot_radius_posydon(:)=0.0d0
         cz_top_mass_posydon(:)=0.0d0
         cz_top_radius_posydon(:)=0.0d0
         n_conv_regions_posydon = 0
         pot_cz_bot_mass_posydon = 0.0d0
         pot_cz_bot_radius_posydon = 0.0d0
         pot_bot_bdy = 0
         pot_n_zones_of_region = 0

         in_convective_region = (s% mixing_type(nz) == convective_mixing)
         if (in_convective_region) then
            pot_cz_bot_mass_posydon = s% M_center
            pot_cz_bot_radius_posydon = 0.0d0
            pot_bot_bdy = nz
         end if

         !write(*,*) 'initial in_convective_region', in_convective_region
         do k=nz-1, 2, -1
            if (in_convective_region) then
               if (s% mixing_type(k) /= convective_mixing) then ! top of convective region
                  pot_top_bdy = k
                  pot_n_zones_of_region = pot_bot_bdy - pot_top_bdy
                  if (pot_n_zones_of_region >= min_zones_for_convective_tides) then
                     if (n_conv_regions_posydon < max_num_mixing_regions) then
                        n_conv_regions_posydon = n_conv_regions_posydon + 1
                     end if
                     cz_top_mass_posydon(n_conv_regions_posydon) = &
                        s% M_center + (s% q(k) - s% cz_bdy_dq(k))*s% xmstar
                     cz_bot_mass_posydon(n_conv_regions_posydon) = pot_cz_bot_mass_posydon
                     cz_top_radius_posydon(n_conv_regions_posydon) = s% r(k)/Rsun
                     cz_bot_radius_posydon(n_conv_regions_posydon) = pot_cz_bot_radius_posydon
                     top_bdy(n_conv_regions_posydon) = pot_top_bdy
                     bot_bdy(n_conv_regions_posydon) = pot_bot_bdy
                     n_zones_of_region(n_conv_regions_posydon) = pot_n_zones_of_region
                  end if
                  in_convective_region = .false.
               end if
            else if (s% mixing_type(k) == convective_mixing) then ! bottom of convective region
               pot_cz_bot_mass_posydon = s% M_center + (s% q(k) - s% cz_bdy_dq(k))*s% xmstar
               pot_cz_bot_radius_posydon = s% r(k)/Rsun
               pot_bot_bdy = k
               in_convective_region = .true.
            end if
         end do
         if (in_convective_region) then
            pot_top_bdy = 1
            pot_n_zones_of_region = pot_bot_bdy - pot_top_bdy
            if (pot_n_zones_of_region >= min_zones_for_convective_tides) then
               if (n_conv_regions_posydon < max_num_mixing_regions) then
                  n_conv_regions_posydon = n_conv_regions_posydon + 1
               end if
               cz_top_mass_posydon(n_conv_regions_posydon) = s% mstar
               cz_top_radius_posydon(n_conv_regions_posydon) = s% r(1)/Rsun
               top_bdy(n_conv_regions_posydon) = 1
               cz_bot_mass_posydon(n_conv_regions_posydon) = pot_cz_bot_mass_posydon
               cz_bot_radius_posydon(n_conv_regions_posydon) = pot_cz_bot_radius_posydon
               bot_bdy(n_conv_regions_posydon) = pot_bot_bdy
               n_zones_of_region(n_conv_regions_posydon) = pot_n_zones_of_region
           end if
         end if

         !write(*,*)
         !write(*,2) 'set_mixing_info n_conv_regions_posydon', n_conv_regions_posydon
         !do j = 1, n_conv_regions_posydon
         !   write(*,2) 'conv region', j, cz_bot_mass_posydon(j)/Msun, cz_top_mass_posydon(j)/Msun
         !   write(*,2) 'conv region', j, cz_bot_radius_posydon(j), cz_top_radius_posydon(j)
         !end do
         !write(*,*)
      end subroutine loop_conv_layers

      real(dp) function k_div_T_posydon(b, s, conv_layer_calculation)
         type(binary_info), pointer :: b
         type(star_info), pointer :: s
         !logical, intent(in) :: has_convective_envelope
         logical, intent(in) :: conv_layer_calculation

         integer :: k, i, l, h1, top_bound_zone, bot_bound_zone
         real(dp) :: osep, qratio, m, r_phot,porb, m_env, Dr_env, Renv_middle, tau_conv, P_tid, f_conv,E2, Xs, m_conv_core
         real(dp) :: k_div_T_posydon_new, conv_mx_top, conv_mx_bot, conv_mx_top_r, conv_mx_bot_r ,omega_conv_region,r_top, r_bottom
         integer :: n_conv_regions_posydon
         integer, dimension (max_num_mixing_regions) :: n_zones_of_region, bot_bdy, top_bdy
         real(dp), dimension (max_num_mixing_regions) :: cz_bot_mass_posydon
         real(dp) :: cz_bot_radius_posydon(max_num_mixing_regions)
         real(dp), dimension (max_num_mixing_regions) :: cz_top_mass_posydon, cz_top_radius_posydon
	      real(dp) :: dm, dmsum, omega_k, omega_l, omega_sum

         ! k/T computed as in Hurley, J., Tout, C., Pols, O. 2002, MNRAS, 329, 897
         ! Kudos to Francesca Valsecchi for help implementing and testing this

         k_div_T_posydon = 0d0

         osep = b% separation
         qratio = b% m(b% a_i) / b% m(b% d_i)
         if (is_donor(b, s)) then
            m = b% m(b% d_i)
            r_phot = b% r(b% d_i)
         else
            qratio = 1.0d0/qratio
            m = b% m(b% a_i)
            r_phot = b% r(b% a_i)
         end if
         porb = b% period

         if (conv_layer_calculation) then
            m_conv_core = mass_conv_core(s)
            ! In POSYDON the calculation is done for the most important convective layer, found below
            n_zones_of_region(:)=0
            bot_bdy(:)=0
            top_bdy(:)=0
            cz_bot_mass_posydon(:)=0.0d0
            cz_bot_radius_posydon(:)=0.0d0
            cz_top_mass_posydon(:)=0.0d0
            cz_top_radius_posydon(:)=0.0d0
            n_conv_regions_posydon = 0

            call loop_conv_layers(s, n_conv_regions_posydon, n_zones_of_region, bot_bdy, top_bdy, &
                                  cz_bot_mass_posydon, cz_bot_radius_posydon, &
                                  cz_top_mass_posydon, cz_top_radius_posydon)

            if (n_conv_regions_posydon > 0) then
               do k=1, n_conv_regions_posydon ! from inside out
                  m_env = 0.0d0
                  Dr_env = 0.0d0
                  Renv_middle = 0.0d0

		            ! calculate mass averaged omega for convective region
                  ! this is a more stable quantity than simply omega within a zone
                  omega_sum = 0
                  dmsum = 0
                  do l = top_bdy(k), bot_bdy(k) - 1
                     dm = s% dm(l)
                     dmsum = dmsum + dm
                     omega_l = 0.5d0*(s% omega(l) + s% omega(l+1))
                     omega_sum = omega_sum + dm*omega_l
                  end do
                  omega_k = omega_sum / dmsum
  
                  if ((cz_bot_mass_posydon(k) / Msun) >=  m_conv_core) then ! if the conv. region is not inside the conv. core
                     m_env = (cz_top_mass_posydon(k) - cz_bot_mass_posydon(k)) / Msun
                     Dr_env = cz_top_radius_posydon(k) - cz_bot_radius_posydon(k)  ! depth of the convective layer, length of the eddie
                     ! Corresponding to the Renv term in eq.31 of Hurley et al. 2002
                     ! and to (R-Renv) term in eq. 4 of Rasio et al. 1996  (different notation)
                     Renv_middle = (cz_top_radius_posydon(k) + cz_bot_radius_posydon(k) )*0.5d0  ! middle of the convective layer
                     ! Corresponding to the (R-0.5d0*Renv) in eq.31 of Hurley et al 2002
                     ! and to the Renv in eq. 4 of Rasio et al. 1996
                     ! where it represented the base of the convective layer (different notation)
                     tau_conv = 0.431_dp*pow_cr(m_env*Dr_env* &
                                 Renv_middle/3d0/s% L_phot,1.0d0/3.0d0) * secyer
                     P_tid = 1d0/abs(1d0/porb-omega_k/(2d0*pi))
                     f_conv = min(1.0d0, pow_cr(P_tid/(2d0*tau_conv), b% tidal_reduction))
                     !write(*,'(g0)') 'porb, p_from_omega, f_conv = ', porb, &
                     !1 / (s% omega(top_bdy(k))/(2d0*pi)), &
                     !1 /(s% omega_avg_surf/(2d0*pi)), f_conv
                     k_div_T_posydon_new = 2d0/21d0*f_conv/tau_conv*m_env/(m/Msun)
                     !write(*,'(g0)') 'tau_conv, K/T = ', tau_conv, k_div_T_posydon_new, m_env, (m/Msun)
                     if (k_div_T_posydon_new >= k_div_T_posydon) then
                        k_div_T_posydon = k_div_T_posydon_new
                        !write(*,'(g0)') 'M_env, DR_env, Renv_middle, omega_conv_region in conv region ', k ,' is ', &
                        ! m_env, Dr_env, Renv_middle, s% omega(top_bdy(k)), 'spanning number of zones = ', &
                        !top_bdy(k) , bot_bdy(k), &
                        !n_zones_of_region(k)
                     end if
                  end if
               end do
            end if
         else ! assuming a radiative star
            ! New fitting E2 (Qin et al. 2018)
            do i = s% nz, 2, -1
               if (s% brunt_N2(i) >= 0d0) exit
            end do
            !write(*,*) i
			   ! if fully convective or no convective core
	         if ( (s% r(i) / r_phot < 0.01 * r_phot) .or. &
                 (s% r(i) >= 0.99 * r_phot) ) then
	     	      E2 = 1d-99
		      ! has convective core and not fully convective
	         else
	     	      h1 = s% net_iso(ih1)
		         Xs = s% xa(h1,1)
	     	      ! E2 is different for H-rich and He stars (Qin et al. 2018)
	     	      if (Xs < 0.4d0) then ! HeStar
		            E2 = exp10_cr(-0.93_dp)*pow_cr(s% r(i)/r_phot, 6.7_dp)! HeStars
	     	      else
		            E2 = exp10_cr(-0.42_dp)*pow_cr(s% r(i)/r_phot, 7.5_dp)! H-rich stars
	     	         !write(*,*) E2, s% r(i)
	     	      end if
	         end if

            if (isnan(E2)) then  !maybe this won't be used.
               k_div_T_posydon = 1d-99
            else
               k_div_T_posydon = sqrt(standard_cgrav*m*r_phot*r_phot/pow5(osep)/(Msun/pow3(Rsun)))
               k_div_T_posydon = k_div_T_posydon*pow_cr(1d0+qratio,5d0/6d0)
               k_div_T_posydon = k_div_T_posydon * E2
            end if
         end if
		  
         ! protect against NaN, this can happen near sync
         if (isnan(k_div_T_posydon)) k_div_T_posydon = 1d-99

      end function k_div_T_posydon

      real(dp) function acc_radius(b, m_acc) !Calculates Sch. radius of compact object (or surface radius in case of NS) in cm
         type(binary_info), pointer :: b
         real(dp) :: m_acc, a
         real(dp) :: r_isco, Z1, Z2, eq_initial_bh_mass

         if (m_acc/Msun <= 2.50d0) then ! NS
            !Radius refernces for NS:
            ! 1) Miller, M. C., Lamb, F. K., Dittmann, A. J., et al. 2019, ApJL, 887, L2
            ! 2) Riley, T. E., Watts, A. L., Bogdanov, S., et al., 2019, ApJL, 887, L21
            ! 3) Landry, P., Essick, R., & Chatziioannou, K. 2020
            ! 4) E.R. Most, L.R. Weih, L. Rezzolla and J. Schaffner-Bielich, 2018, Phys. Rev. Lett. 120, 261103
            ! 5) Abbott, B. P., Abbott, R., Abbott, T. D., et al. 2020, ApJL, 892, L3
            acc_radius = 12.5E5_dp !* 10 ** 5 !in cm
         else ! Event horizon for Kerr-BH
            ! this part is only relevant for BH accretors
            if (b% initial_bh_spin < 0d0) then
               b% initial_bh_spin = 0d0
               write(*,*) "initial_bh_spin is smaller than zero. It has been set to zero."
            else if (b% initial_bh_spin > 1d0) then
               b% initial_bh_spin = 1d0
               write(*,*) "initial_bh_spin is larger than one. It has been set to one."
            end if
            ! compute isco radius from eq. 2.21 of Bardeen et al. (1972), ApJ, 178, 347
            Z1 = 1d0 + pow_cr(1d0 - pow2(b% initial_bh_spin),one_third) * &
                    (pow_cr(1d0 + b% initial_bh_spin,one_third) + &
                     pow_cr(1d0 - b% initial_bh_spin,one_third))
            Z2 = sqrt(3d0*pow2(b% initial_bh_spin) + pow2(Z1))
            r_isco = 3d0 + Z2 - sqrt((3d0 - Z1)*(3d0 + Z1 + 2d0*Z2))
            ! compute equivalent mass at zero spin from eq. (3+1/2) (ie. the equation between (3) and (4))
            ! of Bardeen (1970), Nature, 226, 65, taking values with subscript zero to correspond to
            ! zero spin (r_isco = sqrt(6)).

            ! If it was already a BH then take the initial mass m2
	         if (initial_mass(2) > 2.5_dp) then
		         eq_initial_bh_mass = b% eq_initial_bh_mass
            ! If it was initially a NS then take 2.5Msun as eq_initial_mass
	         else if (initial_mass(2) <= 2.5_dp) then 
	            eq_initial_bh_mass = 2.5_dp * Msun * sqrt(r_isco/6d0)
	         end if

            a = sqrt(two_thirds) &
                  * (eq_initial_bh_mass/min(b% m(b% point_mass_i), sqrt(6d0)* eq_initial_bh_mass)) &
                  * (4._dp - sqrt(18._dp*pow2(eq_initial_bh_mass / &
                  min(b% m(b% point_mass_i),sqrt(6d0)* eq_initial_bh_mass)) - 2._dp))
            ! Podsiadlowski et al. (2003) assuming a initially non-rotating BH
            acc_radius = (1.0_dp + sqrt(1.0_dp - a*a)) * b% s_donor% cgrav(1) * m_acc/pow2(clight)
         end if
      end function acc_radius

      !! Eddington accreton limits for NS and BH
      subroutine my_mdot_edd(binary_id, mdot_edd, ierr)
         use const_def, only: dp
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot_edd
         integer, intent(out) :: ierr
         real(dp) :: mdot_edd_eta
         real(dp) :: r_isco, Z1, Z2, eq_initial_bh_mass
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         !call my_nova_flag()
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         if (b% m(2)/Msun > 2.50_dp) then ! M2 > 2.5 Msol for BHs
            ! this part is only relevant for BH accretors
            if (b% initial_bh_spin < 0d0) then
               b% initial_bh_spin = 0d0
               write(*,*) "initial_bh_spin is smaller than zero. It has been set to zero."
            else if (b% initial_bh_spin > 1d0) then
               b% initial_bh_spin = 1d0
               write(*,*) "initial_bh_spin is larger than one. It has been set to one."
            end if
            ! compute isco radius from eq. 2.21 of Bardeen et al. (1972), ApJ, 178, 347
            Z1 = 1d0 + pow_cr(1d0 - pow2(b% initial_bh_spin),one_third) * &
                 (pow_cr(1d0 + b% initial_bh_spin,one_third) + &
                  pow_cr(1d0 - b% initial_bh_spin,one_third))
            Z2 = sqrt(3d0*pow2(b% initial_bh_spin) + pow2(Z1))
            r_isco = 3d0 + Z2 - sqrt((3d0 - Z1)*(3d0 + Z1 + 2d0*Z2))
            ! compute equivalent mass at zero spin from eq. (3+1/2) (ie. the equation between (3) and (4))
            ! of Bardeen (1970), Nature, 226, 65, taking values with subscript zero to correspond to
            ! zero spin (r_isco = sqrt(6)).

            if (initial_mass(2) > 2.5_dp) then ! If it was already a BH then take the initial mass m2
               eq_initial_bh_mass = b% eq_initial_bh_mass
            else if (initial_mass(2) <= 2.5_dp) then! If it was initially a NS then take 2.5 as eq_initial_mass
               eq_initial_bh_mass = 2.5_dp * Msun * sqrt(r_isco/6d0)
            end if

            !! mdot_edd_eta for BH following Podsiadlowski, Rappaport & Han (2003), MNRAS, 341, 385
            mdot_edd_eta = 1d0 - sqrt(1d0 - &
                  pow2(min(b% m(b% a_i),sqrt(6d0)*eq_initial_bh_mass)/(3d0*eq_initial_bh_mass)))
         else ! NS
             !! mdot_edd_eta for NS accretors
             mdot_edd_eta = b% s_donor% cgrav(1) * b% m(2) / (pow2(clight) * acc_radius(b, b% m(2)))
         end if
         mdot_edd = 4d0*pi*b% s_donor% cgrav(1)*b% m(b% a_i) &
                  /(clight*0.2d0*(1d0+b% s_donor% surface_h1)*mdot_edd_eta)
         !b% s1% x_ctrl(1) used to adjust the Eddington limit in inlist1
         mdot_edd = mdot_edd * b% s1% x_ctrl(1)
      !if (my_nova_flag() = 1) then 
         !call my_dmHee
         !call my_accretion_efficiency
         !if (my_dmHee = ____) then
            !call my_m_loss
         
      end subroutine my_mdot_edd

!!**********************************************!!

      !!Here begins Novae Perscriptions from Dean Kousiounelos and Seth Gossage
      !!my_m_loss
      subroutine my_m_loss(binary_id, mloss, ierr)
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mloss !mass ejected per nova cycle [Msun]
         integer, intent(out) :: ierr
         real(dp) :: dmHe, eta
         type (binary_info), pointer :: b

         ierr = 0
         mloss = 0d0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /=0) then
         write(*,*) 'failed in binary_ptr in my_m_loss'
         return
      end if

         call my_dm_ig(binary_id, dmHe, ierr) !ignition mass
         if (ierr /= 0) return

         call my_eta_accretion(binary_id, eta, ierr) !retention efficiency
         if (ierr /= 0) return

         ! Ejected mass per nova: Delta_Mej = Delta_Mig * (1 - eta_H).
         mloss = dmHe * (1d0 - eta) !mass ejected per nova cycle

      end subroutine my_m_loss



      !my_eta_accretion
      subroutine my_eta_accretion(binary_id, etaH, ierr)
         ! Returns helium-nova retention efficiency eta_He based on Wu+2017 fits.
         ! Uses raw mtransfer_rate (no xfer_fraction) — the total He flux onto
         ! the WD surface sets eta.

         ! Inputs: WD mass (Msun), log10 accretion rate (Msun/yr). Output: eta_He in [0,1].
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: etaH
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: m_wd, log10mdot
         real(dp) :: fitlowereta, fithighereta
         real(dp) :: M_WDlower, M_WDhigher

         ierr         = 0
         etaH         = 0d0
         fitlowereta  = 0d0
         fithighereta = 0d0
         M_WDlower    = 0d0
         M_WDhigher   = 0d0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr in my_eta_accretion'
            return
         end if

         m_wd      = b% m(b% point_mass_i) / Msun
         ! Wu+2017 fits are clamped to the nova regime WD mass and Mdot limits.
         m_wd      = min(max(m_wd, nova_m_wd_min), nova_m_wd_max)
         log10mdot = log10_cr(max(10d0**nova_log10mdot_min, abs(b% mtransfer_rate) * secyer / Msun))

         if ((m_wd >= 0.2d0) .and. (m_wd <= 0.7d0)) then

            ! 0.2 lower curve
            if (log10mdot > -7.421d0) then
               fitlowereta = 128.69653622d0*pow4(log10mdot) + 3755.62839546d0*pow3(log10mdot) &
                           + 41097.77848236d0*pow2(log10mdot) + 199876.15570843d0*log10mdot + 364523.39925814d0
            else
               fitlowereta = 0.0284d0*log10mdot + 0.4944d0
            end if

            ! 0.7 higher curve
            if (log10mdot > -6.921d0) then
               fithighereta = 128.69653622d0*pow4(log10mdot) + 3498.23532301d0*pow3(log10mdot) &
                            + 35657.38069353d0*pow2(log10mdot) + 161530.75025466d0*log10mdot + 274398.30600877d0
            else
               fithighereta = 0.0271d0*log10mdot + 0.4214d0
            end if

            M_WDlower  = 0.2d0
            M_WDhigher = 0.7d0

         else if ((m_wd > 0.7d0) .and. (m_wd <= 0.8d0)) then

            ! 0.7 lower curve
            if (log10mdot > -6.921d0) then
               fitlowereta = 128.69653622d0*pow4(log10mdot) + 3498.23532301d0*pow3(log10mdot) &
                           + 35657.38069353d0*pow2(log10mdot) + 161530.75025466d0*log10mdot + 274398.30600877d0
            else
               fitlowereta = 0.0271d0*log10mdot + 0.4214d0
            end if

            ! 0.8 higher curve
            if (log10mdot > -6.495d0) then
               fithighereta = -2.81727160d0*pow4(log10mdot) + 842.92305137d0*pow3(log10mdot) &
                            + 17025.57779911d0*pow2(log10mdot) + 111404.51221858d0*log10mdot + 241315.71447017d0
            else
               fithighereta = 0.0376d0*log10mdot + 0.5774d0
            end if

            M_WDlower  = 0.7d0
            M_WDhigher = 0.8d0

         else if ((m_wd > 0.8d0) .and. (m_wd <= 0.9d0)) then

            ! 0.8 lower curve
            if (log10mdot > -6.495d0) then
               fitlowereta = -2.81727160d0*pow4(log10mdot) + 842.92305137d0*pow3(log10mdot) &
                           + 17025.57779911d0*pow2(log10mdot) + 111404.51221858d0*log10mdot + 241315.71447017d0
            else
               fitlowereta = 0.0376d0*log10mdot + 0.5774d0
            end if

            ! 0.9 higher curve
            if (log10mdot > -6.585d0) then
               fithighereta = 371.49410540d0*pow4(log10mdot) + 9601.91634884d0*pow3(log10mdot) &
                            + 93061.48887580d0*pow2(log10mdot) + 400844.20977850d0*log10mdot + 647425.88293061d0
            else
               fithighereta = 0.0328d0*log10mdot + 0.4598d0
            end if

            M_WDlower  = 0.8d0
            M_WDhigher = 0.9d0

         else if ((m_wd > 0.9d0) .and. (m_wd <= 1.0d0)) then

            ! 0.9 lower curve
            if (log10mdot > -6.585d0) then
               fitlowereta = 371.49410540d0*pow4(log10mdot) + 9601.91634884d0*pow3(log10mdot) &
                           + 93061.48887580d0*pow2(log10mdot) + 400844.20977850d0*log10mdot + 647425.88293061d0
            else
               fitlowereta = 0.0328d0*log10mdot + 0.4598d0
            end if

            ! 1.0 higher curve
            if (log10mdot > -6.155d0) then
               fithighereta = 8.3115d0*log10mdot + 51.674d0
            else if (log10mdot > -7.046d0) then
               fithighereta = 0.42950525d0*pow2(log10mdot) + 6.08106227d0*log10mdot + 21.66204702d0
            else
               fithighereta = 0.0274d0*log10mdot + 0.331d0
            end if

            M_WDlower  = 0.9d0
            M_WDhigher = 1.0d0

         else if ((m_wd > 1.0d0) .and. (m_wd <= 1.1d0)) then

            ! 1.0 lower curve
            if (log10mdot > -6.155d0) then
               fitlowereta = 8.3115d0*log10mdot + 51.674d0
            else if (log10mdot > -7.046d0) then
               fitlowereta = 0.42950525d0*pow2(log10mdot) + 6.08106227d0*log10mdot + 21.66204702d0
            else
               fitlowereta = 0.0274d0*log10mdot + 0.331d0
            end if

            ! 1.1 higher curve
            if (log10mdot > -5.921d0) then
               fithighereta = 8.3856d0*log10mdot + 50.358d0
            else if (log10mdot > -7.097d0) then
               fithighereta = -0.19075473d0*pow4(log10mdot) - 4.92869779d0*pow3(log10mdot) &
                            - 47.22636728d0*pow2(log10mdot) - 198.28033421d0*log10mdot - 306.30153647d0
            else
               fithighereta = 0.0275d0*log10mdot + 0.2987d0
            end if

            M_WDlower  = 1.0d0
            M_WDhigher = 1.1d0

         else if ((m_wd > 1.1d0) .and. (m_wd <= 1.2d0)) then

            ! 1.1 lower curve
            if (log10mdot > -5.921d0) then
               fitlowereta = 8.3856d0*log10mdot + 50.358d0
            else if (log10mdot > -7.097d0) then
               fitlowereta = -0.19075473d0*pow4(log10mdot) - 4.92869779d0*pow3(log10mdot) &
                           - 47.22636728d0*pow2(log10mdot) - 198.28033421d0*log10mdot - 306.30153647d0
            else
               fitlowereta = 0.0275d0*log10mdot + 0.2987d0
            end if

            ! 1.2 higher curve
            if (log10mdot > -5.854d0) then
               fithighereta = 6.7216d0*log10mdot + 40.146d0
            else if (log10mdot > -7.222d0) then
               fithighereta = 0.61101079d0*pow4(log10mdot) + 16.16260323d0*pow3(log10mdot) &
                            + 160.39740390d0*pow2(log10mdot) + 708.15701427d0*log10mdot + 1174.46712731d0
            else
               fithighereta = 0.0278d0*log10mdot + 0.2994d0
            end if

            M_WDlower  = 1.1d0
            M_WDhigher = 1.2d0

         else if ((m_wd > 1.2d0) .and. (m_wd <= 1.25d0)) then

            ! 1.2 lower curve
            if (log10mdot > -5.854d0) then
               fitlowereta = 6.7216d0*log10mdot + 40.146d0
            else if (log10mdot > -7.222d0) then
               fitlowereta = 0.61101079d0*pow4(log10mdot) + 16.16260323d0*pow3(log10mdot) &
                           + 160.39740390d0*pow2(log10mdot) + 708.15701427d0*log10mdot + 1174.46712731d0
            else
               fitlowereta = 0.0278d0*log10mdot + 0.2994d0
            end if

            ! 1.25 higher curve
            if (log10mdot > -5.824d0) then
               fithighereta = 7.3403d0*log10mdot + 43.604d0
            else if (log10mdot > -7.222d0) then
               fithighereta = 0.70445443d0*pow4(log10mdot) + 18.61070131d0*pow3(log10mdot) &
                            + 184.34484390d0*pow2(log10mdot) + 811.82891151d0*log10mdot + 1342.09332478d0
            else
               fithighereta = 0.0278d0*log10mdot + 0.2924d0
            end if

            M_WDlower  = 1.2d0
            M_WDhigher = 1.25d0

         else if ((m_wd > 1.25d0) .and. (m_wd <= 1.3d0)) then

            ! 1.25 lower curve
            if (log10mdot > -5.824d0) then
               fitlowereta = 7.3403d0*log10mdot + 43.604d0
            else if (log10mdot > -7.222d0) then
               fitlowereta = 0.70445443d0*pow4(log10mdot) + 18.61070131d0*pow3(log10mdot) &
                           + 184.34484390d0*pow2(log10mdot) + 811.82891151d0*log10mdot + 1342.09332478d0
            else
               fitlowereta = 0.0278d0*log10mdot + 0.2924d0
            end if

            ! 1.3 higher curve
            if (log10mdot > -5.824d0) then
               fithighereta = 5.2446d0*log10mdot + 31.397d0
            else if (log10mdot > -7.222d0) then
               fithighereta = 0.35511452d0*pow4(log10mdot) + 9.43073062d0*pow3(log10mdot) &
                            + 93.97568200d0*pow2(log10mdot) + 416.89622203d0*log10mdot + 695.71187519d0
            else
               fithighereta = 0.0278d0*log10mdot + 0.2754d0
            end if

            M_WDlower  = 1.25d0
            M_WDhigher = 1.3d0

         else if ((m_wd > 1.3d0) .and. (m_wd <= 1.35d0)) then

            ! 1.3 lower curve
            if (log10mdot > -5.824d0) then
               fitlowereta = 5.2446d0*log10mdot + 31.397d0
            else if (log10mdot > -7.222d0) then
               fitlowereta = 0.35511452d0*pow4(log10mdot) + 9.43073062d0*pow3(log10mdot) &
                           + 93.97568200d0*pow2(log10mdot) + 416.89622203d0*log10mdot + 695.71187519d0
            else
               fitlowereta = 0.0278d0*log10mdot + 0.2754d0
            end if

            ! 1.35 higher curve
            if (log10mdot > -5.824d0) then
               fithighereta = 4.742d0*log10mdot + 28.497d0
            else if (log10mdot > -7.301d0) then
               fithighereta = 0.17534165d0*pow2(log10mdot) + 2.81395280d0*log10mdot + 11.31042978d0
            else
               fithighereta = 0.028d0*log10mdot + 0.315d0
            end if

            M_WDlower  = 1.3d0
            M_WDhigher = 1.35d0

         else
            ! m_wd > 1.35 up to the clamped ceiling of 1.44 Msun

            ! 1.35 lower curve
            if (log10mdot > -5.824d0) then
               fitlowereta = 4.742d0*log10mdot + 28.497d0
            else if (log10mdot > -7.301d0) then
               fitlowereta = 0.17534165d0*pow2(log10mdot) + 2.81395280d0*log10mdot + 11.31042978d0
            else
               fitlowereta = 0.028d0*log10mdot + 0.315d0
            end if

            ! 1.44 higher curve
            if (log10mdot > -5.824d0) then
               fithighereta = 4.742d0*log10mdot + 28.597d0
            else if (log10mdot > -7.301d0) then
               fithighereta = 0.17534165d0*pow2(log10mdot) + 2.81395280d0*log10mdot + 11.41042978d0
            else
               fithighereta = 0.028d0*log10mdot + 0.415d0
            end if

            M_WDlower  = 1.35d0
            M_WDhigher = 1.44d0

         end if

         ! Linear interp: equivalent to np.interp(m_wd, [M_WDlower, M_WDhigher], [fitlowereta, fithighereta])
         etaH = fitlowereta + (m_wd - M_WDlower) * (fithighereta - fitlowereta) / (M_WDhigher - M_WDlower)

         if (etaH < 0d0) etaH = 0d0
         if (etaH > 1d0) etaH = 1d0
      end subroutine my_eta_accretion

!**********************************************!

      subroutine my_dmHe(binary_id, dmHe, ierr)
      !from formulae in Kato et al. 2018 (M_WD>=1.0) and Piersanti et al. 2014
      !(M_WD<1.0), as tabulated in binary_c's Henova_dMHe.c (KEMP_NOVAE module)
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: dmHe !ignition mass [Msun]
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: m_wd, log10mdot
         real(dp) :: fitlower, fithigher
         real(dp) :: M_WDlower, M_WDhigher

         ierr       = 0
         dmHe       = 0d0
         fitlower   = 0d0
         fithigher  = 0d0
         M_WDlower  = 0d0
         M_WDhigher = 0d0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /=0) then
            write(*,*) 'failed in binary_ptr in my_dmHe'
            return

         end if

         m_wd = b% m(b% point_mass_i) / Msun
         ! Kato/Piersanti fits are clamped to the nova regime WD mass and Mdot limits.
         m_wd = min(max(m_wd, nova_m_wd_min), nova_m_wd_max)
         log10mdot = log10_cr(max(10d0**nova_log10mdot_min, abs(b% mtransfer_rate) * secyer / Msun)) !convert to log10

         if ((m_wd >= 0.2d0) .and. (m_wd <= 0.6d0)) then

            ! 0.2 lower curve: returns the constant 0.2 verbatim (NOT in
            ! log10 space, unlike every other node here) - kept exactly as
            ! in the reference; after exp10_cr below this node maps to
            ! ~1.585 Msun rather than 0.2 Msun.
            fitlower = 0.2d0

            ! 0.6 higher curve
            if (log10mdot > -7.39d0) then
               fithigher = -1.2658d0*log10mdot - 10.79d0
            else
               fithigher = -2.3314d0*log10mdot - 18.668d0
            end if

            M_WDlower  = 0.2d0
            M_WDhigher = 0.6d0

         else if ((m_wd > 0.6d0) .and. (m_wd <= 0.7d0)) then

            ! 0.6 lower curve
            if (log10mdot > -7.39d0) then
               fitlower = -1.2658d0*log10mdot - 10.79d0
            else
               fitlower = -2.3314d0*log10mdot - 18.668d0
            end if

            ! 0.7 higher curve
            if (log10mdot > -7.39d0) then
               fithigher = -0.9928d0*log10mdot - 8.8399d0
            else
               fithigher = -2.8393d0*log10mdot - 22.535d0
            end if

            M_WDlower  = 0.6d0
            M_WDhigher = 0.7d0

         else if ((m_wd > 0.7d0) .and. (m_wd <= 0.81d0)) then

            ! 0.7 lower curve
            if (log10mdot > -7.39d0) then
               fitlower = -0.9928d0*log10mdot - 8.8399d0
            else
               fitlower = -2.8393d0*log10mdot - 22.535d0
            end if

            ! 0.81 higher curve
            if (log10mdot > -7.40d0) then
               fithigher = -1.2462d0*log10mdot - 10.729d0
            else
               fithigher = -2.7234d0*log10mdot - 21.695d0
            end if

            M_WDlower  = 0.7d0
            M_WDhigher = 0.81d0

         else if ((m_wd > 0.81d0) .and. (m_wd <= 0.92d0)) then

            ! 0.81 lower curve
            if (log10mdot > -7.40d0) then
               fitlower = -1.2462d0*log10mdot - 10.729d0
            else
               fitlower = -2.7234d0*log10mdot - 21.695d0
            end if

            ! 0.92 higher curve
            if (log10mdot > -7.30d0) then
               fithigher = -0.8228d0*log10mdot - 7.8311d0
            else
               fithigher = -2.0714d0*log10mdot - 16.929d0
            end if

            M_WDlower  = 0.81d0
            M_WDhigher = 0.92d0

         else if ((m_wd > 0.92d0) .and. (m_wd <= 1.02d0)) then

            ! 0.92 lower curve
            if (log10mdot > -7.30d0) then
               fitlower = -0.8228d0*log10mdot - 7.8311d0
            else
               fitlower = -2.0714d0*log10mdot - 16.929d0
            end if

            ! 1.02 higher curve
            if (log10mdot > -7.09d0) then
               fithigher = -0.9181d0*log10mdot - 8.6514d0
            else
               fithigher = -2.0406d0*log10mdot - 16.597d0
            end if

            M_WDlower  = 0.92d0
            M_WDhigher = 1.02d0

         else if ((m_wd > 1.02d0) .and. (m_wd <= 1.10d0)) then

            ! 1.02 lower curve
            if (log10mdot > -7.09d0) then
               fitlower = -0.9181d0*log10mdot - 8.6514d0
            else
               fitlower = -2.0406d0*log10mdot - 16.597d0
            end if

            ! 1.10 higher curve
            if (log10mdot < -6.0d0) then
               fithigher = 1.72662751d0*pow4(log10mdot) + 45.88762486d0*pow3(log10mdot) &
                         + 456.63345197d0*pow2(log10mdot) + 2015.69281770d0*log10mdot + 3325.89047580d0
            else
               fithigher = -0.0937d0*log10mdot - 4.0218d0
            end if

            M_WDlower  = 1.02d0
            M_WDhigher = 1.10d0

         else if ((m_wd > 1.10d0) .and. (m_wd <= 1.20d0)) then

            ! 1.10 lower curve
            if (log10mdot < -6.0d0) then
               fitlower = 1.72662751d0*pow4(log10mdot) + 45.88762486d0*pow3(log10mdot) &
                        + 456.63345197d0*pow2(log10mdot) + 2015.69281770d0*log10mdot + 3325.89047580d0
            else
               fitlower = -0.0937d0*log10mdot - 4.0218d0
            end if

            ! 1.20 higher curve
            if (log10mdot < -6.0d0) then
               fithigher = 0.91755524d0*pow4(log10mdot) + 24.10100929d0*pow3(log10mdot) &
                         + 236.98730163d0*pow2(log10mdot) + 1033.14600808d0*log10mdot + 1680.17435446d0
            else
               fithigher = -0.0542d0*log10mdot - 4.1937d0
            end if

            M_WDlower  = 1.10d0
            M_WDhigher = 1.20d0

         else if ((m_wd > 1.20d0) .and. (m_wd <= 1.25d0)) then

            ! 1.20 lower curve
            if (log10mdot < -6.0d0) then
               fitlower = 0.91755524d0*pow4(log10mdot) + 24.10100929d0*pow3(log10mdot) &
                        + 236.98730163d0*pow2(log10mdot) + 1033.14600808d0*log10mdot + 1680.17435446d0
            else
               fitlower = -0.0542d0*log10mdot - 4.1937d0
            end if

            ! 1.25 higher curve
            if (log10mdot < -6.0d0) then
               fithigher = 1.68244304d0*pow4(log10mdot) + 45.15739606d0*pow3(log10mdot) &
                         + 453.76725009d0*pow2(log10mdot) + 2022.35971770d0*log10mdot + 3368.09376231d0
            else
               fithigher = -0.1274d0*log10mdot - 4.7732d0
            end if

            M_WDlower  = 1.20d0
            M_WDhigher = 1.25d0

         else if ((m_wd > 1.25d0) .and. (m_wd <= 1.30d0)) then

            ! 1.25 lower curve
            if (log10mdot < -6.0d0) then
               fitlower = 1.68244304d0*pow4(log10mdot) + 45.15739606d0*pow3(log10mdot) &
                        + 453.76725009d0*pow2(log10mdot) + 2022.35971770d0*log10mdot + 3368.09376231d0
            else
               fitlower = -0.1274d0*log10mdot - 4.7732d0
            end if

            ! 1.30 higher curve
            if (log10mdot < -6.0d0) then
               fithigher = 1.16410162d0*pow4(log10mdot) + 30.86542279d0*pow3(log10mdot) &
                         + 306.21947654d0*pow2(log10mdot) + 1346.49353679d0*log10mdot + 2208.95849065d0
            else
               fithigher = -0.1574d0*log10mdot - 5.3083d0
            end if

            M_WDlower  = 1.25d0
            M_WDhigher = 1.30d0

         else if ((m_wd > 1.30d0) .and. (m_wd <= 1.35d0)) then

            ! 1.30 lower curve
            if (log10mdot < -6.0d0) then
               fitlower = 1.16410162d0*pow4(log10mdot) + 30.86542279d0*pow3(log10mdot) &
                        + 306.21947654d0*pow2(log10mdot) + 1346.49353679d0*log10mdot + 2208.95849065d0
            else
               fitlower = -0.1574d0*log10mdot - 5.3083d0
            end if

            ! 1.35 higher curve
            if (log10mdot < -6.0d0) then
               fithigher = 0.47596972d0*pow4(log10mdot) + 12.61002433d0*pow3(log10mdot) &
                         + 125.08958464d0*pow2(log10mdot) + 549.88858261d0*log10mdot + 898.29715596d0
            else
               fithigher = -0.0934d0*log10mdot - 5.3203d0
            end if

            M_WDlower  = 1.30d0
            M_WDhigher = 1.35d0

         else if ((m_wd > 1.35d0) .and. (m_wd <= 1.38d0)) then

            ! 1.35 lower curve
            if (log10mdot < -6.0d0) then
               fitlower = 0.47596972d0*pow4(log10mdot) + 12.61002433d0*pow3(log10mdot) &
                        + 125.08958464d0*pow2(log10mdot) + 549.88858261d0*log10mdot + 898.29715596d0
            else
               fitlower = -0.0934d0*log10mdot - 5.3203d0
            end if

            ! 1.38 higher curve
            if (log10mdot < -6.0d0) then
               fithigher = 0.54437094d0*pow4(log10mdot) + 14.68324636d0*pow3(log10mdot) &
                         + 148.24381872d0*pow2(log10mdot) + 663.22092056d0*log10mdot + 1103.57417218d0
            else
               fithigher = -0.1056d0*log10mdot - 5.6752d0
            end if

            M_WDlower  = 1.35d0
            M_WDhigher = 1.38d0

         else
            ! m_wd > 1.38 up to the clamped ceiling of 1.44 Msun

            ! 1.38 lower curve
            if (log10mdot < -6.0d0) then
               fitlower = 0.54437094d0*pow4(log10mdot) + 14.68324636d0*pow3(log10mdot) &
                        + 148.24381872d0*pow2(log10mdot) + 663.22092056d0*log10mdot + 1103.57417218d0
            else
               fitlower = -0.1056d0*log10mdot - 5.6752d0
            end if

            ! 1.44 higher curve: M_Heig_chand = log10(8e-7) (previously a
            ! natural-log bug in the reference, fixed upstream to log10)
            fithigher = log10_cr(8d-7)

            M_WDlower  = 1.38d0
            M_WDhigher = 1.44d0

         end if

         ! Linear interp of log10(Delta_M_ig) in mass, then back to linear
         ! space (matches the reference's interp_lin(...) followed by exp10()).
         dmHe = exp10_cr(fitlower + (m_wd - M_WDlower) * (fithigher - fitlower) / (M_WDhigher - M_WDlower))

      end subroutine my_dmHe

!**************************************************!

      subroutine my_Hedetonation(binary_id, dmHedet, ierr)
      !He ignition mass in the low-flash/detonation regime, from formulae in
      !Kato et al. 2018 (M_WD>=1.0) and Piersanti et al. 2014 (M_WD<1.0), as
      !tabulated in binary_c's Hedet_dMHe.c (KEMP_NOVAE module).
      !Shares the same curves as my_dmHe in a band of Mdot for each WD mass;
      !only differs below that band (where it flattens to ~1e-3-0.7 Msun).
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: dmHedet !ignition mass [Msun]
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: m_wd, log10mdot
         real(dp) :: fitlower, fithigher
         real(dp) :: M_WDlower, M_WDhigher

         ierr       = 0
         dmHedet    = 0d0
         fitlower   = 0d0
         fithigher  = 0d0
         M_WDlower  = 0d0
         M_WDhigher = 0d0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /=0) then
            write(*,*) 'failed in binary_ptr in my_Hedetonation'
            return

         end if

         m_wd = b% m(b% point_mass_i) / Msun
         ! Kato/Piersanti fits are clamped to the nova regime WD mass and Mdot limits
         ! (equivalent to np.interp holding endpoint values in the reference).
         m_wd = min(max(m_wd, nova_m_wd_min), nova_m_wd_max)
         log10mdot = log10_cr(max(10d0**nova_log10mdot_min, abs(b% mtransfer_rate) * secyer / Msun)) !convert to log10

         if ((m_wd >= 0.2d0) .and. (m_wd <= 0.6d0)) then

            ! 0.2 lower curve: constant 0.2 verbatim (phantom point, same as my_dmHe)
            fitlower = 0.2d0

            ! 0.6 higher curve
            if (log10mdot >= -7.82d0) then
               fithigher = -2.3314d0*log10mdot - 18.668d0
            else if (log10mdot >= -9d0) then
               fithigher = -0.13105692d0*pow3(log10mdot) - 3.35859429d0*pow2(log10mdot) &
                         - 28.84546767d0*log10mdot - 83.28360749d0
            else
               fithigher = -0.0303d0*log10mdot - 0.4544d0
            end if

            M_WDlower  = 0.2d0
            M_WDhigher = 0.6d0

         else if ((m_wd > 0.6d0) .and. (m_wd <= 0.7d0)) then

            ! 0.6 lower curve
            if (log10mdot >= -7.82d0) then
               fitlower = -2.3314d0*log10mdot - 18.668d0
            else if (log10mdot >= -9d0) then
               fitlower = -0.13105692d0*pow3(log10mdot) - 3.35859429d0*pow2(log10mdot) &
                        - 28.84546767d0*log10mdot - 83.28360749d0
            else
               fitlower = -0.0303d0*log10mdot - 0.4544d0
            end if

            ! 0.7 higher curve
            if (log10mdot >= -7.70d0) then
               fithigher = -2.8393d0*log10mdot - 22.535d0
            else if (log10mdot >= -9d0) then
               fithigher = -0.68502223d0*pow4(log10mdot) - 23.27090657d0*pow3(log10mdot) &
                         - 296.28435839d0*pow2(log10mdot) - 1675.83253716d0*log10mdot - 3553.76035429d0
            else
               fithigher = -0.0251d0*log10mdot - 0.4658d0
            end if

            M_WDlower  = 0.6d0
            M_WDhigher = 0.7d0

         else if ((m_wd > 0.7d0) .and. (m_wd <= 0.81d0)) then

            ! 0.7 lower curve
            if (log10mdot >= -7.70d0) then
               fitlower = -2.8393d0*log10mdot - 22.535d0
            else if (log10mdot >= -9d0) then
               fitlower = -0.68502223d0*pow4(log10mdot) - 23.27090657d0*pow3(log10mdot) &
                        - 296.28435839d0*pow2(log10mdot) - 1675.83253716d0*log10mdot - 3553.76035429d0
            else
               fitlower = -0.0251d0*log10mdot - 0.4658d0
            end if

            ! 0.81 higher curve
            if (log10mdot >= -7.60d0) then
               fithigher = -2.7234d0*log10mdot - 21.695d0
            else if (log10mdot >= -8.82d0) then
               fithigher = -0.89021703d0*pow3(log10mdot) - 22.38723555d0*pow2(log10mdot) &
                         - 187.76315351d0*log10mdot - 525.64755404d0
            else
               fithigher = -0.0265d0*log10mdot - 0.5781d0
            end if

            M_WDlower  = 0.7d0
            M_WDhigher = 0.81d0

         else if ((m_wd > 0.81d0) .and. (m_wd <= 0.92d0)) then

            ! 0.81 lower curve
            if (log10mdot >= -7.60d0) then
               fitlower = -2.7234d0*log10mdot - 21.695d0
            else if (log10mdot >= -8.82d0) then
               fitlower = -0.89021703d0*pow3(log10mdot) - 22.38723555d0*pow2(log10mdot) &
                        - 187.76315351d0*log10mdot - 525.64755404d0
            else
               fitlower = -0.0265d0*log10mdot - 0.5781d0
            end if

            ! 0.92 higher curve
            if (log10mdot >= -7.52d0) then
               fithigher = -2.0714d0*log10mdot - 16.929d0
            else if (log10mdot >= -8.82d0) then
               fithigher = -2.22457174d0*pow4(log10mdot) - 74.06815881d0*pow3(log10mdot) &
                         - 924.25967441d0*pow2(log10mdot) - 5123.20522589d0*log10mdot - 10644.68164309d0
            else
               fithigher = -0.0238d0*log10mdot - 0.6265d0
            end if

            M_WDlower  = 0.81d0
            M_WDhigher = 0.92d0

         else if ((m_wd > 0.92d0) .and. (m_wd <= 1.02d0)) then

            ! 0.92 lower curve
            if (log10mdot >= -7.52d0) then
               fitlower = -2.0714d0*log10mdot - 16.929d0
            else if (log10mdot >= -8.82d0) then
               fitlower = -2.22457174d0*pow4(log10mdot) - 74.06815881d0*pow3(log10mdot) &
                        - 924.25967441d0*pow2(log10mdot) - 5123.20522589d0*log10mdot - 10644.68164309d0
            else
               fitlower = -0.0238d0*log10mdot - 0.6265d0
            end if

            ! 1.02 higher curve
            if (log10mdot >= -7.30d0) then
               fithigher = -2.0406d0*log10mdot - 16.597d0
            else if (log10mdot >= -8.82d0) then
               fithigher = -0.52114281d0*pow3(log10mdot) - 13.26176747d0*pow2(log10mdot) &
                         - 112.71207541d0*log10mdot - 320.54975268d0
            else
               fithigher = -0.0282d0*log10mdot - 0.7734d0
            end if

            M_WDlower  = 0.92d0
            M_WDhigher = 1.02d0

         else if ((m_wd > 1.02d0) .and. (m_wd < 1.10d0) .and. (log10mdot >= -6.0d0)) then
            ! The reference's 1.02-1.10 bridge only covers log10mdot < -6; above
            ! that no branch matches and it falls through to 10**0 = 1 Msun.
            ! Kept verbatim so results match the reference exactly.

            fitlower   = 0d0
            fithigher  = 0d0
            M_WDlower  = 1.02d0
            M_WDhigher = 1.10d0

         else if ((m_wd > 1.02d0) .and. (m_wd <= 1.10d0)) then
            ! Bridge: 1.10 from Kato 2018, 1.02 from Piersanti 2014

            ! 1.02 lower curve
            if (log10mdot >= -7.30d0) then
               fitlower = -2.0406d0*log10mdot - 16.597d0
            else if (log10mdot >= -8.82d0) then
               fitlower = -0.52114281d0*pow3(log10mdot) - 13.26176747d0*pow2(log10mdot) &
                        - 112.71207541d0*log10mdot - 320.54975268d0
            else
               fitlower = -0.0282d0*log10mdot - 0.7734d0
            end if

            ! 1.10 higher curve
            if (log10mdot >= -7.56d0) then
               fithigher = 1.72662751d0*pow4(log10mdot) + 45.88762486d0*pow3(log10mdot) &
                         + 456.63345197d0*pow2(log10mdot) + 2015.69281770d0*log10mdot + 3325.89047580d0
            else
               fithigher = -0.0731d0*log10mdot - 2.1163d0
            end if

            M_WDlower  = 1.02d0
            M_WDhigher = 1.10d0

         else if ((m_wd > 1.10d0) .and. (m_wd <= 1.20d0)) then

            ! 1.10 lower curve
            if (log10mdot >= -7.56d0) then
               fitlower = 1.72662751d0*pow4(log10mdot) + 45.88762486d0*pow3(log10mdot) &
                        + 456.63345197d0*pow2(log10mdot) + 2015.69281770d0*log10mdot + 3325.89047580d0
            else
               fitlower = -0.0731d0*log10mdot - 2.1163d0
            end if

            ! 1.20 higher curve
            if (log10mdot >= -7.67d0) then
               fithigher = 0.91755524d0*pow4(log10mdot) + 24.10100929d0*pow3(log10mdot) &
                         + 236.98730163d0*pow2(log10mdot) + 1033.14600808d0*log10mdot + 1680.17435446d0
            else
               fithigher = -0.0646d0*log10mdot - 2.1626d0
            end if

            M_WDlower  = 1.10d0
            M_WDhigher = 1.20d0

         else if ((m_wd > 1.20d0) .and. (m_wd <= 1.25d0)) then

            ! 1.20 lower curve
            if (log10mdot >= -7.67d0) then
               fitlower = 0.91755524d0*pow4(log10mdot) + 24.10100929d0*pow3(log10mdot) &
                        + 236.98730163d0*pow2(log10mdot) + 1033.14600808d0*log10mdot + 1680.17435446d0
            else
               fitlower = -0.0646d0*log10mdot - 2.1626d0
            end if

            ! 1.25 higher curve
            if (log10mdot >= -7.62d0) then
               fithigher = 1.68244304d0*pow4(log10mdot) + 45.15739606d0*pow3(log10mdot) &
                         + 453.76725009d0*pow2(log10mdot) + 2022.35971770d0*log10mdot + 3368.09376231d0
            else
               fithigher = -0.0957d0*log10mdot - 2.8734d0
            end if

            M_WDlower  = 1.20d0
            M_WDhigher = 1.25d0

         else if ((m_wd > 1.25d0) .and. (m_wd <= 1.30d0)) then

            ! 1.25 lower curve
            if (log10mdot >= -7.62d0) then
               fitlower = 1.68244304d0*pow4(log10mdot) + 45.15739606d0*pow3(log10mdot) &
                        + 453.76725009d0*pow2(log10mdot) + 2022.35971770d0*log10mdot + 3368.09376231d0
            else
               fitlower = -0.0957d0*log10mdot - 2.8734d0
            end if

            ! 1.30 higher curve
            if (log10mdot >= -7.6920763d0) then
               fithigher = 1.16410162d0*pow4(log10mdot) + 30.86542279d0*pow3(log10mdot) &
                         + 306.21947654d0*pow2(log10mdot) + 1346.49353679d0*log10mdot + 2208.95849065d0
            else
               fithigher = -0.0888d0*log10mdot - 2.8977d0
            end if

            M_WDlower  = 1.25d0
            M_WDhigher = 1.30d0

         else if ((m_wd > 1.30d0) .and. (m_wd <= 1.35d0)) then

            ! 1.30 lower curve
            if (log10mdot >= -7.6920763d0) then
               fitlower = 1.16410162d0*pow4(log10mdot) + 30.86542279d0*pow3(log10mdot) &
                        + 306.21947654d0*pow2(log10mdot) + 1346.49353679d0*log10mdot + 2208.95849065d0
            else
               fitlower = -0.0888d0*log10mdot - 2.8977d0
            end if

            ! 1.35 higher curve
            if (log10mdot >= -7.88d0) then
               fithigher = 0.47596972d0*pow4(log10mdot) + 12.61002433d0*pow3(log10mdot) &
                         + 125.08958464d0*pow2(log10mdot) + 549.88858261d0*log10mdot + 898.29715596d0
            else
               fithigher = -0.0999d0*log10mdot - 3.1982d0
            end if

            M_WDlower  = 1.30d0
            M_WDhigher = 1.35d0

         else if ((m_wd > 1.35d0) .and. (m_wd <= 1.38d0)) then

            ! 1.35 lower curve
            if (log10mdot >= -7.88d0) then
               fitlower = 0.47596972d0*pow4(log10mdot) + 12.61002433d0*pow3(log10mdot) &
                        + 125.08958464d0*pow2(log10mdot) + 549.88858261d0*log10mdot + 898.29715596d0
            else
               fitlower = -0.0999d0*log10mdot - 3.1982d0
            end if

            ! 1.38 higher curve (note strict >, as in the reference)
            if (log10mdot > -7.94d0) then
               fithigher = 0.54437094d0*pow4(log10mdot) + 14.68324636d0*pow3(log10mdot) &
                         + 148.24381872d0*pow2(log10mdot) + 663.22092056d0*log10mdot + 1103.57417218d0
            else
               fithigher = -0.1433d0*log10mdot - 4.079d0
            end if

            M_WDlower  = 1.35d0
            M_WDhigher = 1.38d0

         else
            ! m_wd > 1.38 up to the clamped ceiling of 1.44 Msun

            ! 1.38 lower curve
            if (log10mdot > -7.94d0) then
               fitlower = 0.54437094d0*pow4(log10mdot) + 14.68324636d0*pow3(log10mdot) &
                        + 148.24381872d0*pow2(log10mdot) + 663.22092056d0*log10mdot + 1103.57417218d0
            else
               fitlower = -0.1433d0*log10mdot - 4.079d0
            end if

            ! 1.44 higher curve: M_Heig_chand = log10(8e-7)
            fithigher = log10_cr(8d-7)

            M_WDlower  = 1.38d0
            M_WDhigher = 1.44d0

         end if

         ! Linear interp of log10(Delta_M_ig) in mass, then back to linear space.
         dmHedet = exp10_cr(fitlower + (m_wd - M_WDlower) * (fithigher - fitlower) / (M_WDhigher - M_WDlower))

      end subroutine my_Hedetonation

!**************************************************!

      subroutine my_dm_ig(binary_id, dm_ig, ierr)
      !He ignition mass [Msun] for the unsteady-burning regime, switching
      !between the strong-flash (nova) fit and the low-flash/detonation fit
      !at the strong-flash/detonation boundary used by binary_c (KEMP_NOVAE,
      !src/binary_star_functions/He_accretion_regime_WD.c, SFDt_log10mdot,
      !"boundary curves from Piersanti et al. 2014"):
      !   log10(Mdot_det) = -8.313 + 1.018*M_WD   [Msun/yr]
      !As in binary_c's limit_accretion_rates.c, Hedet_dMHe replaces
      !Henova_dMHe below this line. The two fits agree to <~0.01 dex on this
      !boundary for 0.6 <= M_WD <= 1.38, and each extrapolates unphysically
      !on the other side (my_dmHe -> up to ~1e100+ Msun at low Mdot), so the
      !switch must be on Mdot, never min/max.
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: dm_ig
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: m_wd, log10mdot

         ierr  = 0
         dm_ig = 0d0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr in my_dm_ig'
            return
         end if

         m_wd      = b% m(b% point_mass_i) / Msun
         log10mdot = log10_cr(max(10d0**nova_log10mdot_min, abs(b% mtransfer_rate) * secyer / Msun))

         if (log10mdot < -8.313d0 + 1.018d0*m_wd) then
            call my_Hedetonation(binary_id, dm_ig, ierr)
         else
            call my_dmHe(binary_id, dm_ig, ierr)
         end if

      end subroutine my_dm_ig

!**************************************************!

      subroutine my_t_recurrence(binary_id, t_rec, ierr)
      !calculates the recurrence time between novae given ignition mass and mdot
      !note that this recurrence time is only calculated when the accreted mass
      !is greater than the ignition mass at the start of the MESA step. If the step
      !is much longer than the recurrence time, then one canot get information about
      !the other nova which would've occured during that timestep.
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: t_rec
         integer, intent(out) :: ierr
         real(dp) :: dmHe, mdot_msun_yr
         type (binary_info), pointer :: b

         ierr  = 0
         t_rec = 0d0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr in my_t_recurrence'
            return
         end if

         if (b% mtransfer_rate == 0d0) return

         mdot_msun_yr = abs(b% mtransfer_rate) * secyer / Msun
         call my_dm_ig(binary_id, dmHe, ierr)
         if (ierr /= 0) return

         ! Recurrence time: tau_rec = Delta_Mig / Mdot.
         t_rec = dmHe / mdot_msun_yr

      end subroutine my_t_recurrence

!*****************************************************!

      subroutine my_nova_flag(binary_id, nova_flag, ierr)
         ! nova_flag = 1 if system is in nova regime (Mdot below stable-burning threshold), else 0.
         ! Threshold (Wang+2018 eq. 4):
         ! Mdot_st = 1.46 * 10^-6 (-Mwd^3 + 3.45*Mwd^2 - 2.60*Mwd + 0.85)  [Msun/yr]
         ! Nova occurs when accretion rate falls below Mdot_st.
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: nova_flag
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: m_wd, mdot, mdot_st

         ierr       = 0
         nova_flag  = 0d0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr in my_nova_flag'
            return
         end if

         m_wd  = b% m(b% point_mass_i) / Msun
         mdot  = abs(b% mtransfer_rate) * secyer / Msun

         ! Nova regime requires Mdot < Mdot_st; above this threshold steady He-burning retains all mass with eta_H = 1.
         !stability bound from Bo Wang 2018 Res. Astron. Astrophys. 18 049
         mdot_st = 1.46d-6 * (-pow3(m_wd) + 3.45d0*pow2(m_wd) - 2.60d0*m_wd + 0.85d0)
         if (mdot_st > 0d0 .and. mdot > 10d0**nova_log10mdot_min .and. mdot < mdot_st) nova_flag = 1d0

      end subroutine my_nova_flag

!****************************************************!
      
      subroutine my_rlo_mdot(binary_id, mdot, ierr) ! Adapted from a routine kindly provided by Anastasios Fragkos
         integer, intent(in) :: binary_id
         real(dp), intent(out) :: mdot
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp):: mdot_normal, mdot_reverse

         include 'formats.inc'

         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         mdot = 0d0
         mdot_normal = 0d0
         mdot_reverse = 0d0

         if (b% mdot_scheme == "Kolb" .and. b% eccentricity <= 0.0) then
            call get_info_for_ritter(b)
            mdot_normal = b% mdot_thin
            call get_info_for_kolb(b)
            mdot_normal = mdot_normal + b% mdot_thick
         else if (b% mdot_scheme == "Kolb" .and. b% eccentricity > 0.0) then
            call get_info_for_ritter_eccentric(b)
            mdot_normal = b% mdot_thin
            call get_info_for_kolb_eccentric(b)
            mdot_normal = mdot_normal + b% mdot_thick
         end if

         mdot = mdot_normal
         !write(*,*) 'mdot_normal, from donor i:', mdot_normal, b% d_i

         if (b% point_mass_i == 0) then
            ! if mdot = 0d0 then ! no RLO from the donor so far.
            if (b% r(b% d_i) < b% rl(b% d_i)) then   ! Better version
               ! donor  reversed temporarily
               if (b% d_i == 2) then
                  b% d_i = 1
                  b% a_i = 2
                  b% s_donor => b% s1
                  b% s_accretor => b% s2
               else
                  b% d_i = 2
                  b% a_i = 1
                  b% s_donor => b% s2
                  b% s_accretor => b% s1
               end if

               if ((b% mdot_scheme == "Kolb") .and. (b% eccentricity<=0.0)) then
                  call get_info_for_ritter(b)
                  mdot_reverse = b% mdot_thin
                  call get_info_for_kolb(b)
                  mdot_reverse = mdot_reverse + b% mdot_thick
               else if ((b% mdot_scheme == "Kolb") .and. (b% eccentricity>0.0)) then
                  call get_info_for_ritter_eccentric(b)
                  mdot_reverse = b% mdot_thin
                  call get_info_for_kolb_eccentric(b)
                  mdot_reverse = mdot_reverse + b% mdot_thick
               end if

               !write(*,*) 'mdot_reverse, from donor i:', mdot_reverse, b% d_i
               if (abs(mdot_reverse) > abs(mdot_normal)) then
                  mdot = mdot_reverse
               else
                  ! switch donor back to the initial one in the step after the Kolb explicit calculation
                  if (b% d_i == 2) then
                     b% d_i = 1
                     b% a_i = 2
                     b% s_donor => b% s1
                     b% s_accretor => b% s2
                  else
                     b% d_i = 2
                     b% a_i = 1
                     b% s_donor => b% s2
                     b% s_accretor => b% s1
                  end if
               end if
               !write(*,*) 'final mdot, from donor i:', mdot,  b% d_i
            end if
         end if

      end subroutine my_rlo_mdot

      subroutine get_info_for_ritter(b)
         type(binary_info), pointer :: b
         real(dp) :: rho_exponent, F1, q, rho, p, grav, hp, v_th, rl3, q_temp
         include 'formats.inc'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% p(1) ! pressure at surface in dynes/cm^2
         grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2 ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))

         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10_cr(q_temp))
         rl3 = (b% rl(b% d_i))*(b% rl(b% d_i))*(b% rl(b% d_i))
         b% mdot_thin0 = (2.0D0*pi/exp_cr(0.5d0)) * v_th*v_th*v_th * &
             rl3/(b% s_donor% cgrav(1)*b% m(b% d_i)) * rho * F1
         ! Once again, do not extrapolate! Eq. (7) of Ritter 1988
         q_temp = min(max(q,0.04d0),20d0)
         if (q_temp < 1.0d0) then
            b% ritter_h = hp/( 0.954D0 + 0.025D0*log10_cr(q_temp) - 0.038D0*(log10_cr(q_temp))**2 )
         else
            b% ritter_h = hp/( 0.954D0 + 0.039D0*log10_cr(q_temp) + 0.114D0*(log10_cr(q_temp))**2 )
         end if

         b% ritter_exponent = (b% r(b% d_i)-b% rl(b% d_i))/b% ritter_h

         if (b% mdot_scheme == "Kolb") then
            if (b% ritter_exponent > 0) then
               b% mdot_thin = -b% mdot_thin0
            else
               b% mdot_thin = -b% mdot_thin0 * exp_cr(b% ritter_exponent)
            end if
         else
            b% mdot_thin = -b% mdot_thin0 * exp_cr(b% ritter_exponent)
         end if

      end subroutine get_info_for_ritter

      real(dp) function calculate_kolb_mdot_thick(b, indexR, rl_d) result(mdot_thick)
         real(dp), intent(in) :: rl_d
         integer, intent(in) :: indexR
         real(dp) :: F1, F3, G1, dP, q, rho, p, grav, hp, v_th, rl3, q_temp
         integer :: i
         type(binary_info), pointer :: b
         include 'formats.inc'

         !--------------------- Optically thin MT rate -----------------------------------------------
         ! As described in Kolb and H. Ritter 1990, A&A 236,385-392

         ! compute integral in Eq. (A17 of Kolb & Ritter 1990)
         mdot_thick = 0d0
         do i=1,indexR-1
            G1 = b% s_donor% gamma1(i)
            F3 = sqrt(G1) * pow_cr(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
            mdot_thick = mdot_thick + F3*sqrt(kerg*b% s_donor% T(i) / &
                         (mp*b% s_donor% mu(i))) * &
                         (b% s_donor% P(i+1)-b% s_donor% P(i))
         end do
         ! only take a fraction of dP for last cell
         G1 = b% s_donor% gamma1(i)
         F3 = sqrt(G1) * pow_cr(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
         dP = (b% s_donor% r(indexR) - rl_d) / &
            (b% s_donor% r(indexR) - b% s_donor% r(indexR+1)) * (b% s_donor% P(i+1)-b% s_donor% P(i))
         mdot_thick = mdot_thick + F3*sqrt(kerg * b% s_donor% T(i) / (mp*b% s_donor% mu(i)))*dP

         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10_cr(q_temp))
         mdot_thick = -2.0D0*pi*F1*rl_d*rl_d*rl_d/(b% s_donor% cgrav(1)*b% m(b% d_i))*mdot_thick

      end function calculate_kolb_mdot_thick

      subroutine get_info_for_kolb(b)
         type(binary_info), pointer :: b
         real(dp) :: F3, FF, G1, x_L1, q, g
         real(dp) :: mdot_thick0,  R_gas, dP, rl, s_div_rl
         integer :: i, indexR
         include 'formats.inc'

         !--------------------- Optically thick MT rate -----------------------------------------------
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         ! First we need to find how deep inside the star the Roche lobe reaches. In other words the mesh point of the star at which R=R_RL
         b% mdot_thick = 0d0
         indexR=-1
         if(b% r(b% d_i)-b% rl(b% d_i) > 0.0d0) then
            i=1
            do while (b% s_donor% r(i) > b% rl(b% d_i))
               i=i+1
            end do

            if (i .eq. 1) then
               b% mdot_thick = 0d0
            else
               b% mdot_thick = calculate_kolb_mdot_thick(b, i-1, b% rl(b% d_i))
            end if
         end if

      end subroutine get_info_for_kolb

      subroutine get_info_for_ritter_eccentric(b)
         type(binary_info), pointer :: b
         integer :: i
         real(dp) :: rho_exponent, F1, q, q_temp, rho, p, grav, hp, v_th, dm
         real(dp), DIMENSION(b% anomaly_steps):: mdot0, mdot, Erit, rl_d
         include 'formats.inc'

         ! Optically thin MT rate adapted for eccentric orbits
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

         rho = b% s_donor% rho(1) ! density at surface in g/cm^3
         p = b% s_donor% p(1) ! pressure at surface in dynes/cm^2
         grav = b% s_donor% cgrav(1)*b% m(b% d_i)/(b% r(b% d_i))**2 ! local gravitational acceleration
         hp = p/(grav*rho) ! pressure scale height
         v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1))) ! kerg = Boltzmann's constant
         ! phase dependant RL radius
         rl_d = b% rl(b% d_i) * (1d0 - b% eccentricity**2) / &
                (1 + b% eccentricity * cos(b% theta_co) )
         q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                       ! (Kolb & Ritter 1990 use the opposite!)
         q_temp = min(max(q,0.5d0),10d0)
         F1 = (1.23d0  + 0.5D0* log10_cr(q_temp))
         mdot0 = (2.0D0*pi/exp_cr(0.5d0)) * pow3(v_th) * rl_d*rl_d*rl_d / &
             (b% s_donor% cgrav(1)*b% m(b% d_i)) * rho * F1
         q_temp = min(max(q,0.04d0),20d0)
         if (q_temp < 1.0d0) then
            b% ritter_h = hp/( 0.954D0 + 0.025D0*log10_cr(q_temp) - 0.038D0*(log10_cr(q_temp))**2 )
         else
            b% ritter_h = hp/( 0.954D0 + 0.039D0*log10_cr(q_temp) + 0.114D0*(log10_cr(q_temp))**2 )
         end if
         Erit = (b% r(b% d_i)- rl_d) / b% ritter_h
         if (b% mdot_scheme == "Kolb") then
            do i = 1,b% anomaly_steps
               if (Erit(i) > 0) then
                  mdot(i) = -1 * mdot0(i)
               else
                  mdot(i) = -1 * mdot0(i) * exp(Erit(i))
               end if
            end do
         else
            mdot = -1 * mdot0 * exp(Erit)
         end if
         b% mdot_donor_theta = mdot
         !integrate to get total massloss
         dm = 0d0
         do i = 2,b% anomaly_steps ! trapezoidal integration
            dm = dm + 0.5d0 * (mdot(i-1) + mdot(i)) * (b% time_co(i) - b% time_co(i-1))
         end do
         b% mdot_thin = dm
      end subroutine get_info_for_ritter_eccentric

      subroutine get_info_for_kolb_eccentric(b)
         type(binary_info), pointer :: b
         real(dp) :: e, dm
         integer :: i, j
         real(dp), DIMENSION(b% anomaly_steps):: rl_d_i, mdot_thick_i
         include 'formats.inc'
         ! Optically thick MT rate adapted for eccentric orbits
         ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392
         b% mdot_thick = 0d0
         e = b% eccentricity
         ! If the radius of the donor is smaller as the smallest RL radius,
         ! there is only atmospheric RLOF, thus return.
         if ( b% r(b% d_i) < b% rl(b% d_i) * (1-e**2)/(1+e) ) then
            return
         end if
         ! phase dependant RL radius
         rl_d_i = b% rl(b% d_i) * (1d0 - b% eccentricity**2) / &
                  (1 + b% eccentricity * cos(b% theta_co) )
         ! For each point in the orbit calculate mdot_thick
         do i = 1,b% anomaly_steps
            ! find how deep in the star we are
            j=1
            do while (b% s_donor% r(j) > rl_d_i(i))
               j=j+1
            end do
            ! calculate mdot_thick
            if (j .eq. 1) then
               mdot_thick_i(i) = 0d0
            else
               mdot_thick_i(i) = calculate_kolb_mdot_thick(b, j-1, rl_d_i(i))
            end if
         end do
         b% mdot_donor_theta = b% mdot_donor_theta + mdot_thick_i
         ! Integrate mdot_thick over the orbit
         dm = 0d0
         do i = 2,b% anomaly_steps ! trapezoidal integration
            dm = dm + 0.5d0 * (mdot_thick_i(i-1) + mdot_thick_i(i)) * &
                              (b% time_co(i) - b% time_co(i-1))
         end do
         b% mdot_thick = dm
      end subroutine get_info_for_kolb_eccentric

      ! This subroutine determines which torque prescription to use for magnetic braking
      ! To use custom magnetic braking prescriptions, need at least these options...
      ! >> In inlist_project:
      !      do_jdot_mb = .true.
      !      use_other_jdot_mb = .true.
      !
      ! >> In inlist1:
      !      x_character_ctrl(1) = <<MB option>>
      ! >> In inlist2:
      !      x_character_ctrl(1) = <<MB option>>
      !
      ! Replace <<MB option>> to select a prescription:
      !    <<MB option>> = 1 (Garraffo et al. 2018)
      !    <<MB option>> = 2 (Matt et al. 2015)
      !    <<MB option>> = 3 (Van & Ivanova 2019 -- CARB)
      subroutine mb_torque_selector(binary_id, ierr)
         use star_lib, only: star_ptr
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         real(dp) :: dJdt
 
         ierr = 0
         ! call star_ptr(id, s, ierr)
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
 
         b% jdot_mb = 0d0
         dJdt = 0d0
         ! turn on Garraffo+ 2018 style braking?
         if (b% s1% x_character_ctrl(1) == 'g18') then
            if (b% model_number == 0) then
               write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
               write(*,*) 'Garraffo+ 2016/18 torque enabled (star 1)', b% d_i
               write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            end if
            call garraffo_torque(binary_id, b% s_donor, dJdt, ierr)
            if (.not. b% do_jdot_ls) then
               b% jdot_mb = b% jdot_mb + dJdt
            end if
 
            ! check if braking should be applied from the accretor as well
            if ((b% point_mass_i == 0) .and. (b% include_accretor_mb)) then
               if (b% model_number == 0) then
                  write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                  write(*,*) 'Garraffo+ 2016/18 torque enabled (star 2)', b% a_i
                  write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
               end if
               call garraffo_torque(binary_id, b% s_accretor, dJdt, ierr)
               if (.not. b% do_jdot_ls) then
                  b% jdot_mb = b% jdot_mb + dJdt
               end if
            end if
 
         ! turn on Matt+ 2015 style braking?
         else if (b% s1% x_character_ctrl(1) == 'm15') then
            if (b% model_number == 0) then
               write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
               write(*,*) 'Matt+ 2015 torque enabled (star 1)'
               write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            end if
            call matt_torque(binary_id, b% s_donor, dJdt, ierr)
            if (.not. b% do_jdot_ls) then
               b% jdot_mb = b% jdot_mb + dJdt
            end if
 
            if ((b% point_mass_i == 0) .and. (b% include_accretor_mb)) then
               if (b% model_number == 0) then
                  write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                  write(*,*) 'Matt+ 2015 torque enabled (star 2)'
                  write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
               end if
               call matt_torque(binary_id, b% s_accretor, dJdt, ierr)
               if (.not. b% do_jdot_ls) then
                  b% jdot_mb = b% jdot_mb + dJdt
               end if
            end if
         
         ! turn on Van & Ivanova 2019 (CARB) style braking?
         else if (b% s1% x_character_ctrl(1) == 'carb') then
            if (b% model_number == 0) then
               write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
               write(*,*) 'Van & Ivanova 2019 (CARB) torque enabled (star 1)'
               write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            end if
            call carb_torque(binary_id, b% s_donor, dJdt, ierr)
            if (.not. b% do_jdot_ls) then
               b% jdot_mb = b% jdot_mb + dJdt
            end if
 
            if ((b% point_mass_i == 0) .and. (b% include_accretor_mb)) then
               if (b% model_number == 0) then
                  write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                  write(*,*) 'Van & Ivanova 2019 (CARB) torque enabled (star 2)'
                  write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
               end if
               call carb_torque(binary_id, b% s_accretor, dJdt, ierr)
               if (.not. b% do_jdot_ls) then
                  b% jdot_mb = b% jdot_mb + dJdt
               end if
            end if
         end if
      end subroutine mb_torque_selector

      ! Matt et al. (2015), ApJ, 799, L23 magnetic braking prescription
      subroutine matt_torque(binary_id, s, dJdt, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         integer :: j, k, mix_reg_extent, mix_reg_bot_k, mix_reg_top_k, &
                 nz, n_conv_bdy, i, k_ocz_bot, k_ocz_top
 
         real(dp) :: Prot, Ro, Rosol, Rosat, K_const, m, p, u, gamma, dJdt, &
                     tau_convective, mixing_length_at_bcz, MOI, Om, rsol, &
                     msol, omega_sol, tau_cz_sol, chi, T0, t_spindown, &
                     ocz_bot_mass
 
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
 
         ! number of convective regions; corresponds to outer CZ as an index
         i = s% n_conv_regions
         Prot = 0d0
         Ro = 0d0
         dJdt = 0d0
         s% extra_omegadot(:) = 0d0
 
         ! This only calculates if a convective envelope exists and not fully conv.
         ! (ocz_bot_mass > 0d0) is used to detect whether fully convective
         if ((s% n_conv_regions > 0)) then

            ocz_bot_mass = s% cz_bot_mass(i)

            if ((s% cz_top_mass(i)/s% mstar > 0.99d0) .and. &
               ((s% cz_top_mass(i)-s% cz_bot_mass(i))/s% mstar > 1d-11) .and. &
               (ocz_bot_mass > 0d0)) then
  
               call calc_tau_convective(binary_id, s, tau_convective, ierr)  
 
               MOI = dot_product(s% dm_bar(1:s% nz), s% i_rot(1:s% nz))
               Om = s% omega_avg_surf !s% total_angular_momentum / MOI
               ! surface rotation period (assuming tidal locking if do_jdot_mb = .true.)
               Prot = b% period
 
               ! rossby number
               Ro = Prot / tau_convective
               Rosol = 2d0 ! Amard+ 2019
               Rosat = 0.14d0 ! Amard+ 2019
               K_const = 1.4e30 ! Solar calibrated as in Gossage et al. 2021, ApJ 912, 65
               m = 0.22d0 ! Solar calibrated '                                    '
               p = 2.6d0 ! Solar calibrated '                                     '
               tau_cz_sol = 12.9d0*86400d0 ! 12.9 days Matt+ 2015 [sec]
               omega_sol = 2.6d-6 !2.6E-6 s^-1 solar solid body ang. rot. rate from measured Prot of Sun...
               u = s% v_div_v_crit_avg_surf 
               msol = 1.99d33 ! g
               rsol = 6.96d10 ! cm
               chi = Rosol / Rosat
 
               gamma = pow_cr(1.0+powi_cr(u / 0.072, 2), 0.5d0)
               T0 = K_const * pow_cr(pow_cr(10d0,s% log_surface_radius),3.1d0) * &
                    pow_cr(s% star_mass,0.5d0) * pow_cr(gamma,-2.0d0 * m)
 
               ! saturated regime
               if (Ro .lt. Rosat) then
                  dJdt = T0 * pow_cr(chi, p) * (Om / omega_sol)
               ! unsaturated regime
               else
                  dJdt = T0 * pow_cr(tau_convective/tau_cz_sol, p) * &
                         pow_cr(Om / omega_sol, p + 1.0d0)
               end if
 
               ! angular momentum change per second.
               dJdt = max(0.0d0, dJdt)
               dJdt = -dJdt
 
               ! Check if spindown timescale is shorter than timestep. Print a warning in case.
               ! In extras_finish_step check that dt < t_spindown. If not, decrease timestep
               t_spindown = abs(s% total_angular_momentum / dJdt) ! Estimate spindown timescale
 
               ! If tidal sync is enforced, remove AM from the orbit
               if (.not. b% do_jdot_ls) then
                  return
               ! If tidal sync is not enforced, remove AM from the indiv. stars
               else
                  do k = s% nz, 1, -1
                     ! angular velocity loss per second. If d(omega)/ dt
                     ! would be too large for current cell:
                     if (s% omega(k) < s% dt * abs(dJdt / MOI) ) then
                        ! use omega(k) / dt to as a cap on the 'max rate of change' for omega
                        s% extra_omegadot(k) = - s% omega(k) / s% dt
                     else
                        ! or else if d(omega)/dt * dt is < current cell's 
                        ! omega, use the calculated value. 
                        s% extra_omegadot(k) = dJdt / MOI 
                     end if
                  end do
                  ! Reset dJdt so it is not also removed from the orbit
                  dJdt = 0d0
                  return
               end if
            else
               ! To avoid decreasing the timestep in extras_finish_step
               t_spindown = 100 * s% dt
            end if
 
         else
            ! To avoid decreasing the timestep in extras_finish_step
            t_spindown = 100 * s% dt
         end if
 
      end subroutine matt_torque

      ! As implemented in Van & Ivanova 2019 (CARB magnetic braking), ApJ, 886, L31
      ! from files hosted on Zenodo: https://zenodo.org/record/3647683#.Y_TfedLMKUk
      ! Slightly modified to avoid INF values.
      subroutine carb_torque(binary_id, s, dJdt, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         integer :: k, nz, i
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: tau_convective, period, mdot, mstar, ocz_bot_mass
         real(dp) :: rsun4, two_pi_div_p3, two_pi_div_p2, K2
         real(dp) :: tt_ratio, tt4
         real(dp) :: rot_ratio, rot4
         real(dp) :: rad4, rsurf
         real(dp) :: v_esc2, v_mod2
         real(dp) :: alfven_no_R, R_alfven
         real(dp) :: dJdt, MOI
         logical :: conv_env_found
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr .ne. 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
 
         nz = s% nz
         i = s% n_conv_regions
         tau_convective = 0d0
         period = 2d0 * pi / s% omega_avg_surf
         MOI = dot_product(s% dm_bar(1:s% nz), s% i_rot(1:s% nz))
         mdot = s% star_mdot * msun / 3.154d7 ! [g/sec]
         rsurf = rsun * pow_cr(10d0, s% log_surface_radius) ! [cm]
         mstar = s% star_mass * msun ! [g]

         s% extra_jdot(:) = 0d0
         s% extra_omegadot(:) = 0d0
 
         if (s% n_conv_regions > 0) then

            ocz_bot_mass = s% cz_bot_mass(i)

            if ((s% cz_top_mass(i) / s% mstar > 0.99d0) .and. &
               ((s% cz_top_mass(i) - s% cz_bot_mass(i)) / s% mstar > 1d-11) .and. &
               (ocz_bot_mass > 0d0)) then

               call calc_tau_convective(binary_id, s, tau_convective, ierr)

               ! MAGNETIC BRAKING CALCULATION
               rsun4 = pow4(rsun)
               two_pi_div_p2 = (2.0 * pi / period) * (2.0 * pi / period)
               K2 = 0.07 * 0.07
               tt_ratio = tau_convective / 2.8d6
               tt4 = pow4(tt_ratio)
               rot_ratio = (2073600d0 / period)
               rot4 = pow4(rot_ratio)
               rad4 = pow4(rsurf)
 
               ! escape speed
               v_esc2 = 2.0 * standard_cgrav * b% m(b% d_i) / b% r(b% d_i)
               ! modified escape speed, e.g., Matt et al. 2012/Reville et al. 2015
               v_mod2 = v_esc2 + 2.0 * two_pi_div_p2 * b% r(b% d_i) * b% r(b% d_i) / K2 
                     
               ! SSG edit to prevent INF values when b% mdot_system_wind(b% d_i) = 0
               if (abs(b% mdot_system_wind(b% d_i)) > 0d0) then
                  alfven_no_R = rad4 * rot4 * tt4 / (b% mdot_system_wind(b% d_i) * &
                                b% mdot_system_wind(b% d_i)) * (1.0 / v_mod2)
               else
                  alfven_no_R = 0d0
               end if
   
               R_alfven = b% r(b% d_i) * alfven_no_R**(1.d0/3.d0)
               dJdt = 1d0 * (2.0/3.0) * (2.0*pi/b% period) * b% mdot_system_wind(b% d_i) * R_alfven * R_alfven
   
               ! If tidal sync is enforced, remove AM from the orbit
               if (.not. b% do_jdot_ls) then
                  return
               ! If tidal sync is not enforced, remove AM from the individual stars
               else
                  do k = s% nz, 1, -1
                     ! angular velocity loss per second. If d(omega)/ dt would be too large for current cell:
                     if (s% omega(k) < s% dt * abs(dJdt / MOI) ) then
                        ! use omega(k) / dt to as a cap on the 'max rate of change' for omega
                        s% extra_omegadot(k) = - s% omega(k) / s% dt
                     else
                        ! or else if d(omega)/dt * dt is < current cell's omega, use the calculated value. 
                        s% extra_omegadot(k) = dJdt / MOI 
                     end if
                  end do
                  ! Reset dJdt so it is not also removed from the orbit
                  dJdt = 0d0
                  return
               end if   
            end if
         end if
 
      end subroutine carb_torque

      ! Garraffo et al. (2018), ApJ, 862, 90 torque prescription
      subroutine garraffo_torque(binary_id, s, dJdt, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         integer :: j, k, mix_reg_extent, mix_reg_bot_k, mix_reg_top_k, &
                    nz, n_conv_bdy, i
 
         real(dp) :: Prot, Ro, n, Qn, dJdt, tau_convective, MOI, Om, &
                     scale_height_at_bcz, c_const, omega_crit, & 
                     residual_jdot, a_constant, b_constant, c_constant, &
                     t_spindown
 
         ierr = 0
 
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
             write(*,*) 'failed in binary_ptr'
            return
         end if
 
         ! output info about the CONV. ENV.: the CZ location, turnover time
         nz = s% nz
         n_conv_bdy = s% num_conv_boundaries
         i = s% n_conv_regions
         tau_convective = 0d0
         tau_convective = 0d0
         Prot = 0d0
         Ro = 0d0
         a_constant = 0.03d0  ! Solar calibrated as in Gossage et al. 2021, ApJ 912, 65
         b_constant = 0.5d0  ! Solar calibrated '                                     '
         c_constant = 3d41 ! Solar calibrated '                                      '
         n = 0d0
         Qn = 0d0
         scale_height_at_bcz = 0d0
         residual_jdot = 0d0
         dJdt = 0d0
         s% extra_omegadot(:) = 0d0
 
         if ((s% n_conv_regions > 0)) then 
            if ((s% cz_top_mass(i)/s% mstar > 0.99d0) .and. &
               ((s% cz_top_mass(i)-s% cz_bot_mass(i))/s% mstar > 1d-11) .and. &
               (s% star_age > s% x_ctrl(2))) then
 
               ! calculate convective turnover time
               call calc_tau_convective(binary_id, s, tau_convective, ierr)
 
               ! spin down according to Garraffo et al. 2018
               MOI = dot_product(s% dm_bar(1:s% nz), s% i_rot(1:s% nz))
               Om = s% omega_avg_surf
               ! surface rotation period (assuming tidal locking if do_jdot_mb = .true.)
               Prot = b% period
               ! Rossby number
               Ro = Prot / tau_convective
               ! calculate n
               ! a = 0.02, b = 2.0
               n = (a_constant/pow_cr(Ro,1d0)) + (b_constant*Ro) + 1d0
               if (n < 1d0) then
                 n = 1d0
               else if (n > 1d1) then
                 n = 1d1
               end if
 
               ! magnetic supression factor
               Qn = 4.05_dp*exp_cr(-1.4_dp*n)
 
               ! angular momentum change per second.
               dJdt = c_constant * powi_cr(Om, 3) * tau_convective * Qn
               dJdt = max(0.0_dp, dJdt)
               dJdt = -dJdt
 
               ! Check if spindown timescale is shorter than timestep. Print a warning in case.
               ! In extras_finish_step check that dt < t_spindown. If not, decrease timestep
               t_spindown = abs(s% total_angular_momentum / dJdt) ! Estimate spindown timescale
 
               ! If tidal sync is enforced, remove AM from the orbit
               if (.not. b% do_jdot_ls) then
                  return

               ! if tidal sync is not enforced, remove AM from the indiv. stars
               else
                  do k = s% nz, 1, -1
                     ! angular velocity loss per second. 
                     ! If d(omega)/ dt would be too large for current cell:
                     if (s% omega(k) < s% dt * abs(dJdt / MOI) ) then
                        ! use omega(k) / dt to as a cap on the 'max rate of change' for omega
                        s% extra_omegadot(k) = - s% omega(k) / s% dt
                     else
                        ! or else if d(omega)/dt * dt is < current cell's omega, use the calculated value. 
                        s% extra_omegadot(k) = dJdt / MOI 
                     end if
                  end do
                  ! Reset dJdt so it is not also removed from the orbit
                  dJdt = 0d0
                  return
               end if
            else
               t_spindown = 100 * s% dt ! To avoid decreasing the timestep in extras_finish_step
            end if
         else
            t_spindown = 100 * s% dt ! To avoid decreasing the timestep in extras_finish_step     
         end if
 
      end subroutine garraffo_torque

      ! As implemented in Van & Ivanova 2019 (CARB magnetic braking), ApJ, 886, L31
      ! from files hosted on Zenodo: https://zenodo.org/record/3647683#.Y_TfedLMKUk
      subroutine check_radiative_core(b)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
             
         real(dp) :: sum_conv, q_loc, sum_div_qloc 
         integer :: i, k, id
 
         include 'formats.inc'
 
         do i=1,2
            if (i == 1) then
               s => b% s_donor
               id = b% d_i
            else if (b% point_mass_i == 0 .and. b% include_accretor_mb) then
               s => b% s_accretor
               id = b% a_i
            else
               exit
            end if
 
            ! calculate how much of inner region is convective
            sum_conv = 0; q_loc = 0
            do k = s% nz, 1, -1
               q_loc = s% q(k)
               if (q_loc > 0.5d0) exit 
               if (s% mixing_type(k) == convective_mixing) &
                  sum_conv = sum_conv + s% dq(k)
            end do
                 
            sum_div_qloc = (b% sum_div_qloc(id) + sum_conv/q_loc)/2
            b% sum_div_qloc(id) = sum_div_qloc   

            if (b% have_radiative_core(id)) then ! check if still have rad core
               if (sum_div_qloc > 0.75d0) then
                  b% have_radiative_core(id) = .false.
                  write(*,*)
                  write(*,*) 'turn off magnetic braking because radiative core has gone away'
                  write(*,*)
                  ! required mdot for the implicit scheme may drop drastically,
                  ! so its neccesary to increase change factor to avoid implicit 
                  ! scheme from getting stuck
                  b% change_factor = b% max_change_factor
               end if
            else if (sum_div_qloc < 0.25d0) then ! check if now have rad core
               if (.not. b% have_radiative_core(id)) then
                  write(*,*)
                  write(*,*) 'turn on magnetic braking'
                  write(*,*)
               end if
               b% have_radiative_core(id) = .true.
            end if
         end do
                 
       end subroutine check_radiative_core

      ! calculate the convective turnover time of one star, half a 
      ! pressure scale height above the bottom of the convection zone
      subroutine calc_tau_convective(binary_id, s, ocz_turnover_time, ierr)
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         integer :: k, nz, n_conv_bdy, i, k_ocz_bot, k_ocz_top
         real(dp) :: ocz_top_mass, ocz_bot_mass, mixing_length_at_bcz, &
                     ocz_turnover_time, ocz_top_radius, ocz_bot_radius
 
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
 
         ! output info about the CONV. ENV.: the CZ location, turnover time
         nz = s% nz
         n_conv_bdy = s% num_conv_boundaries
         i = s% n_conv_regions
         k_ocz_bot = 0
         k_ocz_top = 0
         ocz_turnover_time = 0d0
         ocz_top_mass = 0d0
         ocz_bot_mass = 0d0
         ocz_top_radius = 0d0
         ocz_bot_radius = 0d0
         mixing_length_at_bcz = 0d0

         ocz_bot_mass = s% cz_bot_mass(i)
         ocz_top_mass = s% cz_top_mass(i)

         !get top radius information
         !start from k=2 (second most outer zone) in order to access k-1
         do k=2,nz
            if (s% m(k) < ocz_top_mass) then
               ocz_top_radius = s% r(k-1)
               k_ocz_top = k-1
               exit
            end if
         end do

         ! get bottom radius information
         if (ocz_bot_mass == 0d0) then
            ocz_bot_radius = s% r(nz)
            k_ocz_bot = nz
         else
            do k=2,nz
               if (s% m(k) < ocz_bot_mass) then
                  ocz_bot_radius = s% r(k-1)
                  k_ocz_bot = k-1
                  exit
               end if
            end do
         end if

         !if the star is fully convective, then the bottom boundary is the center
         if ((k_ocz_bot == 0) .and. (k_ocz_top > 0)) then
            k_ocz_bot = nz
         end if

         mixing_length_at_bcz = s% mlt_mixing_length(k_ocz_bot)
         !scale_height_at_bcz = s% scale_height(k_ocz_bot)
         !compute the "local" turnover time a scale height above the BCZ
         do k=k_ocz_top,k_ocz_bot
            if (s% r(k) < (s% r(k_ocz_bot) + 0.5d0*s% scale_height(k)) ) then
               ocz_turnover_time = s% mixing_length_alpha * s% scale_height(k) / s% conv_vel(k)
               exit
            end if
         end do   

      end subroutine calc_tau_convective

      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 27
      end function how_many_extra_binary_history_columns

      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         use const_def, only: dp
         type (binary_info), pointer :: b
		   type (star_info), pointer :: s
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         integer:: i_don, i_acc
         real(dp) :: beta, trap_rad, mdot_edd, accretor_radius, MOI, prot
         real(dp) :: tau_convective, Ro, n_complexity, a_constant, b_constant
         real(dp) :: mdot_msun_yr, m_wd, log10mdot, t_rec, nova_flag, dt_yr, ratio

         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         call my_mdot_edd(binary_id,mdot_edd,ierr)

         if (b% point_mass_i == 0) then ! if there is no compact object then trappping radius is 0
            trap_rad = 0.0_dp
            accretor_radius = 0.0_dp
         else ! Begelman 1997 and King & Begelman 1999 eq. 1: accretor is star 2
            trap_rad = 0.5_dp*abs(b% mtransfer_rate) * acc_radius(b, b% m(2)) / mdot_edd
            accretor_radius = acc_radius(b, b% m(2))
         end if

         names(1) = 'trap_radius'
         vals(1) = trap_rad/Rsun ! in Rsun units
         names(2) = 'acc_radius'
         vals(2) = accretor_radius ! in cm units

         names(3) = 't_sync_rad_1'
         names(4) = 't_sync_conv_1'
         names(5) = 't_sync_rad_2'
         names(6) = 't_sync_conv_2'
         if (b% point_mass_i /= 1) then
            vals(3) = b% s1% xtra1
            vals(4) = b% s1% xtra2
         else
            vals(3) = -1.0d0
            vals(4) = -1.0d0
         end if
         if (b% point_mass_i /= 2) then
            vals(5) = b% s2% xtra1
            vals(6) = b% s2% xtra2
         else
            vals(5) = -1.0d0
            vals(6) = -1.0d0
         end if

         names(7) = 'star_mb_jdot_1'
         if (b% point_mass_i /= 1) then
            MOI = dot_product(b% s1% dm_bar(1:b% s1% nz),b% s1% i_rot(1:b% s1% nz))
            vals(7) = b% s1% extra_omegadot(1) * MOI
         else
            vals(7) = 0d0
         end if
		
		   names(8) = 'star_mb_jdot_2'
         if (b% point_mass_i /= 2) then
            MOI = dot_product(b% s2% dm_bar(1:b% s2% nz),b% s2% i_rot(1:b% s2% nz))
            vals(8) = b% s2% extra_omegadot(1) * MOI
         else
            vals(8) = 0d0
         end if

         ! calc rossby, garraffo magnetic complexity param, conv turnover time:
         a_constant = 0.002d0  ! Solar calibrated as in Gossage et al. 2021, ApJ 912, 65
         b_constant = 0.5d0  ! Solar calibrated '                                     '
         tau_convective = 0d0
         Ro = 0d0
         n_complexity = 0d0

         if (b% point_mass_i /= 1) then
            call star_ptr(1, s, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in star_ptr'
               return
            end if
            
            if (b% s1% n_conv_regions > 0) then
               call calc_tau_convective(binary_id, s, tau_convective, ierr)
            else
               tau_convective = 1d-99
            end if
            prot = 2d0 * pi / b% s1% omega_avg_surf
            Ro = prot / tau_convective

            n_complexity = (a_constant/pow_cr(Ro,1d0)) + (b_constant*Ro) + 1d0
            if (n_complexity < 1d0) then
               n_complexity = 1d0
            else if (n_complexity > 1d99) then
               n_complexity = 1d99
            end if
         end if

         names(9) = "conv_tau_sec_1"
         vals(9) = tau_convective ! [s]

         names(10) = "rossby_number_1"
         vals(10) = Ro

         names(11) = "n_complexity_1"
         vals(11) = n_complexity

         ! calc values for star 2
         tau_convective = 0d0
         Ro = 0d0
         n_complexity = 0d0

         if (b% point_mass_i /= 2) then

            call star_ptr(2, s, ierr)
            if (ierr /= 0) then
               write(*,*) 'failed in star_ptr'
               return
            end if

            if (b% s2% n_conv_regions > 0) then
               call calc_tau_convective(binary_id, s, tau_convective, ierr)
            else
               tau_convective = 1d-99
            end if
            prot = 2d0 * pi / b% s2% omega_avg_surf
            Ro = prot / tau_convective

            n_complexity = (a_constant/pow_cr(Ro,1d0)) + (b_constant*Ro) + 1d0
            if (n_complexity < 1d0) then
               n_complexity = 1d0
            else if (n_complexity > 1d99) then
               n_complexity = 1d99
            end if
         end if

         names(12) = "conv_tau_sec_2"
         vals(12) = tau_convective ! [s]

         names(13) = "rossby_number_2"
         vals(13) = Ro

         names(14) = "n_complexity_2"
         vals(14) = n_complexity
		
         !write(*,*) "synchr timescales: ", b% s1% xtra1, b% s1% xtra2, b% s2% xtra1, b% s2% xtra2

!**************************!extra columns from novae physics***********!
         names(15) = 'eta_retention'
         vals(15)  = b% s_donor% xtra3
         names(16) = 'm_ignition'
         vals(16)  = b% s_donor% xtra4
         names(17) = 't_recurrence'
         names(18) = 'nova_validity'
         names(19) = 'mass_transfer_beta'
         names(20) = 'nova_flag'
         names(21) = 'n_novae'
         names(22) = 'ratio_trec_dt'
         names(23)  = 'nova_mass_lost'
         vals(23)   = b% s_donor% xtra5
         names(24) = 'wd_mass_ref_nova'
         vals(24)  = b% s_donor% xtra6 / Msun  ! WD mass [Msun] at last nova reset
         names(25) = 'ejected_mass_cum'
         vals(25)  = b% s_donor% xtra8
         names(26) = 'accreted_mass_cum'
         vals(26)  = b% s_donor% xtra9
         names(27) = 'ejected_fraction'

! 15 eta_retention     xtra3  He retention efficiency (Wu+2017); holds last value outside nova regime
! 16 m_ignition        xtra4  ignition mass dm_ig [Msun] from my_dm_ig; holds last value outside nova regime
! 17 t_recurrence      live   dm_ig / Mdot [yr]; 0 unless nova_flag = 1 and beta > 0
! 18 nova_validity     live   1 if WD mass and Mdot are inside the trusted fit range, else 0
! 19 mass_transfer_beta       beta = 1 - eta in nova regime, 0 otherwise
! 20 nova_flag         live   1 if Mdot < Mdot_st (Wang+2018), else 0
! 21 n_novae           xtra7  dt / t_rec, novae this step
! 22 ratio_trec_dt     live   t_rec / dt (step-size diagnostic)
! 23 nova_mass_lost    xtra5  beta * |Mdot| * dt [Msun]
! 24 wd_mass_ref_nova  xtra6  WD mass [Msun] at last beta > 0 step (unused by physics)
! 25 ejected_mass_cum  xtra8  cumulative ejected mass [Msun]
! 26 accreted_mass_cum xtra9  cumulative retained mass [Msun]
! 27 ejected_fraction  live   100 * ejected / (ejected + accreted)

!******************************************************************!

         if (b% s_donor% xtra9 + b% s_donor% xtra8 > 0d0) then
            vals(27) = b% s_donor% xtra8 / (b% s_donor% xtra9 + b% s_donor% xtra8) * 100d0
         else
            vals(27) = 0d0
         end if
         if (b% mtransfer_rate /= 0d0) then
            mdot_msun_yr = abs(b% mtransfer_rate) * secyer / Msun
            m_wd         = b% m(b% point_mass_i) / Msun
            log10mdot    = log10_cr(mdot_msun_yr)
            if (in_nova_regime(m_wd, log10mdot)) then
               vals(18) = 1d0
            else
               vals(18) = 0d0
            end if
            call my_nova_flag(binary_id, nova_flag, ierr)
            if (ierr /= 0) return
            vals(20) = nova_flag
            ! Only log t_recurrence (and the derived ratio) on nova-active
            ! steps (mass_transfer_beta > 0)
            if (nova_flag == 1d0 .and. b% mass_transfer_beta > 0d0) then
               call my_t_recurrence(binary_id, t_rec, ierr)
               if (ierr /= 0) return
               vals(17) = t_rec
               dt_yr   = b% s_donor% dt / secyer
               ratio   = t_rec / dt_yr
               vals(21) = b% s_donor% xtra7
               vals(22) = ratio
            else
               vals(17) = 0d0
               vals(21) = 0d0
               vals(22) = 0d0
            end if
         else
            vals(17) = 0d0
            vals(18) = 0d0
            vals(20) = 0d0
            vals(21) = 0d0
            vals(22) = 0d0
         end if
         vals(19) = b% mass_transfer_beta
      end subroutine data_for_extra_binary_history_columns


      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         real(dp) :: rl_gap_1
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         if (.not. restart) then
            ! extras are used to store the two tidal sychronization timescales (rad/conv) for each star.
            ! -1 if they are point masses
            if (b% point_mass_i /= 1) then
               b% s1% xtra1 = -1.0d0 ! t_sync_rad_1
               b% s1% xtra2 = -1.0d0 ! t_sync_conv_1
            end if
            if (b% point_mass_i /= 2) then
               b% s2% xtra1 = -1.0d0 ! t_sync_rad_2
               b% s2% xtra2 = -1.0d0 ! t_sync_conv_2
            end if
         end if
         extras_binary_startup = keep_going

         if (b% doing_first_model_of_run .and. b% terminate_if_initial_overflow &
            .and. (.not. b% ignore_rlof_flag .or. b% model_twins_flag)) then

            rl_gap_1 = (b% s1% photosphere_r - b% rl(1)/Rsun * (1 - b% eccentricity) )/(b% rl(1)/Rsun)
            if (b% rl_relative_gap(b% d_i) >= 0.0d0 &
                .or. (b% point_mass_i == 0 .and. b% rl_relative_gap(b% a_i) >= 0.0d0) &
                .or. rl_gap_1 >= 0.0 ) then

               extras_binary_startup = terminate
               write(*,'(g0)') "termination code: Terminate because of overflowing initial model"
            end if
         end if

         ! This is if we want to modulate orbital AML via LS coupling
         if (b% do_jdot_ls .and. b% use_other_jdot_ls) then
            ! ensure these options are set b/c LS is incompatible 
            ! with do_jdot_mb, which assumes tidal locking
            b% do_jdot_mb = .false.
            b% do_tidal_sync = .true.

            write(*,*) '++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) 'use_other_jdot_ls (jdot_ls_with_mb) = ', b% use_other_jdot_ls
            write(*,*) 'set do_jdot_mb = ', b% do_jdot_mb
            write(*,*) 'set do_tidal_sync = ', b% do_tidal_sync
            write(*,*) '++++++++++++++++++++++++++++++++++++++++++'

            b% s1% use_other_torque = .true.
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) 'Magnetic braking torque enabled for star 1'
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

            b% s2% use_other_torque = .true.
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) 'Magnetic braking torque enabled for star 2'
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

         end if

      end function  extras_binary_startup

!****************************************************************!

      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         real(dp) :: nova_flag, eta, m_wd
         real(dp) :: t_recurr, mdot_msun_yr

         extras_binary_start_step = keep_going
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) return

         mdot_msun_yr = abs(b% mtransfer_rate) * secyer / Msun

         if (b% point_mass_i /= 0) then
          if (b% mtransfer_rate /= 0d0 .and. &
            b% m(b% point_mass_i) < nova_m_wd_max * Msun .and. &
            log10_cr(mdot_msun_yr) >= nova_log10mdot_min) then

            call my_nova_flag(binary_id, nova_flag, ierr)
            if (ierr /= 0) return

            if (nova_flag == 1d0) then
               call my_t_recurrence(binary_id, t_recurr, ierr)
               if (ierr /= 0) return
            endif

            m_wd = b% m(b% point_mass_i)

            if (nova_flag == 1d0) then
               ! SIMPLE-BETA VERSION: time-averaged nova retention on every
               ! nova-regime step, beta = 1 - eta. No ignition-mass gate (the
               ! gated version only set beta when the mass accreted since the
               ! last nova reached dmHe, which leaves the WD almost fully
               ! conservative whenever dt < t_rec). 
               if (in_nova_regime(m_wd/Msun, log10_cr(mdot_msun_yr))) then
                  call my_eta_accretion(binary_id, eta, ierr)
                  if (ierr /= 0) return
                  b% mass_transfer_beta = 1d0 - eta
               else
                  ! Outside the trusted fit range: fully conservative.
                  b% mass_transfer_beta = 0d0
               end if
            else
               ! Stable burning: fully conservative
               b% mass_transfer_beta = 0d0
            end if

          else
            ! No active mass transfer this step, or WD mass outside the
            ! range the nova physics applies to: nova_flag/nova_validity
            ! will read 0, so beta must not retain a stale nonzero value
            ! from an earlier nova-active step.
            b% mass_transfer_beta = 0d0
          end if
         end if

      end function  extras_binary_start_step

!***************************************************************!

      !Return either rety,backup,keep_going or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer:: i_don, i_acc
         real(dp) :: q
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         extras_binary_check_model = keep_going


!         if (b% point_mass_i/=0 .and. ((b% rl_relative_gap(1) .ge. 0.d0) &
!             .or. (abs(b% mtransfer_rate/(Msun/secyer)) .ge. 1.0d-10))) then
          if (.true.) then     
            if (b% point_mass_i/=1) then
               i_don = 1
               b% s_donor => b% s1
            else
               i_don = 2
               b% s_donor => b% s2
            end if
            ! Turning back on binary orbital evolution
!            if ((.not. b% s_donor% x_logical_ctrl(6)) .and. (.not. b% use_other_jdot_ls)) then
!               b% do_jdot_mb = .true. ! turn on magnetic braking for RLOFing HMS stars only
!           end if
            b% do_jdot_gr = .true.
            b% do_jdot_ml = .true.
            b% do_jdot_ls = .true.
            b% do_jdot_missing_wind = .true.
            b% do_j_accretion = .true.
         end if

      end function extras_binary_check_model


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer:: i_don, i_acc
	      real(dp) :: r_l2, d_l2
         integer :: ierr, star_id, i
         real(dp) :: q, mdot_limit_low, mdot_limit_high, &
                     center_h1, center_h1_old, center_he4, center_he4_old, &
                     rl23,rl2_1,trap_rad, mdot_edd, tau_kh, tau_macc, nova_flag, &
                     etaH, m_ig, t_rec, dt_yr, nova_mass_lost, mass_transferred, n_novae, &
                     mdot_msun_yr, m_wd_msun, log10mdot, nova_validity
         logical :: is_ne_biggest, superthermal_accretion
         real(dp) :: gamma1_integral, integral_norm, Pdm_over_rho

         extras_binary_finish_step = keep_going

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if

         if (b% point_mass_i == 0) then
            ! Check for simultaneous RLOF from both stars after TAMS of one star
            if (b% s2% center_h1 < 1.0d-6 .or. b% s1% center_h1 < 1.0d-6) then
                if (b% rl_relative_gap(1) > 0.0_dp .and. b% rl_relative_gap(2) > 0.0_dp) then
                  extras_binary_finish_step = terminate
                  write(*,'(g0)') "termination code: Both stars fill their Roche Lobe and at least one of them is off MS"
		            return
                end if
            end if
         end if


         !remove gradL_composition term after MS, it can cause the convective helium core to recede
         if (b% point_mass_i /= 1 .and. b% s1% center_h1 < 1.0d-6) then
            b% s1% num_cells_for_smooth_gradL_composition_term = 0
         end if
         if (b% point_mass_i /= 2 .and. b% s2% center_h1 < 1.0d-6) then
            b% s2% num_cells_for_smooth_gradL_composition_term = 0
         end if

         !check if mass transfer rate reached maximun, assume unstable regime if it happens
         if (abs(b% mtransfer_rate/(Msun/secyer)) >= 1d-1) then            !stop when larger than 0.1 Msun/yr
            extras_binary_finish_step = terminate
            write(*,'(g0)') "termination code: Reached maximum mass transfer rate: 1d-1"
	         return
         end if

         ! check trapping radius only for runs with a compact object
         if (b% point_mass_i == 2) then
            call my_mdot_edd(binary_id,mdot_edd,ierr)

            ! Begelman 1997 and King & Begelman 1999 eq. 1: accretor is star 2
            trap_rad = 0.5_dp*abs(b% mtransfer_rate) * acc_radius(b, b% m(2)) / mdot_edd

            !check if mass transfer rate reached maximun, assume unstable regime if it happens
            if (trap_rad >= b% rl(2)) then                                     !stop when trapping radius larger than rl(2)
               !if (abs(b% mtransfer_rate/(Msun/secyer)) >= 1d-1) then            !stop when larger than 0.1 Msun/yr
               extras_binary_finish_step = terminate
               write(*,'(g0)') "termination code: Reached maximum mass transfer rate: Exceeded photon trapping radius"
	            return
            end if
         end if

         ! check for termination due to carbon depletion
         if (b% point_mass_i /= 1) then
            if (b% s1% center_c12 < 1.0d-2 .and. b% s1% center_he4 < 1.0d-6) then
               write(*,'(g0)') "termination code: Primary has depleted central carbon"
               extras_binary_finish_step = terminate
               return
            !else
            !   ! check if neon is by far greatest source of energy
            !   is_ne_biggest = .true.
            !   do i=1, num_categories
            !      if(i /= i_burn_ne .and. b% s1% L_by_category(i_burn_ne) < 10*b% s1% L_by_category(i)) then
            !         is_ne_biggest = .false.
            !         exit
            !      end if
            !   end do
            !   if (is_ne_biggest .and. b% s1% max_eps_z_m/b% s1% xmstar > 0.01) then
            !         write(*,'(g0)') "offcenter neon ignition for primary at q=",  b% s1% max_eps_z_m/b% s1% xmstar, &
            !            b% s1% max_eps_z_m
            !         extras_binary_finish_step = terminate
            !         write(*,'(g0)') "termination code: offcenter neon ignition for primary"
            !   end if
            end if
         end if

         ! check for termination due to carbon depletion
         if (b% point_mass_i /= 2) then
            if (b% s2% center_c12 < 1.0d-2 .and. b% s2% center_he4 < 1.0d-6) then
               write(*,'(g0)') "termination code: Secondary has depleted central carbon"
               extras_binary_finish_step = terminate
               return
            !else
            !   ! check if neon is by far greatest source of energy
            !   is_ne_biggest = .true.
            !   do i=1, num_categories
            !      if(i /= i_burn_ne .and. b% s2% L_by_category(i_burn_ne) < 10._dp*b% s2% L_by_category(i)) then
            !         is_ne_biggest = .false.
            !         exit
            !      end if
            !   end do
            !   if (is_ne_biggest .and. b% s2% max_eps_z_m/b% s2% xmstar > 0.01_dp) then
            !         write(*,'(g0)') "offcenter neon ignition for secondary at q=",  b% s2% max_eps_z_m/b% s2% xmstar, &
            !            b% s2% max_eps_z_m
            !         extras_binary_finish_step = terminate
            !         write(*,'(g0)') "termination code: offcenter neon ignition for secondary"
            !   end if
            end if
         end if

         ! check for L2 overflow after ZAMS, but before TAMS as in Marchant et al. 2016
         if(.not. b% ignore_rlof_flag &
            .and. extras_binary_finish_step /= terminate &
            .and. (b% point_mass_i == 0)) then ! only when we evolve both stars in MS

            if (b% s1% center_h1 > 1d-6 .and. b% s2% center_h1 > 1d-6) then
               if (b% m(1) > b% m(2)) then
                  q = b% m(2) / b% m(1)
                  star_id = 2
               else
                  q = b% m(1) / b% m(2)
                  star_id = 1
               end if
               if (b% rl_relative_gap(star_id) > 0.29858997d0*atan_cr(1.83530121d0*pow_cr(q,0.39661426d0))) then
                  write(*,'(g0)') "termination code: Terminate due to L2 overflow during case A"
                  extras_binary_finish_step = terminate
		            return
               end if
            end if
         end if



         if (b% point_mass_i /= 1) then !Check for L2 overflow for primary when not in MS
            ! Misra et al. 2020 L2 overflow check starts only after TAMS of one of the two stars.
            ! Before we use Marchant et al. 2016 L2 overflow check implemented already in MESA
            if (b% s1% center_h1 < 1.0d-6) then  
               i_don = 1
               i_acc = 2
               if (b% m(i_don) .gt. b% m(i_acc)) then !mdon>macc, q<1
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don) * (0.784_dp * pow_cr(q,1.05_dp) * exp_cr(-0.188_dp*q) + 1.004_dp)
                  d_l2 = b% rl(i_don) * (3.334_dp * pow_cr(q, 0.514_dp) * exp_cr(-0.052_dp*q) + 1.308_dp)
                  !Condition to stop when star overflows L2
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 1'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 1'
                     return
                  end if
               ! mdonor < maccretor, Condition to stop when mass loss from L2 (previously it was L3) q>1
               else             
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don) * (0.29066811_dp * pow_cr(q, 0.82788069_dp) * exp_cr(-0.01572339_dp*q) + 1.36176161_dp)
                  d_l2 = b% rl(i_don) * (-0.04029713_dp * pow_cr(q, 0.862143_dp) * exp_cr(-0.04049814_dp*q) + 1.88325644_dp)
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 1'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 1'
                     return
                  end if
               end if
            end if
         end if

         if (b% point_mass_i /= 2) then  !Check for L2 overflow for primary when not in MS
            ! Misra et al. 2020 L2 overflow check starts only after TAMS of one of the two stars.
            ! Before we use Marchant et al. 2016 L2 overflow check implemented already in MESA
            if (b% s2% center_h1 < 1.0d-6) then  
               i_don = 2
               i_acc = 1
               ! mdon > macc, q < 1
               if (b% m(i_don) .gt. b% m(i_acc)) then 
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don)*(0.784_dp*pow_cr(q, 1.05_dp)*exp_cr(-0.188_dp*q) + 1.004_dp)
                  d_l2 = b% rl(i_don)*(3.334_dp*pow_cr(q,  0.514_dp)*exp_cr(-0.052_dp*q) + 1.308_dp)
                  !Condition to stop when star overflows L2
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)<1, donor is star 2'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)<1, donor is star 2'
                     return
                  end if
               ! mdonor<maccretor, Condition to stop when mass loss from L2 (previously it was L3) q>1
               else             
                  q = b% m(i_acc) / b% m(i_don)
                  r_l2 = b% rl(i_don)*(0.29066811_dp*pow_cr(q, 0.82788069_dp)*exp_cr(-0.01572339_dp*q) + 1.36176161_dp)
                  d_l2 = b% rl(i_don)*(-0.04029713_dp*pow_cr(q, 0.862143_dp)*exp_cr(-0.04049814_dp*q) + 1.88325644_dp)
                  if (b% r(i_don) .ge. (r_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (R_L2) surface for q(=Macc/Mdon)>1, donor is star 2'
                     return
                  end if
                  if (b% r(i_don) .ge. (d_l2)) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: overflow from L2 (D_L2) distance for q(=Macc/Mdon)>1, donor is star 2'
                     return
                  end if
               end if
            end if
         end if

         ! check for termination due to pair-instability in primary
         if (b% point_mass_i /= 1) then
            ! calculate volumetric pressure-weighted average adiabatic index -4/3, following Renzo et al. 2020
            integral_norm = 0.0d0
            gamma1_integral = 0.0d0
            do i=1,b% s1% nz
               Pdm_over_rho = b% s1% P(i)*b% s1% dm(i)/b% s1% rho(i)
               integral_norm = integral_norm + Pdm_over_rho
               gamma1_integral = gamma1_integral + &
                                 (b% s1% gamma1(i)-4.0d0/3.0d0)*Pdm_over_rho
            end do
            gamma1_integral = gamma1_integral/max(1.0d-99,integral_norm)
            if (gamma1_integral < 0.0d0) then
               ! check central value of adiabatic index to differentiate between full and pulsational pair-instability
               if (b% s1% gamma1(b% s1% nz)-4.0d0/3.0d0 < 0.0d0) then
                  write(*,'(g0)') "termination code: Primary enters pair-instability regime"
                  extras_binary_finish_step = terminate
                  return
               else
                  write(*,'(g0)') "termination code: Primary enters pulsational pair-instability regime"
                  extras_binary_finish_step = terminate
                  return
               end if
            end if
         end if

         ! check for termination due to pair-instability in secondary
         if (b% point_mass_i /= 2) then
            ! calculate volumetric pressure-weighted average adiabatic index -4/3, following Renzo et al. 2020
            integral_norm = 0.0d0
            gamma1_integral = 0.0d0
            do i=1,b% s2% nz
               Pdm_over_rho = b% s2% P(i)*b% s2% dm(i)/b% s2% rho(i)
               integral_norm = integral_norm + Pdm_over_rho
               gamma1_integral = gamma1_integral + &
                                 (b% s2% gamma1(i)-4.0d0/3.0d0)*Pdm_over_rho
            end do
            gamma1_integral = gamma1_integral/max(1.0d-99,integral_norm)
            if (gamma1_integral < 0.0d0) then
               ! check central value of adiabatic index to differentiate between full and pulsational pair-instability
               if (b% s2% gamma1(b% s2% nz)-4.0d0/3.0d0 < 0.0d0) then
                  write(*,'(g0)') "termination code: Secondary enters pair-instability regime"
                  extras_binary_finish_step = terminate
                  return
               else
                  write(*,'(g0)') "termination code: Secondary enters pulsational pair-instability regime"
                  extras_binary_finish_step = terminate
                  return
               end if
            end if
         end if

         if (extras_binary_finish_step == terminate) then
            !write(*,*) "saving final profilesA"
            !call star_write_profile_info(b% s1% id, "LOGS1/prof_9FINAL.data", b% s1% id, ierr)
            !if (ierr /= 0) return ! failure in profile
            !call star_write_profile_info(b% s2% id, "LOGS2/prof_9FINAL.data", b% s2% id, ierr)
            !if (ierr /= 0) return ! failure in profile
         else
            if (b% point_mass_i /= 1) then
               if (b% s1% center_h1 < 1d-6 .and. b% mdot_scheme .ne. "Kolb") then ! Changing from 'contact' scheme to Kolb if one star reaches TAMS
                  b% mdot_scheme = "Kolb"
                  write(*,*) "Primary reached TAMS, changing mdot_scheme to ", b% mdot_scheme, &
                             " and changing L2 overflow check according to Misra et al. 2020"
                  b% terminate_if_L2_overflow = .false.
               end if
            end if
            if (b% point_mass_i /= 2) then
               if (b% s2% center_h1 < 1d-6 .and. b% mdot_scheme .ne. "Kolb") then
                  b% mdot_scheme = "Kolb"
                  write(*,*) "Secondary reached TAMS, changing mdot_scheme to", b% mdot_scheme, &
                             " and changing L2 overflow check according to Misra et al. 2020"
                  b% terminate_if_L2_overflow = .false.
               end if
            end if
            !write(*,*) "still using: ", b% mdot_scheme

            !if (b% model_number == 1 ) then ! Saving initial_profiles
            !   write(*,*) "saving initial profiles"
            !   if (b% point_mass_i /= 1) then
            !        call star_write_profile_info(b% s1% id, "LOGS1/initial_profile.data", b% s1% id, ierr)
            !   end if
            !   if (ierr /= 0) return ! failure in profile
            !   if (b% point_mass_i /= 2) then
            !        call star_write_profile_info(b% s2% id, "LOGS2/initial_profile.data", b% s2% id, ierr)
            !   end if
            !   if (ierr /= 0) return ! failure in profile
            !end if
            if (b% model_number == 1 ) then ! Saving initial_models
               write(*,*) "saving initial models"
               if (b% point_mass_i /= 1) then
                  call star_write_model(b% s1% id, "initial_star1.mod", ierr)
               end if
               if (ierr /= 0) return ! failure
               if (b% point_mass_i /= 2) then
                  call star_write_model(b% s2% id, "initial_star2.mod", ierr)
               end if
               if (ierr /= 0) return ! failure
            end if
         end if

	      if (b% point_mass_i == 0) then
            if (b% s_accretor% x_logical_ctrl(4)) then
               if (b% s_accretor% w_div_w_crit_avg_surf >= 0.97d0 .and. b% d_i == 2) then
	               b% mass_transfer_beta = 1.0d0
                  b% s_accretor% max_wind = 1d-12
	            end if
	            if (b% mass_transfer_beta == 1.0d0 .and. abs(b% mtransfer_rate/(Msun/secyer)) <= 1d-7) then
	               b% mass_transfer_beta = 0d0
	               b% s_accretor% max_wind = 0d0
	            end if
            end if
	      end if

         ! conditions to check for termination in the case of superthermal accretion w/ contact and
         ! critical rotation + accretion disk
         if (b% point_mass_i /= b% a_i) then
            tau_macc = b% s_accretor% star_mass/abs(b% s_accretor% mstar_dot/Msun*secyer)
            tau_kh = b% s_accretor% kh_timescale
            ! (via experimentation, we can follow evolution up to about when this ratio is 1/50)
            ! Pols & Marinus 1994, A&A, 288, 475 cite effects kicking in as early as a ratio of 1/10
            if (tau_macc / tau_kh < 0.02d0) then
               superthermal_accretion = .true.
            else
               superthermal_accretion = .false.
            end if

            if (superthermal_accretion) then
               ! condition to check for a contact binary
               if (b% point_mass_i /= b% d_i) then
                  if ((b% r(b% d_i) .ge. b% rl(b% d_i)) .and. (b% r(b% a_i) .ge. b% rl(b% a_i))) then
                     extras_binary_finish_step = terminate
                     write(*,'(g0)') 'termination code: Both stars fill their Roche Lobe and t_kh > t_acc'
                     return
                  end if
               end if

               ! check if accretor is accreting at a superthermal rate and critically rotating 
               ! (decretion + expansion due to rapid accretion, here we are assuming this leads to L2 overflow)
               if (superthermal_accretion .and. &
                  (b% s_accretor% w_div_w_crit_avg_surf >= 0.99d0*b% s_accretor% surf_w_div_w_crit_limit)) then

                  ! terminate as a case of L2 overflow
                  extras_binary_finish_step = terminate
                  if (b% d_i == 1) then
                     write(*,'(g0)') 'termination code: overflow from L2, t_kh > t_acc and w > w_crit_lim, donor is star 1'
                  else
                     write(*,'(g0)') 'termination code: overflow from L2, t_kh > t_acc and w > w_crit_lim, donor is star 2'
                  end if
		  
                  return
		  
               end if
            end if
         end if 
!**********************************************************************!        
!        Nova finish_step
         if (b% point_mass_i /= 0 .and. b% mtransfer_rate /= 0d0 .and. &
            b% m(b% point_mass_i) < nova_m_wd_max * Msun) then

            b% s_donor% xtra7 = 0d0

            call my_nova_flag(binary_id, nova_flag, ierr)
            if (ierr /= 0) return

            ! Same in_nova_regime check used for the 'nova_validity' history column.
            mdot_msun_yr = abs(b% mtransfer_rate) * secyer / Msun
            m_wd_msun    = b% m(b% point_mass_i) / Msun
            log10mdot    = log10_cr(mdot_msun_yr)
            if (in_nova_regime(m_wd_msun, log10mdot)) then
               nova_validity = 1d0
            else
               nova_validity = 0d0
            end if

            if (nova_flag == 1d0 .and. nova_validity == 1d0) then
               call my_eta_accretion(binary_id, etaH, ierr)
               if (ierr /= 0) return
               call my_dm_ig(binary_id, m_ig, ierr)
               if (ierr /= 0) return
               b% s_donor% xtra3 = etaH
               b% s_donor% xtra4 = m_ig

               if (b% mass_transfer_beta > 0d0) then
                  call my_t_recurrence(binary_id, t_rec, ierr)
                  if (ierr /= 0) return
                  dt_yr = b% s_donor% dt / secyer
                  if (t_rec > 0d0) then
                     n_novae = dt_yr / t_rec
                     b% s_donor% xtra7 = n_novae
                  end if
               end if

               ! Reset WD mass reference after nova fires, deferred here so history captures beta > 0.
               if (b% mass_transfer_beta > 0d0) &
                  b% s_donor% xtra6 = b% m(b% point_mass_i)
            end if
            ! Note: eta_retention (xtra3) / m_ignition (xtra4) are ONLY
            ! written above, when nova_flag = 1 and nova_validity = 1. Any
            ! other case (flag/validity 0, or beta = 0 this step) leaves
            ! them untouched so they keep the last physically computed
            ! value instead of being zeroed.

            ! Accumulate mass budget [Msun] — applies in nova and
            ! stable-burning regimes alike. nova_mass_lost = beta *
            ! mass_transferred is 0 exactly when beta = 0, which is always
            ! the case when nova_flag = 0 and nova_validity = 0 (see
            ! extras_binary_start_step); nothing else forces it to 0.
            mass_transferred  = abs(b% mtransfer_rate) * b% s_donor% dt / Msun
            nova_mass_lost    = b% mass_transfer_beta * mass_transferred
            b% s_donor% xtra5 = nova_mass_lost
            b% s_donor% xtra8 = b% s_donor% xtra8 + nova_mass_lost
            b% s_donor% xtra9 = b% s_donor% xtra9 + &
                                 (1d0 - b% mass_transfer_beta) * mass_transferred
         else
            ! No active point-mass mass transfer this step: mass_transfer_beta
            ! is 0 here (see extras_binary_start_step), so nova_mass_lost
            ! correctly reads 0. eta_retention/m_ignition (xtra3, xtra4) are
            ! left untouched so they keep their last computed value instead
            ! of being zeroed.
            b% s_donor% xtra5 = 0d0
            b% s_donor% xtra7 = 0d0
         end if
!************************************TERMNATION CONDITONS************************************        
         ! Stop once the WD reaches the Chandrasekhar mass (1.38 Msun, as
         ! in Kemp+2021). Checked after the nova bookkeeping so the final
         ! history row is consistent.

         !Here we need to add a termination condition for Sub-Chandraskehar detonation
         if (b% point_mass_i /= 0) then
            if (b% m(b% point_mass_i) >= 1.38d0 * Msun) then
               extras_binary_finish_step = terminate
               write(*,'(g0)') 'termination code: WD reached Chandrasekhar mass (1.38 Msun)'
               return
            end if
         end if



!***************************************************************************!
      end function extras_binary_finish_step

      real(dp) function eval_rlobe(m1, m2, a) result(rlobe)
         real(dp), intent(in) :: m1, m2, a
         real(dp) :: q
         q = pow_cr(m1/m2,one_third)
         ! Roche lobe size for star of mass m1 with a
         ! companion of mass m2 at separation a, according to
         ! the approximation of Eggleton 1983, apj 268:368-369
         rlobe = a*0.49d0*q*q/(0.6d0*q*q + log1p_cr(q))
      end function eval_rlobe

      subroutine extras_binary_after_evolve(binary_id, ierr)
         use run_star_support
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         type(star_Info), pointer :: s
         integer :: iounit
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) return
         !if (b% point_mass_i /= 1) then
         call star_write_profile_info(b% s1% id, "LOGS1/final_profile.data", b% s1% id, ierr)
         !end if
         if (ierr /= 0) return ! failure in profile
         if (b% point_mass_i /= 2) then
              call star_write_profile_info(b% s2% id, "LOGS2/final_profile.data", b% s2% id, ierr)
         end if
         if (ierr /= 0) return ! failure in profile

      end subroutine extras_binary_after_evolve

      end module run_binary_extras




