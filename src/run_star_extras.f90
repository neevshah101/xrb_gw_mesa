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

module run_star_extras

  use star_lib
  use star_def
  use const_def
  use math_lib
  use chem_def
  use num_lib
  use binary_def
  use ionization_def

  implicit none

  ! these routines are called by the standard run_star check_model
contains


  subroutine extras_controls(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! this is the place to set any procedure pointers you want to change
    ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


    ! the extras functions in this file will not be called
    ! unless you set their function pointers as done below.
    ! otherwise we use a null_ version which does nothing (except warn).

    s% extras_startup => extras_startup
    s% extras_start_step => extras_start_step
    s% extras_check_model => extras_check_model
    s% extras_finish_step => extras_finish_step
    s% extras_after_evolve => extras_after_evolve
    s% how_many_extra_history_columns => how_many_extra_history_columns
    s% data_for_extra_history_columns => data_for_extra_history_columns
    s% how_many_extra_profile_columns => how_many_extra_profile_columns
    s% data_for_extra_profile_columns => data_for_extra_profile_columns

    s% how_many_extra_history_header_items => how_many_extra_history_header_items
    s% data_for_extra_history_header_items => data_for_extra_history_header_items
    s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
    s% data_for_extra_profile_header_items => data_for_extra_profile_header_items
    s% other_wind => my_wind

  end subroutine extras_controls


  subroutine extras_startup(id, restart, ierr)
    integer, intent(in) :: id
    logical, intent(in) :: restart
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! use 11-16 to avoid overlap with lxtra(1) and lxtra(2) used in run_binary_extras.f90
    if (.not. restart) then
       s% lxtra(11) = .false.
       s% lxtra(12) = .false.
       s% lxtra(13) = .false.
       s% lxtra(14) = .false.
       s% lxtra(15) = .false.
       s% lxtra(16) = .false.
       s% lxtra(17) = .false.
       s% lxtra(18) = .false.
       s% lxtra(19) = .false.
       s% lxtra(20) = .false.
       s% lxtra(21) = .false.
    end if
  end subroutine extras_startup


  integer function extras_start_step(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_start_step = 0

  end function extras_start_step


  ! returns either keep_going, retry, backup, or terminate.
  integer function extras_check_model(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    real(dp) :: error, atol, rtol
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going

    ! check we don't overshoot too much the radii of interest
    if ((s% lxtra(11) .eqv. .false.) .and. (s% r(1)/Rsun >= 20.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 20.0d0)
        
        if (error > 0.1d0) then
            extras_check_model = retry
        end if
        
    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(12) .eqv. .false.) .and. (s% r(1)/Rsun >= 25.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 25.0d0)
        
        if (error > 1d0) then
            print*, "error", error
            extras_check_model = retry
        end if
        
    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(13) .eqv. .false.) .and. (s% r(1)/Rsun >= 30.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 30.0d0)
        
        if (error > 1d0) then
            print*, "error", error
            extras_check_model = retry
            
        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(14) .eqv. .false.) .and. (s% r(1)/Rsun >= 50.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 50.0d0)
        
        if (error > 1d0) then
            print*, "error", error
            extras_check_model = retry
            
        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(15) .eqv. .false.) .and. (s% r(1)/Rsun >= 100.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 100.0d0)
        
        if (error > 5d0) then
            print*, "error", error
            extras_check_model = retry
            
        end if

    end if

    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(16) .eqv. .false.) .and. (s% r(1)/Rsun >= 500.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 500.0d0)
        
        if (error > 10d0) then
            print*, "error", error
            extras_check_model = retry
            
        end if

    end if
    
    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(17) .eqv. .false.) .and. (s% r(1)/Rsun >= 1000.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 1000.0d0)
        
        if (error > 20d0) then
            print*, "error", error
            extras_check_model = retry
            
        end if

    end if
    
    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(18) .eqv. .false.) .and. (s% r(1)/Rsun >= 1200.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 1200.0d0)
        
        if (error > 20d0) then
            print*, "error", error
            extras_check_model = retry
            
        end if

    end if
    
    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(19) .eqv. .false.) .and. (s% r(1)/Rsun >= 1300.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 1300.0d0)
        
        if (error > 20d0) then
            print*, "error", error
            extras_check_model = retry
            
        end if

    end if
    
    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(20) .eqv. .false.) .and. (s% r(1)/Rsun >= 1400.0d0)) then

        error = 1d0*(s%r(1)/Rsun - 1400.0d0)
        
        if (error > 20d0) then
            print*, "error", error
            extras_check_model = retry
            
        end if

    end if
    
    !check we don't overshoot too much the radii of interest
    if ((s% lxtra(21) .eqv. .false.) .and. (s% r(1)/Rsun >= 1500.0d0)) then


        error = 1d0*(s%r(1)/Rsun - 1500.0d0)
        
        if (error > 20d0) then
            print*, "error", error
            extras_check_model = retry
            
        end if
        
    end if

    ! by default, indicate where (in the code) MESA terminated
    if (extras_check_model == terminate) s% termination_code = t_extras_check_model
  end function extras_check_model


  integer function how_many_extra_history_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_columns = 0
  end function how_many_extra_history_columns


  subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
    use chem_def, only: chem_isos
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    integer, intent(out) :: ierr
    ! -------------------------------------
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

  end subroutine data_for_extra_history_columns


  integer function how_many_extra_profile_columns(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_columns = 0
  end function how_many_extra_profile_columns


  subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
    integer, intent(in) :: id, n, nz
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(nz,n)
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    integer :: k
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! note: do NOT add the extra names to profile_columns.list
    ! the profile_columns.list is only for the built-in profile column options.
    ! it must not include the new column names you are adding here.

    ! here is an example for adding a profile column
    !if (n /= 1) stop 'data_for_extra_profile_columns'
    !names(1) = 'beta'
    !do k = 1, nz
    !   vals(k,1) = s% Pgas(k)/s% P(k)
    !end do

  end subroutine data_for_extra_profile_columns


  integer function how_many_extra_history_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_history_header_items = 3
  end function how_many_extra_history_header_items


  subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_history_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    integer :: i
    real(dp) :: Initial_X, Initial_Y, Initial_Z, initial_m

    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return

  end subroutine data_for_extra_history_header_items


  integer function how_many_extra_profile_header_items(id)
    integer, intent(in) :: id
    integer :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    how_many_extra_profile_header_items = 0
  end function how_many_extra_profile_header_items


  subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
    integer, intent(in) :: id, n
    character (len=maxlen_profile_column_name) :: names(n)
    real(dp) :: vals(n)
    type(star_info), pointer :: s
    integer, intent(out) :: ierr
    ierr = 0
    call star_ptr(id,s,ierr)
    if(ierr/=0) return

    ! here is an example for adding an extra profile header item
    ! also set how_many_extra_profile_header_items
    ! names(1) = 'mixing_length_alpha'
    ! vals(1) = s% mixing_length_alpha

  end subroutine data_for_extra_profile_header_items


  ! returns either keep_going or terminate.
  ! note: cannot request retry or backup; extras_check_model can do that.
  integer function extras_finish_step(id)
    integer, intent(in) :: id
    integer :: ierr, i
    real(dp) :: envelope_mass_fraction, L_He, L_tot, min_center_h1_for_diff, &
         critmass, feh, rot_full_off, rot_full_on, frac2, TAMS_h1_treshold
    real(dp), parameter :: huge_dt_limit = 3.15d16 ! ~1 Gyr
    real(dp), parameter :: new_varcontrol_target = 1d-3
    real(dp), parameter :: Zsol = 0.0142_dp
    logical :: diff_test1, diff_test2, diff_test3, is_ne_biggest
    !character (len=strlen) :: photoname, fname
    character (len=30) :: photoname, fname
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_finish_step = keep_going

    if (s%lxtra(11) .eqv. .false.) then
       ! save profile for R=20Rsun
       if (s% r(1)/Rsun >= 20) then
          s% lxtra(11) = .true.
          write(fname, fmt="(a11)") '20Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a10)") '20Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a6)") '20Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(12) .eqv. .false.) then
       ! save profile for R=25Rsun
       if (s% r(1)/Rsun >= 25) then
          s% lxtra(12) = .true.
          write(fname, fmt="(a12)") '25Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '25Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '25Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(13) .eqv. .false.) then
       ! save profile for R=30Rsun
       if (s% r(1)/Rsun >= 30) then
          s% lxtra(13) = .true.
          write(fname, fmt="(a12)") '30Rsun.data'
          print*, "saving profile "// fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '30Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '30Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(14) .eqv. .false.) then
       ! save profile for R=50Rsun
       if (s% r(1)/Rsun >= 50) then
          s% lxtra(14) = .true.
          write(fname, fmt="(a12)") '50Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '50Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '50Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(15) .eqv. .false.) then
       ! save profile for R=100Rsun
       if (s% r(1)/Rsun >= 100) then
          s% lxtra(15) = .true.
          write(fname, fmt="(a12)") '100Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a11)") '100Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a7)") '100Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (s%lxtra(16) .eqv. .false.) then
       ! save profile for R=500Rsun
       if (s% r(1)/Rsun >= 500) then
          s% lxtra(16) = .true.
          write(fname, fmt="(a13)") '500Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '500Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '500Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if
    
    if (s%lxtra(17) .eqv. .false.) then
       ! save profile for R=1000Rsun
       if (s% r(1)/Rsun >= 1000) then
          s% lxtra(17) = .true.
          write(fname, fmt="(a13)") '1000Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1000Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1000Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if
    
    if (s%lxtra(18) .eqv. .false.) then
       ! save profile for R=1200Rsun
       if (s% r(1)/Rsun >= 1200) then
          s% lxtra(18) = .true.
          write(fname, fmt="(a13)") '1200Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1200Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1200Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if
    
    if (s%lxtra(19) .eqv. .false.) then
       ! save profile for R=1300Rsun
       if (s% r(1)/Rsun >= 1300) then
          s% lxtra(19) = .true.
          write(fname, fmt="(a13)") '1300Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1300Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1300Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if
    
    if (s%lxtra(20) .eqv. .false.) then
       ! save profile for R=1400Rsun
       if (s% r(1)/Rsun >= 1400) then
          s% lxtra(20) = .true.
          write(fname, fmt="(a13)") '1400Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1400Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1400Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if
    
    if (s%lxtra(21) .eqv. .false.) then
       ! save profile for R=1500Rsun
       if (s% r(1)/Rsun >= 1500) then
          s% lxtra(21) = .true.
          write(fname, fmt="(a13)") '1500Rsun.data'
          print*, "saving profile "//fname
          call star_write_profile_info(id, trim(s%log_directory)//'/'//trim(fname), ierr)
          write(fname, fmt="(a12)") '1500Rsun.mod'
          call star_write_model(id, trim(fname), ierr)
          write(fname, fmt="(a8)") '1500Rsun'
          call star_save_for_restart(id, trim(s%photo_directory)//'/'//trim(fname), ierr)
       end if
    end if

    if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
  end function extras_finish_step


  subroutine extras_after_evolve(id, ierr)
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type (star_info), pointer :: s
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
  end subroutine extras_after_evolve
  

    subroutine my_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
   use star_def
   type (star_info), pointer :: s
   integer, intent(in) :: id
   real(dp), parameter :: Zbjork = 0.014d0, T_high = 11000, T_low = 10000
   real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
   real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
   integer, intent(out) :: ierr
   real(dp) :: M1, L1, T1, R1, w1, w2, alfa
   w = 0
   ierr = 0
   call star_ptr(id, s, ierr)

   M1 = Msurf 
   L1 = Lsurf
   T1 = Tsurf

   if (T1 <= T_low) then
      call eval_lowT(w)
      !print*, "Decin"

   else if (T1 >= T_high) then 
        call eval_highT(w)

   else ! transition
      call eval_lowT(w1)
      call eval_highT(w2)
      alfa = (T1 - T_low) / (T_high - T_low)
      w = (1-alfa)*w1 + alfa*w2
   end if

contains

    subroutine eval_de_Jager_wind(w4)
       real(dp), intent(out) :: w4
       real(dp) :: log10w

       log10w = 1.769d0*log10(L1/Lsun) - 1.676d0*log10(T1) - 8.158d0
       w4 = exp10(log10w)
   end subroutine eval_de_Jager_wind

   subroutine eval_Decin_wind(w)
      ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
      real(dp), intent(out) :: w
      real(dp) :: log10w
      include 'formats'
      !log10w = -17.44 &
      !         +2.80*log10(L1/Lsun) &
      !         -0.13*s% initial_mass! &
      log10w = -20.71 &
               +3.50*log10(L1/Lsun) &
               -0.17*s% initial_mass+0.49! &
      !write(*,*) "Check Decin wind", log10w 
      w = 10**(log10w)
   end subroutine eval_Decin_wind

  subroutine eval_lowT(w)
        real(dp), intent(out) :: w
        if (s% x_character_ctrl(2) == 'Decin') then
                call eval_Decin_wind(w)
        end if

        if (s% x_character_ctrl(2) == 'nowind') then
                call eval_Decin_wind(w1)
                w = min(w1, 1d-6)
        end if
        
   end subroutine eval_lowT

   subroutine eval_highT(w5)
      real(dp), intent(out) :: w5

      if (X < 0.4d0) then
         w5 = 1d-11 * pow(L1/Lsun, 1.29d0) * pow(Y, 1.7d0) * sqrt(Z)
         !print*, "Nugis & Lamers"
      else if (s% x_character_ctrl(1) == 'Bjorklund') then
         call eval_Bjorklund_wind(w5)
      else if (s% x_character_ctrl(1) == 'Vink') then
         call eval_Vink_wind(w5)
      else if (s% x_character_ctrl(1) == 'Sabhahit') then
        call eval_Sabhahit_wind(w5)
      end if
   end subroutine eval_highT

   subroutine eval_Bjorklund_wind(w6)
      real(dp), intent(out) :: w6
      real(dp) :: logw, Meff

      Meff = M1 * (1d0 - 0.34d0 * L1 / (pi4 * clight * s% cgrav(1) * M1)) ! effective mass

      logw = -5.52d0 &
           + 2.39d0 * log10(L1 / (1d6 * Lsun)) &
           - 1.48d0 * log10(Meff / (4.5d1 * Msun)) &
           + 2.12d0 * log10(T1 / 4.5d4) &
           + (0.75d0 - 1.87d0 * log10(T1 / 4.5d4)) * log10(Z / Zbjork)
      w6 = exp10(logw)
   end subroutine eval_Bjorklund_wind

    subroutine eval_Vink_wind(w)
    real(dp), intent(inout) :: w
    real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc
    real(dp), parameter :: Zsolar = 0.019d0 ! for Vink et al formula

    ! alfa = 1 for hot side, = 0 for cool side
    if (T1 > 27500d0) then
       alfa = 1
    else if (T1 < 22500d0) then
       alfa = 0
    else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
       Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
       dT = 100d0
       if (T1 > Teff_jump + dT) then
          alfa = 1
       else if (T1 < Teff_jump - dT) then
          alfa = 0
       else
          alfa = 0.5d0*(T1 - (Teff_jump - dT)) / dT
       end if
    end if

    if (alfa > 0) then ! eval hot side wind (eqn 24)
       vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
       vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
       logMdot = &
          - 6.697d0 &
          + 2.194d0*log10(L1/Lsun/1d5) &
          - 1.313d0*log10(M1/Msun/30d0) &
          - 1.226d0*log10(vinf_div_vesc/2d0) &
          + 0.933d0*log10(T1/4d4) &
          - 10.92d0*pow2(log10(T1/4d4)) &
          + 0.85d0*log10(Z/Zsolar)
       w1 = exp10(logMdot)
    else
       w1 = 0
    end if

    if (alfa < 1) then ! eval cool side wind (eqn 25)
       vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
       vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
       logMdot = &
          - 6.688d0 &
          + 2.210d0*log10(L1/Lsun/1d5) &
          - 1.339d0*log10(M1/Msun/30d0) &
          - 1.601d0*log10(vinf_div_vesc/2d0) &
          + 1.07d0*log10(T1/2d4) &
          + 0.85d0*log10(Z/Zsolar)
       w2 = exp10(logMdot)
    else
       w2 = 0
    end if

    w = alfa*w1 + (1 - alfa)*w2

   ! if (dbg) write(*,*) 'vink wind', w

 end subroutine eval_Vink_wind

subroutine eval_Sabhahit_wind(w)
            real(dp), intent(inout) :: w
            real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc
            real(dp), parameter :: Zsolar = 0.019d0 ! for Vink et al formula
            ! alfa = 1 for hot side, = 0 for cool side
            if (T1 > 27500d0) then
               alfa = 1
            else if (T1 < 22500d0) then
               alfa = 0
            else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
               Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
               dT = 100d0
               if (T1 > Teff_jump + dT) then
                  alfa = 1
               else if (T1 < Teff_jump - dT) then
                  alfa = 0
               else
                  alfa = 0.5d0*(T1 - (Teff_jump - dT)) / dT
               end if
            end if
            if (alfa > 0) then ! eval hot side wind (eqn 10 in Sabhahit)
               vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
               logMdot = &
                  - 6.697d0 &
                  + 2.194d0*log10(L1/Lsun/1d5) &
                  - 1.313d0*log10(M1/Msun/30d0) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.933d0*log10(T1/4d4) &
                  - 10.92d0*pow2(log10(T1/4d4)) &
                  + 0.85d0*log10(Z/Zsolar)
               w1 = exp10(logMdot)
            else
               w1 = 0
            end if
            if (alfa < 1) then ! eval cool side wind (eqn 25)
               vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
               vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
               logMdot = &
                  - 9.552 &
                  + 4.77d0*log10(1+X) &
                  + 4.77d0*log10(L1/Lsun/1d5) &
                  - 3.99d0*log10(M1/Msun/30d0) &
                  - 1.226d0*log10(vinf_div_vesc/2d0) &
                  + 0.5d0*log10(Z/Zsolar)
               w2 = exp10(logMdot)
            else
               w2 = 0
            end if
            w = alfa*w1 + (1 - alfa)*w2
      !      if (dbg) write(*,*) 'Sabhahit wind', w
         end subroutine eval_Sabhahit_wind

end subroutine my_wind

      


  !include 'POSYDON_single_stars.inc'

end module run_star_extras
