program dft
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 300)

    ! Constants
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp), two_pi = 2.0_dp * pi

    ! Signal parameters
    integer, parameter :: N = 1000          ! Number of samples
    real(dp), parameter :: dt = 0.01_dp      ! Time step
    real(dp), parameter :: amp = 1.0_dp      ! Signal amplitude

    ! Window parameters
    integer :: window_type = 0               ! Window type
    real(dp), parameter :: gaussian_width = 0.3_dp ! Width for the Gaussian window

    ! Leakage control parameters
    integer :: cycle_numerator = 7, cycle_denominator = 2
    real(dp) :: period = 0.0_dp             ! Sine wave period
    real(dp) :: fractional_cycles = 0.0_dp  ! Number of cycles in window
    logical :: auto_period = .true.          ! Flag: auto calculate period if true

    integer :: signal_type = 0
    real(dp) :: pulse_width = 0.0_dp      ! For rectangle and triangle pulses
    real(dp) :: sigma_width = 0.0_dp      ! For Gaussian pulse

    ! Arrays for signal processing
    complex(dp), dimension(0:N-1) :: signal         ! Time domain signal
    complex(dp), dimension(0:N-1) :: dft_result     ! Frequency domain DFT result
    real(dp), dimension(0:N-1)    :: window_function  ! Window function values

    integer :: ios, leakage_choice       ! I/O status, user choice
    character(len=100) :: output_prefix     ! Prefix for output file names

    ! Get Signal Type
    write(*,*) "Select signal type:"
    write(*,*) " 0 - Sine wave (uses period/leakage settings)"
    write(*,*) " 1 - Rectangular pulse"
    write(*,*) " 2 - Triangle pulse"
    write(*,*) " 3 - Gaussian pulse"
    read(*,*, iostat=ios) signal_type
    ! Input Validation
    if (ios /= 0) stop "Error: Invalid numeric input for signal type."
    if (signal_type < 0 .or. signal_type > 3) then
        write(*,*) "Error: Signal type must be between 0 and 3."
        stop 1
    end if

    ! Branch based on signal type
    if (signal_type == 0) then
        ! For sine wave only ask about windowing and period
        write(*,*) "Window types: 0 - none, 1 - triangle, 2 - cosine(Hann), 3 - gaussian"
        write(*,*) "Enter window type (0-3): "
        read(*,*, iostat=ios) window_type
        ! Input Validation
        if (ios /= 0) stop "Error: Invalid numeric input for window type."
        if (window_type < 0 .or. window_type > 3) then
            write(*,*) "Error: Window type must be between 0 and 3."
            stop 2
        end if

        ! Get Leakage Configuration
        write(*,*) "Sine Wave Spectral leakage configurations:"
        write(*,*) " 1 - Maximally incommensurate (e.g., 7.5 cycles)"
        write(*,*) " 2 - Moderately incommensurate (e.g., 7.25 cycles)"
        write(*,*) " 3 - Commensurate (e.g., 7.0 cycles)"
        write(*,*) " 4 - Custom period value"
        write(*,*) "Enter configuration (1-4): "
        read(*,*, iostat=ios) leakage_choice
        ! Input Validation
        if (ios /= 0) stop "Error: Invalid numeric input for leakage configuration."
        if (leakage_choice < 1 .or. leakage_choice > 4) then
            write(*,*) "Error: Leakage configuration choice must be between 1 and 4."
            stop 3
        end if

        select case(leakage_choice)
        case(1)
            cycle_numerator = 15 ; cycle_denominator = 2  ! 7.5 cycles
        case(2)
            cycle_numerator = 29 ; cycle_denominator = 4  ! 7.25 cycles
        case(3)
            cycle_numerator = 7  ; cycle_denominator = 1  ! 7.0 cycles
        case(4)
            auto_period = .false.
            write(*,*) "Enter custom period value (must be > 0): "
            read(*,*, iostat=ios) period
            ! Input Validation
            if (ios /= 0) stop "Error: Invalid numeric input for period."
            if (period <= 0.0_dp) then
                write(*,*) "Error: Custom period must be positive."
                stop 4
            end if
        end select

        ! Calculate period automatically if needed
        if (auto_period) then
            fractional_cycles = real(cycle_numerator, dp) / real(cycle_denominator, dp)
            period = real(N, dp) * dt / fractional_cycles
            write(*,*) "Configured for", fractional_cycles, "cycles in window"
            write(*,*) "Calculated period:", period
            if (period <= 0.0_dp) stop "Error: Calculated period is not positive. Check N, dt, cycles."
        end if
    else
        ! For pulse signals no window is needed
        window_type = 0

        ! Get parameters for pulse types
        select case(signal_type)
        case(1,2) ! Rectangular or Triangle pulse
            write(*,*) "Enter pulse width L (must be > 0): "
            read(*,*, iostat=ios) pulse_width
            ! Input Validation
            if (ios /= 0) stop "Error: Invalid input for pulse width."
            if (pulse_width <= 0.0_dp) then
                write(*,*) "Error: Pulse width L must be positive."
                stop 5
            end if
        case(3) ! Gaussian pulse
            write(*,*) "Enter Gaussian pulse width sigma (must be > 0): "
            read(*,*, iostat=ios) sigma_width
            ! Input Validation
            if (ios /= 0) stop "Error: Invalid numeric input for width."
            if (sigma_width <= 0.0_dp) then
                write(*,*) "Error: Gaussian sigma width must be positive."
                stop 6
            end if
        end select
    end if

    ! Processing Steps
    call update_output_prefix(output_prefix, signal_type, window_type, &
                              auto_period, cycle_numerator, cycle_denominator, &
                              period, pulse_width, sigma_width)

    call create_window_function(window_function, N, window_type, output_prefix)

    call generate_signal(signal, N, dt, signal_type, period, pulse_width, sigma_width, output_prefix)

    ! Apply window only for sine wave
    if (signal_type == 0) then
         call apply_window(signal, window_function, N)
    end if

    call calculate_dft(dft_result, signal, N)

    call output_results(dft_result, N, dt, output_prefix)

    write(*,*) "Processing complete. Files saved with prefix: ", trim(output_prefix)

contains

    ! Generates the filename prefix based on configuration
    subroutine update_output_prefix(output_prefix, signal_type, window_type, &
                                    auto_period, cycle_numerator, cycle_denominator, &
                                    period, pulse_width, sigma_width)
        character(len=*), intent(out) :: output_prefix
        integer, intent(in) :: signal_type, window_type
        logical, intent(in) :: auto_period
        integer, intent(in) :: cycle_numerator, cycle_denominator
        real(dp), intent(in) :: period, pulse_width, sigma_width
        character(len=20) :: window_str, cycle_str, signal_str, param_str

        select case(signal_type)
        case(0); signal_str = "sine"
        case(1); signal_str = "rect"
        case(2); signal_str = "tri"
        case(3); signal_str = "gauss_pulse"
        case default; signal_str = "unknown" ! Should not be reached
        end select

        select case(window_type)
        case(0); window_str = "none"
        case(1); window_str = "triangle"
        case(2); window_str = "cosine"
        case(3); window_str = "gaussian_win"
        case default; window_str = "unknown" ! Should not be reached
        end select

        if (signal_type == 0) then
            if (auto_period) then
                write(cycle_str, '(I0,"_",I0)') cycle_numerator, cycle_denominator
                param_str = "c" // trim(adjustl(cycle_str))
            else
                write(cycle_str, '(G10.3E2)') period ! Use scientific notation for consistency
                param_str = "p" // trim(adjustl(cycle_str))
            end if
            output_prefix = trim(signal_str) // "_" // trim(param_str) // "_w" // trim(window_str)
        else
            if (signal_type == 1 .or. signal_type == 2) then
                write(cycle_str, '(G10.3E2)') pulse_width ! Use scientific notation
                param_str = "_L" // trim(adjustl(cycle_str))
            else if (signal_type == 3) then
                write(cycle_str, '(G10.3E2)') sigma_width ! Use scientific notation
                param_str = "_sigma" // trim(adjustl(cycle_str))
            else
                 param_str = "" ! Should not happen
            end if
            output_prefix = trim(signal_str) // trim(param_str) // "_w" // trim(window_str)
        end if
    end subroutine update_output_prefix

    ! Computes the specified window function values
    subroutine create_window_function(window_function, n_samples, window_type, output_prefix)
        integer, intent(in) :: n_samples, window_type
        real(dp), dimension(0:n_samples-1), intent(out) :: window_function
        character(len=*), intent(in) :: output_prefix
        integer :: i, file_unit, ios
        real(dp) :: x, center_x
        character(len=120) :: filename

        do i = 0, n_samples-1
            if (n_samples > 1) then
                 x = real(i, dp) / real(n_samples-1, dp)
            else
                 x = 0.5_dp ! Handle n_samples=1 case
            end if

            select case(window_type)
            case(0) ! Rectangular (None)
                window_function(i) = 1.0_dp
            case(1) ! Triangle
                window_function(i) = 1.0_dp - 2.0_dp * abs(x - 0.5_dp)
            case(2) ! Cosine (Hann)
                window_function(i) = 0.5_dp * (1.0_dp - cos(two_pi * x))
            case(3) ! Gaussian
                center_x = 0.5_dp
                window_function(i) = exp(-((x - center_x)**2) / (2.0_dp * gaussian_width**2))
            case default
                ! Should not be reached due to input validation
                window_function(i) = 1.0_dp
            end select
        end do

        ! Output window function to file (optional for pulses, but done for consistency)
        filename = trim(output_prefix)//"_window.csv"
        open(newunit=file_unit, file=filename, status="replace", action="write", iostat=ios)
        if (ios /= 0) stop "Error opening window file."
        write(file_unit, '(A)') "index,x,window_value"
        do i = 0, n_samples-1
             if (n_samples > 1) then
                 x = real(i, dp) / real(n_samples-1, dp)
             else
                 x = 0.5_dp
             end if
            write(file_unit, '(I0,",",F12.5,",",F12.5)') i, x, window_function(i)
        end do
        close(file_unit)
    end subroutine create_window_function

    ! Generates the specified time domain signal without window
    subroutine generate_signal(signal, n_samples, dt_step, signal_type, &
                               period, pulse_width, sigma_width, output_prefix)
        integer, intent(in) :: n_samples, signal_type
        complex(dp), dimension(0:n_samples-1), intent(out) :: signal
        real(dp), intent(in) :: dt_step, period, pulse_width, sigma_width
        character(len=*), intent(in) :: output_prefix
        integer :: i, file_unit, ios
        real(dp) :: t, center_time, value
        character(len=120) :: filename

        center_time = (real(n_samples-1, dp) * dt_step) / 2.0_dp
        signal = cmplx(0.0_dp, 0.0_dp, kind=dp) ! Initialize signal

        select case(signal_type)
        case(0) ! Sine wave
             if (period <= 0.0_dp) stop "Error in generate_signal: Sine wave period must be positive."
            do i = 0, n_samples-1
                t = real(i, dp) * dt_step
                signal(i) = cmplx(amp * sin(two_pi * t / period), 0.0_dp, kind=dp)
            end do
        case(1) ! Rectangular pulse
            if (pulse_width <= 0.0_dp) stop "Error in generate_signal: Rectangular pulse width must be positive."
            do i = 0, n_samples-1
                t = real(i, dp) * dt_step
                if (abs(t - center_time) <= pulse_width / 2.0_dp) then
                    signal(i) = cmplx(amp, 0.0_dp, kind=dp)
                ! else signal remains 0
                end if
            end do
        case(2) ! Symmetric triangle pulse
            if (pulse_width <= 0.0_dp) stop "Error in generate_signal: Triangle pulse width must be positive."
            do i = 0, n_samples-1
                t = real(i, dp) * dt_step
                if (abs(t - center_time) <= pulse_width / 2.0_dp) then
                    value = amp * (1.0_dp - 2.0_dp * abs(t - center_time) / pulse_width)
                    signal(i) = cmplx(value, 0.0_dp, kind=dp)
                ! else signal remains 0
                end if
            end do
        case(3) ! Gaussian pulse
             if (sigma_width <= 0.0_dp) stop "Error in generate_signal: Gaussian sigma must be positive."
            do i = 0, n_samples-1
                t = real(i, dp) * dt_step
                value = amp * exp(-((t - center_time)**2) / (2.0_dp * sigma_width**2))
                signal(i) = cmplx(value, 0.0_dp, kind=dp)
            end do
        case default
             ! This case should not be reached due to input validation
             write(*,*) "Warning: Unknown signal type encountered in generate_signal."
             ! Optionally stop or handle error, but signal remains zero
        end select

        ! Output the generated signal (before windowing)
        filename = trim(output_prefix)//"_signal_raw.csv"
        open(newunit=file_unit, file=filename, status="replace", action="write", iostat=ios)
        if (ios /= 0) stop "Error opening raw signal file."
        write(file_unit, '(A)') "time,real,imaginary,signal_magnitude"
        do i = 0, n_samples-1
            t = real(i, dp) * dt_step
            write(file_unit, '(F15.8,",",ES18.8,",",ES18.8,",",ES18.8)') &
                t, real(signal(i)), aimag(signal(i)), abs(signal(i))
        end do
        close(file_unit)
    end subroutine generate_signal

    ! Applies the window function using element-wise multiplication
    subroutine apply_window(signal_data, window_vals, n_samples)
        integer, intent(in) :: n_samples
        complex(dp), dimension(0:n_samples-1), intent(inout) :: signal_data
        real(dp), dimension(0:n_samples-1), intent(in) :: window_vals

        ! Apply window using efficient array syntax
        signal_data = signal_data * window_vals
    end subroutine apply_window

    ! Calculates the DFT using the direct O(N^2) formula
    subroutine calculate_dft(dft_result, signal, n_samples)
        integer, intent(in) :: n_samples
        complex(dp), dimension(0:n_samples-1), intent(out) :: dft_result
        complex(dp), dimension(0:n_samples-1), intent(in) :: signal
        complex(dp) :: W_base, W_power_k_n  ! Renamed for clarity
        integer :: n, k  ! Renamed indices (k for time, n for frequency)

        if (n_samples <= 0) then
           dft_result = cmplx(0.0_dp, 0.0_dp, kind=dp)
           return
        end if

        ! DFT Definition: H_n = sum_{k=0}^{N-1} h_k * exp(2*pi*i*k*n / N)
        ! W_base = exp(2*pi*i / N)
        W_base = exp(cmplx(0.0_dp, two_pi / real(n_samples, dp), kind=dp))
        dft_result = cmplx(0.0_dp, 0.0_dp, kind=dp) ! Initialize result array

        do n = 0, n_samples-1 ! Loop over frequency components H_n
            do k = 0, n_samples-1 ! Sum over time samples h_k
                 W_power_k_n = W_base**(k*n) ! W^(kn) term
                 dft_result(n) = dft_result(n) + signal(k) * W_power_k_n
            end do
        end do

    end subroutine calculate_dft

    ! Outputs the DFT results and magnitude spectrum to CSV files
    subroutine output_results(dft_result, n_samples, dt_step, output_prefix)
        integer, intent(in) :: n_samples
        complex(dp), dimension(0:n_samples-1), intent(in) :: dft_result
        real(dp), intent(in) :: dt_step
        character(len=*), intent(in) :: output_prefix
        integer :: n, dft_unit, mag_unit, ios ! Use 'n' for frequency index
        real(dp) :: freq, magnitude, fs
        character(len=120) :: filename_dft, filename_mag

        if (n_samples <= 0) return
        if (dt_step <= 0.0_dp) stop "Error in output_results: dt_step must be positive."

        fs = 1.0_dp / dt_step ! Sampling frequency

        ! Output Full DFT Results
        filename_dft = trim(output_prefix)//"_dft.csv"
        open(newunit=dft_unit, file=filename_dft, status="replace", action="write", iostat=ios)
        if (ios /= 0) stop "Error opening DFT file: "//trim(filename_dft)
        write(dft_unit, '(A)') "index,frequency,magnitude,real,imaginary"
        do n = 0, n_samples-1
            ! Calculate frequency corresponding to index n
            if (n <= n_samples/2) then
                freq = real(n, dp) * fs / real(n_samples, dp)
            else
                ! Negative frequencies for n > N/2
                freq = real(n - n_samples, dp) * fs / real(n_samples, dp)
            end if
            magnitude = abs(dft_result(n))
            write(dft_unit, '(I0,",",F15.8,",",ES18.8E3,",",ES18.8E3,",",ES18.8E3)') &
                n, freq, magnitude, real(dft_result(n)), aimag(dft_result(n))
        end do
        close(dft_unit)

        ! Output Magnitude Spectrum Only 
        filename_mag = trim(output_prefix)//"_magnitude.csv"
        open(newunit=mag_unit, file=filename_mag, status="replace", action="write", iostat=ios)
        if (ios /= 0) stop "Error opening magnitude file: "//trim(filename_mag)
        write(mag_unit, '(A)') "frequency,magnitude"
        do n = 0, n_samples-1
             ! Calculate frequency corresponding to index n 
            if (n <= n_samples/2) then
                freq = real(n, dp) * fs / real(n_samples, dp)
            else
                freq = real(n - n_samples, dp) * fs / real(n_samples, dp)
            end if
            magnitude = abs(dft_result(n))
            ! Write frequency and magnitude
            write(mag_unit, '(F15.8,",",ES18.8E3)') freq, magnitude
        end do
        close(mag_unit)
    end subroutine output_results

end program dft