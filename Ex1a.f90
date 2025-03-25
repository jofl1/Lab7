program dft
   implicit none
   
   ! Define precision for real numbers 
   integer, parameter :: dp = selected_real_kind(15, 300)
   ! Constants
   real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp), two_pi = 2.0_dp * pi
   ! Signal parameters
   integer, parameter :: N = 1000 ! Number of samples
   real(dp), parameter :: period = 5.0_dp, amp = 1.0_dp  ! Wave properties
   
   integer :: ios  ! For I/O status
   character(len=100) :: output_prefix = "sine_wave"
   
   ! Arrays for time domain signal and frequency domain result
   complex(dp), dimension(0:N-1) :: signal, dft_result
   
   call generate_signal() ! Create sine wave
   call calculate_dft()  ! Perform DFT
   call output_results()  ! Write results to file

contains

   subroutine generate_signal()
       real(dp) :: dt = 5.0_dp * period / real(N, dp)  ! Time step
       integer :: i
       real(dp) :: t
       character(len=100) :: filename
       
       ! Generate sine wave samples
       do i = 0, N-1
           t = real(i, dp) * dt
           signal(i) = cmplx(amp * sin(two_pi * t / period), 0.0_dp, kind=dp)
       end do
       
       ! Write time domain signal to CSV file
       filename = trim(output_prefix) // "_data.csv"
       open(unit=10, file=filename, status="replace", action="write", iostat=ios)
       if (ios /= 0) stop "Error: Could not open file for writing."
       
       write(10, '(A)') "time,real,imaginary"  ! CSV header
       do i = 0, N-1
           t = real(i, dp) * dt
           write(10, '(F12.5,",",ES12.5,",",ES12.5)') t, real(signal(i)), aimag(signal(i))
       end do
       
       close(10)
   end subroutine generate_signal
   
   subroutine calculate_dft()
       complex(dp) :: W                        ! Nth root of unity
       integer :: i, j
       
       ! Define the principal Nth root of unity for DFT
       W = exp(cmplx(0.0_dp, two_pi/real(N, dp), kind=dp))
       
       ! Direct implementation of DFT formula
       dft_result = 0.0_dp
       do j = 0, N-1  ! For each frequency bin
           do i = 0, N-1  ! For each time sample
               dft_result(j) = dft_result(j) + signal(i) * W ** (i * j)
           end do
       end do
   end subroutine calculate_dft
   
   subroutine output_results()
       integer :: j
       character(len=100) :: filename
       
       ! Write DFT results to CSV file
       filename = trim(output_prefix) // "_dft.csv"
       open(unit=20, file=filename, status="replace", action="write", iostat=ios)
       if (ios /= 0) stop "Error: Could not open file for writing."
       
       write(20, '(A)') "index,real,imaginary"  ! CSV header
       do j = 0, N-1  ! Write each frequency component
           write(20, '(I6,",",ES12.5,",",ES12.5)') j, real(dft_result(j)), aimag(dft_result(j))
       end do
       
       close(20)
   end subroutine output_results

end program dft