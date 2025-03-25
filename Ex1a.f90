program dft_main
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 300)
   
   ! Parameters for the DFT and sine wave
   integer, parameter :: N = 1000
   real(dp), parameter :: period = 5.0_dp, amp = 1.0_dp
   
   real(dp) :: dt, t, two_pi
   integer :: i, j, ios
   
   ! Arrays for the time domain signal and its DFT result
   complex(dp), dimension(0:N-1) :: signal, dft_result
   ! DFT matrix: each element is the Nth root of unity raised to (j*k)
   complex(dp), dimension(0:N-1,0:N-1) :: dftMatrix
   complex(dp) :: W

   ! Open CSV files for output with error handling
   open(unit=10, file="data.csv", status="replace", action="write", iostat=ios)
   if (ios /= 0) then
       print *, "Error: Could not open data.csv for writing."
       error stop
   end if

   open(unit=20, file="dft.csv", status="replace", action="write", iostat=ios)
   if (ios /= 0) then
       print *, "Error: Could not open dft.csv for writing."
       error stop
   end if

   ! Write header lines for CSV files
   write(10, '(A)') "time,real,imaginary"
   write(20, '(A)') "index,real,imaginary"

   ! Calculate 2*pi 
   two_pi = 2.0_dp * (4.0_dp * atan(1.0_dp))

   ! Calculate time step so that data spans 5 periods
   dt = 5.0_dp * period / real(N, dp)

   ! Generate the sine wave data
   do i = 0, N-1
       t = real(i, dp) * dt
       signal(i) = cmplx(amp * sin(two_pi * t / period), 0.0_dp, kind=dp)
       write(10, '(F12.5, ",", ES12.5, ",", ES12.5)') t, real(signal(i)), aimag(signal(i))
   end do

   ! Define the primitive Nth root of unity for the forward DFT
   W = exp(cmplx(0.0_dp, two_pi/real(N, dp), kind=dp))

   ! Build the DFT matrix
   do j = 0, N-1
       do i = 0, N-1
           dftMatrix(j, i) = W ** (j * i)
       end do
   end do

   ! Compute the DFT 
   dft_result = matmul(dftMatrix, signal)

   ! Write the DFT output (index, real, and imaginary parts) to dft.csv
   do j = 0, N-1
       write(20, '(I6, ",", ES12.5, ",", ES12.5)') j, real(dft_result(j)), aimag(dft_result(j))
   end do

   close(10)
   close(20)

   print *, "Output written to data.csv and dft.csv"
end program dft_main


 
 
  
