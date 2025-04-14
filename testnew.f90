program edge_detection
    use pgm_utils 
    implicit none
    
    ! External declarations for FFTW library
    integer, parameter :: i64 = selected_int_kind(18) 
    
    ! Interface block for external FFTW procedures
    interface
        subroutine dfftw_plan_dft_r2c_2d(plan, n0, n1, in, out, flags)
            import i64, dp 
            integer(kind=i64) :: plan 
            integer :: n0, n1, flags
            real(dp) :: in(*)      
            complex(dp) :: out(*)   
        end subroutine
        
        subroutine dfftw_plan_dft_c2r_2d(plan, n0, n1, in, out, flags)
            import i64, dp
            integer(kind=i64) :: plan 
            integer :: n0, n1, flags
            complex(dp) :: in(*)    
            real(dp) :: out(*)      
        end subroutine
        
        subroutine dfftw_execute(plan)
            import i64
            integer(kind=i64) :: plan 
        end subroutine
        
        subroutine dfftw_destroy_plan(plan)
            import i64
            integer(kind=i64) :: plan 
        end subroutine
    end interface
    
    ! Parameters
    character(len=*), parameter :: input_file = "clown.pgm"     
    character(len=*), parameter :: output_file = "clown_edges.pgm"
    character(len=*), parameter :: blurred_file = "blurred.pgm" 
    character(len=*), parameter :: deblurred_file = "deblurred.pgm"
    
    ! Variables
    integer, allocatable :: input_image(:,:)   
    real(dp), allocatable :: result_image(:,:) 
    integer :: width, height, image_maxval     
    integer :: method            
    integer :: mask_radius = 20  
    character(len=100) :: user_input 
    logical :: valid_input       
    integer :: ios               
    integer :: blur_width = 10   
    
    ! Main Program Execution 
    call read_pgm(input_file, input_image, width, height, image_maxval)
    
    ! User Method Selection 
    valid_input = .false.
    do while (.not. valid_input)
        print *, "Please select the image processing method:"
        print *, "  1: Sobel edge detection (spatial domain)"
        print *, "  2: High-pass filter edge detection (frequency domain using FFT)"
        print *, "  3: Deblur image (deconvolution using FFT)"
        write (*, '(A)', advance='no') "Enter choice (1, 2, or 3): " 
        read (*, '(A)') user_input 
        read(user_input, *, iostat=ios) method
        if (ios == 0) then 
            if (method >= 1 .and. method <= 3) then
                valid_input = .true.; print *, "Method selected: ", method 
            else
                print *, "-> Invalid choice. Please enter 1, 2, or 3."
            end if
        else
            print *, "-> Invalid input. Please enter a number (1, 2, or 3)."
        end if
    end do
    
    ! Get Method specific parameters 
    if (method == 2) then
        valid_input = .false.
        do while (.not. valid_input)
            print *, ""
            write (*, '(A, I0, A)', advance='no') "Enter mask radius for high-pass filter: " 
            read (*, '(A)') user_input
            if (len_trim(user_input) == 0) then 
                mask_radius = 20; valid_input = .true.; print *, "Using default mask radius: ", mask_radius 
            else
                read(user_input, *, iostat=ios) mask_radius
                if (ios == 0 .and. mask_radius > 0) then 
                    valid_input = .true.; print *, "Mask radius set to: ", mask_radius 
                else if (ios /= 0) then
                     print *, "-> Invalid input. Please enter a positive number or press Enter for default."
                else 
                     print *, "-> Invalid radius. Please enter a positive number."
                end if
            end if
        end do
    end if
    
    if (method == 3) then
        valid_input = .false.
        do while (.not. valid_input)
            print *, ""
            write (*, '(A, I0, A)', advance='no') "Enter width of the horizontal motion blur: "
            read (*, '(A)') user_input
            if (len_trim(user_input) == 0) then 
                blur_width = 10; valid_input = .true.; print *, "Using default blur width: ", blur_width 
            else
                read(user_input, *, iostat=ios) blur_width
                if (ios == 0 .and. blur_width > 0 .and. blur_width <= width) then 
                    valid_input = .true.; print *, "Blur width set to: ", blur_width 
                else if (ios /= 0) then
                     print *, "-> Invalid input. Please enter a positive number or press Enter for default."
                else if (blur_width <= 0) then
                     print *, "-> Invalid width. Please enter a positive number."
                else 
                     print *, "-> Invalid width. Blur width cannot exceed image width (", width, ")."
                end if
            end if
        end do
        
        if (allocated(input_image)) deallocate(input_image)
        print *, ""; print *, "Reading blurred image: ", blurred_file
        call read_pgm(blurred_file, input_image, width, height, image_maxval) 
        print *, "Blurred image dimensions:", width, "x", height, " Max value:", image_maxval
    end if
    
    ! Allocate result image, initialising to zero 
    allocate(result_image(height, width), stat=ios, source=0.0_dp) 
    if (ios /= 0) then
        print *, "Error: Could not allocate memory for result image."; stop 1 
    end if
    
    ! Call selected processing method
    select case(method)
        case(1) 
            print *, "Using Sobel operator method"
            call sobel_method(input_image, result_image, width, height, image_maxval)
            print *, "Writing edge-detected image to:", output_file
            call write_pgm(output_file, result_image, width, height, image_maxval)
        case(2) 
            print *, "Using FFT-based high-pass filter method with mask radius:", mask_radius
            call fft_method(input_image, result_image, width, height, image_maxval, mask_radius)
            print *, "Writing edge-detected image to:", output_file
            call write_pgm(output_file, result_image, width, height, image_maxval)
        case(3) 
            print *, "Using deblurring method with blur width:", blur_width
            call deblur_method(input_image, result_image, width, height, image_maxval, blur_width)
            print *, "Writing deblurred image to:", deblurred_file
            call write_pgm(deblurred_file, result_image, width, height, image_maxval)
        case default
            print *, "Error: Invalid method selected:", method; stop 2 
    end select
    
    ! Cleanup
    if (allocated(input_image)) deallocate(input_image)
    if (allocated(result_image)) deallocate(result_image)

contains 

    ! Sobel edge detection subroutine
    subroutine sobel_method(input_image, gradient, width, height, image_maxval)
        integer, intent(in) :: input_image(:,:), width, height, image_maxval 
        real(dp), intent(out) :: gradient(:,:) 
        
        real(dp), allocatable :: grad_x(:,:), grad_y(:,:) 
        integer :: i, j
        integer, parameter :: Gx(3,3) = reshape((/ -1, 0, 1, -2, 0, 2, -1, 0, 1 /), shape(Gx))
        integer, parameter :: Gy(3,3) = reshape((/  1, 2, 1,  0, 0, 0, -1,-2,-1 /), shape(Gy))
        
        ! Allocate and initialise gradient components to zero 
        allocate(grad_x(height, width), grad_y(height, width), source=0.0_dp)
        
        ! Apply Sobel operators
        do j = 2, height-1 
            do i = 2, width-1
                grad_x(j,i) = sum( real(Gx(:,:), dp) * real(input_image(j-1:j+1, i-1:i+1), dp) )
                grad_y(j,i) = sum( real(Gy(:,:), dp) * real(input_image(j-1:j+1, i-1:i+1), dp) )
            end do
        end do
        
        gradient = sqrt(grad_x**2 + grad_y**2)
        
        ! Clip result
        where (gradient < 0.0_dp) gradient = 0.0_dp
        where (gradient > real(image_maxval, dp)) gradient = real(image_maxval, dp) 
        
        deallocate(grad_x, grad_y)
    end subroutine sobel_method

    ! FFT High pass filter edge detection subroutine
    subroutine fft_method(input_image, gradient, width, height, image_maxval, mask_radius)
        implicit none
    
        integer, intent(in)     :: input_image(:,:), width, height, image_maxval, mask_radius 
        real(dp), intent(out) :: gradient(:,:) 
    
        integer, parameter :: FFTW_ESTIMATE = 64 
        real(dp),    allocatable :: image_real(:,:)    
        complex(dp), allocatable :: image_complex(:,:) 
        integer(i64) :: plan_forward, plan_backward
        integer :: i, j, nx, ny 
        real(dp) :: dist_sq, max_abs_value 
        integer(i64) :: zeroed_freqs, total_freqs 
        real(dp) :: grad_val 
    
        ! Allocate and initialise FFT arrays 
        allocate(image_real(height, width), source=0.0_dp)
        allocate(image_complex(height, (width/2)+1), source=cmplx(0.0_dp, 0.0_dp, kind=dp))
    
        image_real = real(input_image, kind=dp)
    
        ! Create FFT Plans
        call dfftw_plan_dft_r2c_2d(plan_forward,  height, width,  image_real, image_complex, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_2d(plan_backward, height, width,  image_complex, image_real, FFTW_ESTIMATE)
        if (plan_forward == 0_i64 .or. plan_backward == 0_i64) then
            print *, "Error: Failed to create FFTW plans in fft_method."; stop 3 
        end if

        ! Execute Forward FFT
        call dfftw_execute(plan_forward)
        
        ! Apply high pass filter in frequency domain
        ! Zero out low-frequency components within the specified radius around the DC component
        zeroed_freqs = 0_i64
        total_freqs  = int(height, kind=i64) * int((width/2)+1, kind=i64)
        do j = 1, height 
            ! Calculate the effective y frequency index (ny) considering wrap around
            ny = j - 1; if (ny > height/2) ny = ny - height 
            do i = 1, (width/2)+1 
                nx = i - 1; dist_sq = real(nx*nx + ny*ny, dp) 
                ! High-pass filter - Zero out frequencies inside the radius
                if (dist_sq <= real(mask_radius, dp)**2) then 
                    image_complex(j,i) = cmplx(0.0_dp, 0.0_dp, kind=dp) 
                    zeroed_freqs = zeroed_freqs + 1_i64
                end if
            end do
        end do
        print *, "FFT filter: Zeroed ", zeroed_freqs, " out of ", total_freqs, " frequency components."

        ! Execute Backward FFT
        call dfftw_execute(plan_backward)

        ! Normalise and Finalise Output
        image_real = image_real / real(width*height, dp) 
        max_abs_value = 0.0_dp 
        do j = 1, height
            do i = 1, width
                grad_val = abs(image_real(j,i))
                ! Check for NaN or Inf manually (value /= value is true for NaN)
                if ((grad_val /= grad_val) .or. (abs(grad_val) >= huge(grad_val))) then
                    gradient(j,i) = 0.0_dp
                else
                    ! Clip valid numbers
                    if (grad_val > max_abs_value) max_abs_value = grad_val
                    if (grad_val > real(image_maxval, dp)) then
                       gradient(j,i) = real(image_maxval, dp) 
                    else if (grad_val < 0.0_dp) then 
                       gradient(j,i) = 0.0_dp
                    else
                       gradient(j,i) = grad_val
                    end if
                end if
            end do
        end do
        print *, "Gradient values clipped to [0, ", image_maxval, "]." 

        call dfftw_destroy_plan(plan_forward); call dfftw_destroy_plan(plan_backward)
        deallocate(image_real); deallocate(image_complex)
    end subroutine fft_method

    ! FFT deblurring subroutine
    subroutine deblur_method(input_image, output_image, width, height, image_maxval, blur_width)
        implicit none
    
        integer, intent(in)     :: input_image(:,:), width, height, image_maxval, blur_width 
        real(dp), intent(out) :: output_image(:,:) 
    
        integer, parameter :: FFTW_ESTIMATE = 64
        real(dp),    allocatable :: image_real(:,:), blur_real(:,:)
        complex(dp), allocatable :: image_complex(:,:), blur_complex(:,:)
        integer(i64) :: plan_forward_image, plan_backward_image, plan_forward_blur
        integer :: i, j
        real(dp) :: max_value_out, min_value_out, magnitude_blur       
        real(dp), parameter :: eps = 1.0e-10_dp 
        complex(dp) :: H 
        real(dp) :: pixel_val

        ! Allocate and initialize FFT arrays 
        allocate(image_real(height, width), source=0.0_dp)
        allocate(image_complex(height, (width/2)+1), source=cmplx(0.0_dp, 0.0_dp, kind=dp))
        allocate(blur_real(height, width), source=0.0_dp)
        allocate(blur_complex(height, (width/2)+1), source=cmplx(0.0_dp, 0.0_dp, kind=dp))
    
        image_real = real(input_image, kind=dp)
    
        ! Create Blur Function - horizontal motion blur at top left
        if (blur_width > 0 .and. blur_width <= width) then
            blur_real(1, 1:blur_width) = 1.0_dp / real(blur_width, dp) 
        else
            print *, "Warning: Invalid blur_width (", blur_width, ") in deblur_method. Using width=1."
            blur_real(1, 1) = 1.0_dp 
        end if
        
        ! Create FFT Plans
        call dfftw_plan_dft_r2c_2d(plan_forward_image, height, width, image_real, image_complex, FFTW_ESTIMATE)
        call dfftw_plan_dft_r2c_2d(plan_forward_blur, height, width, blur_real, blur_complex, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_2d(plan_backward_image, height, width, image_complex, image_real, FFTW_ESTIMATE)
        if (plan_forward_image == 0_i64 .or. plan_forward_blur == 0_i64 .or. plan_backward_image == 0_i64) then
             print *, "Error: Failed to create FFTW plans in deblur_method."; stop 4 
        end if

        ! Execute Forward FFTs
        call dfftw_execute(plan_forward_image)
        call dfftw_execute(plan_forward_blur)

        ! Deconvolve in frequency domain deblurred = Blurred / Blur_PSF
        ! Includes check for near zero magnitude in Blur_PSF to avoid noise amplification
        do j = 1, height
            do i = 1, (width/2)+1
                H = blur_complex(j,i) 
                ! Reverted back to manual calculation for complex magnitude
                magnitude_blur = sqrt(real(H, dp)**2 + aimag(H)**2) 
                
                if (magnitude_blur > eps) then
                    image_complex(j,i) = image_complex(j,i) / H
                else
                    ! Avoid division by zero/small number
                    image_complex(j,i) = cmplx(0.0_dp, 0.0_dp, kind=dp) 
                end if
            end do
        end do
    
        ! Execute Backward FFT
        call dfftw_execute(plan_backward_image)

        ! Normalise and Finalise Output
        image_real = image_real / real(width*height, dp) 
        min_value_out = minval(image_real); max_value_out = maxval(image_real)
        print *, "After inverse FFT (raw deblurred) - Min:", min_value_out, "Max:", max_value_out
        print *, "Clipping deblurred image to [0, ", image_maxval, "]" 
        do j = 1, height
            do i = 1, width
                pixel_val = image_real(j,i) 
                ! Check for NaN or Inf manually (value /= value is true for NaN)
                if ((pixel_val /= pixel_val) .or. (abs(pixel_val) >= huge(pixel_val))) then
                    output_image(j,i) = 0.0_dp 
                else 
                    ! Clip valid numbers
                    if (pixel_val < 0.0_dp) then 
                        output_image(j,i) = 0.0_dp 
                    else if (pixel_val > real(image_maxval, dp)) then 
                        output_image(j,i) = real(image_maxval, dp) 
                    else
                        output_image(j,i) = pixel_val 
                    end if
                end if
            end do
        end do
    
        call dfftw_destroy_plan(plan_forward_image)
        call dfftw_destroy_plan(plan_forward_blur)
        call dfftw_destroy_plan(plan_backward_image)
        deallocate(image_real, image_complex, blur_real, blur_complex)
    end subroutine deblur_method
    
end program edge_detection
