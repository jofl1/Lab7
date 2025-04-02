program sobel_edge_detection
    use pgm_utils
    implicit none
    
    ! Parameters
    character(len=*), parameter :: input_file = "clown.pgm"
    character(len=*), parameter :: output_file = "clown_edges.pgm"
    
    ! Variables for image
    integer, allocatable :: input_image(:,:)
    real(dp), allocatable :: grad_x(:,:), grad_y(:,:), gradient(:,:)
    integer :: width, height, maxval, i, j
    
    ! Sobel operators (3x3 convolution kernels)
    integer :: Gx(3,3) = reshape((/ -1, 0, 1, -2, 0, 2, -1, 0, 1 /), shape(Gx))
    integer :: Gy(3,3) = reshape((/ 1, 2, 1, 0, 0, 0, -1, -2, -1 /), shape(Gy))
    
    ! Read the PGM image file
    print *, "Reading image file:", input_file
    call read_pgm(input_file, input_image, width, height, maxval)
    print *, "Image dimensions:", width, "x", height, "with max value:", maxval
    
    ! Allocate memory for gradient images
    allocate(grad_x(height, width), grad_y(height, width), gradient(height, width))
    grad_x = 0.0_dp
    grad_y = 0.0_dp
    gradient = 0.0_dp
    
    ! Apply Sobel operator in x-direction
    print *, "Applying Sobel filter in x-direction..."
    do j = 2, height-1
        do i = 2, width-1
            grad_x(j,i) = &
                Gx(1,1) * input_image(j-1,i-1) + Gx(1,2) * input_image(j-1,i) + Gx(1,3) * input_image(j-1,i+1) + &
                Gx(2,1) * input_image(j,i-1)   + Gx(2,2) * input_image(j,i)   + Gx(2,3) * input_image(j,i+1)   + &
                Gx(3,1) * input_image(j+1,i-1) + Gx(3,2) * input_image(j+1,i) + Gx(3,3) * input_image(j+1,i+1)
        end do
    end do
    
    ! Apply Sobel operator in y-direction
    print *, "Applying Sobel filter in y-direction..."
    do j = 2, height-1
        do i = 2, width-1
            grad_y(j,i) = &
                Gy(1,1) * input_image(j-1,i-1) + Gy(1,2) * input_image(j-1,i) + Gy(1,3) * input_image(j-1,i+1) + &
                Gy(2,1) * input_image(j,i-1)   + Gy(2,2) * input_image(j,i)   + Gy(2,3) * input_image(j,i+1)   + &
                Gy(3,1) * input_image(j+1,i-1) + Gy(3,2) * input_image(j+1,i) + Gy(3,3) * input_image(j+1,i+1)
        end do
    end do
    
    ! Calculate the gradient magnitude
    print *, "Calculating gradient magnitude..."
    do j = 1, height
        do i = 1, width
            ! Use the square root of the sum of squares for the gradient magnitude
            gradient(j,i) = sqrt(grad_x(j,i)**2 + grad_y(j,i)**2)
        end do
    end do
    
    ! Apply clipping to the gradient values (instead of normalization)
    print *, "Clipping gradient values..."
    do j = 1, height
        do i = 1, width
            ! Clip the gradient values to the range [0, maxval]
            if (gradient(j,i) > maxval) then
                gradient(j,i) = maxval
            else if (gradient(j,i) < 0) then
                gradient(j,i) = 0
            end if
        end do
    end do
    
    ! Write the gradient image to a PGM file
    print *, "Writing edge-detected image to:", output_file
    call write_pgm(output_file, gradient, width, height, maxval)
    
    ! Clean up
    deallocate(input_image, grad_x, grad_y, gradient)
    print *, "Edge detection completed successfully."

end program sobel_edge_detection