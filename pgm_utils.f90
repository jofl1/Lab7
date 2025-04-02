module pgm_utils
    implicit none
    ! Define double precision kind for real numbers
    integer, parameter :: dp = selected_real_kind(15, 300)
    ! Make subroutines publicly accessible 
    public :: read_pgm, write_pgm
    
  contains
  
    ! Subroutine to read an ASCII PGM file 
    subroutine read_pgm(filename, image, width, height, maxval)
      implicit none
      character(len=*), intent(in) :: filename
      integer, allocatable, intent(out) :: image(:,:)
      integer, intent(out) :: width, height, maxval
      integer :: unit, i, j, k, ios
      character(len=2) :: magic
  
      unit = 10
      open(unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print*, "Error opening file ", filename
         stop
      end if
  
      ! Read magic number
      read(unit, *) magic
      if (magic /= 'P2') then
         print*, "Error - file should be PGM P2 format"
         close(unit)
         stop
      end if
  
      ! Read dimensions
      read(unit, *) width, height
      if (width < 1 .or. height < 1) then
         print*, "Error - invalid grid size in PGM file"
         close(unit)
         stop
      end if
  
      ! Read maximum gray value
      read(unit, *) maxval
      if (maxval <= 0 .or. maxval > 255) then
         print*, "Error - invalid maxval in PGM file"
         close(unit)
         stop
      end if
  
      ! Allocate memory for the image
      allocate(image(height, width), stat=ios)
      if (ios /= 0) then
         print*, "Error in allocating memory for image"
         close(unit)
         stop
      end if
  
      ! Read pixel values in blocks of 17 
      do j = 1, height
        do i = 1, width, 17
           read(unit, *, iostat=ios) (image(j, i+k-1), k=1, min(17, width-i+1))
           if (ios /= 0) then
              print*, "Error reading image data at line", j, "position", i
              close(unit)
              stop
           end if
        end do
     end do
     
  
      close(unit)
    end subroutine read_pgm
  
    ! Subroutine to write a PGM file in ASCII (P2 format)
    subroutine write_pgm(filename, image, width, height, maxval)
      implicit none
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: image(:,:)
      integer, intent(in) :: width, height, maxval
      integer :: unit, i, j, k, ios
  
      unit = 20
      open(unit, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         print*, "Error opening output file ", filename
         stop
      end if
  
      ! Write header
      write(unit, '(A2)') 'P2'
      write(unit, '(I3,1X,I3)') width, height
      write(unit, '(I5)') maxval
  
      ! Write pixel values in blocks of 17 per line
      do j = 1, height
        do i = 1, width, 17
           write(unit, '(17(I3.3,1X))') (nint(image(j, i+k-1)), k=1, min(17, width-i+1))
        end do
     end do
     
      close(unit)
    end subroutine write_pgm
  
  end module pgm_utils