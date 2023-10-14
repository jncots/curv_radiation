module util_test

  implicit none
  private
  save

public :: cwd_parent_dir

 contains

!===========================================================================
! FUNCTION: cwd_parent_dir
!
! DESCRIPTION:
!   This function returns the parent directory path of the current working
!   directory up to a specified level.
!
! INPUT:
!   level - An integer specifying the number of levels to move up from
!           the current working directory.
!
! OUTPUT:
!   cwd_parent_dir - A character string representing the parent directory
!                    path up to the specified level.
!
! USAGE:
!   use util_test, only : cwd_parent_dir
!   write(*,*) trim(cwd_parent_dir(i))
!
!===========================================================================

 function cwd_parent_dir(level)
    implicit none

    ! Input Parameters
    integer, intent(in) :: level

    ! Output
    character(300) :: cwd_parent_dir
    character(300) :: cwd
    integer :: i, idx

    ! Get the current working directory
    call getcwd(cwd)
    cwd = trim(cwd)

    ! Find the index of the last character in the path
    idx = len(cwd)

    ! Move up the directory hierarchy to the specified level
    do i = 1, level
        idx = index(cwd(:idx), '/', .true.) - 1
        if (idx < 1) then
            cwd = "/"  ! Reached the root directory
            exit
        end if
        cwd = cwd(:idx)
    end do

    ! Set the output variable to the parent directory path
    cwd_parent_dir = cwd
end function cwd_parent_dir


end module util_test




