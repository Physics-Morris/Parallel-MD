! A derived type for storing data.
module particles_data
  implicit none

  private
  public :: part_data
  public :: part_ptr

  ! Data is stored in part_data
  type :: part_data
     real :: x
  end type part_data

  ! A trick to allow us to store pointers in the list
  type :: part_ptr
     type(part_data), pointer :: p
  end type part_ptr
end module particles_data