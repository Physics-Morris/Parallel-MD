! A simple generic linked list test program
program part_listest
  use particles_list
  implicit none

  type(part_list), pointer :: local_part => null()
  type(part_data), target :: part_a
  type(part_data), target :: part_b
  type(part_data), target :: part_c
  type(part_ptr) :: ptr

  type(particle), dimension(:), allocatable :: local

  integer  ::  i

  ! Initialize two data objects
  part_a%pos_x = 1.d0
  part_a%pos_y = 2.d0

  part_b%pos_x = 3.0d0
  part_b%pos_y = 4.0d0

  part_c%pos_x = 5.0d0
  part_c%pos_y = 6.0d0


  !> Initialize the list with part_a
  allocate(ptr%p)
  ptr%p%pos_x = 1.d0
  ptr%p%pos_y = 2.d0
  ! ptr%p => part_a
  call part_list_init(self=local_part, DATA=transfer(ptr, mold=part_list_data))
  print *, 'Initializing list with data:', ptr%p%pos_x
  print *, 'Initializing list with data:', ptr%p%pos_y

  !> Insert part_b into the list
  allocate(ptr%p)
  ptr%p%pos_x = 3.d0
  ptr%p%pos_y = 4.d0
  ! ptr%p => part_b
  call part_list_insert(self=local_part, DATA=transfer(ptr, mold=part_list_data))
  print *, 'Inserting node with data:', ptr%p%pos_x
  print *, 'Inserting node with data:', ptr%p%pos_y

  !> Insert part_c into the list
  allocate(ptr%p)
  ptr%p%pos_x = 5.d0
  ptr%p%pos_y = 6.d0
  ! ptr%p => part_c
  call part_list_insert(self=part_list_next(local_part), DATA=transfer(ptr, mold=part_list_data))
  print *, 'Inserting node with data:', ptr%p%pos_x
  print *, 'Inserting node with data:', ptr%p%pos_y

  !> Insert part_c into the list
  allocate(ptr%p)
  ptr%p%pos_x = 7.d0
  ptr%p%pos_y = 8.d0
  ! ptr%p => part_c
  call part_list_insert(self=part_list_next(part_list_next(local_part)), DATA=transfer(ptr, mold=part_list_data))
  print *, 'Inserting node with data:', ptr%p%pos_x
  print *, 'Inserting node with data:', ptr%p%pos_y


  !> Get the head node
  ptr = transfer(part_list_get_current(local_part), ptr)
  print *, 'Head node data:', ptr%p%pos_x
  print *, 'Head node data:', ptr%p%pos_y

  !> Get the next node
  ptr = transfer(part_list_get_current(part_list_next(local_part)), ptr)
  print *, 'Second node data:', ptr%p%pos_x
  print *, 'Second node data:', ptr%p%pos_y

  !> Get the next node
  ptr = transfer(part_list_get_current(part_list_next(part_list_next(local_part))), ptr)
  print *, 'Third node data:', ptr%p%pos_x
  print *, 'Third node data:', ptr%p%pos_y

  !> Get the next node
  ptr = transfer(part_list_get_current(part_list_next(part_list_next(part_list_next(local_part)))), ptr)
  print *, 'Forth node data:', ptr%p%pos_x
  print *, 'Forth node data:', ptr%p%pos_y

  !> Get the specific node
  ! i = 1
  ! ptr = transfer(part_list_get_index(local_part, index=i), ptr)
  ! print *, i, 'th node data:', ptr%p%pos_x
  ! print *, i, 'th node data:', ptr%p%pos_y
  do i = 1, 3
    ptr = transfer(part_list_get_index(local_part, index=i), ptr)
    print *, i, 'th node data:', ptr%p%pos_x
    print *, i, 'th node data:', ptr%p%pos_y
  end do
  i = 2
  ptr = transfer(part_list_get_index(local_part, index=i), ptr)
  print *, i, 'th node data:', ptr%p%pos_x
  print *, i, 'th node data:', ptr%p%pos_y

  !> Free the list
  call part_list_free(local_part)

  !> try another method
  allocate(local(5))

end program part_listest