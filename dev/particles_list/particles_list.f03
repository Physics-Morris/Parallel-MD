module particles_list
    implicit none

    private

    !> store pointer in a list
    public :: part_list
    !> public variable to use as a MOLD for transfer function
    public :: part_list_data
    !> initilize particle list
    public :: part_list_init
    !> free particle list
    public :: part_list_free
    !> insert part list data into the particle list
    public :: part_list_insert
    public :: part_list_put
    !> get current data at current NODE
    public :: part_list_get_current
    !> get data at ith index NODE
    public :: part_list_get_index
    !> go to next in part_list
    public :: part_list_next

    !> particle data type
    public :: part_data
    !> pointer to particle data
    public :: part_ptr

    !> A public variable to use as a MOLD for transfer()
    integer, dimension(:), allocatable :: part_list_data


    !> Linked list node data type
    type :: part_list
        private
        integer, dimension(:), pointer :: data => null()
        type(part_list), pointer :: next => null()
    end type part_list

    !> Data is stored in part_data
    type :: part_data
        integer          :: id, species_id
        double precision :: mass, charge
        double precision :: pos_x, pos_y, pos_z
        double precision :: vel_x, vel_y, vel_z
    end type part_data

    !> A trick to allow us to store pointers in the list
    type :: part_ptr
        type(part_data), pointer :: p
    end type part_ptr


    contains


        !> Initialize a head node SELF and optionally store the provided DATA.
        subroutine part_list_init(self, data)
            type(part_list), pointer :: self
            integer, dimension(:), intent(in), optional :: data

            allocate(self)
            nullify(self%next)

            if (present(data)) then
                allocate(self%data(size(data)))
                self%data = data
            else
                nullify(self%data)
            end if
        end subroutine part_list_init


        !> Free the entire list and all data, beginning at SELF
        subroutine part_list_free(self)
            type(part_list), pointer :: self
            type(part_list), pointer :: current
            type(part_list), pointer :: next

            current => self
            do while (associated(current))
                next => current%next
                if (associated(current%data)) then
                    deallocate(current%data)
                    nullify(current%data)
                end if
                deallocate(current)
                nullify(current)
                current => next
            end do
        end subroutine part_list_free


        ! Return the next node after SELF
        function part_list_next(self) result(next)
            type(part_list), pointer :: self
            type(part_list), pointer :: next
            next => self%next
        end function part_list_next


        !> Insert a list node after SELF containing DATA (optional)
        subroutine part_list_insert(self, data)
            type(part_list), pointer :: self
            integer, dimension(:), intent(in), optional :: data
            type(part_list), pointer :: next

            allocate(next)

            if (present(data)) then
                allocate(next%data(size(data)))
                next%data = data
            else
                nullify(next%data)
            end if

            next%next => self%next
            self%next => next
        end subroutine part_list_insert


        !> Store the encoded DATA in list node SELF
        subroutine part_list_put(self, data)
            type(part_list), pointer :: self
            integer, dimension(:), intent(in) :: data

            if (associated(self%data)) then
                deallocate(self%data)
                nullify(self%data)
            end if
            self%data = data
        end subroutine part_list_put


        !> Return the DATA stored in the node SELF
        function part_list_get_current(self) result(data)
            type(part_list), pointer :: self
            integer, dimension(:), pointer :: data
            data => self%data
        end function part_list_get_current


        !> get DATA stored in ith index
        function part_list_get_index(self, index) result(data)
            type(part_list), pointer :: self
            type(part_list), pointer :: indexNode
            integer, dimension(:), pointer :: data
            integer :: index, i

            nullify(indexNode)
            if (index == 1) then
                data => self%data
            else
                indexNode => self%next
                do i = 1, index-2
                    indexNode = indexNode%next
                end do
                data => indexNode%data
            end if
        end function


end module particles_list