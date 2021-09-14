program test_input
    implicit none
    !> simulation domain related shared data
    integer                     :: sim_dimension
    double precision            :: x_min, x_max, y_min, y_max, z_min, z_max

    !> particle related shared data
    integer                     :: total_particles
    double precision            :: species_mass, species_charge
    character(len=60)           :: species_distribution

    integer                     :: number_snapshots

    integer                     :: ierr

    namelist / basics_block    / sim_dimension, x_min, x_max, y_min, y_max, z_min, z_max
    namelist / species_block   / total_particles
    namelist / output_block    / number_snapshots

    open(10, file='../../inp/default')
    read(10, nml=basics_block)
    read(10, nml=species_block)
    read(10, nml=output_block)
    if (ierr .eq. 0) then
        ! read(10, nml=species)
        ! read(10, nml=output)
        close(10)
    else
        write(*, *) ierr
        stop
    end if
    print*, species_distribution
end program test_input