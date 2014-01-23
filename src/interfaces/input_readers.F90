!------------------------------------------------------------------------------
!> @brief Integral interface for the FMO method.
!!
!! @author Casper Steinmann
module mb_input_readers

    use mb_utilities, only : ntotal_calculations
    use mb_precision
    use mb_variables

    implicit none

    private

#if defined (PRG_DALTON)
    public :: dalton_input
#endif

    contains

#if defined (PRG_DALTON)
!------------------------------------------------------------------------------
!> @brief reads input section from DALTON
!!
!! @author Casper Steinmann
!!
!! @details the input is read TWICE to allow for memory allocation to occur
!!          before actual readin of parameters
!!
!! @param[in] word
!! @param[in] luinp
!! @param[in] lupri
subroutine dalton_input(word, luinp, lupri)

    use mb_io, only : change_case
    use mb_utilities, only : fill_fragment_binomials

    character(len=7), intent(inout) :: word
    integer, intent(in) :: luinp
    integer, intent(in) :: lupri

    character(len=7) :: option
    integer, dimension(:), allocatable :: tempread

    integer :: iat, idx, ifrag
    integer :: ncalcs

    if (mb_initialized) return

    luout = lupri

    ! we start out by reading only what is needed for memory allocation,
    ! that is: read the following keywords
    ! NFRAG
    ! NBODY
    ! INDAT (only first record for total number of atoms)

    do
        read(luinp,'(a7)') option
        call change_case(option)

        ! check to see if we really are running a many-body calculation
        if (trim(option(2:)) == 'MB') then
            mbrun = .true.

        ! read number of fragments. No default.
        else if (trim(option(2:)) == 'NFRAG') then
            read(luinp,*) nfrag
            cycle
        ! read number of atoms from INDAT.
        else if (trim(option(2:)) == 'INDAT') then
            read(luinp,*) natmb
            cycle
        ! read n-body level. No default.
        else if (trim(option(2:)) == 'NBODY') then
            read(luinp,*) nbody
        else if (option(1:1) == '*') then
            word = option
            exit
        else if (option(1:1) == '!' .or. option(1:1) == '#') then
            cycle
        endif
    end do

    ! some sanity checks on the above variables
    if (nfrag <= 0) then
        write(luout,9000)
        stop
    endif

    ! compatibility checks
    if (nbody <= 0) then
        write(luout,9001)
        stop
    endif

    if (nfrag < nbody) then
        write(luout,9002) nfrag, nbody
        stop
    endif

    ! allocate memory blocks here based on the above variables.
    allocate(fragment_charges(nfrag))
    fragment_charges = 0

    allocate(fragment_atom_counts(nfrag))
    fragment_atom_counts = 0

    allocate(fragment_aosizes(nfrag))
    fragment_aosizes = 0

    allocate(atom_in_fragment_ids(natmb))
    atom_in_fragment_ids = 0

    allocate(nuclear_coordinates(3,natmb))
    nuclear_coordinates = 0.0d0

    allocate(nuclear_charges(natmb))
    nuclear_charges = 0.0d0

    allocate(fragment_atomic_charges(nbody,natmb))
    fragment_atomic_charges = 0.0d0

    call ntotal_calculations(ncalcs)
    allocate(fragment_energies(ncalcs))
    fragment_energies = 0.0d0

    allocate(fragment_embed_energies(nfrag))
    fragment_embed_energies = 0.0d0

    allocate(fragment_dipoles(3,ncalcs))
    fragment_dipoles = 0.0d0

    allocate(fragment_binomials(nbody,nbody))
    call fill_fragment_binomials(nbody,fragment_binomials)

    ! nullify pointers, just to be sure
    nullify(p_mb_nmer)
    nullify(p_mb_atoms)
    nullify(p_mb_fock)
    nullify(p_nmer_fock)
    nullify(p_kmer_dens)
    nullify(p_densities)

    ! start over by rewinding and read again until we hit the .MB keyword, then parse
    mbrun = .false.
    rewind(luinp)


    ! now parse the .dal file properly
    do
        read(luinp,'(a7)') option
        call change_case(option)

        ! check to see if we really are running a many-body calculation
        if (trim(option(2:)) == 'MB') then
            mbrun = .true.

        else if (.not. mbrun) then
            cycle

        ! do dummy reads of the above memory-allocating variables
        else if (trim(option(2:)) == 'NFRAG') then
            read(luinp,*) option

        ! read n-body level. No default.
        else if (trim(option(2:)) == 'NBODY') then
            read(luinp,*) option

        ! READ REAL OPTIONS

        ! read in fragment charges
        else if (trim(option(2:)) == 'CHARGE') then
            ! for now, just read a single line, but
            ! this has to be extended later on to
            ! read up to 10 items per line until
            ! nfrag items has been read
            read(luinp,*) fragment_charges

        ! read atom indices in fragments. dummy read the number of atoms
        else if (trim(option(2:)) == 'INDAT' ) then
            read(luinp,*) option
            ! read all nat fragment indices. each line is terminated with a
            ! zero (0) so read a total nat+nfrag entries
            ! find the 0 which terminates the fragment and assign to fragments
            allocate(tempread(natmb+nfrag))
            read(luinp,*) tempread
            ifrag = 1
            do idx = 1, natmb+nfrag
                if (tempread(idx) == 0) then
                    ifrag = ifrag +1
                    if (ifrag > nfrag) exit
                    cycle
                endif
                iat = tempread(idx)
                atom_in_fragment_ids(iat) = ifrag
                fragment_atom_counts(ifrag) = fragment_atom_counts(ifrag) +1
            enddo
            deallocate(tempread)

        ! disable the electrostatic potential
        else if (trim(option(2:)) == 'NOESP') then
            doesp = .false.

        ! whether or not to use atomic charges for embedding
        else if (trim(option(2:)) == 'DOATMQ') then
            use_atomic_charges = .true.

        ! MB 1-BODY CONVERGENCE THRESHOLDS
        ! maximum number of SCC iterations
        else if (trim(option(2:)) == 'MXSCCI') then
            read(luinp,*) max_scc_iterations

        else if (trim(option(2:)) == 'SCCDTO') then
            read(luinp,*) scc_d_convergence

        ! force a basis set to be used throughout
        else if (trim(option(2:)) == 'BASIS') then
            read(luinp,*) mb_basis_name

        ! verbose output
        else if (trim(option(2:)) == 'VERBOS') then
            mb_verbose = .true.

        ! debug output
        else if (trim(option(2:)) == 'DEBUG') then
            mb_debug = .true.
            mb_verbose = .true.
        else if (option(1:1) == '*') then
            word = option
            exit
        else if (option(1:1) == '!' .or. option(1:1) == '#') then
            cycle
        else
            write(luout,*) 'Unknown option:', option, ' in *MB.'
        end if
    end do

    if (.not.doesp) max_scc_iterations = 1

    ! if there is no basis set, then there is no calculation
    ! this is somewhat dangerous and I would rather read it
    ! from the .mol file. But that is a blocker right now
    if (len(mb_basis_name) == 0) then
        write(luout,9004)
        stop
    endif

9000 format(/1x,'ERROR: NFRAG must be larger than 0.')
9001 format(/1x,'ERROR: NBODY must be larger than 0.')
9002 format(/1x,'ERROR: NFRAG should be larger than or equal to the NBODY level.')
9004 format(/1x,'ERROR: You have to specify a basis set through the .BASIS keyword. in the *MB group.')

end subroutine dalton_input

#endif

end module mb_input_readers
