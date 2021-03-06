!> @brief A generalized many-body calculation module
!!
!! @author Casper Steinmann
!!
!! @details This module is developed in the hope that it will
!!          be useful to people, either in terms of inspiration
!!          or for further development
!! @bug DFT-calculations do multiple read-ins of coefficients from DALTON
module mb

    use mb_precision
    use mb_variables

    implicit none

    private

    intrinsic :: allocated, present, min, minval, max, maxval, size, cpu_time

    public :: mb_init            ! initializes 
    public :: mb_input           ! reads and parses input
    public :: mb_executer        ! executes a many-body calculation.
    public :: mb_fock            ! adds a possibe contribution to the fock matrix
    public :: mb_store_nmer_property ! storage of properties
    public :: mb_output          ! handles many-body output
    public :: mb_twoints

contains

!------------------------------------------------------------------------------
!> @brief controls and executes a many-body calculation
!!
!! @author Casper Steinmann
!!
!! @param[in,out] dalwrk dalton memory pointer
!! @param[in] iwsize dalton memory size
!! @param[in] wctrl dalton memory control variable
subroutine mb_executer(dalwrk,iwsize,wctrl)

    ! dimension dalwrk properly here or you will spend many hours tracing a bug
    integer, intent(in) :: iwsize
    real(dp), dimension(0:iwsize+1), target, intent(inout) :: dalwrk
    real(dp), intent(in) :: wctrl

    integer :: kbody

    work => dalwrk

    ! pass memory pointers and checks along to internal storage
    idalton_work_size = iwsize
    dalton_work_control = wctrl


    ! loop over all nbody calculations we have to do.
    ! The imb_current_calculation variable keeps
    ! track for storing _ALL_ properties whereas
    ! i_mb_kbody is for storing the current n-body level.
    imb_current_calculation = 0
    do kbody=1,nbody
        i_mb_kbody = kbody
        if (i_mb_kbody == 1) then
            do i_scc_iteration=1,max_scc_iterations

                imb_current_calculation = 0
                max_density_difference = 0.0d0

                ! evaluate the energy and properties for the current
                ! kbody level (here, explicitly 1-body)
                call kbody_calculation(kbody)
                write(luout,9000) i_scc_iteration, max_density_difference, scc_d_convergence

                if ((max_density_difference < scc_d_convergence .and. &
                    i_scc_iteration > 1).or. &
                    (i_scc_iteration == max_scc_iterations)) then
                    write(luout,9005) i_scc_iteration
                    exit
                endif

            enddo
        else
            ! evaluate the energy and properties for the current
            ! kbody level (here kbody > 1)
            call kbody_calculation(kbody)
        endif

    enddo

    nullify(work)

9000 format(72(1H-),/3x,"MB one-body iter",i3," max(dD) =",2f12.9,/72(1H-))
9005 format(/72(1H-),///3x,"SCC converged in",i3," iterations!",///72(1H-))

end subroutine

!------------------------------------------------------------------------------
!> @brief calculate all k-body energies (and properties) of the system in the presence of all the other fragments
!!
!! @author Casper Steinmann
!!
!! @detals calculates all the individual terms for a given @f$ k @f$-body.
!!         All 1-body densities are stored in memory and written to disk when all
!!         calculations for this SCC iteration is finished.
!!
!! @param[in] kbody the current k-body level
subroutine kbody_calculation(kbody)

    use mb_utilities, only : calc_n_kbodies

    integer, intent(in) :: kbody

    integer :: icalc, ncalcs, total_nao
    integer, dimension(:,:), allocatable :: fragment_list
    real(dp), dimension(:), allocatable, target :: densities

    ! obtain the numer of calculations for a given kbody
    call calc_n_kbodies(kbody, ncalcs)
    allocate(fragment_list(ncalcs,kbody))

    ! build list of fragments on which to calculate, look in the subroutine
    ! for the definition of the format
    call build_fragment_list(fragment_list)

    write(luout,8000) kbody, ncalcs

    ! ready the storage of densities. this is only possible after the first i_scc_iteration
    if (i_scc_iteration > 1.and. kbody == 1) then
        total_nao = sum(fragment_aosizes)
        allocate(densities(total_nao*(total_nao+1)/2))
        densities = 0.0d0
        p_densities => densities
    endif

    ! perform the ncalcs calculations for the current ibody level
    do icalc=1,ncalcs
        imb_current_calculation = imb_current_calculation + 1
        call calculate_nmer(fragment_list(icalc,:))
        write(luout,*)
    enddo

    ! store the densities to disk here and free the memory
    if (i_scc_iteration > 1.and. kbody == 1) then
        if (allocated(densities)) then
            call store_nmer_densities(ncalcs)
            nullify(p_densities)
            deallocate(densities)
        endif
    endif

    deallocate(fragment_list)

    write(luout,9000) kbody

8000 format(/1x,'----',i3,'-body calculation has',i5,' calculations ----')
9000 format(/1x,'----',i3,'-body calculation complete ----',/)

end subroutine kbody_calculation

!> @brief stores all densities after they have been calculated
!! @details Stores all densities currently held in memory to disk after each SCC iteration.
!!          This is needed because the previous SCC iteration densities are located on disk.
!! @note this is only supposed to be executed for 1-body calculations
!! @info uses the pointer p_densitites
!! @param[in] ncalcs the number of calculations at the 1-body level
subroutine store_nmer_densities(ncalcs)

    use mb_utilities, only : get_nmer_aosize, &
                              set_nmer_density

    integer, intent(in) :: ncalcs

    integer :: icalc, nao_offset, nao_imer, nao
    integer, dimension(1) :: imer

    ! if we are here and the densities pointer is not associated, then
    ! something went horribly wrong
    if (.not. associated(p_densities)) then
        write(luout,*) "ERROR: memory-pointer to densities is not associated"
        stop
    endif

    nao_offset = 0
    do icalc=1,ncalcs
        imer(1) = icalc
        call get_nmer_aosize(imer, nao)
        nao_imer = nao*(nao+1)/2
        call set_nmer_density(icalc, 'old', p_densities(1+nao_offset:nao_imer+nao_offset))
        nao_offset = nao_offset + nao_imer
    enddo

end subroutine

!------------------------------------------------------------------------------
!> @brief Builds a list of fragments for a given n-body calculation
!!
!! @author Casper Steinmann
!!
!! @details storage for 4 fragments in a three-body calculation will
!!          calculate the following table of calculations for a given
!!          k-body calculation
!!          1-body
!!          1, (1)
!!          2, (2)
!!          3, (3)
!!          4, (4)
!!
!!          2-body
!!          (1,2)
!!          (1,3)
!!          (1,4)
!!          (2,3)
!!          (2,4)
!!          (3,4)
!!
!!          3-body
!!          (1,2,3)
!!          (1,2,4)
!!          (1,3,4)
!!          (2,3,3)
!!
!! @param[out] fragment_list list of fragments and their constituion monomers
!!                           for a given n-mer level
subroutine build_fragment_list(fragment_list)

    integer, dimension(:,:), intent(out) :: fragment_list

    integer :: ibody, icalc, imer
    integer :: ncalc

    fragment_list = 0
    ncalc = size(fragment_list,1)
    ibody = size(fragment_list,2)

    icalc = 1
    do imer=1,nfrag
        ! avoid errornous icalc counting for ibody > 1
        if(icalc<=ncalc) fragment_list(icalc:,1) = imer
        call get_related_nmers(fragment_list, icalc, ibody, imer)
    enddo

end subroutine build_fragment_list

!------------------------------------------------------------------------------
!> @brief fills fragment_list with fragment indices for a given n-mer level
!!
!! @author Casper Steinmann
!!
!! @param fragment_list[in,out] list of fragments and their constituion monomers
!!                              for a given n-mer level
!! @param[in,out] icalc the icalc'th calculation
!! @param[in] ilevel how far down the recursion we are.
!! @param[in] ifrg the fragment index, offset for next recursion
recursive subroutine get_related_nmers(fragment_list,icalc,ilevel,ifrg)

    integer, intent(in) :: ilevel, ifrg
    integer, intent(inout) :: icalc
    integer, dimension(:,:), intent(inout) :: fragment_list
    integer :: ii, idlevel, icount, ncalc, kbody

    ncalc = size(fragment_list,1)
    kbody = size(fragment_list,2)

    if(ilevel == 1) then
        icalc = icalc +1
        return
    else
        icount = 1
        do ii=ifrg+1,nfrag
            idlevel = kbody - ilevel + 2
            if(icalc<= ncalc) fragment_list(icalc:, idlevel) = ii
            call get_related_nmers(fragment_list, icalc, ilevel-1, ii)
        enddo
    endif

end subroutine get_related_nmers

!------------------------------------------------------------------------------
!> @brief Evaluates the energy (and other properties) of a given nmer
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer the nmer to evaluate
subroutine calculate_nmer(nmer)

    use mb_utilities, only : get_nmer_atoms

    integer, dimension(:), intent(in), target :: nmer

    integer, dimension(:), allocatable, target :: nmer_atoms
    real(dp), dimension(:), allocatable, target :: nmer_fock

    integer :: nmer_size, imer, nmer_charge, nmer_aosize
    integer :: nmer_nat
    integer :: nao

    ! nmer_size=1 is a monomer,
    ! nmer_size=2 is a dimer and so on
    nmer_size = size(nmer)

    ! get number of atoms in n-mer for memory.
    ! get total nmer charge too
    nmer_nat = 0
    nmer_charge = 0
    nmer_aosize = 0
    do imer=1,nmer_size
        nmer_nat = nmer_nat + fragment_atom_counts(nmer(imer))
        nmer_charge = nmer_charge + fragment_charges(nmer(imer))
        nmer_aosize = nmer_aosize + fragment_aosizes(nmer(imer))
    enddo
    allocate(nmer_atoms(nmer_nat))

    ! extract the atoms for the current n-mer
    call get_nmer_atoms(nmer,nmer_atoms)

    ! set the pointers
    p_mb_atoms => nmer_atoms
    p_mb_nmer => nmer

    nao = nmer_aosize*(nmer_aosize+1)/2
    if (nmer_aosize > 0 .and. doesp .and. &
        i_scc_iteration > 1) then
        allocate(nmer_fock(nao))
        p_mb_fock => nmer_fock

        ! calculate the ESP of the surroundings onto the current nmer
        call environmental_esp(nmer,nmer_charge,nmer_atoms,nmer_fock)
    endif

    ! Evaluate fragment SCF
    ! MB_FOCK is called from sirius to add environment
    ! MB_DENSITY is called from sirius to store density
    ! MB_STORE_NMER_PROPERTY is called from sirius, abacus
    !                         and (eventually) the response module

#if defined (PRG_DALTON)
    call execute_dalton(nmer_charge, nmer_nat, nmer_atoms)
#endif

    nullify(p_mb_atoms)
    nullify(p_mb_nmer)
    nullify(p_mb_fock)

    if (nmer_aosize > 0 .and. doesp .and. &
        i_scc_iteration > 1) then
        deallocate(nmer_fock)
    endif
    deallocate(nmer_atoms)

end subroutine calculate_nmer

!> @brief Evaluates the controbutions from the environmental esp onto the
!!        current nmer
!!
!! @details The potential is given as a sum of 1e- and 2e- contributions
!! of fragments @f$K@f$ is used
!! @f[
!!   V^X_{\mu\nu} = \sum_{K\not = X} u^K_{\mu\nu} + \upsilon^K_{\mu\nu}
!! @f]
!! where @f$ X @f$ is a generalized nmer
!! @todo clean up 1e- integrals so it looks more like 2e- part
subroutine environmental_esp(nmer, nmer_charge, nmer_atoms, nmer_fock)

    use mb_integral_interface, only : one_electron_integrals, &
                                       two_electron_integrals
    use mb_variables, only : nuclear_coordinates, nuclear_charges
    use mb_utilities, only : has_element, get_nmer_atoms, calc_n_kbodies

    integer, intent(in), dimension(:) :: nmer
    integer, intent(in) :: nmer_charge
    integer, intent(in), dimension(:) :: nmer_atoms
    real(dp), intent(out), dimension(:) :: nmer_fock
    ! ------
    integer :: ncalcs, jmer, nat, iat, nmer_ao
    integer :: nmer_nat
    integer, dimension(:,:), allocatable :: fragment_list
    integer, dimension(:), allocatable :: atoms
    real(dp), dimension(:,:), allocatable :: coords
    real(dp), dimension(:), allocatable :: charges

    ! 1e- integral contributions to the fock matrix
    real(dp), dimension(:), allocatable :: one_integrals
    ! 2e- integral contributions to the fock matrix
    real(dp), dimension(:), allocatable :: two_integrals

    nmer_ao = size(nmer_fock)
    nmer_nat = size(nmer_atoms)
    nmer_fock = 0.0d0

    allocate(one_integrals(nmer_ao))
    allocate(two_integrals(nmer_ao))

    ! we have to obtain a list of individual nmers
    call calc_n_kbodies(1, ncalcs)
    allocate(fragment_list(ncalcs,1))

    ! build list of fragments on which to calculate, look in the subroutine
    ! for the definition of the format
    call build_fragment_list(fragment_list)

#if defined (PRG_DALTON)
    ! ready the integrals in gen1int
    call execute_hermit(nmer_charge, nmer_nat, nmer_atoms)
#endif

    ! here we do 1e- contribution
    do jmer=1,size(fragment_list)
        if (has_element(fragment_list(jmer,1), nmer)) cycle

        nat = fragment_atom_counts(fragment_list(jmer,1))
        if(mb_debug) &
            write(luout,6000) nmer_ao,fragment_list(jmer,1),nat

        allocate(atoms(nat))
        allocate(coords(3,nat))
        allocate(charges(nat))

        ! extract the atoms for the current n-mer
        call get_nmer_atoms(fragment_list(jmer,:),atoms)

        ! now we have the atom indices, lets fetch the coordinates
        ! and the charges
        do iat=1,nat
            coords(:,iat) = nuclear_coordinates(:,atoms(iat))*aa2au

            ! if we do the true potential (1e- and 2e-) then use
            ! the true nuclear charges
            if (.not.use_atomic_charges) then
                charges(iat) = real(nuclear_charges(atoms(iat)))
            else
                charges(iat) = fragment_atomic_charges(1,atoms(iat))
            endif

            ! else use the partial atomic charges (mulliken)
            if(mb_debug) &
                write(luout,6005) charges(iat), coords(:,iat)
        enddo

        ! reset integrals before we compute them for the other nmer
        one_integrals = 0.0d0
        call one_electron_integrals(charges, coords, one_integrals)

        nmer_fock = nmer_fock + one_integrals

        deallocate(atoms)
        deallocate(coords)
        deallocate(charges)
    enddo
    deallocate(one_integrals)

    ! here we do 2e- coulomb contribution
    do jmer=1,size(fragment_list)
        if (has_element(fragment_list(jmer,1), nmer)) cycle

        ! reset integrals before we compute them
        two_integrals = 0.0d0
        if(.not.use_atomic_charges) &
            call two_electron_integrals(nmer, fragment_list(jmer,1), two_integrals)

        nmer_fock = nmer_fock + two_integrals
    enddo
    deallocate(two_integrals)

    deallocate(fragment_list)

6000 format(/2x,'ifg-nao:',i4,' jfg=',i3,' natj=',i3)
6005 format('Z=',f10.6,' x,y,z=',3f9.4)

end subroutine environmental_esp

!------------------------------------------------------------------------------
!> @brief Adds Fock matrix contributions to a fragment from the environment
!!
!! @author Casper Steinmann
!
!! @param[out] fock_matrix contribution to the Fock matrix
!!
!! @details The environmental contribution to the fock matrix has been calculated before
!!          so here it is merely added.
!!
!! @bug If one specifies ESP, but forgets the MULLIKEN part of the input
!!      then it does a gas phase calculation. Perhaps do inform the user what
!!      is going on and enable the mulliken charges unless .NOESP is selected.
subroutine mb_fock(fock_matrix)

    !----- INPUT -----
    real(dp), dimension(:), intent(out) :: fock_matrix

    fock_matrix = 0.0d0

    if (associated(p_mb_fock)) then
        if (size(p_mb_fock) /= size(fock_matrix)) then
            stop 'ERROR: fock matrix size inconsistency'
        endif
        fock_matrix = fock_matrix + p_mb_fock
    endif

end subroutine mb_fock

!------------------------------------------------------------------------------
!> @brief stores the density of the current nmer. only execute for monomers.
!!
!! @author Casper Steinmann
!!
!! @param[in] density density of the current nmer
!!
!! @note the densities for k-bodies larger than 1 are @b not stored.
subroutine store_density(density)

    use mb_utilities, only : openfile, get_max_difference, &
                              set_nmer_density, get_nmer_density
    !----- INPUT -----
    real(dp), dimension(:), intent(in) :: density
    !----- LOCAL -----
    real(dp), dimension(:), allocatable :: odensity
    integer :: lden, nao, imer, nao_offset, kmer
    real(dp) :: ddiff
    !character(len=99) :: cl
    !character(len=3) :: tcl


    imer = 0
    if (i_mb_kbody > 1) return
    imer = imb_current_calculation
    lden = size(density)
    nao = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * lden) - 1.0d0))

    ! get the filename of the density we want to store to
    !write(cl,*) imer
    !tcl = adjustl(cl)
    !write(cl,'(a,a,a)') 'fdens_',trim(tcl),'.bin'


    ! if we haven't been here before, open a new file
    if(fragment_aosizes(imer) == 0) then
        fragment_aosizes(imer) = nao
        !call openfile(trim(cl),lu,'new','unformatted')
        call set_nmer_density(imer,'new',density)
    else
        ! otherwise read the old density file to get the difference
        ! and check for convergence
        if(fragment_aosizes(imer) /= nao) then
            write(luout,9000) nao,imer,fragment_aosizes(imer)
            stop
        endif

        ! calculate the difference between the two density matrices
        ! to see if we have convergence
        allocate(odensity(lden))
        call get_nmer_density(imer,odensity)
        ddiff = get_max_difference(density,odensity)
        max_density_difference = max(ddiff, max_density_difference)
        deallocate(odensity)
        if (mb_debug) write(luout,8000) imer, ddiff

        !write(luout,*) "here1"
        ! find the place in the global density array and store the new density
        nao_offset = 0
        do kmer=1,imer
            nao = fragment_aosizes(kmer)
            nao_offset = nao_offset+nao*(nao+1)/2
        enddo
        ! It will be stored later in store_nmer_densities
        !write(luout,*) "here2", associated(p_densities), size(density), nao_offset, lden
        if (associated(p_densities)) then
            p_densities(1+nao_offset-lden:nao_offset) = density
        endif
        !write(luout,*) "here3"

    endif

8000 format(/1x,"MB DEBUG: nmer-",i4," density difference =",f10.6)
9000 format(/1x,"ERROR: attempt to write",I5," records for",/1x, &
     & "density",I4," but found that",I5,"was previous written")

end subroutine store_density

!------------------------------------------------------------------------------
!> @brief stores an nmer property using predefined actions for each property.
!!
!! @author Casper Steinmann
!!
!! @param[in] property property type, 'energy' etc.
!! @param[in] property_value value to store
subroutine mb_store_nmer_property(property, property_value)

    character(*), intent(in) :: property
    real(dp), dimension(:), intent(in) :: property_value

    integer :: nva, nat, i

    ! currently, the following properties are 'allowed'
    !
    ! Energy
    ! Mulliken charges
    ! Density
    ! Dipole

    ! requires: global pointer to atoms

    ! the energies are stored per nmer calculations
    if (property == 'energy') &
        & fragment_energies(imb_current_calculation) = property_value(1)

    ! store correlated energies similarly
    if (property == 'correnergy') &
        & fragment_corr_energies(imb_current_calculation) = property_value(1)

    ! only store the embedding energy for 1-mers to correct double
    ! counting the 1-body term
    if (property == 'embedding_energy' .and. i_mb_kbody == 1) &
        & fragment_embed_energies(imb_current_calculation) = property_value(1)

    ! mulliken charges are stored on an atom basis
    if (property == 'mulliken') then
        nat = size(p_mb_atoms)
        nva = size(property_value)
        if (nat /= nva) then
            write(luout,*) "ERROR: Wrong size of properties"
            stop
        endif
        do i=1,nat

        ! reset one-body property every iteration to update the charges
        if (i_mb_kbody==1) &
          & fragment_atomic_charges(i_mb_kbody,p_mb_atoms(i)) = 0.0d0

            fragment_atomic_charges(i_mb_kbody,p_mb_atoms(i)) = &
            fragment_atomic_charges(i_mb_kbody,p_mb_atoms(i)) + &
                                 & property_value(i)
        enddo
        has_atomic_charges = .true.
    endif

    if (property == 'density') call store_density(property_value)

    if (property == 'dipole') then
        fragment_dipoles(:,imb_current_calculation) = property_value(:)
        has_dipole = .true.
    endif

    if (property == 'molgrad') then
        ! there are always 3 x nat elements we need and they are
        ! stored as x1, y1, z1, ..., xn, yn, zn
        nat = size(p_mb_atoms)
        do i=1,nat

        if (i_mb_kbody == 1) &
            & fragment_molecular_gradient(i_mb_kbody,p_mb_atoms(i), :) = 0.0d0

        fragment_molecular_gradient(i_mb_kbody,p_mb_atoms(i), :) = &
        fragment_molecular_gradient(i_mb_kbody,p_mb_atoms(i), :) + &
                                      & property_value(1+3*(i-1):3*i)

        enddo
        has_molecular_gradient = .true.
    endif

end subroutine mb_store_nmer_property

!------------------------------------------------------------------------------
!> @brief initializes mb variables and prints out information on the run
!!
!! @author Casper Steinmann
!!
!! @param[in] coords nuclear coordinates
!! @param[in] charges nuclear charges
subroutine mb_init(coords, charges)

    real(dp), dimension(:), intent(in), optional :: charges
    real(dp), dimension(:,:), intent(in), optional :: coords

    if (mb_initialized) return

    ! a last series of checks before we let the many-body run loose
    if (size(charges) /= natmb) then
        write(luout,9000) size(charges), natmb
        stop
    endif

    ! store charges and coords in mb-arrays
    nuclear_coordinates = coords
    nuclear_charges = charges

    write(luout,'(//2x,a)')  'Generalized Embedded Many-Body Method'
    write(luout,'(2x,a)')    '---------------------------------------'
    write(luout,'(2x,a,i4)') 'nbody        = ', nbody
    write(luout,'(2x,a,i4)') 'nfrag        = ', nfrag
    write(luout,'(2x,a,i4)') 'total charge = ', sum(fragment_charges)
    write(luout,'(2x,a,a)')  'basis set    = ', mb_basis_name
    write(luout,'(2x,a,i4)') 'scc iters    = ', max_scc_iterations
    write(luout,'(2x,a,es12.3)') 'scc D tol    = ', scc_d_convergence
    if(.not. doesp) then
        write(luout,'(/2x,a)')    'Electrostatic potential will be ignored.'
    else
        write(luout,'(/2x,a)')    'Electrostatic potential will be included:'
        if (use_atomic_charges) then
            write(luout,'(6x,a)')    'Using atomic (Mulliken) charges.'
        else
            write(luout,'(6x,a)')    'Using full 1e- and 2e- potential.'
        endif
    endif

    mb_initialized = .true.

9000 format('ERROR: Number of atoms in the .mol file (',I4,') and ', &
    & 'in .INDAT in the .dal file (',I4,') must match!')

end subroutine mb_init

!------------------------------------------------------------------------------
!> @brief Reads the input from the file luinp
!!
!! @author Casper Steinmann
!!
!! @details the input is read TWICE to allow for memory allocation to occur
!!          before actual readin of parameters
!!
!! @param[in] word
!! @param[in] luinp
!! @param[in] lupri
subroutine mb_input(word, luinp, lupri)

    use mb_input_readers

    character(len=7), intent(inout) :: word
    integer, intent(in) :: luinp
    integer, intent(in) :: lupri

#if defined (PRG_DALTON)
    call dalton_input(word, luinp, lupri)
#endif

end subroutine mb_input


!------------------------------------------------------------------------------
!> @brief sums up and prints the energies (and properties) for the many-body calculation
!!
!! @author Casper Steinmann
subroutine mb_output

    use mb_utilities, only : calc_n_kbodies

    real(dp), external :: ddot

    real(dp) :: ibody_energy, ibody_corr_energy

    integer :: ibody, icalc, ioff, k, idip
    integer :: ncalcs
    integer, dimension(:,:), allocatable :: fragment_list

    real(dp), dimension(3) :: dipdip

    ! format statement for energies
    integer :: neg_labels = 1
    character(len=80) :: label, header, separator
    character(len=1), dimension(5) :: id_labels
    character(len=16), dimension(2) :: eg_labels

    ! storage for the i-body'th energy
    real(dp), dimension(:), allocatable :: energy
    real(dp), dimension(:), allocatable :: energy_corr
    ! storage for the i-body'th dipole
    real(dp), dimension(:,:), allocatable :: dipole
    ! storage for the i-body'th atomic charges
    real(dp), dimension(:,:), allocatable :: atomic_charges
    ! storage for the i-body'th molecular gradient
    real(dp), dimension(:,:,:), allocatable :: molecular_gradient

    id_labels = (/"I", "J", "K", "L", "M"/)
    eg_labels = (/"      Euncorr   ", "      Ecorr     "/)

    write(luout,'(//2x,a)')  'Generalized Embedded Many-Body Method'
    write(luout,'(2x,a)')    '-------------------------------------'

    allocate(energy(nbody))
    allocate(energy_corr(nbody))
    allocate(dipole(3,nbody))
    allocate(atomic_charges(nbody,natmb))
    allocate(molecular_gradient(nbody,natmb,3))

    ! if the sum of correlated energies is not zero, then we calculated them
    if (sum(fragment_corr_energies) .ne. 0.0d0) &
        neg_labels = 2

    ioff = 0
    do ibody = 1,nbody
        ! setup the correct accounting information for properties
        call calc_n_kbodies(ibody, ncalcs)
        allocate(fragment_list(ncalcs,ibody))
        call build_fragment_list(fragment_list)

        ! write the correct labels and headers depending on what is calculated
        write(label,'("(",i2,"i4,",i2,"f20.9)")') ibody, neg_labels
        write(header,'("(",i2,"a4,",i2,"a20)")') ibody, neg_labels
        write(separator, '("(",i2,"(1H-),",i2,"(1H-))")') 4*ibody, 20*neg_labels

        ! write the i-body header
        write(luout,'(//2x,i3,a)')  ibody,'-body result(s)'
        write(luout,'(2x,a,/)') '===================='

        write(luout, header) id_labels(1:ibody), eg_labels(1:neg_labels)
        write(luout, separator)
        do icalc=1,ncalcs
            if (neg_labels == 1) then
                write(luout,label) (fragment_list(icalc,k),k=1,ibody), fragment_energies(icalc+ioff)
            else
                write(luout,label) (fragment_list(icalc,k),k=1,ibody), &
                    & fragment_energies(icalc+ioff), &
                    & fragment_corr_energies(icalc+ioff)
            endif
        enddo

        !------------
        ! total i-body energy
        energy(ibody) = sum(fragment_energies(1+ioff:ncalcs+ioff))
        energy_corr(ibody) = sum(fragment_corr_energies(1+ioff:ncalcs+ioff))

        ibody_energy = ddot(ibody, energy(1:ibody), 1, fragment_binomials(ibody,1:ibody), 1)
        ibody_corr_energy = ddot(ibody, energy_corr(1:ibody), 1, fragment_binomials(ibody,1:ibody), 1)

        ! remove the embedding energy (if present) from the monomers to
        ! avoid double counting
        if (ibody==1) then
            ibody_energy = ibody_energy - sum(fragment_embed_energies)
            ibody_corr_energy = ibody_corr_energy - sum(fragment_embed_energies)
        endif
        write(luout,9000) ibody, ibody_energy
        if (neg_labels == 2) write(luout,9001) ibody, ibody_corr_energy
        !------------

        !------------
        ! dipole
        if (has_dipole) then
            do idip=1,3
                dipole(idip,ibody) = sum(fragment_dipoles(idip,1+ioff:ncalcs+ioff))
            enddo
            do idip=1,3
               dipdip(idip) = ddot(ibody, dipole(idip,1:ibody),1,fragment_binomials(ibody,1:ibody),1)
            enddo
            write(luout,9005) ibody, dipdip, sqrt(ddot(3,dipdip,1,dipdip,1))
        endif
        !------------

        !------------
        ! write the atomic charges
        if (has_atomic_charges) then
            write(luout,'(//2x,i3,a)') ibody,'-body mulliken charges'
            write(luout,'(2x,a)') '---------------------------'
            atomic_charges = 0.0d0
            do k=1,ibody
               atomic_charges(k,:) = fragment_atomic_charges(k,:) * &
                        & fragment_binomials(ibody,k)
            enddo
            do k=1,natmb
                write(luout,9010) k,sum(atomic_charges(:,k))
            enddo
        endif
        !------------

        !------------
        ! write the molecular gradient
        if (has_molecular_gradient) then
            write(luout,'(//2x,i3,a)') ibody,'-body molecular gradient'
            write(luout,'(2x,a)') '-----------------------------'
            molecular_gradient = 0.0d0
            do k=1,ibody
               molecular_gradient(k,:,:) = fragment_molecular_gradient(k,:,:) * &
                        & fragment_binomials(ibody,k)
            enddo
            do k=1,natmb
                write(luout,9015) k, sum(molecular_gradient(:,k,1)), &
                                     sum(molecular_gradient(:,k,2)), &
                                     sum(molecular_gradient(:,k,3))
            enddo
        endif
        !------------

        ioff = ioff+ncalcs
        deallocate(fragment_list)
    enddo

    deallocate(energy)
    deallocate(dipole)
    deallocate(atomic_charges)

    write(luout,*)

9000 format(/3x,'Eu(',i2,') =',f20.12)
9001 format(/3x,'Ec(',i2,') =',f20.12)
9005 format(/3x,'D(',i4,') =',3f12.5,' |D| =',f12.5)
9010 format(4x,'Q(',i2,') =',f12.5)
9015 format(4x,'x(',i2,') =',3f20.12)

end subroutine

!> @brief wrapper routine for the calculation of two-electron integrals
!!        needed for the exact coulomb embedding potential.
!! @param[in] nbas number of atomic orbitals
!! @param[in,out] wrk work memory
!! @param[out] iconv whether or not the calculation is converged.
!! @details The iconv flag is used to signal whether any calling subroutine
!!          should exit or not. In the case of DALTON, the SCF driver is called
!!          but signalling iconv=1 makes it exit gracefully from that calculation
!!          before the SCF is progressing.
!!          Similar behavior is expected in other programs
subroutine mb_twoints(nbas, wrk, iconv)

    integer, intent(in) :: nbas
    real(dp), intent(inout), dimension(:) :: wrk
    integer, intent(out) :: iconv

    iconv = 0
#if defined (PRG_DALTON)
    call dalton_twoints(nbas, wrk, iconv)
#endif

end subroutine mb_twoints


end module mb
