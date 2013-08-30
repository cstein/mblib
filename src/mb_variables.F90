!> @file
!! contains the MBLIB variables module with global variables
!! for the many-body calculations.

!> @brief module that holds all MBLIB-specific variables
!!
!! @author Casper Steinmann
module mb_variables

    use mb_precision

    implicit none

    !> whether or not we are running MBLIB.
    logical, save :: mbrun = .false.
    !> whether or not to be more verbose with output from
    !! the MBLIB module.
    logical, save :: mb_verbose = .false.
    !> whether or not to print debug output.
    logical, save :: mb_debug = .false.
    !> whether or not mb has been initialized.
    logical, save :: mb_initialized = .false.

    !
    ! Run-time settings to be read in from the .dal file
    !
    !> n-body level of the calculation
    integer, save :: nbody = 0
    !> number of fragments
    integer, save :: nfrag = 0
    !> total number of atoms. This value is given in INDAT
    !! and must match the number of atoms in the input file
    integer, save :: natmb = 0
    !> whether or not to calculate the ESP
    logical, save :: doesp = .true.
    !> maximum number of self-consistent charge (monomer) iterations
    !! @see i_scc_iterations
    integer, save :: max_scc_iterations = 30
    !> ! convergence threshold for change in density for scc iterations
    real(dp), save :: scc_d_convergence = 1.0d-6
    !> the integer charge of the fragments
    integer, dimension(:), save, allocatable :: fragment_charges
    !> basis set name to be used for all fragments
    character(len=80), save :: mb_basis_name = ''
    !> whether or not to use atomic charges instead of the true 2e- potential
    !! when evaluating the coulomb potential.
    !! @note this can later be changed into some distance criterium (RESPPC)
    logical, save :: use_atomic_charges = .false.

    !
    ! properties that are calculated based on the input
    !
    !> holds the fragment indices of each atom
    integer, dimension(:), save, allocatable :: atom_in_fragment_ids
    !> the number of atoms in each fragment
    integer, dimension(:), save, allocatable :: fragment_atom_counts
    !> binomial coefficients to calculate the n-body energy
    real(dp), dimension(:,:), save, allocatable :: fragment_binomials
    !> nuclear charges of \b all atoms in the calculation
    real(dp), dimension(:), save, allocatable :: nuclear_charges
    !> nuclear coordinates of \b all atoms in the calculation
    real(dp), dimension(:,:), save, allocatable :: nuclear_coordinates

    !
    ! properties that are calculated during or when the calculation is finished
    !
    ! fragment-related properties such as energies
    !> calculated fragment energies for @b all fragments (1-mers, 2-mers and so on ...)
    real(dp), dimension(:), save, allocatable :: fragment_energies
    !> calculated fragment embedding energies for @b 1-body fragments
    real(dp), dimension(:), save, allocatable :: fragment_embed_energies
    !> whether calculated dipoles are available
    logical, save :: has_dipole = .false.
    !> calculated fragment dipoles for @b all fragments
    real(dp), dimension(:,:), save, allocatable :: fragment_dipoles
    !> number of basis functions per fragment
    integer, dimension(:), save, allocatable :: fragment_aosizes
    !
    ! atom-related properties such as mulliken charges
    !> whether calculated atomic charges are available
    logical, save :: has_atomic_charges = .false.
    !> atomic charges (mulliken or others)
    real(dp), dimension(:,:), save, allocatable :: fragment_atomic_charges


    !> unit to write to for output (default is stdout)
    integer, save :: luout = 6
    !> unit for the storage of densitities
    integer, save :: luden


    ! constants
    !> @f$ \pi @f$
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: aa2au = 1.0_dp / 0.5291772109217_dp

    ! pointers for mbody calculations.
    ! NB! These will most likely fuck stuff up during parallel runs
    !     so think this trough in the future.
    !     or perhaps it will not.
    !> pointer to the current nmer
    integer, pointer :: p_mb_nmer(:)
    !> pointer to the current nmer atoms
    integer, pointer :: p_mb_atoms(:)
    !> pointer to the current nmer fock matrix esp
    real(dp), pointer :: p_mb_fock(:)
    real(dp), pointer :: p_nmer_fock(:)
    !> pointer to the array of densities
    real(dp), pointer :: p_densities(:)

    !> pointer to the current kmer density. Yes kmer. used during the 2e- integrals
    real(dp), pointer :: p_kmer_dens(:)

    ! global values to be stored
    real(dp), save :: max_density_difference = 0.0d0

    ! memory pointers in DALTON
    !> work memory array
    real(dp), dimension(:), pointer :: work
    !> work memory size
    integer, save :: idalton_work_size
    !> work memory control value
    real(dp), save :: dalton_work_control

    !
    ! global counters
    !
    !> current nmer calculation
    !! @warning this will likely not work well in parallel because it is used for storage
    integer, save :: imb_current_calculation
    !> signalling the current k-body level.
    integer, save :: i_mb_kbody
    !> the current self-consistent charge (monomer) iteration
    !! @see max_scc_iteration
    integer, save :: i_scc_iteration

end module mb_variables
