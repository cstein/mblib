!------------------------------------------------------------------------------
!> @brief Integral interface for the FMO method.
!!
!! @author Casper Steinmann
module mb_integral_interface

    use mb_precision
    use mb_variables

    implicit none

    private

    public :: one_electron_integrals
    public :: two_electron_integrals

    contains

!------------------------------------------------------------------------------
!> @brief evaluates the one-electron contribution to the fock matrix.
!!
!! @author Casper Steinmann
!!
!! @param[in] charges Charges of the nuclei of fragment @f$ K @f$
!! @param[in] coords Coordinates of the nuclei of fragment @f$ K @f$
!! @param[out] integrals the one-electron contribution to the fock matrix
!!
!! @note lifted and adapted from the polarizable embedding library
!!
!! @details If the two-electron part is included explcitly, the nuclear charges
!! of fragment @f$K@f$ is used
!! @f[
!!   u^{K}_{\mu\nu} = -\sum_{A\in K} \big \langle \mu \big| \frac{Z_A}{|\mathbf{r}-\mathbf{R}_A|} \big| \nu \big\rangle
!! @f]
!!
!! An approximation is to change the two-electron contribution into a
!! one-electron part using mulliken charges yielding the following
!! expression for the one-electron terms
!! @f[
!!   u^{K}_{\mu\nu} = -\sum_{A\in K} \big \langle \mu \big| \frac{Z_A-Q_A}{|\mathbf{r}-\mathbf{R}_A|} \big| \nu \big\rangle
!! @f]
subroutine one_electron_integrals(charges, coords, integrals)

    real(dp), intent(in), dimension(:) :: charges
    real(dp), intent(in), dimension(:,:) :: coords
    real(dp), dimension(:), intent(out) :: integrals

    integrals = 0.0d0

#if defined (BUILD_GEN1INT)
    call one_e_integrals_gen1int(charges, coords, integrals)
#endif

end subroutine one_electron_integrals

#if defined (BUILD_GEN1INT)

subroutine one_e_integrals_gen1int(charges, coords, integrals)

    use gen1int_api

    real(dp), intent(in), dimension(:) :: charges
    real(dp), intent(in), dimension(:,:) :: coords
    real(dp), dimension(:), intent(out) :: integrals

    ! -- GEN1INT VARIABLES --
    integer :: ierr
    integer :: prop_sym
    integer :: num_prop
    integer :: num_ao
    integer :: num_geo_bra
    integer :: num_geo_ket
    integer :: num_geo_total
    integer :: io
    integer :: printlvl = 0
    logical :: symmetric
    logical :: triangular
    type(one_prop_t) :: prop_operator
    type(nary_tree_t) :: nary_tree_bra
    type(nary_tree_t) :: nary_tree_ket
    type(nary_tree_t) :: nary_tree_total
    type(matrix), dimension(:), allocatable :: intmats

    ! non-zero components for the operator, the first dimension is for bra and
    ! ket sub-shells, the last is the number of non-zero components, which should
    ! be 1 for non-relativistic calcualtions
    integer nnz_comp(2,1)

    ! -- SUBROUTINE VARIABLES --
    integer :: nat, nnbas
    integer, dimension(:), allocatable :: nuclei
    real(dp), dimension(:), allocatable :: real_charges

    integrals = 0.0d0
    nnbas = size(integrals)

    ! treat the nuclei as non-nuclei because we want to also
    ! be able to evaluate other types of atomic charges
    nat = size(charges)
    allocate(nuclei(nat))
    allocate(real_charges(nat))
    nuclei = -1

    ! sets the non-zero components for the one-electron operator, here we only
    ! consider the (large,large) part
    nnz_comp(1,1) = 1
    nnz_comp(2,1) = 1

    ! do the correct sign according to the operator
    real_charges = -1.0d0 * charges

    call OnePropCreate(prop_name=INT_POT_ENERGY, &
                       one_prop=prop_operator,   &
                       info_prop=ierr,           &
                       idx_nuclei=nuclei,        &
                       coord_nuclei=coords,      &
                       charge_nuclei=real_charges,     &
                       order_geo_pot=0)

    if (ierr /= 0) stop 'ERROR: Failed to create property operator'

    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeCreate(max_num_cent=0,      &
                                  order_geo=0,         &
                                  num_geo_atoms=0,     &
                                  idx_geo_atoms=(/0/), &
                                  nary_tree=nary_tree_total)

    ! gets the number of property integrals and their symmetry
    call OnePropGetNumProp(one_prop=prop_operator, &
                           num_prop=num_prop)
    call OnePropGetSymmetry(one_prop=prop_operator, &
                            prop_sym=prop_sym)

    ! gets the number of geometric derivatives
    call NaryTreeGetNumGeo(nary_tree=nary_tree_bra, num_unique_geo=num_geo_bra)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_ket, num_unique_geo=num_geo_ket)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_total, num_unique_geo=num_geo_total)

    call Gen1IntAPIGetNumAO(num_ao=num_ao)
    !write(luout,*) "gen1int: num_props", num_prop, num_ao

    if (num_prop /= 1) stop 'ERROR: Integral property failed.'

    triangular = .true.
    symmetric = (prop_sym == SYMM_INT_MAT)

    allocate(intmats(num_prop), stat=ierr)

    call MatAssociate(work_alpha=integrals(:), &
                      num_row=num_ao,         &
                      A=intmats(num_prop),     &
                      info_mat=ierr,           &
                      triangular=triangular,   &
                      symmetric=symmetric)

    if (ierr /= 0) stop 'ERROR: Failed to associate matrices.'

    call Gen1IntOnePropGetIntExpt(nnz_comp=nnz_comp, &
                                  one_prop=prop_operator, &
                                  nary_tree_bra=nary_tree_bra, &
                                  nary_tree_ket=nary_tree_ket, &
                                  nary_tree_total=nary_tree_total, &
                                  num_ints=num_prop, &
                                  val_ints=intmats, &
                                  num_dens=1,       &
                                  io_viewer=io,       &
                                  level_print=printlvl)

    call OnePropDestroy(one_prop=prop_operator)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)

    call MatNullify(A=intmats(num_prop))
    deallocate(intmats)
    deallocate(nuclei)
    deallocate(real_charges)

end subroutine one_e_integrals_gen1int

#endif

!------
!> @brief evaluates the two-electron contribution to the fock matrix.
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer nmer that is embedded in ESP from kmer @f$ K @f$
!! @param[in] kmer the kmer that the nmer feels
!! @param[out] integrals The two-electron contribution to the fock matrix
!!
!! @note lifted and adapted from the polarizable embedding library
!!
!! @details The two-electron part (if included explicitly) from fragment @f$K@f$ is
!! @f[
!!   \upsilon^{K}_{\mu\nu} = \sum_{\sigma\lambda\in K} D^K_{\lambda\sigma} \big( \mu \nu \big| \lambda \sigma \big)
!! @f]
subroutine two_electron_integrals(nmer, kmer, integrals)

    use mb_utilities, only : get_nmer_charge, &
                              get_nmer_atoms,  &
                              get_nmer_aosize, &
                              get_nmer_density,&
                              get_nmer_atom_count

    integer, dimension(:), intent(in) :: nmer
    integer, intent(in) :: kmer
    real(dp), dimension(:), intent(out) :: integrals
    ! -----
    integer, dimension(:), allocatable :: mmer
    real(dp), dimension(:), allocatable, target :: kmer_dens
    real(dp), dimension(:), allocatable, target :: nmer_fock
    integer :: nmer_size, imer, kmer_nao, nmer_nao
    integer, dimension(1) :: kmer_array

    ! K-Nmer fragment variables
    integer, allocatable, dimension(:) :: atoms
    integer :: charge, nat

    integrals = 0.0d0
    nmer_size = size(nmer)

    ! build fragment K-N and extract the overall charge, number of
    ! atoms and the atoms of the K-Nmer
    allocate(mmer(nmer_size+1))
    mmer(1) = kmer
    do imer=2,nmer_size+1
        mmer(imer) = nmer(imer-1)
    enddo
    call get_nmer_atom_count(mmer, nat)
    call get_nmer_charge(mmer, charge)
    allocate(atoms(nat))
    call get_nmer_atoms(mmer,atoms)

    ! load the density of the kmer
    kmer_array(1) = kmer
    call get_nmer_aosize(kmer_array, kmer_nao)
    allocate(kmer_dens(kmer_nao*(kmer_nao+1)/2))
    call get_nmer_density(kmer, kmer_dens)
    p_kmer_dens => kmer_dens

    ! ready storage for the fock contribution to the nmer
    call get_nmer_aosize(nmer, nmer_nao)
    allocate(nmer_fock(nmer_nao*(nmer_nao+1)/2))
    p_nmer_fock => nmer_fock

    ! contributions that are calculated below are added
    ! to the p_nmer_fock pointer
#if defined (PRG_DALTON)
    call execute_integral_dalton(charge, nat, atoms)
#endif

    integrals = nmer_fock

    nullify(p_nmer_fock)
    nullify(p_kmer_dens)
    deallocate(nmer_fock)
    deallocate(kmer_dens)
    deallocate(atoms)
    deallocate(mmer)

end subroutine two_electron_integrals
end module mb_integral_interface
