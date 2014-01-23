!> @brief Utilities for the MBLIB library
!!
!! @author Casper Steinmann
!!
!! @note Some of these subroutines come from other people who willingly
!!       donated their work
module mb_utilities

    use mb_precision

    implicit none

    contains

!------------------------------------------------------------------------------
!> @brief Calculates the factorial of n
!!
!! @author Jógvan Magnus Olsen
!!
!! @see http://www.star.le.ac.uk/~cgp/f90course/f90.html
recursive function factorial(n) result(nfact)

    integer, intent(in) :: n
    integer :: nfact

    if (n > 0) then
        nfact = n * factorial(n-1)
    else
        nfact = 1
    end if

end function factorial

!------------------------------------------------------------------------------
!> @brief Opens a file. Returns the first available unit in lunit
!!
!! @author Jógvan Magnus Olsen
!!
!! @param[in] filename the filename to use for the file
!! @param[out] lunit file handle
!! @param[in] stat status of the file. 'OLD', 'NEW' and 'SCRATCH' are supported
!! @param[in] frmt format
subroutine openfile(filename, lunit, stat, frmt)

    character(*), intent(in) :: filename, stat, frmt
    integer, intent(out) :: lunit
    integer :: i
    logical :: lexist, lopen

    if (stat == 'old') then
        inquire(file=filename, exist=lexist)

        if (.not. lexist) then
            print *, filename, ' not found!'
            stop
        end if
    end if

    do i = 21, 99
        inquire(unit=i, opened=lopen)
        if (lopen) then
            cycle
        else
            lunit = i
            open(unit=lunit, file=filename, status=stat, form=frmt)
            exit
        end if
    end do

    return

end subroutine openfile

!------------------------------------------------------------------------------
!> @brief Calculate binomial coefficients using an efficient recursive formula
!!
!! @author Jógvan Magnus Olsen
!!
!! @see S. Rettrup, R. Pauncz, Int. J. Quant. Chem. 60, 91 (1996)
!!
!! @param[in] n
!! @param[in] k
!!
!! @details It is provided in courtesy of Sten Rettrup et al.
function binomial_coefficient(n,k) result(binom)

    implicit none
    integer, intent(in) :: n, k
    integer :: i, j, nk
    integer*8 :: binom
    integer*8, dimension(0:n) :: w

    if ((k .eq. 0) .or. (n .eq. k)) then
        binom = 1
        return
    end if
    nk = n - k

    w(0:nk) = 1
    !forall(j = 0:nk)
    !  w(j) = 1
    !end forall

    do i = 1, k
        binom = 0
        do j = 0, nk
            binom = binom + w(j)
            w(j) = binom
        end do
     end do

end function

!------------------------------------------------------------------------------
!> @brief returns the maximum \em absolute difference between two vectors
!!
!! @author Casper Steinmann
!!
!! @param[in] v1 the first vector to use for comparison
!! @param[in] v2 the second vector to use for comparison
function get_max_difference(v1,v2) result(diff)

    real(dp), intent(in), dimension(:) :: v1, v2
    real(dp) :: diff

    real(dp) :: adiff
    integer :: n1,n2,i

    n1 = size(v1)
    n2 = size(v2)

    if(n1 /= n2) then
        write(*,*) 'ERROR: get_max_difference argument lengths ', &
                 & 'do not match.'
        stop
    endif

    diff = 0.0d0
    do i=1,n2
        adiff = abs(v1(i)-v2(i))
        if(adiff > diff) diff = adiff
    enddo

end function get_max_difference

!------------------------------------------------------------------------------
!> @brief returns true if an element is in the list
!!
!! @author Casper Steinmann
!!
!! @param[in] element value to search for in the list
!! @param[in] list list of values to search for element
function has_element(element, list)

    integer, intent(in) :: element
    integer, dimension(:), intent(in) :: list

    logical :: has_element
    integer :: ielem

    has_element = .false.
    do ielem=1,size(list)
        if (list(ielem) == element) then
            has_element = .true.
            exit
        endif
    enddo

end function has_element


!------------------------------------------------------------------------------
!> @brief unfolds a packed matrix into full storage
!!
!! @author Casper Steinmann
!!
!! @param[in] folded packed matrix in upper triangular storage
!! @param[out] unfolded matrix in full storage
subroutine unfold_matrix(folded, unfolded)
 
    real(dp), intent(in), dimension(:) :: folded
    real(dp), intent(out), dimension(:,:) :: unfolded

    integer :: nao
    integer :: i,j,l

    nao = size(unfolded,1)

    if ((nao*(nao+1)/2) /= size(folded)) then
        stop 'ERROR: could not unfold matrix.'
    endif

    l=1
    do i=1,nao
        do j=1,i
            if(j == i) then
                unfolded(i,i) = folded(l)
            else
                unfolded(i,j) = 0.5d0 * folded(l)
                unfolded(j,i) = 0.5d0 * folded(l)
            endif
            l = l+1
        enddo
    enddo

end subroutine

!------------------------------------------------------------------------------
!> @brief Extracts the coordinates given a set of atoms.
!!
!! @author Casper Steinmann
!!
!! @param[in] atoms the atom ids to use for the coordinates
!! @param[out] coords the coordinates of the currently selected atoms
!! @note Coordinates are returned in Ångstrom
!! @bug if the coordinates are read from the .mol file in Bohr, then they
!! are returned in bohr as well.
subroutine get_nmer_coords(atoms, coords)

    use mb_variables, only : nuclear_coordinates

    integer, dimension(:), intent(in) :: atoms
    real(dp), dimension(:,:), intent(out) :: coords

    integer :: iat

    do iat=1,size(atoms)
        coords(:,iat) = nuclear_coordinates(:,atoms(iat))
    enddo

end subroutine

!------------------------------------------------------------------------------
!> @brief Calculates the integer charge for an nmer
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer the nmer of which to extract the integer charge
!! @param[out] charge the integer charge of the nmer
subroutine get_nmer_charge(nmer, charge)

    use mb_variables, only : fragment_charges

    integer, dimension(:), intent(in) :: nmer
    integer, intent(out) :: charge

    integer :: imer

    charge = 0
    do imer=1,size(nmer)
        charge = charge + fragment_charges(nmer(imer))
    enddo

end subroutine

!------------------------------------------------------------------------------
!> @brief Calculates the number of atomic orbitals (along one dimension) for an nmer
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer the nmer for which to calculate the number of atomic orbitals
!! @param[out] aosize the number of atomic orbitals for the nmer
subroutine get_nmer_aosize(nmer, aosize)

    use mb_variables, only : fragment_aosizes

    integer, dimension(:), intent(in) :: nmer
    integer, intent(out) :: aosize

    integer :: imer

    aosize = 0
    do imer=1,size(nmer)
        aosize = aosize + fragment_aosizes(nmer(imer))
    enddo

end subroutine

!------------------------------------------------------------------------------
!> @brief Writes the density of an nmer to disk.
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer the nmer for which to store the density.
!! @param[in] stat 'OLD' or 'NEW' depending on status
!! @param[out] density the density of the nmer.
!! @info only to be called by 1-body n-mers
subroutine set_nmer_density(nmer,stat,density)

    !use mb_utilities, only: openfile
    use mb_variables, only : luout, mb_debug

    integer, intent(in) :: nmer
    character(*), intent(in) :: stat
    real(dp), intent(in), dimension(:) :: density
    !---------------------
    character(len=99) :: cl
    character(len=3) :: tcl
    integer :: lu

    ! get the filename of the density we want to write to
    write(cl,*) nmer
    tcl = adjustl(cl)
    write(cl,'(a,a,a)') 'fdens_',trim(tcl),'.bin'

    if(mb_debug) write(luout,9000) size(density), trim(cl)

    call openfile(trim(cl),lu,stat,'unformatted')
    rewind(lu)
    write(lu) density
    close(lu)

9000 format(/1x,'DEBUG: Writing',I5,' records to "',a,'"')

end subroutine

!------------------------------------------------------------------------------
!> @brief Reads the density of an nmer from disk.
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer the nmer for which to extract the density.
!! @param[out] density the density of the nmer.
!! @note The density is read directly from disk. Size is known indirectly from
!!       the size of the density matrix used for storage.
subroutine get_nmer_density(nmer,density)

    !use mb_utilities, only: openfile
    use mb_variables, only : luout, mb_debug

    integer, intent(in) :: nmer
    real(dp), intent(out), dimension(:) :: density
    !---------------------
    character(len=99) :: cl
    character(len=3) :: tcl
    integer :: lu

    ! get the filename of the density we want to read from
    write(cl,*) nmer
    tcl = adjustl(cl)
    write(cl,'(a,a,a)') 'fdens_',trim(tcl),'.bin'

    if(mb_debug) write(luout,9000) size(density), trim(cl)

    call openfile(trim(cl),lu,'old','unformatted')
    rewind(lu)
    read(lu) density
    close(lu)

9000 format(/1x,'DEBUG: Reading',I5,' records from "',a,'"')

end subroutine

!------------------------------------------------------------------------------
!> @brief calculates the total number of atoms in an nmer
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer The nmer for which to calculate the number of atoms
!! @param[out] nat The number of atoms in the nmer
subroutine get_nmer_atom_count(nmer, nat)

    use mb_variables, only : fragment_atom_counts, luout

    integer, dimension(:), intent(in) :: nmer
    integer, intent(out) :: nat

    integer :: imer

    nat = 0
    do imer=1,size(nmer)
        nat = nat + fragment_atom_counts(nmer(imer))
    enddo

    if (nat == 0) write(luout,*) "WARNING: nmer does not contain any atoms!"

end subroutine get_nmer_atom_count

!------------------------------------------------------------------------------
!> @brief fills out fragment_binomials based on the n-body level of the calculation
!!
!! @author Casper Steinmann
!!
!! @param[in] nbody The n-body level at which we evaluate the energy
!! @param[out] fragment_binomials lookup table with coefficient binomials
subroutine fill_fragment_binomials(nbody,fragment_binomials)

    !use mb_utilities, only : binomial_coefficient
    use mb_variables, only : nfrag

    integer, intent(in) :: nbody
    real(dp), intent(out), dimension(nbody,nbody) :: fragment_binomials

    integer :: ibody, kbody
    real(dp) :: bc, sig

    fragment_binomials = 0.0d0
    do kbody = 1,nbody
        do ibody = 1,kbody
            sig = (-1.0d0)**(ibody+kbody)
            bc = binomial_coefficient(nfrag-ibody-1, kbody-ibody)
            fragment_binomials(kbody,ibody) = sig*bc
        enddo
    enddo

end subroutine fill_fragment_binomials

!------------------------------------------------------------------------------
!> @brief obtains the atoms of the nmer given as input
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer the nmer of which to extract the atoms
!! @param[out] atoms the atoms for the nmer
subroutine get_nmer_atoms(nmer, atoms)

    use mb_variables, only : natmb, atom_in_fragment_ids

    integer, intent(in), dimension(:) :: nmer
    integer, intent(out), dimension(:) :: atoms

    integer :: imer, iat, id, nmer_size

    nmer_size = size(nmer)
    id = 1
    do imer=1,nmer_size
        do iat=1,natmb
            if(atom_in_fragment_ids(iat)==nmer(imer)) then
                atoms(id) = iat
                id = id +1
            endif
        enddo
    enddo

end subroutine get_nmer_atoms



!------------------------------------------------------------------------------
!> @brief calculates the number of calculations required for a kbody calculation
!!
!! @author Casper Steinmann
!!
!! @param[in] kbody the k-body level
!! @param[out] ncalcs the number of calculations
subroutine calc_n_kbodies(kbody,ncalcs)

    use mb_variables, only : nfrag

    integer, intent(in) :: kbody
    integer, intent(out) :: ncalcs

    integer :: i

    ncalcs = nfrag
    do i=2,kbody
        ncalcs = ncalcs * (nfrag-i+1)
    enddo
    ncalcs = ncalcs / factorial(kbody)

end subroutine calc_n_kbodies

!------------------------------------------------------------------------------
!> @brief calculate the total number of calculations for the entire calculation
!!
!! @author Casper Steinmann
!!
!! @param[out] ncalcs the total number of calculations
subroutine ntotal_calculations(ncalcs)

    use mb_variables, only : nbody

    integer, intent(out) :: ncalcs

    integer :: i
    integer :: tcalcs

    ncalcs = 0
    do i=1,nbody
        call calc_n_kbodies(i,tcalcs)
        ncalcs = ncalcs + tcalcs
    enddo

end subroutine ntotal_calculations

end module mb_utilities
