#if defined (PRG_DALTON)

!------------------------------------------------------------------------------
!> @brief Executes only the HERMIT module of DALTON for an nmer
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer_charge total charge on the nmer
!! @param[in] nat the number of atoms to write
!! @param[in] atoms the atoms to write
      subroutine execute_hermit(nmer_charge, nat, atoms)

      use mb_variables, only : work, 
     *                          idalton_work_size,
     *                          dalton_work_control

      integer nmer_charge
      integer nat
      integer atoms
      dimension atoms(nat)

      ! write the input file
      call write_input(nmer_charge, nat, atoms)

      call reset_for_hermit
      call exeher(work(0),
     *            idalton_work_size,
     *            dalton_work_control)

      end subroutine execute_hermit
!------------------------------------------------------------------------------
!> @brief Executes the SCF of an nmer trough dalton
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer_charge total charge on the nmer
!! @param[in] nat the number of atoms to write
!! @param[in] atoms the atoms to write
      subroutine execute_dalton(nmer_charge, nat, atoms)
 
      use mb_variables, only : work, 
     *                          idalton_work_size,
     *                          dalton_work_control

#include "implicit.h"
#include "maxorb.h"
#include "gnrinf.h"

      integer nmer_charge
      integer nat
      integer atoms
      dimension atoms(nat)

      ! write the input file
      call write_input(nmer_charge, nat, atoms)

      ! execute DALTON
      call reset_for_hermit
      call reset_for_sirius
      call exeher(work(0),
     *            idalton_work_size,
     *            dalton_work_control)

      call reset_for_sirius
      if(rnsiri) call exesir(work(0),
     *                       idalton_work_size,
     *                       dalton_work_control)

      if(rnabac) then
          call abainp('**PROPE',work(1),idalton_work_size)
          call exeaba(work(0),
     *                idalton_work_size,
     *                dalton_work_control)
      endif

      end subroutine
!------------------------------------------------------------------------------
!> @brief Executes the code needed for the 2e- integrals of dalton
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer_charge total charge on the nmer
!! @param[in] nat the number of atoms to write
!! @param[in] atoms the atoms to write
      subroutine execute_integral_dalton(nmer_charge, nat, atoms)
 
      use mb_variables, only : work, 
     *                          idalton_work_size,
     *                          dalton_work_control

#include "implicit.h"
#include "maxorb.h"
#include "gnrinf.h"

      integer nmer_charge
      integer nat
      integer atoms
      dimension atoms(nat)

      ! write the input file
      call write_input(nmer_charge, nat, atoms)

      ! execute DALTON
      call reset_for_hermit
      call reset_for_sirius
      call exeher(work(0),
     *            idalton_work_size,
     *            dalton_work_control)

      if(rnsiri) call exesir(work(0),
     *                       idalton_work_size,
     *                       dalton_work_control)

      end subroutine

!------------------------------------------------------------------------------
!> @brief Writes the current atoms to a DALTON .mol file
!!
!! @author Casper Steinmann
!!
!! @param[in] nmer_charge total charge on the nmer
!! @param[in] nat the number of atoms to write
!! @param[in] atoms the atoms to write
      subroutine write_input(nmer_charge, nat, atoms)

      use mb_utilities, only : openfile
      use mb_variables, only : nuclear_coordinates,
     *                          nuclear_charges,
     *                          mb_basis_name

      implicit none

      integer nmer_charge
      integer nat
      integer atoms
      dimension atoms(nat)


      integer :: lu
      integer :: iat
      integer :: zat, k
      character(len=99) :: cl
      character(len=8) :: cl1

      character(len=2), dimension(112) :: elements

      elements = (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ',
     *              'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ',
     *              'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr',
     *              'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
     *              'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',
     *              'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
     *              'In', 'Sn', 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba',
     *              'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
     *              'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
     *              'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     *              'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',
     *              'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm',
     *              'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf',
     *              'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn' /)

      write(cl,*) nat

      write(cl1,'(A7,I1)') 'Charge=',nmer_charge

      call openfile('MOLECULE.INP', lu, 'old', 'formatted')
      rewind(lu)
      write(lu, '(a)') 'BASIS'
      write(lu, '(a)') mb_basis_name
      write(lu, '(a)') ' '
      write(lu, '(a)') ' '
      write(lu, '(5a)') 'AtomTypes=',trim(cl),
     *                    ' NoSymmetry ',trim(cl1),' Angstrom'
C     we write one atomtype per coordinate
      do iat=1,nat
         zat = int(nuclear_charges(atoms(iat)))
         write(cl,*) zat
         write(lu,'(a,a,a)') 'Charge=',trim(cl),'.0 Atoms=1'
         write(lu,'(a,f25.10,2f16.10)') trim(elements(zat)),
     *                 (nuclear_coordinates(k,atoms(iat)),k=1,3)
      enddo
      close(lu)

      end subroutine write_input

!------------------------------------------------------------------------------
C> @brief resets variables for HERMIT
C!
C! @author Casper Steinmann
      subroutine reset_for_hermit

#include "implicit.h"
#include "maxorb.h"
#include "gnrinf.h"
!#include "siripc.h"
!#include "huckel.h"
!include "cbisol.h"

      ! flag that we have not read in the current geometry
      rdinpc = .false.

      end subroutine reset_for_hermit

!------------------------------------------------------------------------------
C> @brief resets variables for SIRIUS
C!
C! @author Casper Steinmann
C!
C! @note this must be called before EXEHER for some
C!       unknown reason
      subroutine reset_for_sirius

#include "implicit.h"
#include "maxorb.h"
#include "siripc.h"
#include "huckel.h"

      ! mark that sirius should re-read options and settings
      inpprc = .false.

      ! always flag a new geometry
      newgeo = .true.

      ! force that we do huckel again
      dohuckel = .true.
      end subroutine

!> @brief evaluates the two-electron contribution of the embedding potential through dalton
!! @param[in] nbas number of atomic orbitals.
!! @param[in] nocc number of occupied orbitals of the entire complex
!! @param[in] norb number of orbitals in total
!! @param[in,out] dalwrk temporary DALTON work
!! @param[out] iconv Whether or not convergence should be forced. If iconv=1
!!             then we are doing a potential calculation, otherwize 0 and exit
!! @note inspired by code from the Polarizable Embedding library
      subroutine dalton_twoints(nbas, dalwrk, iconv)

      use mb_precision
      use mb_utilities, only : unfold_matrix
      use mb_variables, only : p_kmer_dens, p_nmer_fock, work,
     *                          idalton_work_size
  
      external :: sirfck
  
      integer, intent(in) :: nbas
      real(dp), intent(inout), dimension(:) :: dalwrk
      integer, intent(out) :: iconv
      ! ----
      integer :: fbas, kmer_nao, kden, nden
      integer :: isymdm, ifctyp
      real(dp), allocatable, dimension(:,:) :: full_densmat
      real(dp), allocatable, dimension(:,:) :: full_fockmat
      real(dp), allocatable, dimension(:,:) :: full_kmer_densmat
      integer i,j,l
  
      ! if there is no K-mer density assigned, we are here during
      ! a regular SCF. Just return.
      iconv = 0
      if (associated(p_kmer_dens)) iconv = 1
      if (iconv == 0) return
  
      ! otherwise, calculate the electronic repulsion through
      ! V = D^K (mu nu | lambda sigma). D^K is the density of the
      ! kmer which we have to unpack into full storage
      kden = size(p_kmer_dens)
      kmer_nao = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * kden) - 1.0d0))
      allocate(full_kmer_densmat(kmer_nao, kmer_nao))
      call unfold_matrix(p_kmer_dens, full_kmer_densmat)
  
      ! in two_int we built a super-fragment like kmer+nmer.
      ! In DALTON, we can obtain V from the following procedure:
  
      ! fill full_densmat matrix with D^K in first block
      ! <----- nbas ---->
      ! <-fbas ->
      ! +-------+-------+
      ! |   K   |       |
      ! |  D    |   0   |
      ! |       |       |
      ! +-------+-------+
      ! |       |       |
      ! |   0   |   0   |
      ! |       |       |
      ! +-------+-------+
      !nden = size(p_mb_fock)
      !nmer_nao = int(0.5d0 * (sqrt(1.0d0 + 8.0d0 * nden) - 1.0d0))
  
      ! insert the unpacked kmer density into the full density matrix
      allocate(full_densmat(nbas,nbas))
      full_densmat = 0.0d0
      full_densmat(1:kmer_nao,1:kmer_nao) = full_kmer_densmat
      deallocate(full_kmer_densmat)
  
      ! 2e- electron contribution is
      ! +-------+-------+
      ! |       |       |
      ! |   Z   |   X   |
      ! |       |       |
      ! +-------+-------+
      ! |       |       |
      ! |   X   |   Z   |
      ! |       |       |
      ! +-------+-------+
      ! sirfck will calculate the product of D^K and the above matrix
  
      ! IFCTYP = +/-XY
      !   X indicates symmetry about diagonal
      !     X = 0 No symmetry
      !     X = 1 Symmetric
      !     X = 2 Anti-symmetric
      !   Y indicates contributions
      !     Y = 0 no contribution !
      !     Y = 1 Coulomb
      !     Y = 2 Exchange
      !     Y = 3 Coulomb + Exchange
      !   + sign: alpha + beta matrix (singlet)
      !   - sign: alpha - beta matrix (triplet)
      ! sirfck(fckmat, denmat, ?, isymdm, ifctyp, direct, work, nwrk)
      isymdm = 1
      ifctyp = 11
      allocate(full_fockmat(nbas,nbas))
      full_fockmat = 0.0d0
      call sirfck(full_fockmat, full_densmat, 1, isymdm, ifctyp,
     *            .true., work(1), idalton_work_size)
  
      deallocate(full_densmat)
  
      ! now add the contributions from the lower part of the full fock matrix
      ! to the current fragment fock matrix because it contains V
      ! +-------+-------+
      ! |       |       |
      ! |   ?   |   ?   |
      ! |       |       |
      ! +-------+-------+
      ! |       |       |
      ! |   ?   |   V   |
      ! |       |       |
      ! +-------+-------+
      l = 1
      p_nmer_fock = 0.0d0
      do j = kmer_nao + 1, nbas
          do i = kmer_nao + 1, j
              p_nmer_fock(l) = full_fockmat(i,j)
              l = l +1
          enddo
      enddo
  
      deallocate(full_fockmat)

      end subroutine

#endif
