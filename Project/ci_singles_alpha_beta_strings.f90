  Program CI_Singles
  !
  ! This program constructs all Sz-conserving singles substitutions.
  !
  ! The singles excitations are formatted as strings using the number of
  ! alpha electrons, the number of beta electrons, and the number of basis
  ! functions. Additionally, the output will be compared with the expected
  ! number of i-->a substituted determinants.
  !
  ! A.Zamani
  !
  ! Last Updated 4/13/19

  implicit none
  integer :: iDet = 2
  integer :: nDet, nOccBeta, nOccAlpha
  integer :: nVirtAlpha, nVirtBeta
  integer, dimension(:,:,:), allocatable :: iStrings
  
  !The following reads command line input into the variables: nAlpha,
  !nBeta, and nBasis.
 
  integer, parameter :: maximum_integer_digits = 999 !input restriction
  integer :: nAlpha, nBeta, nBasis      
  character(maximum_integer_digits) :: nA, nB, nBas !dummy input strings
        
  call get_command_argument(1, nA) !separate each by a space in CLI
  call get_command_argument(2, nB)
  call get_command_argument(3, nBas)

  read (nA, *) nAlpha
  read (nB, *) nBeta
  read (nBas, *) nBasis

  !Define O/V orbitals and total number of singles determinants.

  nVirtAlpha = nBasis - nAlpha
  nVirtBeta = nBasis - nBeta
  nOccAlpha = nAlpha
  nOccBeta = nBeta

  nDet = ((nOccAlpha*nVirtAlpha)+(nOccBeta*nVirtBeta) + 1)

  !
  !Below is a conditional that requires the number of basis
  !functions to be greater than the number of alpha or beta
  !electrons
  !
  if (nBasis .le. max(nAlpha, nBeta)) then
    write(*,*) ' The number of basis functions must be greater ', &
               'than the number of alpha or beta electrons.'
    stop
  endif
 
!
!     Allocate memory.
!
      allocate(IStrings(NBasis,2,NDet)) ! Spin index is 2 long
!
!     We pass the program's declared variables into the subroutine called
!     stringAlphaBeta.
!
      call stringAlphaBeta(nAlpha,nBeta,nBasis,nOccAlpha,nOccBeta,  &
        iDet,nDet)
!
      end program CI_Singles


      subroutine stringAlphaBeta(NAlpha,NBeta,NBasis,NOccAlpha,  &
        NOccBeta,IDet,NDet,Strings)
!
!     This subroutine builds a list of strings corresponding to singles
!     substitutions. The reference string is taken to be in IStrings(1,1,1).
!
      implicit none
      integer :: NAlpha, NBeta, NBasis, II, IDet
      integer :: NOccBeta, NOccAlpha, IA, NDet
      integer, dimension(:,:,:), allocatable :: IStrings
!
!     Format statements.
!
 1000 format(/,1x,(i1),(i1)) !New line after each 2*NBasis-long string.
!
!
!     Redefine variables for the number of alpha and beta occupied MOs. (I'm not
!     sure this needs to be done here.?)
!
      NOccAlpha = NAlpha !redefine number of alpha occupied MOs
      NOccBeta = NBeta !redefine number of beta occupied MOs
!
!     Below, we provide a reference for the alpha and beta string
!     determinants. Our loops for changing the orbital occupation number will
!     then start after IDet = 1.
!
!     Give reference for alpha string
      do II = 1, NOccAlpha
        do IA = NOccAlpha + 1, NBasis
          IStrings(II,1,1) = 1
          IStrings(IA,1,1) = 0
        endDo
      endDo
!
!     Give reference for beta string
      do II = 1, NOccBeta
        do IA = NOccBeta + 1, NBasis
          IStrings(II,2,1) = 1
          IStrings(IA,2,1) = 0
        endDo
      endDo
!
!     IDet starts at 2 because we skip the reference state
!
      do II = 1, NOccAlpha !starting from orbital 1-2
        do IA = NOccAlpha + 1, NBasis !starting from orbital 3-4
          IStrings(:,:,IDet) = IStrings(:,:,1) !open array indices
          IStrings(II,1,IDet) = 0
          IStrings(IA,1,IDet) = 1
          IDet = IDet + 1 !We loop over states starting at IDet = 1
        endDo
      endDo
      do II = 1, NOccBeta !starting from orbital 1-2
        do IA = NOccBeta + 1, NBasis !starting from orbital 3-4
          IStrings(:,:,IDet) = IStrings(:,:,1) !open array indices
          IStrings(II,2,IDet) = 0
          IStrings(IA,2,IDet) = 1
          IDet = IDet + 1 !We loop over states starting at IDet = 1
        enddo
      enddo
!
!     Format the output by placing every array containing the alpha and beta
!     strings on a new line with a space between them. This format should be
!     recoded for generality.   
!
      write(*,*)
      do IDet = 1, NDet
        write(*,1000,advance='no') IStrings(:,:,IDet)
      endDo
      write(*,*)
      end subroutine stringAlphaBeta
!
!     Can we use a variable as a format parameter? Truncate alpha digits by
!     NBasis? This code can be modified to push the singles array into a doubles
!     subroutine as a reference. All doubles substitutions can then be printed
!     out from that second subroutine. This requires careful tracking of our
!     declared variables.
!
!
