  Program CI_Singles_Doubles
  !
  ! This program constructs all Sz-conserving singles and doubles substitutions.
  !
  ! The singles and doubles excitations are formatted as strings using the number of
  ! alpha electrons, the number of beta electrons, and the number of basis
  ! functions. Additionally, the output will be compared with the expected
  ! number of i-->a and ij-->ab substituted determinants.
  !
  ! A.Zamani
  !
  ! Last Updated 5/8/19

  implicit none
  integer :: iDet = 2
  integer :: nDetD, nDetA, nDetB, nDetS
  integer :: nOccBeta, nOccAlpha 
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


  !Read in inputs
  read (nA, *) nAlpha
  read (nB, *) nBeta
  read (nBas, *) nBasis

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


  !Define O/V orbitals and total number of singles determinants.

  nVirtAlpha = nBasis - nAlpha
  nVirtBeta = nBasis - nBeta
  nOccAlpha = nAlpha
  nOccBeta = nBeta

  nDetS = ((nOccAlpha*nVirtAlpha)+(nOccBeta*nVirtBeta) + 1)
  write(*,*)' nDetS: ',nDetS !print test nDetS     

  !The following combinatorial formulas generate the number of possible
  !subsitutions, represented by nDet, for both the alpha and beta strings.

  nDetA = (choose(nOccAlpha, 2) * choose(nVirtAlpha, 2))
  nDetB = (choose(nOccBeta, 2) * choose(nVirtBeta, 2))

  nDetD = nDetA + nDetB + 1
  write(*,*)' nDetD: ',nDetD !print test nDetD

  ! 
  !We pass the program's declared variables into the subroutine
  !called stringAlphaBetaSingles and stringAlphaBetaDoubles.
  !

  call stringAlphaBetaSingles(nAlpha, nBeta, nBasis, nOccAlpha, nOccBeta, &
  iDet, nDetS)

  !The number of determinants including the reference.
  write(*,*)' Total i-->a Determinants + Reference: ',nDetS 

  call stringAlphaBetaDoubles(nAlpha, nBeta, nBasis, nOccAlpha, nOccBeta, &
  iDet, nDetD, nVirtAlpha, nVirtBeta)

  !The number of determinants including the reference.
  write(*,*)' Total ij-->ab Determinants + Reference: ',nDetD 
  
  !
  !n choose k functions for computing the number of determinants for doubles.
  !
  
  contains

  function factorial (n) result (res) !factorial

    implicit none
    integer, intent (in) :: n
    integer*8 :: res !turn into double
    integer :: i

    res = product ((/(i, i = 1, n)/))

  end function factorial

  function choose (n, k) result (res) !n choose k

    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: k
    integer :: res

    res = factorial (n) / (factorial (k) * factorial (n - k))

  end function choose
  !
  !
  !
  end program CI_Singles_Doubles
  !
  !
  !

  !
  !This subroutine generates all the singles substitutions.
  !

  subroutine stringAlphaBetaSingles(NAlpha, NBeta, NBasis, NOccAlpha, &
  NOccBeta, IDet, NDetS)

  implicit none
  integer :: NAlpha, NBeta, NBasis, II, IDet
  integer :: NOccBeta, NOccAlpha, IA, NDetS
  integer, dimension(:,:,:), allocatable :: IStrings

  NOccAlpha = NAlpha !redefine number of alpha occupied MOs
  NOccBeta = NBeta !redefine number of beta occupied MOs

  allocate(IStrings(NBasis,2,NDetS)) ! Spin index is 2 long

  !
  !Below, we provide a reference for the alpha and beta string
  !determinants. Our loops for changing the orbital occupation
  !number will then start after IDet = 1.
  !

  !Give reference for alpha string
  do II = 1, NOccAlpha
    do IA = NOccAlpha + 1, NBasis
       IStrings(II,1,1) = 1
       IStrings(IA,1,1) = 0
    enddo
  enddo

  !Give reference for beta string
  do II = 1, NOccBeta
    do IA = NOccBeta + 1, NBasis
       IStrings(II,2,1) = 1
       IStrings(IA,2,1) = 0
    enddo
  enddo
  !
  !IDet starts at 2 because we skip the reference state
  !
  !Now perform excitations.
  !
  do II = 1, NOccAlpha !starting from orbital 1-2
    do IA = NOccAlpha + 1, NBasis !starting from orbital 3-4
      IStrings(:,:,IDet) = IStrings(:,:,1) !open array indices
      IStrings(II,1,IDet) = 0
      IStrings(IA,1,IDet) = 1

      IDet = IDet + 1 !We loop over states starting at IDet = 1

    enddo
  enddo

  do II = 1, NOccBeta !starting from orbital 1-2
    do IA = NOccBeta + 1, NBasis !starting from orbital 3-4
      IStrings(:,:,IDet) = IStrings(:,:,1) !open array indices
      IStrings(II,2,IDet) = 0
      IStrings(IA,2,IDet) = 1

      IDet = IDet + 1 !We loop over states starting at IDet = 1

    enddo
  enddo
  !
  !Format the output by placing every array containing the alpha and
  !beta strings on a new line with a space between them. This format should be
  !recoded for generality.   
  !
 1000 format(/,1x,(i1),(i1)) !New line after each 2*NBasis-long string.
  !
  write(*,*)
  write(*,*)' Singles:'
  do IDet = 1, NDetS
    write(*,1000,advance='no') IStrings(:,:,IDet)
  enddo
  write(*,*)
  end subroutine stringAlphaBetaSingles
  !
  !
  !This subroutine generates all the doubles substitutions.
  subroutine stringAlphaBetaDoubles(NAlpha, NBeta, NBasis, NOccAlpha, &
  NOccBeta, IDet, NDetD,NVirtA,NVirtB)

  implicit none
  integer :: NAlpha, NBeta, NBasis, II, IDet, NVirtA, NVirtB
  integer :: NOccBeta, NOccAlpha, IA, NDetD, i, a
  integer, dimension(:,:,:), allocatable :: IStrings

  NOccAlpha = NAlpha !redefine number of alpha occupied MOs
  NOccBeta = NBeta !redefine number of beta occupied MOs
  IDet = 2 !Reinitialize after Singles run
  allocate(IStrings(NBasis,2,NDetD)) ! Spin index is 2 long
  !
  !Below, we provide a reference for the alpha and beta string
  !determinants. Our loops for changing the orbital occupation
  !number will then start after IDet = 1.
  !
  !Give reference for alpha string
  do II = 1, NOccAlpha
    do IA = NOccAlpha + 1, NBasis
       IStrings(II,1,1) = 1
       IStrings(IA,1,1) = 0
    enddo
  enddo

  !Give reference for beta string
  do II = 1, NOccBeta
    do IA = NOccBeta + 1, NBasis
       IStrings(II,2,1) = 1
       IStrings(IA,2,1) = 0
    enddo
  enddo
  !
  !Now perform excitations.
  !
  !alpha doubles looping
  do II = 1, NOccAlpha-1
    do IA = NOccAlpha + 1, NBasis - 1
      do i = 1, NOccAlpha - II
        do a = 1, NBasis - IA
          IStrings(:,:,IDet) = IStrings (:,:,1)
          IStrings(II + i,1,IDet) = 0
          IStrings(IA,1,IDet) = 1
          IStrings(II,1,IDet) = 0
          IStrings(IA + a,1,IDet) = 1

          IDet = IDet + 1

        enddo
      enddo
    enddo
  enddo
  !
  !
  !beta doubles looping
  do II = 1, NOccBeta-1
    do IA = NOccBeta + 1, NBasis - 1
      do i = 1, NOccBeta - II
        do a = 1, NBasis - IA
          IStrings(:,:,IDet) = IStrings (:,:,1)
          IStrings(II + i,2,IDet) = 0
          IStrings(IA,2,IDet) = 1
          IStrings(II,2,IDet) = 0
          IStrings(IA + a,2,IDet) = 1

          IDet = IDet + 1

        enddo
      enddo
    enddo
  enddo

  !
  !Format the output by placing every array containing the alpha and
  !beta strings on a new line with a space between them. This format should be
  !recoded for generality.
  !

 1001 format(/,1x,(i1),(i1)) !New line after each 2*NBasis-long string.
  
  write(*,*)
  write(*,*)' Doubles:'
  do IDet = 1, NDetD
    write(*,1001,advance='no') IStrings(:,:,IDet)
  enddo
  write(*,*)
  end subroutine stringAlphaBetaDoubles
  !
