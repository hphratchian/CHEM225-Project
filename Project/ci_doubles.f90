  Program CI_Doubles
  !
  ! This program constructs all Sz-conserving doubles substitutions.
  !
  ! The doubles excitations are formatted as strings using the number of
  ! alpha electrons, the number of beta electrons, and the number of basis
  ! functions. Additionally, the output will be compared with the expected
  ! number of ij-->ab substituted determinants.
  !
  ! A.Zamani
  !
  ! Last Updated 5/7/19

  implicit none
  integer :: iDet = 2 !Start after reference
  integer :: nDet, nOccBeta, nOccAlpha, nDetA, nDetB
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

  !The following combinatorial formulas generate the number of possible
  !subsitutions, represented by nDet, for both the alpha and beta strings.

  nDetA = (choose(nOccAlpha, 2) * choose(nVirtAlpha, 2))
  nDetB = (choose(nOccBeta, 2) * choose(nVirtBeta, 2))
  
  nDet = nDetA + nDetB + 1

  !
  !Below is a conditional that requires the number of basis
  !functions to be greater than the number of alpha or beta
  !electrons.
  !
  if (nBasis .le. max(nAlpha, nBeta)) then
    write(*,*) ' The number of basis functions must be greater ', &
               'than the number of alpha or beta electrons.'
    stop
  endif
  !
  !We then pass the program's declared variables into the subroutine
  !called stringAlphaBeta.
  !
  call stringAlphaBetaDoubles(nAlpha, nBeta, nBasis, nOccAlpha, nOccBeta, &
  iDet, nDet, nVirtAlpha, nVirtBeta)

  write(*,*)' Total Doubles: ',nDet !The number of determinants minus the reference.
  write(*,*)' Total Determinants: ',iDet !The number of determinants including the reference.

  !n choose k functions
  contains
       
  function factorial (n) result (res) !factorial
 
    implicit none
    integer, intent (in) :: n
    integer :: res
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
  end program CI_Doubles
        
  

  subroutine stringAlphaBetaDoubles(NAlpha, NBeta, NBasis, NOccAlpha, &
  NOccBeta, IDet, NDet,NVirtA,NVirtB)

  implicit none
  integer :: NAlpha, NBeta, NBasis, II, IDet, NVirtA, NVirtB
  integer :: NOccBeta, NOccAlpha, IA, NDet, i, a
  integer, dimension(:,:,:), allocatable :: IStrings

  NOccAlpha = NAlpha !redefine number of alpha occupied MOs
  NOccBeta = NBeta !redefine number of beta occupied MOs

  allocate(IStrings(NBasis,2,NDet)) ! Spin index is 2 long
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
  !
  !alpha
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
  !beta 
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
  !
  write(*,*)
  do IDet = 1, NDet
    write(*,1001,advance='no') IStrings(:,:,IDet)
  enddo
  write(*,*)
  end subroutine stringAlphaBetaDoubles
  !
  !
