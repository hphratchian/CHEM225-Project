      program cis_total_spin_squared_matrix

! This program constructs all Sz-conserving singles substitutions and
! uses them to define the Hilbert space for constructing the S^2 matrix.
!
! The data structures for defining the determinants in an occupation
! number representation are formulated as alpha and beta strings.
!
! The reference and the singles excitations are formatted as strings using 
! the number of alpha electrons, the number of beta electrons, and the number
! of basis functions. These values are found in an unformatted .mat file 
! produced by Gaussian.
!
! The string determinants are used to build the individual elements of
! the S^2 matrix. The total spin-squared expectation value, <S^2>, shall
! provide some indication of the level of spin contamination in the
! molecular system of interest.
!
! The ultimate goal of this project is to develop a post-analysis tool to
! generate spin-adapted configuration state functions and to discern what spin
! states are mixing into the wavefunction describing a molecule.
!
! A.Zamani
!
! Last Updated 6/21/19



! loading MQC tools
! USE connections for accessing modules in MQCPack     
      use mqc_general
      use mqc_gaussian
      use mqc_algebra2
      use iso_fortran_env

!      
! Global Variable Declarations
!

      implicit none      
! input variable declarations
      integer :: nAlpha, nBeta, nBasis
! variable declarations for mqc interface
      character(len = 512) :: matrixFilename
      type(mqc_gaussian_unformatted_matrix_file) :: GMatrixFile
      type(MQC_Variable) :: overlapAO_MQC, cAlpha_MQC, cBeta_MQC 
! redefining MQC matrices as Fortran arrays
      real, dimension(:,:), allocatable :: SMatrixAO, cAlpha, cBeta
! defining O/V orbital count and corresponding temporary variables      
      integer :: nOccAlpha, nOccBeta, nVirtAlpha, nVirtBeta
      integer :: nOccAlphaTemp, nOccBetaTemp      
! defining # of singles determinants and member-index of that set 
      integer :: nDet, iDet = 2
! defining string determinants
      integer, dimension(:,:,:), allocatable :: iStrings
! defining dummy indices
      integer :: i, j
! defining # creation/annihilation operator pairs and spin flips
      integer :: nCAPairs, nCAPairsAlpha, nCAPairsBeta, nSpinFlip
! define other variables used to determine <S^2>
      real :: SzTemp, SSquareSum
      real, dimension(:,:), allocatable :: SSquared
! the temporary matrix below stores the contents of the alpha/beta MO Overlap
      real, dimension(:,:), allocatable :: Temp_SMatrixOccAB
      real, dimension(:,:), allocatable :: Temp_SMatrixOccAB_2
! define format parameter
      integer, parameter :: iOut = 6

! a new definition of the overlap sum term in <S^2>
! remember to allocate and deallocate for other SC rule blocks
      real :: OverlapSum
! passed out of stringsComparison for use in general_contraction      
      integer, dimension(:,:), allocatable :: COP, AOP 




! START: general format statements
 1000 Format(1x,'Enter SSquared.')
 1010 Format(3x,'Matrix File: ',A,/)
 1100 Format(1x,'nAlpha=',I4,'  nBeta=',I4,'  nMOs=',I5)
 1110 Format(1x,'nOccAlpha =,'I4,'  nOccBeta =',I4,/,  &
        1x,'nVirtAlpha=',I5,'  nVirtBeta=',I5,/,       &
        1x,'nDet=',I5)
 8000 Format(1x,'Number of ALPHA creation/annihilation pairs: ',I5,/,  &
        1x,'Number of BETA  creation/annihilation pairs: ',I5,/,  &
        1x,'Number of TOTAL creation/annihilation pairs: ',I5,/,  &
        1x,'Number of SPIN FLIP pairs                  : ',I5)
! END: general format statements

! statement informing the user that we're entering this program
      write(iOut,1000)
     
! open the Gaussian matrix file and load nAlpha, nBeta, and nBasis
! print these variables as formatted statements
      call get_command_argument(1,matrixFilename)
      call GMatrixFile%load(matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      nAlpha = GMatrixFile%getVal('nAlpha')
      nBeta  = GMatrixFile%getVal('nBeta')
      nBasis = GMatrixFile%getVal('nBasisUse')
      write(IOut,1100) nAlpha,nBeta,nBasis

! define expressions
! define O/V orbitals and total number of singles determinants
! print these variables as formatted statements
      nOccAlpha = nAlpha
      nOccBeta = nBeta
      nVirtAlpha = nBasis - nAlpha
      nVirtBeta = nBasis - nBeta
      nDet = ((nOccAlpha * nVirtAlpha)+(nOccBeta * nVirtBeta) + 1)
      write(IOut,1110) nOccAlpha, nOccBeta, nVirtAlpha, nVirtBeta, nDet

! pull AO overlap matrix and alpha/beta MO coefficients from the Gaussian 
! matrix file and load them into mqc_algebra2 objects
      call GMatrixFile%getArray('OVERLAP',mqcVarOut=overlapAO_MQC)
      call overlapAO_MQC%print(header='Overlap Matrix')
      call GMatrixFile%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=cAlpha_MQC)
      call GMatrixFile%getArray('BETA MO COEFFICIENTS',mqcVarOut=cBeta_MQC)
      call cAlpha_MQC%print(header='MO Coefficients, ALPHA')
      call cBeta_MQC%print(header='MO Coefficients, BETA')
     
! copy MQC matrices into intrinsic fortran arrays.
      SMatrixAO = overlapAO_MQC
      cAlpha = cAlpha_MQC
      cBeta = cBeta_MQC


!error print
      SMatrixAO(1,3) = float(360420)
      print*,'printing SMatrixAO from mqc  ',SMatrixAO
      print*,'printing SMatrixAO full REAL '
      call Print_Matrix_Full_Real(IOut,SMatrixAO,nBasis,nBasis)
!error print
      cAlpha(1,3) = float(360420)
      print*,'printing cAlpha from mqc  ',cAlpha
      print*,'printing cAlpha full REAL '
      call Print_Matrix_Full_Real(IOut,cAlpha,nBasis,nBasis)
!error print
      cBeta(1,3) = float(360420)
      print*,'printing cBeta from mqc  ',cBeta
      print*,'printing cBeta full REAL '
      call Print_Matrix_Full_Real(IOut,cBeta,nBasis,nBasis)
!is it different? how is it packed? do we need to transpose?




! START: conditional that requires the # of basis fxns to be
! greater than the # of alpha or beta electrons.
      if (nBasis .le. max(nAlpha, nBeta)) then
        write(IOut,*) ' The number of basis functions must be greater ', &
          'than the number of alpha or beta electrons.'
        stop
      endif 
! END: conditional


! allocate memory for iStrings: spin index is 2 (alpha or beta)
      allocate(iStrings(nBasis,2,nDet))
! allocate memory for SSquared: nDet by nDet in dimensionality
      allocate(SSquared(nDet,nDet))

! pass the declared variables into the subroutine called stringAlphaBeta
      call stringAlphaBeta(iOut,nOccAlpha,nOccBeta,nBasis,iDet,nDet, &
        iStrings)

! START: testing the subroutine, stringsComparison, used to
! characterize the  differences in orbital occupancy between two strings. 
! for example, let's compare strings 3 and 7 and print a formatted
! statement for each value. 
      Write(IOut,*)
      Write(IOut,*)' Sending test code strings 3 and 7...'

! error diagnostics 
! Something's wrong in the subroutine stringsComparison, gives seg-fault
! allocate them in the program??? doesnt seem to work.
      Allocate(COP(nBasis,2),AOP(nBasis,2))
! it's not entering stringsComparison...

      call stringsComparison(IOut,NBasis,IStrings(:,:,1),  &
        IStrings(:,:,3),nCAPairs,nCAPairsAlpha,nCAPairsBeta, &
        nSpinFlip, COP, AOP)
      
      write(IOut,8000) nCAPairsAlpha,nCAPairsBeta,nCAPairs,nSpinFlip
! END: testing


!
! we shall now build the S^2 matrix
!

! remember to allocate and deallocate the temp matrices!
!error print 
      print*,'is the loop working?'

      do i = 1,NDet
        do j = 1,NDet
          call stringsComparison(IOut,NBasis,IStrings(:,:,i),  &
            IStrings(:,:,j),nCAPairs,nCAPairsAlpha,nCAPairsBeta,  &
            nSpinFlip, COP, AOP)
!error print
        print*,'exiting stringsComparison '
!error print
        print*,'printing nCAPairs ',nCAPairs
        print*,'printing nSpinFlip ',nSpinFlip
        print*,'printing i ',i
        print*,'printing j ',j        



          if(nSpinFlip.gt.0.or.nCAPairs.gt.2) then
            SSquared(i,j) = float(0)
! SC Rule 1:
          elseIf(nCAPairs.eq.0) then
            if(i.ne.j) call die('Should be diagonal element, &
                but i.ne.j.?')
            call stringNOcc(NBasis,IStrings(:,:,i),nOccAlphaTemp,  &
              nOccBetaTemp,SzTemp)
! error print
        print*,'istrings from stringnocc: ',IStrings(:,:,i) 

            allocate(Temp_SMatrixOccAB(nOccAlphaTemp,nOccBetaTemp))
!forgot to allocate 2nd temp matrix
        allocate(Temp_SMatrixOccAB_2(nOccAlphaTemp,nOccBetaTemp))

            call Form_AlphaBeta_Occ_Overlap(NBasis,nOccAlphaTemp,  &
              nOccBetaTemp,IStrings(:,:,i),SMatrixAO,CAlpha,CBeta,  &
              Temp_SMatrixOccAB, Temp_SMatrixOccAB_2,SSquareSum)
!error print
        print*,'after calling form alphabetaoccoverlap...'


            SSquared(i,j) = SzTemp*(SzTemp+float(1)) +  &
              float(nOccBetaTemp) - SSquareSum
! placeholder          SSquared(i,j) = float(10)
            deallocate(Temp_SMatrixOccAB)
!deallocate 2nd temp matrix
            deallocate(Temp_SMatrixOccAB_2)



!SC Rule 2:
          elseIf(nCAPairs.eq.1) then
            if(i.eq.j) call die('Should be off-diagonal element, &
                 but i.eq.j. ?')
! placeholder            SSquared(i,j) = float(100)
! AZ coding....

!error print
        print*,'Before entering my contraction routine... '           

!remember to allocate and deallocate the temp MO overlaps...

        allocate(Temp_SMatrixOccAB(nOccAlphaTemp,nOccBetaTemp))
!forgot to allocate 2nd temp matrix
        allocate(Temp_SMatrixOccAB_2(nOccAlphaTemp,nOccBetaTemp))

!try swapping i and j, s^2 should be hermitian... 
!swapping them changes the <s^2> value, check indexing in
!general_contraction indexing is wrong, do by hand and check the
!indexing for SAO and CMOs pulled from MQC. Remember, the SAO 
!matrix is initially constructed in a row-wise vector. Fix it!!!


            call general_contraction(COP, AOP, Temp_SMatrixOccAB, &
                Temp_SMatrixOccAB_2, OverlapSum, nBasis, & 
                IStrings(:,:,i), IStrings(:,:,j), &
                nOccAlpha, nOccBeta) 


!error print
        print*,'printing i after calling gen_contract ',i       
        print*,'printing j after calling gen_contract ',j

        SSquared(i,j) = SzTemp*(SzTemp+float(1)) +  &
              float(nOccBetaTemp) - OverlapSum

!error print
        print*,'printing ssquared(i,j)',SSquared(i,j)



            deallocate(Temp_SMatrixOccAB)
!deallocate 2nd temp matrix
            deallocate(Temp_SMatrixOccAB_2)

! SC Rule 3:
          elseIf(nCAPairs.eq.2) then
            if(i.eq.j) call die('Should be off-diagonal element, &
                 but i.eq.j. ?')
! placeholder            SSquared(i,j) = float(200)
          else
            call die('Confused filling SSquared.')
          endIf
        endDo
      endDo
      write(IOut,*)' Done building the S^2 matrix...'

!
! print S^2 matrix
!
      call Print_Matrix_Full_Real(IOut,SSquared,NDet,NDet)
!
!
!

      end program cis_total_spin_squared_matrix
!
!
! START: subroutines 
!
!

! this subroutine builds a list of alpha and beta strings for every
! singles substitution. The reference determinant/string is taken to be
! iStrings(1,:,1). iStrings is allocated to retain the following
! dimensionality: iStrings(NBasis,alpha(1) or beta(2),iDet).   

      subroutine stringAlphaBeta(IOUT,NOCCALPHA,NOCCBETA,NBASIS,IDET, &
        NDET,ISTRINGS)

      implicit none
      ! variables to pass into the subroutine
      integer, intent(in) :: IOUT, NBASIS, NOCCALPHA, NOCCBETA, NDET
      ! passing in and out of the subroutine
      integer, intent(inout) :: IDET
      ! passing ISTRINGS in and out of the subroutine
      integer, dimension(NBASIS,2,*), intent(inout) :: ISTRINGS
      ! dummy indices 
      integer :: IA, II

      ! below, we provide a reference for the alpha and beta strings
      
      ! alpha reference
      do II = 1, NOCCALPHA
        do IA = NOCCALPHA + 1, NBASIS
          ISTRINGS(II,1,1) = 1
          ISTRINGS(IA,1,1) = 0
        endDo
      endDo

      ! beta reference
      do II = 1, NOCCBETA
        do IA = NOCCBETA + 1, NBASIS
          ISTRINGS(II,2,1) = 1
          ISTRINGS(IA,2,1) = 0
        endDo
      endDo

      ! start singles excitations
      ! IDET starts at 2 because we're skipping the reference

      do II = 1, NOCCALPHA !starting from orbital 1-2
        do IA = NOCCALPHA + 1, NBASIS !starting from orbital 3-4
          ISTRINGS(:,:,IDET) = ISTRINGS(:,:,1) !open array indices
          ISTRINGS(II,1,IDET) = 0
          IStrings(IA,1,IDET) = 1
          IDET = IDET + 1 !We loop over states starting at IDet = 1
        endDo
      endDo
      do II = 1, NOCCBETA !starting from orbital 1-2
        do IA = NOCCBETA + 1, NBASIS !starting from orbital 3-4
          ISTRINGS(:,:,IDET) = ISTRINGS(:,:,1) !open array indices
          ISTRINGS(II,2,IDET) = 0
          IStrings(IA,2,IDET) = 1
          IDET = IDET + 1 !We loop over states starting at IDet = 1
        enddo
      enddo

      ! format and print the singles determinants
      write(IOUT,*)
      do IDET = 1, NDET
      ! using the formatting subroutine 
        call printUnrestrictedString(IOUT,NBASIS,IDET, &
          ISTRINGS(:,:,IDET))
      endDo
      write(IOUT,*)

      return
      end subroutine stringAlphaBeta

! 
!
      subroutine printUnrestrictedString(IOUT, NBASIS,IDET, &
        ISTRING)
!
! this subroutine is used to print a string defining a
! spin-unrestricted determinant. If IDet>=0, the determinant number label 
! is printed with the string.
!
      integer, intent(in) :: IOUT, NBASIS, IDET
      integer, dimension(NBASIS,2), intent(in) :: ISTRING

!
! format statements.
!

! kets formatted and printed in occupation number representation 
 1000 format(' |',i3,' > :  ',(i1))
 1010 format(1x,(i1))
! demarcate alpha and beta strings with a pipe |
 1020 format(' | ',(i1))

!
! print the strings.
!
      if(IDET.ge.0) then
        write(IOUT,1000,advance='no') IDET ,ISTRING(:,1)
      else
        write(IOUT,1010,advance='no') ISTRING(:,1)
      endIf
      write(IOUT,1020,advance='no') ISTRING(:,2)
      write(IOUT,*)
!
      return
      end subroutine printUnrestrictedString

! killswitch (HPH)

      subroutine die(IOUT,message)
!
! this subroutine is called to kill the program after printing our
! an error message.
!
      integer,intent(in) :: IOUT
      character(Len=512),intent(in) :: message
!
 1000 Format(/,1x,A,/,/,1x,'The program has FAILED!')
!
      write(IOUT,*)
      write(IOUT,*)' In DIE!'
      write(IOUT,*)
      call flush(IOUT)
      write(IOUT,1000) 'abc'
      call exit(999)
!
      return
      end subroutine die
!
!


      subroutine stringsComparison(IOUT,NBASIS,ISTRING1,ISTRING2, &
        nCAPairs,nCAPairsAlpha,nCAPairsBeta,nSpinFlip, creationOp, &
        annihilationOp)

!
! This routine is used to compare two strings and determine how they
! relate to one another in terms of creation-annihilation pairs, 
! including spin-flip pairs. Creation and annihilation are defined with 
! the operators operating on the determinant represented by ISTRING2 to 
! give the determinant represented by ISTRING1. (HPH)
!
      implicit none
      integer, intent(in)::IOUT,NBASIS
      integer, dimension(NBASIS,2),intent(in) :: ISTRING1,ISTRING2
      integer, intent(out) :: nCAPairs,nCAPairsAlpha,nCAPairsBeta, &
        nSpinFlip
      integer :: i,nElChange,nCreation,nAnnihilation,nCreationAlpha, &
        nAnnihilationAlpha,nCreationBeta,nAnnihilationBeta
      integer, dimension(:,:), allocatable :: IStringDiff
      !changed out to inout below
      !changed allocatable (:,:)
      integer,dimension(NBASIS,2),intent(inout) :: creationOp
      integer,dimension(NBASIS,2),intent(inout) :: annihilationOp
      logical::DEBUG=.false.
!
 1000 Format(1x,'Enter test code...',/,3x,'Here are the two strings & 
        sent here:' ) 
 1100 Format(1x,'Here is the creation operator...')
 1110 Format(1x,'Here is the annihilation operator...')
 1200 Format(1x,'Number of ALPHA creation/annihilation pairs: ',I5,/,  &
        1x,'Number of BETA  creation/annihilation pairs: ',I5,/,  &
        1x,'Number of TOTAL creation/annihilation pairs: ',I5)
!
! report what we're testing in this particular call to the test
! subroutine and then figure out some critical info about how these two
! determinants relate to one another.
!

!error print
      print*,'Looking for an error in stringsComparison '

      if(DEBUG) then
        write(IOUT,1000)
        call printUnrestrictedString(IOUT,NBASIS,-1,ISTRING1)
        call printUnrestrictedString(IOUT,NBASIS,-1,ISTRING2)
      endIf
!
! build IStringDiff and print it.
!

!error print
      print*,'Looking for an error in stringsComparison2 '


!dont need to allocate AOP and COP
!      Allocate(IStringDiff(Nbasis,2),creationOp(NBASIS,2),  &
!       annihilationOp(NBASIS,2))
!try allocating just istringdiff lol
      Allocate(IStringDiff(NBASIS,2))
!error print
      print*,'Looking for an error in stringsComparison3 '

!error print
      print*,'printing istring2 ',ISTRING2(1,1)
!error print
      print*,'printing istring1 ',ISTRING1(1,2)
!HA it looks like ISTRING 1 and 2 are both the reference determinant
!I need to pass all possible determinants(singles here) to get the cases

      IStringDiff = ISTRING2-ISTRING1
!error print
      print*,'printing istringdiff=istring2-istring1 ',IStringDiff

!error print
      print*,'Looking for an error in stringsComparison4 '

      annihilationOp = 0
      print*,'annihilationOp is the problem!!!! '
      creationOp = 0
!error print
      print*,'Looking for an error in stringsComparison5x4 '
     ! annihilationOp = 0
!error print
      print*,'Looking for an error in stringsComparison6x4 '
!error print
      print*,'printing istringdiffbefore ', IStringDiff
!error print, how many nbasis?
      write(*,*),' THIS IS NBASIS: ',NBASIS


      do i = 1,NBASIS
!error print
      print*,'printing istringdiffinsideDO ', IStringDiff
      print*,'printing istringdiffinsideDO ', IStringDiff(i,1)
      print*,'printing i',i

        select case(IStringDiff(i,1))
!error print
      print*,'printing istringdiffafter ', IStringDiff
!error print
      print*,'Looking for an error in stringsComparison7x4 '



        case(-1)
!error print
      print*,'looking for error in case(-1) '


          annihilationOp(i,1) = 1
        case(1)
!error print
      print*,'looking for error in case(1) '


          creationOp(i,1) = 1
        end select
      endDo
!beta loops
      do i = 1,NBASIS
!error print
        print*,'entering beta comparisons... '
!it was saving to the alpha index, it's fixed.     
        select case(IStringDiff(i,2))
        case(-1)
          annihilationOp(i,2) = 1
        case(1)
          creationOp(i,2) = 1
        end select
      endDo
!error print
        print*,'exiting beta comparisons... '
! THE ERROR IS DOWN HERE VVVVVVV
      if(DEBUG) then
        write(IOUT,1100)
        call printUnrestrictedString(IOUT,NBASIS,-1,creationOp)
        write(IOUT,1110)
        call printUnrestrictedString(IOUT,NBASIS,-1,annihilationOp)
      endIf
      nCreationAlpha = SUM(creationOp(:,1))
      nCreationBeta = SUM(creationOp(:,2))
      nAnnihilationAlpha = SUM(annihilationOp(:,1))
      nAnnihilationBeta = SUM(annihilationOp(:,2))
      nCreation = nCreationAlpha + nCreationBeta
      nAnnihilation = nAnnihilationAlpha + nAnnihilationBeta
      nSpinFlip = 0
      if(nCreation.ne.nAnnihilation) then
        call die(IOUT,'nCreation.ne.nAnnihilation')
        write(IOUT,*)' Hrant - I am back 1!'
      endIf
      if(nCreationAlpha.eq.nAnnihilationAlpha) then
        nCAPairsAlpha = nCreationAlpha
      else
        nSpinFlip = ABS(nCreationAlpha-nAnnihilationAlpha)
        nCAPairsAlpha = MAX(nCreationAlpha,nAnnihilationAlpha)-nSpinFlip
      endIf
      if(nCreationBeta.eq.nAnnihilationBeta) then
        nCAPairsBeta = nCreationBeta
      elseIf(nSpinFlip.ne.(ABS(nCreationBeta-nAnnihilationBeta))) then
        call die(IOUT,'Inconsistent spin flip detected!')
        write(IOUT,*)' Hrant - I am back 2!'
      else
        nCAPairsBeta = MAX(nCreationBeta,nAnnihilationBeta)-nSpinFlip
      endIf
      nCAPairs = nCAPairsAlpha + nCAPairsBeta
      if(DEBUG) write(IOUT,1200) nCAPairsAlpha, nCAPairsBeta, nCAPairs
!
      return
      end subroutine stringsComparison
!

      subroutine stringNOcc(NBASIS,ISTRING,NOCCA,NOCCB,SzTemp)
!
! This routine is used to compute the number of occupied alpha and 
! beta MOs in the determinant defined by IString.
!
      implicit none
      integer,intent(in)::NBASIS
      integer,dimension(NBASIS,2),intent(in)::ISTRING
      integer,intent(out)::NOCCA,NOCCB
      real,intent(out)::SzTemp
!
      NOCCA = SUM(ISTRING(:,1))
      NOCCB = SUM(ISTRING(:,2))
      SzTemp = float(NOCCA-NOCCB)/float(2)
!
      return
      end subroutine stringNOcc


! The subroutine below constructs the overlap sum for the diagonals


      Subroutine Form_AlphaBeta_Occ_Overlap(NBASIS,NOCCA,NOCCB, &
        ISTRING,SMatrixAO,CALPHA,CBETA,SMatrixAlphaBeta, &
        SMatrixAlphaBeta_2, SSquareSum)
!
! This subroutine forms an alpha-beta overlap between occupied MOs for the
! determinant defined by IString.
!
!
! Variable Declarations
!
      implicit none
      integer, intent(in) :: NBASIS, NOCCA, NOCCB
      integer, dimension(NBASIS,2), intent(in) :: ISTRING
      real, dimension(NBASIS,NBASIS), intent(in) :: SMatrixAO, CALPHA, &
        CBETA
      real, dimension(NOCCA,NOCCB), intent(out) :: SMatrixAlphaBeta, &
        SMatrixAlphaBeta_2
      real, intent(out) :: SSquareSum
      integer :: i,ii
      real, dimension(NBASIS,NOCCA) :: TempCAlphaOcc
      real, dimension(NBasis,NOCCB) :: TempCBetaOcc
!
! The starting point is building temporary MO coefficient matrices that run
! only over the occupied MOs of the determinant defined by IString.
!
      ii = 0
      do i = 1,NBASIS
        if(ISTRING(i,1).eq.1) then
          ii = ii+1
          TempCAlphaOcc(:,ii) = CALPHA(:,i)
        endIf
      endDo
      ii = 0
      do i = 1,NBASIS
        if(ISTRING(i,2).eq.1) then
          ii = ii+1
          TempCBetaOcc(:,ii) = CBETA(:,i)
        endIf
      endDo
!
! Using the temp CALPHA and CBETA arrays, form SMatrixAlphaBeta using
! MatMul.
!
      SMatrixAlphaBeta = MatMul(  &
        MatMul(Transpose(TempCAlphaOcc),SMatrixAO),TempCBetaOcc)
      
      ! Second version of matrix to pass out 
      SMatrixAlphaBeta_2 = MatMul(  &
        MatMul(Transpose(TempCAlphaOcc),SMatrixAO),TempCBetaOcc)


      SSquareSum = dot_product(  &
        Reshape(SMatrixAlphaBeta,[NOCCA*NOCCB]),  &
        Reshape(SMatrixAlphaBeta,[NOCCA*NOCCB]))
!
      return
      end subroutine Form_AlphaBeta_Occ_Overlap

! A subroutine to print a formatted matrix (HPH)

      subroutine Print_Matrix_Full_Real(IOUT,AMat,M,N)
!
! This subroutine prints a real matrix that is fully dimension - i.e.,
! not stored in packed form. AMat is the matrix, which is dimensioned
! (M,N).

!
! Variable Declarations
!
      implicit none
      integer,intent(in)::IOUT,M,N
      real,dimension(M,N),intent(in)::AMat
!
! Local variables
      integer,parameter::NColumns=5
      integer::i,j,IFirst,ILast
!
 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)
!
      Do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOUT,2000) (i,i=IFirst,ILast)
        Do i = 1,M
          write(IOUT,2010) i,(AMat(i,j),j=IFirst,ILast)
        endDo
      endDo
!
      return
      end subroutine Print_Matrix_Full_Real

! END: subroutines 
!
!
      subroutine general_contraction(cOp, aOp, overlapMO_1, &
        overlapMO_2, overlapSum, NBASIS, string1, string2, &
        NOCCA, NOCCB) 

      implicit none
!forgot to add dimensions for overlapmo: cant allocate when fixed
      integer, intent(in) :: NOCCA, NOCCB
      integer :: numCreate, numAnnihilate
      integer, intent(in) :: NBASIS
      integer, dimension(NBASIS,2), intent(in) :: cOp
      integer, dimension(NBASIS,2), intent(in) :: aOp
      real, intent(out) :: overlapSum
      integer :: r, s     
!change from :, : to nocca, noccb
      real, dimension(NOCCA,NOCCB), intent(in) :: overlapMO_1, & 
        overlapMO_2
      
! I gotta pass the SMatrixAlphaBeta into here...Need 2 versions to
! replace overlapMO_1 and overlapMO_2 
      
      integer, dimension(NBASIS,2),intent(in) :: string1, string2
      integer :: positionAlpha_1, positionAlpha_2, positionBeta_1, &
        positionBeta_2
 
!error print
        print*, 'before defining numcreate, numan, and overlapsum'

!error print
        print*, 'cOp passed in ',cOp
!error print 
        print*, 'aOp passed in ',aOp

      numCreate = sum(cOp)

!error print
        print*,'numcreate',numcreate
        
      numAnnihilate = sum(aOp)

!error print
        print*,'numAnnihilate',numAnnihilate

!error print 
        print*,'overlapsumbefore ',overlapSum
       
      overlapSum = 0.0

!error print
        print*, 'print string1 before loop ',string1
        print*,'print string2 before loop ',string2

!try swapping r and s (original r=1, NBASIS...) [6/24/19]

      do r = 1, NBASIS
        positionAlpha_1 = string1(r,1)
        positionBeta_1 = string1(r,2)


!error print
        print*, 'inside first do loop with r'
!error print
        print*, 'print string1 in loop r ',string1
        print*,'print string2 in loop r ',string2


        if (positionAlpha_1.eq.1 .or. positionBeta_1.eq.1) then
          do s = 1, NBASIS
!error print
        print*, 'inside first do loop with s'

            positionAlpha_2 = string2(s,1)
            positionBeta_2 = string2(s,2)
            if (positionAlpha_2.eq.1 .or. positionBeta_2.eq.1) then
              overlapSum = overlapMO_1(numCreate,r) * &
                 overlapMO_2(numAnnihilate,s) + overlapSum 
            endif
          enddo
        endif
      enddo

!error print
        print*, 'at the end of do loops, exit gen_contract'


      end subroutine general_contraction 



!CONTINUE RECODING HERE



