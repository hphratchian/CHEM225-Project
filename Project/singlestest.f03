      Program CI_Singles
!
!     This program builds the S^2 Matrix using all Sz-conserving
!     singles substitutions. Each matrix element is formed based on the
!     Slater-Condon Rules. ***Expand Description*** 
!            
!
!     The singles excitations are formatted as strings using the number of alpha
!     electrons, the number of beta electrons, and the number of basis
!     functions. Additionally, the output will be compared with the expected
!     number of i-->a substituted determinants.
!
!     A.Zamani
!
!     Last Updated 6/6/19
!
!
!     MQC Utilities 
!     USE Connections
!
      use mqc_general
      use mqc_gaussian
      use mqc_algebra2
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
      integer::iDet=2
      integer::i,j,nDet,nOccAlpha,nOccBeta,nOccAlphaTemp,nOccBetaTemp
      integer::nVirtAlpha,nVirtBeta,nCAPairs,nCAPairsAlpha,nCAPairsBeta,  &
        nSpinFlip
      integer,dimension(:,:,:),allocatable::IStrings
      integer,parameter::maximum_integer_digits = 999 !input restriction
      integer::nAlpha,nBeta,nBasis
      real::SzTemp,SSquareSum
      real,dimension(:,:),allocatable::SSquare,SMatrixAO,CAlpha,CBeta,  &
        Temp_SMatrixOccAB, Temp_SMatrixOccAB_2
      character(maximum_integer_digits)::nA,nB,nBas !dummy input strings
      character(len=512)::matrixFilename
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile
      type(MQC_Variable)::SMatrixAO_mqc,CAlpha_mqc,CBeta_mqc
!
!     Format Statements
!
 1000 Format(1x,'Enter SSquare.')
 1010 Format(3x,'Matrix File: ',A,/)
 1100 Format(1x,'nAlpha=',I4,'  nBeta=',I4,'  nMOs=',I5)
 1110 Format(1x,'nOccAlpha =,'I4,'  nOccBeta =',I4,/,  &
        1x,'nVirtAlpha=',I5,'  nVirtBeta=',I5,/,       &
        1x,'nDet=',I5)
 8000 Format(1x,'Number of ALPHA creation/annihilation pairs: ',I5,/,  &
        1x,'Number of BETA  creation/annihilation pairs: ',I5,/,  &
        1x,'Number of TOTAL creation/annihilation pairs: ',I5,/,  &
        1x,'Number of SPIN FLIP pairs                  : ',I5)
!
!
      write(IOut,1000)
!
!     Open the Gaussian matrix file and load the number of alpha electrons
!     (nAlpha), number of beta electrons (nBeta), and number of MOs (nBasis).
!
      call get_command_argument(1,matrixFilename)
      call GMatrixFile%load(matrixFilename)
      write(IOut,1010) TRIM(matrixFilename)
      nAlpha = GMatrixFile%getVal('nAlpha')
      nBeta  = GMatrixFile%getVal('nBeta')
      nBasis = GMatrixFile%getVal('nBasisUse')
      write(IOut,1100) nAlpha,nBeta,nBasis
!
!     Define O/V orbitals and total number of singles determinants.
!
      nOccAlpha = nAlpha
      nOccBeta = nBeta
      nVirtAlpha = nBasis - nAlpha
      nVirtBeta = nBasis - nBeta
      nDet = ((nOccAlpha*nVirtAlpha)+(nOccBeta*nVirtBeta) + 1)
      write(IOut,1110) nOccAlpha,nOccBeta,nVirtAlpha,nVirtBeta,nDet
!
!     To evaluate matrix elements below, we need atomic orbital overlap matrix
!     elements and MO coefficient matrices. Pull those form the Gaussian matrix
!     file and load them into MQC_Algebra2 objects.
!
      call GMatrixFile%getArray('OVERLAP',mqcVarOut=SMatrixAO_mqc)
      call SMatrixAO_mqc%print(header='Overlap Matrix')
      call GMatrixFile%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=CAlpha_mqc)
      call GMatrixFile%getArray('BETA MO COEFFICIENTS',mqcVarOut=CBeta_mqc)
      call CAlpha_mqc%print(header='MO Coefficients, ALPHA')
      call Cbeta_mqc%print(header='MO Coefficients, BETA')
!
!     For now, we'll take the MQC matrices and copy them into intrinsic fortran
!     arrays.
!
      SMatrixAO = SMatrixAO_mqc
      CAlpha = CAlpha_mqc
      CBeta= CBeta_mqc
!
!     Below is a conditional that requires the number of basis functions to be
!     greater than the number of alpha or beta electrons
!
      if (nBasis .le. max(nAlpha, nBeta)) then
        write(IOut,*) ' The number of basis functions must be greater ', &
          'than the number of alpha or beta electrons.'
        stop
      endif
!
!     Allocate memory.
!
      allocate(IStrings(NBasis,2,NDet)) ! Spin index is 2 long
      allocate(SSquare(NDet,NDet))
!
!     We pass the program's declared variables into the subroutine called
!     stringAlphaBeta.
!
      call stringAlphaBeta(IOut,nOccAlpha,nOccBeta,nBasis,iDet,nDet,  &
        IStrings)
!
!     Try out the test code to build routines that characterize the
!     difference(s) between pairs of strings.
!
      Write(IOut,*)
      Write(IOut,*)' Sending test code strings 3 and 7...'
      call stringsComparison(IOut,NBasis,IStrings(:,:,1),  &
        IStrings(:,:,3),nCAPairs,nCAPairsAlpha,nCAPairsBeta,nSpinFlip)
      write(IOut,8000) nCAPairsAlpha,nCAPairsBeta,nCAPairs,nSpinFlip
!
!     Build the S-Squared matrix.
!
      do i = 1,NDet
        do j = 1,NDet
          call stringsComparison(IOut,NBasis,IStrings(:,:,i),  &
            IStrings(:,:,j),nCAPairs,nCAPairsAlpha,nCAPairsBeta,  &
            nSpinFlip)
          if(nSpinFlip.gt.0.or.nCAPairs.gt.2) then
            SSquare(i,j) = float(0)
          elseIf(nCAPairs.eq.0) then
            if(i.ne.j) call die('Should be diagonal element, &
                 but i.ne.j. ?')
            call stringNOcc(NBasis,IStrings(:,:,i),nOccAlphaTemp,  &
              nOccBetaTemp,SzTemp)
            allocate(Temp_SMatrixOccAB(nOccAlphaTemp,nOccBetaTemp))
            call Form_AlphaBeta_Occ_Overlap(NBasis,nOccAlphaTemp,  &
              nOccBetaTemp,IStrings(:,:,i),SMatrixAO,CAlpha,CBeta,  &
              Temp_SMatrixOccAB,SSquareSum)
            SSquare(i,j) = SzTemp*(SzTemp+float(1)) +  &
              float(nOccBetaTemp) - SSquareSum
!hph            SSquare(i,j) = float(10)
            deallocate(Temp_SMatrixOccAB)
          elseIf(nCAPairs.eq.1) then
!            SSquare(i,j) = float(100)
! SC Rule: 1 Diff
          if(i.eq.j) call die('Should be off-diagonal element, &
                 but i.eq.j. ?') 
!Print these variables before calling subroutine
        write(*,*), 'nOccAlpha: ',nOccAlphaTemp
        write(*,*), 'nOccAlpha: ',nOccBetaTemp
        write(*,*), 'sZTemp: ',SzTemp !Print to see if Sz is correct
          call stringNOcc(NBasis,IStrings(:,:,i),nOccAlphaTemp,  &
            nOccBetaTemp,SzTemp)

!Print these variables after calling subroutine
        write(*,*), 'nOccAlpha: ',nOccAlphaTemp
        write(*,*), 'nOccAlpha: ',nOccBetaTemp
        write(*,*), 'sZTemp: ',SzTemp !Print to see if Sz is correct


        write(*,*), 'write j ',j !check if j is changing
        

!Need 2 temp arrays for each
            allocate(Temp_SMatrixOccAB(nOccAlphaTemp,nOccBetaTemp))
            allocate(Temp_SMatrixOccAB_2(nOccAlphaTemp,nOccBetaTemp))

!Calling twice for i and j: loops over both work ok
          call Form_AlphaBeta_Occ_Overlap(NBasis,nOccAlphaTemp,  &
            nOccBetaTemp,IStrings(:,:,i),SMatrixAO,CAlpha,CBeta,  &
            Temp_SMatrixOccAB,SSquareSum)
!Call again to fill Temp_SMatrixOccAB_2 & SMatrixAlphaBeta_2         

          call Form_AlphaBeta_Occ_Overlap(NBasis,nOccAlphaTemp,  &
            nOccBetaTemp,IStrings(:,:,i),SMatrixAO,CAlpha,CBeta,  &
            Temp_SMatrixOccAB_2,SSquareSum)

!
!***New routine for forming full S^2 matrix and pulling off-diagonals***
!
          call Form_AlphaBeta_Occ_Overlap_Off_Diagonal(nOccAlpha,nOccBeta, &
            Temp_SMatrixOccAB, Temp_SMatrixOccAB_2,SzTemp,nOccBetaTemp,SSquare) 



!Overlap not changing and SzTemp is not reinitialized
            deallocate(Temp_SMatrixOccAB)
!This TEST version of the overlap temp matrix includes beta 'j' strings
            deallocate(Temp_SMatrixOccAB_2)

!SC Rule: 2 Diff
          elseIf(nCAPairs.eq.2) then
            SSquare(i,j) = float(200)
          else
            call die('Confused filling SSquare.')
          endIf
        endDo
      endDo
      write(IOut,*)' Done building the S^2 matrix...'
      call Print_Matrix_Full_Real(IOut,SSquare,NDet,NDet)
      
!
      end program CI_Singles


      subroutine stringAlphaBeta(IOut,NOccAlpha,NOccBeta,NBasis,IDet,  &
        NDet,IStrings)
!
!     This subroutine builds a list of strings corresponding to singles
!     substitutions. The reference string is taken to be in IStrings(1,1,1).
!
      implicit none
      integer,intent(in)::IOut,NBasis,NOccAlpha,NOccBeta,NDet
      integer,intent(inOut)::IDet
      integer,dimension(NBasis,2,*),intent(inOut)::IStrings
      integer::IA,II
!
!     Format statements.
!
 1000 format(/,' |',i3,' > :  ',(i1))
 1010 format(' | ',(i1))
!
!
!     Below, we provide a reference for the alpha and beta string determinants.
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
      write(IOut,*)
      do IDet = 1, NDet
!hph        write(*,1000,advance='no') IDet,IStrings(:,1,IDet)
!hph        write(*,1010,advance='no') IStrings(:,2,IDet)
        call printUnrestrictedString(IOut,NBasis,IDet,IStrings(:,:,IDet))
      endDo
      write(IOut,*)
!
      return
      end subroutine stringAlphaBeta



!
!     Can we use a variable as a format parameter? Truncate alpha digits by
!     NBasis? This code can be modified to push the singles array into a doubles
!     subroutine as a reference. All doubles substitutions can then be printed
!     out from that second subroutine. This requires careful tracking of our
!     declared variables.
!
!

      subroutine printUnrestrictedString(iOut,NBasis,IDet,IString)
!
!     This subroutine is used to print a string defining a spin-unrestricted
!     determinant. If IDet>=0, the determinant number label is printed with the
!     string.
!
      integer,intent(in)::iOut,NBasis,IDet
      integer,dimension(NBasis,2),intent(in)::IString
!
!     Format statements.
!
 1000 format(' |',i3,' > :  ',(i1))
 1010 format(1x,(i1))
 1020 format(' | ',(i1))
!
!     Print the string.
!
      if(IDet.ge.0) then
        write(iOut,1000,advance='no') IDet,IString(:,1)
      else
        write(iOut,1010,advance='no') IString(:,1)
      endIf
      write(iOut,1020,advance='no') IString(:,2)
      write(iOut,*)
!
      return
      end subroutine printUnrestrictedString


      subroutine die(IOut,message)
!
!     This subroutine is called to kill the program after printing our an error
!     message.
!
      integer,intent(in)::IOut
      character(Len=512),intent(in)::message
!
 1000 Format(/,1x,A,/,/,1x,'The program has FAILED!')
!
      write(IOut,*)
      write(IOut,*)' In DIE!'
      write(IOut,*)
      call flush(IOut)
      write(IOut,1000) 'abc'
      call exit(999)
!
      return
      end subroutine die


      subroutine stringsComparison(IOut,NBasis,IString1,IString2,nCAPairs,  &
        nCAPairsAlpha,nCAPairsBeta,nSpinFlip)
!
!     This routine is used to compare two strings and determine how they relate
!     to one another in terms of creation-annihilation pairs, including
!     spin-flip pairs. Creation and annihilation are defined with the operators
!     operating on the determinant represented by IString2 to give the
!     determinant represented by IString1.
!
      implicit none
      integer,intent(in)::IOut,NBasis
      integer,dimension(NBasis,2),intent(in)::IString1,IString2
      integer,intent(out)::nCAPairs,nCAPairsAlpha,nCAPairsBeta,nSpinFlip
      integer::i,nElChange,nCreation,nAnnihilation,nCreationAlpha,  &
        nAnnihilationAlpha,nCreationBeta,nAnnihilationBeta
      integer,dimension(:,:),allocatable::IStringDiff,creationOp,  &
        annihilationOp
      logical::DEBUG=.false.
!
 1000 Format(1x,'Enter test code...',/,3x,'Here are the two strings sent here:')
 1100 Format(1x,'Here is the creation operator...')
 1110 Format(1x,'Here is the annihilation operator...')
 1200 Format(1x,'Number of ALPHA creation/annihilation pairs: ',I5,/,  &
        1x,'Number of BETA  creation/annihilation pairs: ',I5,/,  &
        1x,'Number of TOTAL creation/annihilation pairs: ',I5)
!
!     Report what we're testing in this particular call to the test subroutine
!     and then figure out some critical info about how these two determinants
!     relate to one another.
!
      if(DEBUG) then
        write(IOut,1000)
        call printUnrestrictedString(IOut,NBasis,-1,IString1)
        call printUnrestrictedString(IOut,NBasis,-1,IString2)
      endIf
!
!     Build IStringDiff and print it.
!
      Allocate(IStringDiff(Nbasis,2),creationOp(Nbasis,2),  &
        annihilationOp(Nbasis,2))
      IStringDiff    = IString2-IString1
      creationOp     = 0
      annihilationOp = 0
      do i = 1,NBasis
        select case(IStringDiff(i,1))
        case(-1)
          annihilationOp(i,1) = 1
        case(1)
          creationOp(i,1) = 1
        end select
      endDo
      do i = 1,NBasis
        select case(IStringDiff(i,2))
        case(-1)
          annihilationOp(i,1) = 1
        case(1)
          creationOp(i,1) = 1
        end select
      endDo
      if(DEBUG) then
        write(IOut,1100)
        call printUnrestrictedString(IOut,NBasis,-1,creationOp)
        write(IOut,1110)
        call printUnrestrictedString(IOut,NBasis,-1,annihilationOp)
      endIf
      nCreationAlpha = SUM(creationOp(:,1))
      nCreationBeta = SUM(creationOp(:,2))
      nAnnihilationAlpha = SUM(annihilationOp(:,1))
      nAnnihilationBeta = SUM(annihilationOp(:,2))
      nCreation = nCreationAlpha + nCreationBeta
      nAnnihilation = nAnnihilationAlpha + nAnnihilationBeta
      nSpinFlip = 0
      if(nCreation.ne.nAnnihilation) then
        call die(IOut,'nCreation.ne.nAnnihilation')
        write(Iout,*)' Hrant - I am back 1!'
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
        call die(IOut,'Inconsistent spin flip detected!')
        write(Iout,*)' Hrant - I am back 2!'
      else
        nCAPairsBeta = MAX(nCreationBeta,nAnnihilationBeta)-nSpinFlip
      endIf
      nCAPairs = nCAPairsAlpha+nCAPairsBeta
      if(DEBUG) write(IOut,1200) nCAPairsAlpha,nCAPairsBeta,nCAPairs
!
      return
      end subroutine stringsComparison


      subroutine stringNOcc(NBasis,IString,NOccA,NOccB,SzTemp)
!
!     This routine is used to compute the number of occupied alpha and beta MOs
!     in the determinant defined by IString.
!
      implicit none
      integer,intent(in)::NBasis
      integer,dimension(NBasis,2),intent(in)::IString
      integer,intent(out)::NOccA,NOccB
      real,intent(out)::SzTemp
!
      NOccA = SUM(IString(:,1))
      NOccB = SUM(IString(:,2))
      SzTemp = float(NOccA-NOccB)/float(2)
!
      return
      end subroutine stringNOcc


      Subroutine Form_AlphaBeta_Occ_Overlap(NBasis,NOccA,NOccB,IString,  &
        SMatrixAO,CAlpha,CBeta,SMatrixAlphaBeta,SSquareSum)
!
!     This subroutine forms an alpha-beta overlap between occupied MOs for the
!     determinant defined by IString.
!
!
!     Variable Declarations
!
      implicit none
      integer,intent(in)::NBasis,NOccA,NOccB
      integer,dimension(NBasis,2),intent(in)::IString
      real,dimension(NBasis,NBasis),intent(in)::SMatrixAO,CAlpha,CBeta
      real,dimension(NOccA,NOccB),intent(out)::SMatrixAlphaBeta
      real,intent(out)::SSquareSum
      integer::i,ii
      real,dimension(NBasis,NOccA)::TempCAlphaOcc
      real,dimension(NBasis,NOccB)::TempCBetaOcc
!
!     The starting point is building temporary MO coefficient matrices that run
!     only over the occupied MOs of the determinant defined by IString.
!
      ii = 0
      do i = 1,NBasis
        if(IString(i,1).eq.1) then
          ii = ii+1
          TempCAlphaOcc(:,ii) = CAlpha(:,i)
        endIf
      endDo
      ii = 0
      do i = 1,NBasis
        if(IString(i,2).eq.1) then
          ii = ii+1
          TempCBetaOcc(:,ii) = CBeta(:,i)
        endIf
      endDo
!
!     Using the temp CAlpha and CBeta arrays, form SMatrixAlphaBeta using
!     MatMul.
!
      SMatrixAlphaBeta = MatMul(  &
        MatMul(Transpose(TempCAlphaOcc),SMatrixAO),TempCBetaOcc)
        print*,'MO OVERLAP MATRIX: ' 
        call Print_Matrix_Full_Real(6,SMatrixAlphaBeta,NOccA,NOccB)

      SSquareSum = dot_product(  &
        Reshape(SMatrixAlphaBeta,[NOccA*NOccB]),  &
        Reshape(SMatrixAlphaBeta,[NOccA*NOccB]))

        
!
      return
      end subroutine Form_AlphaBeta_Occ_Overlap


      Subroutine Print_Matrix_Full_Real(IOut,AMat,M,N)
!
!     This subroutine prints a real matrix that is fully dimension - i.e.,
!     not stored in packed form. AMat is the matrix, which is dimensioned
!     (M,N).
!
!
!     Variable Declarations
!
      implicit none
      integer,intent(in)::IOut,M,N
      real,dimension(M,N),intent(in)::AMat
!
!     Local variables
      integer,parameter::NColumns=5
      integer::i,j,IFirst,ILast
!
 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)
!
      Do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOut,2000) (i,i=IFirst,ILast)
        Do i = 1,M
          write(IOut,2010) i,(AMat(i,j),j=IFirst,ILast)
        endDo
      endDo
!
      return
      end subroutine Print_Matrix_Full_Real



!call Form_AlphaBeta_Occ_Overlap_Off_Diagonal(nOccAlpha,nOccBeta, &
!        Temp_SMatrixOccAB,Temp_SMatrixOccAB_2,SzTemp,nOccBetaTemp,SSquare) 

      subroutine Form_AlphaBeta_Occ_Overlap_Off_Diagonal(NOccA,NOccB, &
        SMatrixAlphaBeta,SMatrixAlphaBeta_2,SzTemp,nOccBetaTemp,SSquare)
       
      implicit none
      integer :: i, j
      integer :: NOccA, NOccB
      integer :: nOccBetaTemp
      real, intent(in) ::  SzTemp
      real,dimension(:,:), allocatable :: SSquare
      real, dimension(NOccA,NOccB) :: OverlapSum
      real, dimension(NOccA,NOccB), intent(in) :: SMatrixAlphaBeta
      real, dimension(NOccA,NOccB), intent(in) :: SMatrixAlphaBeta_2
      ! ^ this comes from Form_AlphaBeta_Occ_Overlap , the 2nd to last term in
      ! the call statement: Temp_SMatrixOccAB ^
      !real, intent(out) SSquare


      OverlapSum = MatMul(  &
      (SMatrixAlphaBeta),  &
      (SMatrixAlphaBeta_2))

      write(*,*), 'Overlap Sum: ', OverlapSum

      do i=1, NoccA
        do j=1, NoccB
          if (i.eq.j) then
            print*,' i != j must be the case '
            elseif (i.ne.j) then
        SSquare(i,j) = SzTemp*(SzTemp+float(1)) +  &
            float(nOccBetaTemp) - OverlapSum(i,j)
          endif
        enddo
      enddo
      

      end subroutine Form_AlphaBeta_Occ_Overlap_Off_Diagonal
