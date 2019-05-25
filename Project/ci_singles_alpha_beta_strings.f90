      Program CI_Singles
!
!     This program constructs all Sz-conserving singles substitutions.
!
!     The singles excitations are formatted as strings using the number of alpha
!     electrons, the number of beta electrons, and the number of basis
!     functions. Additionally, the output will be compared with the expected
!     number of i-->a substituted determinants.
!
!     A.Zamani
!
!     Last Updated 4/13/19

      implicit none
      integer,parameter::IOut=6
      integer::iDet = 2
      integer::nDet,nOccBeta,nOccAlpha
      integer::nVirtAlpha,nVirtBeta,nCAPairs,nCAPairsAlpha,nCAPairsBeta,  &
        nSpinFlip
      integer,dimension(:,:,:),allocatable::IStrings
      integer,parameter::maximum_integer_digits = 999 !input restriction
      integer::nAlpha,nBeta,nBasis      
      character(maximum_integer_digits)::nA,nB,nBas !dummy input strings
!
 8000 Format(1x,'Number of ALPHA creation/annihilation pairs: ',I5,/,  &
        1x,'Number of BETA  creation/annihilation pairs: ',I5,/,  &
        1x,'Number of TOTAL creation/annihilation pairs: ',I5,/,  &
        1x,'Number of SPIN FLIP pairs                  : ',I5)
!
!
!     The following reads command line input into the variables: nAlpha, nBeta,
!     and nBasis.
! 
!        
      call get_command_argument(1, nA) !separate each by a space in CLI
      call get_command_argument(2, nB)
      call get_command_argument(3, nBas)
      read (nA, *) nAlpha
      read (nB, *) nBeta
      read (nBas, *) nBasis
!
!     Define O/V orbitals and total number of singles determinants.
!
      nVirtAlpha = nBasis - nAlpha
      nVirtBeta = nBasis - nBeta
      nOccAlpha = nAlpha
      nOccBeta = nBeta
      nDet = ((nOccAlpha*nVirtAlpha)+(nOccBeta*nVirtBeta) + 1)

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
      call hphTest(IOut,NBasis,IStrings(:,:,3),IStrings(:,:,7),nCAPairs,  &
        nCAPairsAlpha,nCAPairsBeta,nSpinFlip)
      write(IOut,8000) nCAPairsAlpha,nCAPairsBeta,nCAPairs,nSpinFlip
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
      write(IOut,1000) TRIM(message)
      STOP
!
      return
      end subroutine die


      subroutine hphTest(IOut,NBasis,IString1,IString2,nCAPairs,  &
        nCAPairsAlpha,nCAPairsBeta,nSpinFlip)
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
      if(nCreation.ne.nAnnihilation) call die(IOut,  &
        'hphTest: nCreation.ne.nAnnihilation')
      if(nCreationAlpha.eq.nAnnihilationAlpha) then
        nCAPairsAlpha = nCreationAlpha
      else
        nSpinFlip = ABS(nCreationAlpha-nAnnihilationAlpha)
        nCAPairsAlpha = MAX(nCreationAlpha,nAnnihilationAlpha)-nSpinFlip
      endIf
      if(nCreationBeta.eq.nAnnihilationBeta) then
        nCAPairsBeta = nCreationBeta
      elseIf(nSpinFlip.ne.(ABS(nCreationBeta-nAnnihilationBeta))) then
        call die(IOut,'hphTest: Inconsistent spin flip detected!')
      else
        nCAPairsBeta = MAX(nCreationBeta,nAnnihilationBeta)-nSpinFlip
      endIf
      nCAPairs = nCAPairsAlpha+nCAPairsBeta
      if(DEBUG) write(IOut,1200) nCAPairsAlpha,nCAPairsBeta,nCAPairs
!
      return
      end subroutine hphTest



