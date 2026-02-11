SUBROUTINE Create_ntID_Arrays()
!
! Create or update (if MutProtPos isn't zero) nucleotide id arrays.
! Note that INTEGER(KIND=4) can have only 9 digits!  INTEGER(KIND=8) can
! hold 17 digits...
!
! A ntID array holds an n-digit integer in place of the sequence.
! A=-1 T=1 C=-3 G=3
!
! ACGTACGTACGTACGT with a RepLen = 8 would be shown as
! ........        ntID_Rep(1) = 12341234
!  ........       ntID_Rep(2) = 23412341
!   ........      ntID_Rep(3) = 34123412
!    ........     ntID_Rep(4) = 41234123
!
! and so on...

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,a1,a2,b1,b2,fin,gc_count,at_count

! Forward strand digit: A=1, C=2, G=3, T=4 (indexed by NUMseq value -3..3)
  INTEGER, SAVE :: nt_fwd(-3:3)
  DATA nt_fwd /2, 0, 1, 0, 4, 0, 3/

! Reverse complement digit: A->T=1, C->G=2, G->C=3, T->A=4 (wait, needs to match
! original mapping: A(-1)->4, C(-3)->3, G(3)->2, T(1)->1)
  INTEGER, SAVE :: nt_rc(-3:3)
  DATA nt_rc /3, 0, 4, 0, 1, 0, 2/

! Precomputed powers of 10 (max index 8 for KIND=4 integers)
  INTEGER, SAVE :: pow10(0:8)
  DATA pow10 /1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000/

  IF (TEST2) PRINT *,'Create_ntID_Arrays'

  IF (MutProtPos.eq.0) THEN
    a1=1
    a2=DNAlen-MPTip+1
    b1=1
    b2=DNAlen-RepLen+1
  ELSE
    a1=(MAX(1,(MutNtPos(1)-MPTip)))
    a2=(MIN((DNAlen-MPTip+1),(MutNtPos(MutNtNum)+1)))
    b1=(MAX(1,(MutNtPos(1)-RepLen)))
    b2=(MIN((DNAlen-RepLen+1),(MutNtPos(MutNtNum)+1)))
  END IF

! update misprime arrays

  fin=MPTip-1
  DO i=a1,a2
    CurrDNA%ntID_Tip(i)=0
    CurrDNA%ntID_TipRC(i)=0
    DO j=0,fin
      CurrDNA%ntID_Tip(i)=CurrDNA%ntID_Tip(i)+nt_fwd(CurrDNA%NUMseq(i+fin-j))*pow10(j)
      CurrDNA%ntID_TipRC(i)=CurrDNA%ntID_TipRC(i)+nt_rc(CurrDNA%NUMseq(i+j))*pow10(j)
    END DO
  END DO

! update repeat arrays

  fin=RepLen-1
  DO i=b1,b2
    CurrDNA%ntID_Rep(i)=0
    CurrDNA%ntID_RepRC(i)=0
    DO j=0,fin
      CurrDNA%ntID_Rep(i)=CurrDNA%ntID_Rep(i)+nt_fwd(CurrDNA%NUMseq(i+fin-j))*pow10(j)
      CurrDNA%ntID_RepRC(i)=CurrDNA%ntID_RepRC(i)+nt_rc(CurrDNA%NUMseq(i+j))*pow10(j)
    END DO
  END DO

! update GC and AT arrays
! When MutProtPos==0 (full scan), use sliding window for O(n) instead of O(n*k)

  fin=RepLen-1
  IF (MutProtPos.eq.0 .and. b2.gt.b1) THEN
    ! Compute initial window for position b1
    gc_count=0
    at_count=0
    DO j=b1, b1+fin
      IF (ABS(CurrDNA%NUMseq(j)).eq.1) gc_count=gc_count+1
      IF (ABS(CurrDNA%NUMseq(j)).eq.3) at_count=at_count+1
    END DO
    CurrDNA%ntID_GC(b1)=gc_count
    CurrDNA%ntID_AT(b1)=at_count
    ! Slide the window
    DO i=b1+1, b2
      ! Remove element leaving the window (position i-1)
      IF (ABS(CurrDNA%NUMseq(i-1)).eq.1) gc_count=gc_count-1
      IF (ABS(CurrDNA%NUMseq(i-1)).eq.3) at_count=at_count-1
      ! Add element entering the window (position i+fin)
      IF (ABS(CurrDNA%NUMseq(i+fin)).eq.1) gc_count=gc_count+1
      IF (ABS(CurrDNA%NUMseq(i+fin)).eq.3) at_count=at_count+1
      CurrDNA%ntID_GC(i)=gc_count
      CurrDNA%ntID_AT(i)=at_count
    END DO
  ELSE
    ! Partial update or single position: use direct counting
    DO i=b1,b2
      CurrDNA%ntID_GC(i)=0
      CurrDNA%ntID_AT(i)=0
      DO j=i, i+fin
        IF (ABS(CurrDNA%NUMseq(j)).eq.1) CurrDNA%ntID_GC(i)=CurrDNA%ntID_GC(i)+1
        IF (ABS(CurrDNA%NUMseq(j)).eq.3) CurrDNA%ntID_AT(i)=CurrDNA%ntID_AT(i)+1
      END DO
    END DO
  END IF

END SUBROUTINE Create_ntID_Arrays
SUBROUTINE Sort_Misprime_Arrays()

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,n

  IF (TEST2) PRINT *,"Sort_Misprime_Arrays" !TEST2

! Heapsort misprime pairs by M1

  n = CurrDNA%MN
  IF (n.le.1) RETURN

  DO i=n/2,1,-1
    CALL SiftDown(i,n)
  END DO

  DO i=n,2,-1
    CALL IntSwap(CurrDNA%M1(1),CurrDNA%M1(i))
    CALL IntSwap(CurrDNA%M2(1),CurrDNA%M2(i))
    CALL IntSwap(CurrDNA%MX(1),CurrDNA%MX(i))
    CALL SiftDown(1,i-1)
  END DO

CONTAINS

  SUBROUTINE SiftDown(root,bottom)
    INTEGER, INTENT(IN) :: root,bottom
    INTEGER :: parent,child

    parent = root
    child = 2*parent
    DO WHILE (child.le.bottom)
      IF (child.lt.bottom) THEN
        IF (CurrDNA%M1(child).lt.CurrDNA%M1(child+1)) child = child+1
      END IF
      IF (CurrDNA%M1(parent).lt.CurrDNA%M1(child)) THEN
        CALL IntSwap(CurrDNA%M1(parent),CurrDNA%M1(child))
        CALL IntSwap(CurrDNA%M2(parent),CurrDNA%M2(child))
        CALL IntSwap(CurrDNA%MX(parent),CurrDNA%MX(child))
        parent = child
        child = 2*parent
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE SiftDown

END SUBROUTINE Sort_Misprime_Arrays
SUBROUTINE Sort_Repeat_Arrays

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,n

  IF (TEST2) PRINT *,"Sort_Repeat_Arrays" !TEST2

! Rearrange repeat pairs so RS1 <= RS2

  DO i=1,CurrDNA%RN
    IF (CurrDNA%RS1(i).gt.CurrDNA%RS2(i)) THEN
      CALL IntSwap(CurrDNA%RS1(i),CurrDNA%RS2(i))
    END IF
  END DO

! Heapsort repeat pairs by RS1

  n = CurrDNA%RN
  IF (n.le.1) RETURN

  DO i=n/2,1,-1
    CALL SiftDown(i,n)
  END DO

  DO i=n,2,-1
    CALL IntSwap(CurrDNA%RS1(1),CurrDNA%RS1(i))
    CALL IntSwap(CurrDNA%RS2(1),CurrDNA%RS2(i))
    CALL IntSwap(CurrDNA%RLn(1),CurrDNA%RLn(i))
    CALL IntSwap(CurrDNA%RX(1),CurrDNA%RX(i))
    CALL SiftDown(1,i-1)
  END DO

CONTAINS

  SUBROUTINE SiftDown(root,bottom)
    INTEGER, INTENT(IN) :: root,bottom
    INTEGER :: parent,child

    parent = root
    child = 2*parent
    DO WHILE (child.le.bottom)
      IF (child.lt.bottom) THEN
        IF (CurrDNA%RS1(child).lt.CurrDNA%RS1(child+1)) child = child+1
      END IF
      IF (CurrDNA%RS1(parent).lt.CurrDNA%RS1(child)) THEN
        CALL IntSwap(CurrDNA%RS1(parent),CurrDNA%RS1(child))
        CALL IntSwap(CurrDNA%RS2(parent),CurrDNA%RS2(child))
        CALL IntSwap(CurrDNA%RLn(parent),CurrDNA%RLn(child))
        CALL IntSwap(CurrDNA%RX(parent),CurrDNA%RX(child))
        parent = child
        child = 2*parent
      ELSE
        EXIT
      END IF
    END DO
  END SUBROUTINE SiftDown

END SUBROUTINE Sort_Repeat_Arrays
SUBROUTINE Translate_Protein
!
! Translate the mutatable protein residues into DNA sequence

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,k,p,d,x
  REAL :: rand
  LOGICAL :: no_codons
  INTEGER,EXTERNAL :: NT2Int
  CHARACTER(LEN=3) :: tempCodonSeq

  IF (TEST0) PRINT *,"Translate_Protein" !TEST0

! Reset MutProtPos

  MutProtPos=0

  IF (.not.SequenceTranslated) THEN  ! avoid the first time
    main: DO i=1,mutPROTnum
      p=mutPROT2prot(i)
      d=prot2nt(p)

! Choose the codon randomly unless the codon is not allowed

      k = 1
      IF (CodonRandom) THEN
        CALL RANDOM_NUMBER(rand)
        k=(INT(rand*(AAT(prot2aa(p))%NumOfActiveCodons)))+1
      END IF
      tempCodonSeq=CFT(AAT(prot2aa(p))%Codon(k))%Seq
  
! Create the codon and insert it into the DNA sequence.

! If the chain is reversed, put in reverse complement

      IF (ChainReverse(prot2chain(p))) CALL RevComplStr(tempCodonSeq)
      CurrDNA%DNAseq(d-1:d+1)=tempCodonSeq
  
! Fill prot2cod array

      CurrDNA%prot2cod(p) = AAT(prot2aa(p))%Codon(k)
      CurrDNA%nt2cod(d) = AAT(prot2aa(p))%Codon(k)

! Fill the numerical sequence array

      CurrDNA%NUMseq(d-1)=NT2Int(CurrDNA%DNAseq(d-1:d-1))
      CurrDNA%NUMseq(d)=NT2Int(CurrDNA%DNAseq(d:d))
      CurrDNA%NUMseq(d+1)=NT2Int(CurrDNA%DNAseq(d+1:d+1))
  
    END DO main
    SequenceTranslated=.TRUE.
  END IF

END SUBROUTINE Translate_Protein
