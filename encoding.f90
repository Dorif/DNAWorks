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

  INTEGER :: i,j,m,n,a1,a2,b1,b2,t1,t2,fin

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
      SELECT CASE(CurrDNA%NUMseq(i+fin-j))
        CASE(-1)
          CurrDNA%ntID_Tip(i)=CurrDNA%ntID_Tip(i)+(1*(10**j))
        CASE(-3)
          CurrDNA%ntID_Tip(i)=CurrDNA%ntID_Tip(i)+(2*(10**j))
        CASE(3)
          CurrDNA%ntID_Tip(i)=CurrDNA%ntID_Tip(i)+(3*(10**j))
        CASE(1)
          CurrDNA%ntID_Tip(i)=CurrDNA%ntID_Tip(i)+(4*(10**j))
      END SELECT
      SELECT CASE(CurrDNA%NUMseq(i+j))
        CASE(-1)
          CurrDNA%ntID_TipRC(i)=CurrDNA%ntID_TipRC(i)+(4*(10**j))
        CASE(-3)
          CurrDNA%ntID_TipRC(i)=CurrDNA%ntID_TipRC(i)+(3*(10**j))
        CASE(3)
          CurrDNA%ntID_TipRC(i)=CurrDNA%ntID_TipRC(i)+(2*(10**j))
        CASE(1)
          CurrDNA%ntID_TipRC(i)=CurrDNA%ntID_TipRC(i)+(1*(10**j))
      END SELECT
    END DO
  END DO

! update repeat arrays

  fin=RepLen-1
  DO i=b1,b2
    CurrDNA%ntID_Rep(i)=0
    CurrDNA%ntID_RepRC(i)=0
    DO j=0,fin
      SELECT CASE(CurrDNA%NUMseq(i+fin-j))
        CASE(-1)
          CurrDNA%ntID_Rep(i)=CurrDNA%ntID_Rep(i)+(1*(10**j))
        CASE(-3)
          CurrDNA%ntID_Rep(i)=CurrDNA%ntID_Rep(i)+(2*(10**j))
        CASE(3)
          CurrDNA%ntID_Rep(i)=CurrDNA%ntID_Rep(i)+(3*(10**j))
        CASE(1)
          CurrDNA%ntID_Rep(i)=CurrDNA%ntID_Rep(i)+(4*(10**j))
      END SELECT
      SELECT CASE(CurrDNA%NUMseq(i+j))
        CASE(-1)
          CurrDNA%ntID_RepRC(i)=CurrDNA%ntID_RepRC(i)+(4*(10**j))
        CASE(-3)
          CurrDNA%ntID_RepRC(i)=CurrDNA%ntID_RepRC(i)+(3*(10**j))
        CASE(3)
          CurrDNA%ntID_RepRC(i)=CurrDNA%ntID_RepRC(i)+(2*(10**j))
        CASE(1)
          CurrDNA%ntID_RepRC(i)=CurrDNA%ntID_RepRC(i)+(1*(10**j))
      END SELECT
    END DO
  END DO

! update GC array

  fin=RepLen-1
  DO i=b1,b2
    CurrDNA%ntID_GC(i)=0
    DO j=(i+0),(i+fin)
      IF (ABS(CurrDNA%NUMseq(j)).eq.1) CurrDNA%ntID_GC(i)=CurrDNA%ntID_GC(i)+1
    END DO
  END DO

! update AT array

  fin=RepLen-1
  DO i=b1,b2
    CurrDNA%ntID_AT(i)=0
    DO j=(i+0),(i+fin)
      IF (ABS(CurrDNA%NUMseq(j)).eq.3) CurrDNA%ntID_AT(i)=CurrDNA%ntID_AT(i)+1
    END DO
  END DO

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
