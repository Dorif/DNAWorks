SUBROUTINE AT_Score
!
! Find all the 8 nt windows of solid AT content and update AScore.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j

  IF (TEST1) PRINT *,"AT_Score" !TEST1

  DO i=1,DNAlen
    CurrDNA%AScore(i) = 0
  END DO

  DO i=1,DNAlen-7
    IF (CurrDNA%ntID_AT(i).eq.0) THEN
      DO j=i,(i+7)
        CurrDNA%AScore(j)=CurrDNA%AScore(j)+1
      END DO
    END IF
  END DO

  CurrDNA%TotalAScore = 0.0     ! Initialize the repeat scores
  DO i=1,DNAlen
    CurrDNA%TotalAScore=CurrDNA%TotalAScore+CurrDNA%AScore(i)
  END DO
  CurrDNA%TotalAScore=CurrDNA%TotalAScore*20/DNAlen

END SUBROUTINE AT_Score
SUBROUTINE Average_Evaluate_Scores()
!
! This subroutine determines the average scores for the current sequence, each time
! changing the degenerate sequences.  It updates
! TScore, CScore, RScore, and PScore arrays, the Total*Score values.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i
  REAL :: dC, dL, dT, dR, dM, dG, dA, dF, dP, dTotal

  IF (TEST1) PRINT *,"Average_Evaluate_Scores" !TEST1

    dC=0
    dL=0
    dT=0
    dR=0
    dM=0
    dG=0
    dA=0
    dF=0
    dP=0
    dTotal=0

    DO i=1,NumDegPos*10
      CALL Fix_Degenerates
      CALL Create_ntID_Arrays
      CALL Temp_Score                 ! TScore(i) based on olaps
      CALL Misprime_Score             ! MScore(i) based on nt
      CALL Length_Score               ! LScore(i) based on nt
      CALL GapFix_Score               ! FScore(i) based on nt
      IF (ScoreCodons) CALL Codon_Score ! CScore(i) based on codons
      CALL Repeat_Score             ! RScore(i) based on nt
      CALL GC_Score                 ! GScore(i) based on nt
      CALL AT_Score                 ! AScore(i) based on nt
      CALL Pattern_Score            ! PScore(i) based on nt
    
      dC=CurrDNA%TotalCScore+dC
      dL=CurrDNA%TotalLScore+dL
      dT=CurrDNA%TotalTScore+dT
      dR=CurrDNA%TotalRScore+dR
      dM=CurrDNA%TotalMScore+dM
      dG=CurrDNA%TotalGScore+dG
      dA=CurrDNA%TotalAScore+dA
      dF=CurrDNA%TotalFScore+dF
      dP=CurrDNA%TotalPScore+dP
      dTotal=dC+dL+dT+dR+dM+dG+dA+dF+dP+dTotal

    END DO

    CurrDNA%TotalCScore=dC/(NumDegPos*10)
    CurrDNA%TotalLScore=dL/(NumDegPos*10)
    CurrDNA%TotalTScore=dT/(NumDegPos*10)
    CurrDNA%TotalRScore=dR/(NumDegPos*10)
    CurrDNA%TotalMScore=dM/(NumDegPos*10)
    CurrDNA%TotalGScore=dG/(NumDegPos*10)
    CurrDNA%TotalAScore=dA/(NumDegPos*10)
    CurrDNA%TotalFScore=dF/(NumDegPos*10)
    CurrDNA%TotalPScore=dP/(NumDegPos*10)

    CurrDNA%OverallScore = (Cwt*CurrDNA%TotalCScore)+&
                         (Lwt*CurrDNA%TotalLScore)+&
                         (Twt*CurrDNA%TotalTScore)+&
                         (Rwt*CurrDNA%TotalRScore)+&
                         (Mwt*CurrDNA%TotalMScore)+&
                         (Gwt*CurrDNA%TotalGScore)+&
                         (Awt*CurrDNA%TotalAScore)+&
                         (Fwt*CurrDNA%TotalFScore)+&
                         (Pwt*CurrDNA%TotalPScore)

    CALL Revert_Degenerates

END SUBROUTINE Average_Evaluate_Scores
SUBROUTINE Codon_Score
!
! This subroutine calculates a global score for codons based on frequency.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i

  IF (TEST1) PRINT *,"Codon_Score" !TEST1

  CurrDNA%TotalCScore=0.0

  IF (MutProtPos.eq.0) THEN
    DO i=1,PROTlen
      CurrDNA%CScore(i)=(1-(CFT(CurrDNA%prot2cod(i))%Freq/AAT(prot2aa(i))%Freq(1)))**4
    END DO
  ELSE
    CurrDNA%CScore(MutProtPos)=(1-(CFT(CurrDNA%prot2cod(MutProtPos))%Freq/AAT(prot2aa(MutProtPos))%Freq(1)))**4
  END IF

  DO i=1,PROTlen
    CurrDNA%TotalCScore=CurrDNA%TotalCScore+CurrDNA%CScore(i)
  END DO

  CurrDNA%TotalCScore = CurrDNA%TotalCScore/DNAlen

END SUBROUTINE Codon_Score
SUBROUTINE Decrement_Misprime_Arrays()
!
! Removes potential misprime pairs within the current mutant range

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,y,ct
  INTEGER :: TempM1(9999)
  INTEGER :: TempM2(9999)
  INTEGER :: TempMX(9999)

  IF (TEST2) PRINT *,"Decrement_Misprime_Arrays" !TEST2

  ct=0
  y=MutNtPos(MutNtNum)+1

  IF (CurrDNA%MN.gt.0) THEN
    loop: DO i=1,CurrDNA%MN
      IF ((((CurrDNA%M1(i)+MPLn).ge.MutNtPos(1)).and.&
            (CurrDNA%M1(i).le.y)).or.&
          (((CurrDNA%M2(i)+MPLn).ge.MutNtPos(1)).and.&
            (CurrDNA%M2(i).le.y))) THEN
        CYCLE loop
      ELSE
        ct=ct+1
        TempM1(ct)=CurrDNA%M1(i)
        TempM2(ct)=CurrDNA%M2(i)
        TempMX(ct)=CurrDNA%MX(i)
      END IF
    END DO loop

    CurrDNA%MN=ct

    DO i=1,CurrDNA%MN
      CurrDNA%M1(i)=TempM1(i)
      CurrDNA%M2(i)=TempM2(i)
      CurrDNA%MX(i)=TempMX(i)
    END DO
  END IF

END SUBROUTINE Decrement_Misprime_Arrays
SUBROUTINE Decrement_Repeat_Arrays()
!
! Remove repeat pairs and erase scores

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,y,ct
  INTEGER :: TempRS1(9999)
  INTEGER :: TempRS2(9999)
  INTEGER :: TempLn(9999)
  INTEGER :: TempRX(9999)

  IF (TEST2) PRINT *,"Decrement_Repeat_Arrays" !TEST2

  ct=0
  y=MutNtPos(MutNtNum)+1

  IF (CurrDNA%RN.gt.0) THEN
    loop: DO i=1,CurrDNA%RN
      IF ((((CurrDNA%RS1(i)+CurrDNA%RLn(i)).ge.MutNtPos(1)).and.&
            (CurrDNA%RS1(i).le.y)).or.&
          (((CurrDNA%RS2(i)+CurrDNA%RLn(i)).ge.MutNtPos(1)).and.&
            (CurrDNA%RS2(i).le.y))) THEN
        DO k=CurrDNA%RS1(i),(CurrDNA%RS1(i)+CurrDNA%RLn(i)-1)
          CurrDNA%RScore(k) = CurrDNA%RScore(k)-1
        END DO
        DO k=CurrDNA%RS2(i),(CurrDNA%RS2(i)+CurrDNA%RLn(i)-1)
          CurrDNA%RScore(k) = CurrDNA%RScore(k)-1
        END DO
        CYCLE loop
      ELSE
        ct=ct+1
        TempRS1(ct)=CurrDNA%RS1(i)
        TempRS2(ct)=CurrDNA%RS2(i)
        TempLn(ct)=CurrDNA%RLn(i)
        TempRX(ct)=CurrDNA%RX(i)
      END IF
    END DO loop

    CurrDNA%RN=ct

    DO i=1,CurrDNA%RN
      CurrDNA%RS1(i)=TempRS1(i)
      CurrDNA%RS2(i)=TempRS2(i)
      CurrDNA%RLn(i)=TempLn(i)
      CurrDNA%RX(i)=TempRX(i)
    END DO
  END IF

END SUBROUTINE Decrement_Repeat_Arrays
LOGICAL FUNCTION DegCmpr(instr,seq)
!
! This function compares a restriction site in degenerate form with a sequence.
! It returns .TRUE. if the site matches, and .FALSE. if it does not.  The two
! strings MUST be the same length.
!
! The function uses the NEB format of nucleotide degeneracy:
!
! B = C or G or T         rev. compl. = V
! D = A or G or T         rev. compl. = H
! H = A or C or T         rev. compl. = D
! K = G or T              rev. compl. = M
! M = A or C              rev. compl. = K
! N = A or C or G or T    rev. compl. = N
! R = A or G              rev. compl. = Y
! S = C or G              rev. compl. = S
! V = A or C or G         rev. compl. = B
! W = A or T              rev. compl. = W
! Y = C or T              rev. compl. = R
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  CHARACTER(LEN=100) :: instr,seq
  INTEGER :: i
  INTEGER :: Slen
  INTEGER :: stot

  IF (TEST3) PRINT *,"DegCmpr" !TEST3

  Slen=LEN_TRIM(instr)
  DegCmpr=.FALSE.

  stot=0

! First, check site in sense orientation

  DO i=1,Slen
    SELECT CASE(instr(i:i))
      CASE('A')
        IF (seq(i:i).EQ.'A') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('C')
        IF (seq(i:i).EQ.'C') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('G')
        IF (seq(i:i).EQ.'G') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('T')
        IF (seq(i:i).EQ.'T') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('B')
        IF (seq(i:i).EQ.'C'.OR.seq(i:i).EQ.'G'.OR.seq(i:i).EQ.'T') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('D')
        IF (seq(i:i).EQ.'A'.OR.seq(i:i).EQ.'G'.OR.seq(i:i).EQ.'T') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('H')
        IF (seq(i:i).EQ.'A'.OR.seq(i:i).EQ.'C'.OR.seq(i:i).EQ.'T') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('K')
        IF (seq(i:i).EQ.'G'.OR.seq(i:i).EQ.'T') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('M')
        IF (seq(i:i).EQ.'A'.OR.seq(i:i).EQ.'C') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('N')
        IF (seq(i:i).EQ.'A'.OR.seq(i:i).EQ.'C'.OR.seq(i:i).EQ.'G'.OR.seq(i:i).EQ.'T') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('R')
        IF (seq(i:i).EQ.'A'.OR.seq(i:i).EQ.'G') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('S')
        IF (seq(i:i).EQ.'C'.OR.seq(i:i).EQ.'G') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('V')
        IF (seq(i:i).EQ.'A'.OR.seq(i:i).EQ.'C'.OR.seq(i:i).EQ.'G') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('W')
        IF (seq(i:i).EQ.'A'.OR.seq(i:i).EQ.'T') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
      CASE('Y')
        IF (seq(i:i).EQ.'C'.OR.seq(i:i).EQ.'T') THEN
          stot=stot+1 ; ELSE ; EXIT ; END IF
    END SELECT
  END DO

  IF (stot.EQ.Slen) DegCmpr=.TRUE.

END FUNCTION DegCmpr
SUBROUTINE Equalize_Scores()
!
! Converts individual scores to codon-based Xscores for mutation rounds.
! XScore is a combination of all scores applied to each codon.  This should
! allow for a more targeted mutation.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j
  REAL :: CodPerOlap(999)   ! number of codons per overlap, for each overlap
  REAL :: TScorePerCod(999) ! average TScore per codon, for each overlap
  REAL :: r1

  IF (TEST1) PRINT *,"Equalize_Scores"

! Initialize XScore

  DO i=1,mutPROTnum  ! only use the mutatable codons
    XScore(i)=0
  END DO

! Initialize values for CodPerOlap

  DO i=1,999
    CodPerOlap(i)=0
  END DO

! Find how many codons are in each overlap, avoiding non-coding regions
! If nt is within a codon, is within an overlap, and is unique protein residue

  DO i=1,DNAlen
    IF ((nt2prot(i).ne.0).and.(nt2overlap(i).ne.0)) THEN
      IF (i.eq.1 .or. nt2prot(i).ne.nt2prot(i-1)) THEN
        CodPerOlap(nt2overlap(i))=CodPerOlap(nt2overlap(i))+1
      END IF
    END IF
  END DO

! The TScore for each codon a fraction of the total TScore(i) for the overlap

  DO j=1,CurrDNA%NumOlaps
    IF (CodPerOlap(j).ne.0) TScorePerCod(j)=Twt*(CurrDNA%TScore(j)/CodPerOlap(j))
  END DO

! Assign the XScore for TScore, CScore, RScore, PScore, GScore, LScore,
! and AScore for each codon

  DO i=1,mutPROTnum

    j=prot2nt(mutPROT2prot(i))  ! the middle nt of the codon

! if the middle nt of a codon is within an overlap, the XScore for that codon 
! is the average TScore per codon for the overlap

    IF (nt2overlap(j).ne.0) XScore(i)=TScorePerCod(nt2overlap(j))

! the CScore contribution is already for the codon

    XScore(i)=XScore(i)+(Cwt*CurrDNA%CScore(i))+&
 (Rwt*REAL(CurrDNA%RScore(j-1)+CurrDNA%RScore(j)+CurrDNA%RScore(j+1)))+&
 (Mwt*REAL(CurrDNA%MScore(j-1)+CurrDNA%MScore(j)+CurrDNA%MScore(j+1)))+&
 (Gwt*REAL(CurrDNA%GScore(j-1)+CurrDNA%GScore(j)+CurrDNA%GScore(j+1)))+&
 (Awt*REAL(CurrDNA%AScore(j-1)+CurrDNA%AScore(j)+CurrDNA%AScore(j+1)))+&
 (Lwt*REAL(CurrDNA%LScore(j-1)+CurrDNA%LScore(j)+CurrDNA%LScore(j+1)))+&
 (Fwt*REAL(CurrDNA%FScore(j-1)+CurrDNA%FScore(j)+CurrDNA%FScore(j+1)))+&
 (Pwt*REAL(CurrDNA%PScore(j-1)+CurrDNA%PScore(j)+CurrDNA%PScore(j+1)))

  END DO

END SUBROUTINE Equalize_Scores
SUBROUTINE Evaluate_Scores()
!
! This subroutine determines the scores for the current sequence.  It updates
! TScore, CScore, RScore, and PScore arrays, the Total*Score values.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j
  REAL :: dC, dL, dT, dR, dM, dG, dA, dF, dP, dTotal

  IF (TEST1) PRINT *,"Evaluate_Scores" !TEST1

! Degenerate sequences: loop several times and get average scores

  IF (NumDegPos.eq.0) THEN

! If the sequence is recently translated, create the ntID arrays

    CALL Create_ntID_Arrays
    CALL Temp_Score                 ! TScore(i) based on olaps
    CALL Misprime_Score             ! MScore(i) based on nt
    CALL Length_Score               ! LScore(i) based on nt
    CALL GapFix_Score               ! FScore(i) based on nt

! The following scores will not change when the overlap positions are moved,
! but only when the sequence is re-translated after a mutation (also does not
! apply for DNA-only runs)

    IF (ScoreCodons) CALL Codon_Score ! CScore(i) based on codons
    CALL Repeat_Score             ! RScore(i) based on nt
    CALL GC_Score                 ! GScore(i) based on nt
    CALL AT_Score                 ! AScore(i) based on nt
    CALL Pattern_Score            ! PScore(i) based on nt

! Update CurrDNA%OverallScore

    CurrDNA%OverallScore = (Cwt*CurrDNA%TotalCScore)+&
                         (Lwt*CurrDNA%TotalLScore)+&
                         (Twt*CurrDNA%TotalTScore)+&
                         (Rwt*CurrDNA%TotalRScore)+&
                         (Mwt*CurrDNA%TotalMScore)+&
                         (Gwt*CurrDNA%TotalGScore)+&
                         (Awt*CurrDNA%TotalAScore)+&
                         (Fwt*CurrDNA%TotalFScore)+&
                         (Pwt*CurrDNA%TotalPScore)
  ELSE
    CALL Average_Evaluate_Scores
  END IF

END SUBROUTINE Evaluate_Scores
SUBROUTINE Find_Actual_Misprimes()
!
! If one of the positions in a misprime pair aligns to the end of an overlap,
! and if the tip of the overlap is identical (direct or inverse), then raise
! the score on the nts (CurrDNA%MScore).

! 1. direct-sense(DS): forward primer mispriming on the sense strand
!
!      -------------->               -------------->
!      |||||||||||||||               ..........|||||
! -------------------------------------------------------
!
! 2. inverse-sense(IS): reverse primer mispriming on the sense strand
! NOTE THAT IF THE FORWARD OLIGO MATCHES M2, MSX = 5, NOT 2
!
!                                    -------------->
!                                    ..........|||||
! -------------------------------------------------------
!      |||||||||||||||
!      <--------------
!
! 3. inverse-antisense(IA): forward primer mispriming on the antisense strand
! NOTE THAT IF THE REVERSE OLIGO MATCHES M2, MSX = 6, NOT 3
!
!      -------------->
!      |||||||||||||||
! -------------------------------------------------------
!                                    |||||..........
!                                    <--------------
!
! 4. direct-antisense(DA): reverse primer mispriming on the antisense strand
!
! -------------------------------------------------------
!      |||||||||||||||               |||||..........
!      <--------------               <--------------
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,mp
  INTEGER :: o1,o2,m1,m2,mx

  IF (TEST2) PRINT *,"Find_Actual_Misprimes" !TEST2

  mp=(CurrDNA%NumOlaps+1)/2

! Initialize actual misprime arrays

  DO i=1,DNAlen
    CurrDNA%MScore(i) = 0
  END DO
  CurrDNA%MSN=0

  mpair: DO i=1,CurrDNA%MN
    m1=CurrDNA%M1(i)
    m2=CurrDNA%M2(i)
    mx=CurrDNA%MX(i)
    olap: DO j=1,CurrDNA%NumOlaps
      o1=CurrDNA%OlapsPos(j,1)
      o2=(CurrDNA%OlapsPos(j,2)-MPLn+1)
      IF (TBIO) THEN
        IF (j.lt.mp) THEN
          SELECT CASE(mx)
            CASE(1)                                          ! direct-sense
              IF (o2.eq.m1) THEN
                CALL Increment_Misprime_Scores(o2,m2,1,j)
              ELSE IF (o2.eq.m2) THEN
                CALL Increment_Misprime_Scores(o2,m1,1,j)
              END IF
            CASE(2)                                         ! inverse-sense
              IF (o2.eq.m2) CALL Increment_Misprime_Scores(m2,m1,5,j)
            CASE(3)                                     ! inverse-antisense
              IF (o2.eq.m1) CALL Increment_Misprime_Scores(o2,m2,3,j)
          END SELECT
        ELSE IF (j.gt.mp) THEN
          SELECT CASE(mx)
            CASE(2)                                         ! inverse-sense
              IF (o1.eq.m1) CALL Increment_Misprime_Scores(o1,m2,2,j)
            CASE(3)                                     ! inverse-antisense
              IF (o1.eq.m2) CALL Increment_Misprime_Scores(o1,m1,6,j)
            CASE(4)                                      ! direct-antisense
              IF (o1.eq.m1) THEN
                CALL Increment_Misprime_Scores(o1,m2,4,j)
              ELSE IF (o1.eq.m2) THEN
                CALL Increment_Misprime_Scores(o1,m1,4,j)
              END IF
          END SELECT
        ELSE
          SELECT CASE(mx)
            CASE(1)                                          ! direct-sense
              IF (o2.eq.m1) THEN
                CALL Increment_Misprime_Scores(o2,m2,1,j)
              ELSE IF (o2.eq.m2) THEN
                CALL Increment_Misprime_Scores(o2,m1,1,j)
              END IF
            CASE(2)                                         ! inverse-sense
              IF (o2.eq.m2) THEN
                CALL Increment_Misprime_Scores(m2,m1,5,j)
              ELSE IF (o1.eq.m1) THEN
                CALL Increment_Misprime_Scores(o1,m2,2,j)
              END IF
            CASE(3)                                     ! inverse-antisense
              IF (o2.eq.m1) THEN
                CALL Increment_Misprime_Scores(o2,m2,3,j)
              ELSE IF (o1.eq.m2) THEN
                CALL Increment_Misprime_Scores(o1,m1,6,j)
              END IF
            CASE(4)                                      ! direct-antisense
              IF (o1.eq.m1) THEN
                CALL Increment_Misprime_Scores(o1,m2,4,j)
              ELSE IF (o1.eq.m2) THEN
                CALL Increment_Misprime_Scores(o1,m1,4,j)
              END IF
          END SELECT
        END IF
      ELSE
        IF (MOD(j,2).eq.0) THEN
          CYCLE olap
        ELSE
          SELECT CASE(mx)
            CASE(1)                                          ! direct-sense
              IF (o2.eq.m1) THEN
                CALL Increment_Misprime_Scores(o2,m2,1,j)
              ELSE IF (o2.eq.m2) THEN
                CALL Increment_Misprime_Scores(o2,m1,1,j)
              END IF
            CASE(2)                                         ! inverse-sense
              IF (o2.eq.m2) THEN
                CALL Increment_Misprime_Scores(m2,m1,5,j)
              ELSE IF (o1.eq.m1) THEN
                CALL Increment_Misprime_Scores(o1,m2,2,j)
              END IF
            CASE(3)                                     ! inverse-antisense
              IF (o2.eq.m1) THEN
                CALL Increment_Misprime_Scores(o2,m2,3,j)
              ELSE IF (o1.eq.m2) THEN
                CALL Increment_Misprime_Scores(o1,m1,6,j)
              END IF
            CASE(4)                                      ! direct-antisense
              IF (o1.eq.m1) THEN
                CALL Increment_Misprime_Scores(o1,m2,4,j)
              ELSE IF (o1.eq.m2) THEN
                CALL Increment_Misprime_Scores(o1,m1,4,j)
              END IF
          END SELECT
        END IF
      END IF
    END DO olap
  END DO mpair

END SUBROUTINE Find_Actual_Misprimes
SUBROUTINE Find_Potential_Misprimes
!
! This subroutine finds all misprimes in the sequence, both direct and
! inverse, regardless of position, equal or longer than MPLn.
! The inverse search allows palindromic misprimes (i=j).
! The number of potential misprimes is in CurrDNA%MN
! It records the positions and sizes in the global arrays CurrDNA%M1,
! CurrDNA%M2, and CurrDNA%MX.
! The actual misprimes are determined by Find_Actual_Misprimes

! 1. direct-sense(DS): forward primer mispriming on the sense strand
!
!      -------------->               -------------->
!      |||||||||||||||               ..........|||||
! -------------------------------------------------------
!
! 2. inverse-sense(IS): reverse primer mispriming on the sense strand
! NOTE THAT IF THE FORWARD OLIGO MATCHES M2, MSX = 5, NOT 2
!
!                                    -------------->
!                                    ..........|||||
! -------------------------------------------------------
!      |||||||||||||||
!      <--------------
!
! 3. inverse-antisense(IA): forward primer mispriming on the antisense strand
! NOTE THAT IF THE REVERSE OLIGO MATCHES M2, MSX = 6, NOT 3
!
!      -------------->
!      |||||||||||||||
! -------------------------------------------------------
!                                    |||||..........
!                                    <--------------
!
! 4. direct-antisense(DA): reverse primer mispriming on the antisense strand
!
! -------------------------------------------------------
!      |||||||||||||||               |||||..........
!      <--------------               <--------------
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j
  LOGICAL,EXTERNAL :: HMatchNum

  IF (TEST2) PRINT *,"Find_Potential_Misprimes" !TEST2

! Initialize the potential misprime arrays

  CurrDNA%MN=0

! Find the potential misprimes

  DO i=1,DNAlen-MPLn+1
    DO j=i,DNAlen-MPLn+1
      IF (HMatchNum(i,j,1)) THEN
        IF (CurrDNA%ntID_Tip(i+MPLn-MPTip).eq.CurrDNA%ntID_Tip(j+MPLn-MPTip)) &
          CALL Increment_Misprime_Arrays(i,j,1)
        IF (CurrDNA%ntID_Tip(i).eq.CurrDNA%ntID_Tip(j)) &
          CALL Increment_Misprime_Arrays(i,j,4)
      END IF
      IF (HMatchNum(i,j,-1)) THEN
        IF (CurrDNA%ntID_Tip(i).eq.CurrDNA%ntID_TipRC(j+MPLn-MPTip)) &
          CALL Increment_Misprime_Arrays(i,j,2)
        IF (CurrDNA%ntID_Tip(i+MPLn-MPTip).eq.CurrDNA%ntID_TipRC(j)) &
          CALL Increment_Misprime_Arrays(i,j,3)
      END IF
    END DO
  END DO

END SUBROUTINE Find_Potential_Misprimes
SUBROUTINE Find_Repeats()
!
! This subroutine finds all repeats in the sequence, both direct and
! inverse, regardless of position, equal or longer than RepLen.
! The inverse search allows palindromic repeats (i=j).
! It records the positions and sizes in the global arrays CurrDNA%RS1,
! CurrDNA%RS2, and CurrDNA%RLn.  Then it overwrites the array CurrDNA%RScore.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j
  LOGICAL,EXTERNAL :: PairWithinKnownRepeat

  IF (TEST2) PRINT *,"Find_Repeats" !TEST2

  DO i=1,DNAlen
    CurrDNA%RScore(i) = 0
  END DO
  CurrDNA%RN=0

  DO i=1,DNAlen-RepLen+1
    DO j=i,DNAlen-RepLen+1
      IF (i.ne.j) THEN
        IF (.not.PairWithinKnownRepeat(i,j,1)) THEN
          IF (CurrDNA%ntID_Rep(i).eq.CurrDNA%ntID_Rep(j)) &
            CALL Increment_Repeat_Arrays(i,j,1)
        END IF
      END IF
      IF (.not.PairWithinKnownRepeat(i,j,-1)) THEN
        IF (CurrDNA%ntID_Rep(i).eq.CurrDNA%ntID_RepRC(j)) &
          CALL Increment_Repeat_Arrays(i,j,-1)
      END IF
    END DO
  END DO

END SUBROUTINE Find_Repeats
SUBROUTINE GC_Score
!
! Find all the 8 nt windows of solid GC content and update GScore.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j

  IF (TEST1) PRINT *,"GC_Score" !TEST1

  DO i=1,DNAlen
    CurrDNA%GScore(i) = 0
  END DO

  DO i=1,DNAlen-7
    IF (CurrDNA%ntID_GC(i).eq.0) THEN
      DO j=i,(i+7)
        CurrDNA%GScore(j)=CurrDNA%GScore(j)+1
      END DO
    END IF
  END DO

  CurrDNA%TotalGScore = 0.0     ! Initialize the repeat scores
  DO i=1,DNAlen
    CurrDNA%TotalGScore=CurrDNA%TotalGScore+CurrDNA%GScore(i)
  END DO
  CurrDNA%TotalGScore=CurrDNA%TotalGScore*20/DNAlen

END SUBROUTINE GC_Score
SUBROUTINE GapFix_Score
!
! This subroutine returns the GapFix scores for each nt position in the
! GapFixPos array.
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,c

  IF (TEST1) PRINT *,"GapFix_Score" !TEST1

! initialize scores

  CurrDNA%TotalFScore=0
  DO i = 1,DNAlen
    CurrDNA%FScore(i)=0
  END DO

  DO i=1,DNAlen

! if the position should be within a gap but is in an overlap instead

    IF (CurrDNA%GapFixPos(i).and.nt2overlap(i).gt.0) THEN
      CurrDNA%FScore(i)=10
    END IF
  END DO

! generate summary of scores

  DO i=1,DNAlen
    CurrDNA%TotalFScore=CurrDNA%TotalFScore+CurrDNA%FScore(i)
  END DO
  CurrDNA%TotalFScore=CurrDNA%TotalFScore*20/DNAlen

END SUBROUTINE GapFix_Score
LOGICAL FUNCTION HMatchNum(pos1,pos2,dir)
!
! If two positions of equal length are homologous (MaxNonId or fewer
! non-identical nts), returns true

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: pos1,pos2,i,a,b,dir,ct

  IF (TEST3) PRINT *,"HMatchNum" !TEST3

  ct=0

  IF (dir.eq.1) THEN
    direct: DO i=1,MPLn
      HMatchNum=.FALSE.
      IF (pos1.eq.pos2) EXIT direct
      a=pos1+i-1
      b=pos2+i-1
      IF (b.gt.DNAlen) EXIT direct
      IF ((CurrDNA%NUMseq(a)-CurrDNA%NUMseq(b)).ne.0) THEN
        ct=ct+1
        IF (ct.gt.MaxNonId) EXIT direct
      END IF
!      PRINT *,dir,pos1,pos2,CurrDNA%NUMseq(a),CurrDNA%NUMseq(b)
      HMatchNum=.TRUE.
    END DO direct
  ELSE
    inverse: DO i=1,MPLn
      HMatchNum=.FALSE.
      a=pos1+i-1
      b=pos2+MPLn-i
      IF (b.gt.DNAlen) EXIT inverse
      IF ((CurrDNA%NUMseq(a)+CurrDNA%NUMseq(b)).ne.0) THEN
        ct=ct+1
        IF (ct.gt.MaxNonId) EXIT inverse
      END IF
      HMatchNum=.TRUE.
    END DO inverse
  END IF

END FUNCTION HMatchNum
SUBROUTINE Increment_Misprime_Arrays(pos1,pos2,dir)
!
! Add another misprime pair to the arrays and update scores

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: pos1,pos2,dir,i
  CHARACTER(LEN=80) :: text

  IF (TEST2) PRINT *,"Increment_Misprime_Arrays" !TEST2

  CurrDNA%MN=CurrDNA%MN+1
  IF (CurrDNA%MN.ge.DNAlen) THEN
    WRITE(text,FMT="('MN = ',i9,' Too many misprimes.')") CurrDNA%MN
!    DO i=1,CurrDNA%MN
!      PRINT *,CurrDNA%M1(i),CurrDNA%M2(i),CurrDNA%MX(i)
!    END DO
    CALL Stop_Program(text)
  END IF
  CurrDNA%M1(CurrDNA%MN)=pos1
  CurrDNA%M2(CurrDNA%MN)=pos2
  CurrDNA%MX(CurrDNA%MN)=dir

END SUBROUTINE Increment_Misprime_Arrays
SUBROUTINE Increment_Misprime_Scores(o,m,t,j)
!
! Increment the scores and update the actual misprime arrays

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: o,m,t,i,j
  LOGICAL :: x

  IF (TEST2) PRINT *,"Increment_Misprime_Scores" !TEST2

  CurrDNA%MSN=CurrDNA%MSN+1
  CurrDNA%MS1(CurrDNA%MSN)=o
  CurrDNA%MS2(CurrDNA%MSN)=m
  CurrDNA%MSX(CurrDNA%MSN)=t
  CurrDNA%MOL(CurrDNA%MSN)=j
  DO i=o,o+MPLn-1
    CurrDNA%MScore(i)=CurrDNA%MScore(i)+1
  END DO
  DO i=m,m+MPLn-1
    CurrDNA%MScore(i)=CurrDNA%MScore(i)+1
  END DO

END SUBROUTINE Increment_Misprime_Scores
SUBROUTINE Increment_Repeat_Arrays(i,j,dir)
!
! Add another repeat pair to the arrays and update scores after expansion

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k
  INTEGER :: pos1,pos2,dir,length
  INTEGER :: last,diff,a

  IF (TEST2) PRINT *,"Increment_Repeat_Arrays" !TEST2

  pos1=i
  pos2=j
  diff=j-i
  length=RepLen

! Expand direct repeats

  IF (dir.eq.1) THEN
    starting: DO pos1=(i-1),1,-1
      pos2=pos1+diff
      IF (CurrDNA%NUMseq(pos1).ne.CurrDNA%NUMseq(pos2)) THEN
        pos2=pos2+1
        EXIT starting
      END IF
    END DO starting
    pos1=pos2-diff
    length=RepLen+(i-pos1)   ! In case pos2 is DNAlen-RepLen-1
    ending: DO last=(j+RepLen),DNAlen+1
      IF (last.eq.(DNAlen+1)) EXIT ending
      IF (CurrDNA%NUMseq(last-diff).ne.CurrDNA%NUMseq(last)) EXIT ending
    END DO ending
    length=last-pos2                ! Final answer

  ELSE

! Expand inverse repeats

    startingRC: DO a=1,DNAlen
      pos1=i-a
      last=j+length-1+a
      IF ((pos1.lt.1).or.(last.gt.DNAlen).or.&
          (CurrDNA%NUMseq(pos1).ne.(-1*(CurrDNA%NUMseq(last))))) THEN
        pos1=pos1+1
        last=last-1
        EXIT startingRC
      END IF
    END DO startingRC
    endingRC: DO a=1,DNAlen
      pos2=j-a
      last=i+length-1+a
      IF ((pos2.lt.1).or.(last.gt.DNAlen).or.&
          (CurrDNA%NUMseq(last).ne.(-1*(CurrDNA%NUMseq(pos2))))) THEN
        pos2=pos2+1
        last=last-1
        EXIT endingRC
      END IF
    END DO endingRC
    length=last-pos1+1                ! Final answer
  END IF

  CurrDNA%RN=CurrDNA%RN+1
  IF (CurrDNA%RN.ge.DNAlen) CALL Stop_Program("Too many repeats.")
  CurrDNA%RS1(CurrDNA%RN)=pos1
  CurrDNA%RS2(CurrDNA%RN)=pos2
  CurrDNA%RLn(CurrDNA%RN)=length
  CurrDNA%RX(CurrDNA%RN)=dir
  DO k=pos1,(pos1+length-1)
    CurrDNA%RScore(k) = CurrDNA%RScore(k)+1
  END DO
  DO k=pos2,(pos2+length-1)
    CurrDNA%RScore(k) = CurrDNA%RScore(k)+1
  END DO

END SUBROUTINE Increment_Repeat_Arrays
SUBROUTINE Length_Score
!
! This subroutine evaluates the length of the oligos and gives a penalty to
! all the nts in the oligo if it exceeds OligoLen (except for the first and
! last oligos, of course).

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j
  INTEGER :: overrun     !the length of the oligo goes past OligoLen

  IF (TEST1) PRINT *,"Length_Score" !TEST1

  DO i=CurrDNA%OlapsPos(1,1),CurrDNA%OlapsPos(CurrDNA%NumOlaps,2)
    CurrDNA%LScore(i)=0
  END DO

!PRINT *,'START'

  DO i=2,CurrDNA%NumOlaps
    overrun=CurrDNA%OlapsPos(i,2)-CurrDNA%OlapsPos((i-1),1)-OligoLen+1

!PRINT *,i,CurrDNA%OlapsPos(i,2),CurrDNA%OlapsPos((i-1),1),OligoLen,overrun,CurrDNA%MeltT(i)

    IF (overrun.gt.0) THEN
      DO j=CurrDNA%OlapsPos((i-1),1),CurrDNA%OlapsPos(i,2)
        CurrDNA%LScore(j)=(overrun+2)**2
      END DO
    END IF
  END DO

!PRINT *,'FINISH'

  CurrDNA%TotalLScore = 0.0
  DO i=CurrDNA%OlapsPos(1,1),CurrDNA%OlapsPos(CurrDNA%NumOlaps,2)
     CurrDNA%TotalLScore = CurrDNA%TotalLScore + CurrDNA%LScore(i)
  END DO
  CurrDNA%TotalLScore=CurrDNA%TotalLScore*20/DNAlen

END SUBROUTINE Length_Score
SUBROUTINE Misprime_Score
!
! Determine the current mispriming score.  If at the beginning of a run
! (MutProtPos=0), then generate the potential misprime arrays, and then
! find the actual misprimes.
!
! During the run, the potential misprime arrays only need to be regenerated
! once after each mutation.  The potential misprime arrays are only modified
! around the site of mutation.  The actual misprimes are then determined after
! every overlap set generation.
!
! Only evaluating the mutation site speeds up the calculation more than 10-fold.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i

  IF (TEST1) PRINT *,"Misprime_Score" !TEST1

  IF (MutProtPos.eq.0) THEN
    IF (SequenceTranslated) CALL Find_Potential_Misprimes
    CALL Find_Actual_Misprimes
  ELSE
    IF (SequenceTranslated) CALL Find_Mut_Pot_Misprimes
    CALL Find_Actual_Misprimes
  END IF

  CurrDNA%TotalMScore = 0.0     ! Initialize the mispriming scores
  DO i=1,DNAlen
    CurrDNA%TotalMScore=CurrDNA%TotalMScore+CurrDNA%MScore(i)
  END DO
  CurrDNA%TotalMScore=CurrDNA%TotalMScore*20/DNAlen

END SUBROUTINE Misprime_Score
LOGICAL FUNCTION PairWithinKnownRepeat(i,j,dir)
!
! Returns true if pair of residues is already part of a repeat pair

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,dir

  IF (TEST3) PRINT *,"PairWithinKnownRepeat" !TEST3

  PairWithinKnownRepeat=.FALSE.

  IF (CurrDNA%RN.gt.0) THEN
    loop2: DO k=1,CurrDNA%RN
      IF ((CurrDNA%RS1(k).le.i).and.&
          ((CurrDNA%RS1(k)+CurrDNA%RLn(k)-RepLen).ge.i).and.&
          (CurrDNA%RS2(k).le.j).and.&
          ((CurrDNA%RS2(k)+CurrDNA%RLn(k)-RepLen).ge.j).and.&
           (CurrDNA%RX(k).eq.dir)) THEN
        PairWithinKnownRepeat=.TRUE.
        EXIT loop2
      END IF
    END DO loop2
  END IF

END FUNCTION PairWithinKnownRepeat
SUBROUTINE Pattern_Score
!
! This subroutine looks through the DNA sequence corresponding to the protein
!   region and identifies sequence patterns, either for restriction sites or
!   user-input sequences.  When it finds the pattern, it increases the score
!   for those nucleotides in the pattern.
!
! There are two situtations -- when degenerate patterns are present, and when
!   they are not.

!
! This scoring evaluation currently uses text-based comparisons, so it will
! be slow...

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,k,n
  CHARACTER(LEN=9999) :: text,ftext,rtext
  LOGICAL,EXTERNAL :: DegCmpr
  INTEGER :: start,curr,finis

  IF (TEST1) PRINT *,"Pattern_Score" !TEST1

  CurrDNA%TotalPScore = 0.0
  DO k=1,DNAlen
    CurrDNA%PScore(k) = 0
  END DO

  main: DO i=1,PTNnum

! skip the site if it is an isoschizomer

    IF (PTN(i)%Isoschiz) CYCLE main

! treat degenerate sites differently

    IF (PTN(i)%Degen) THEN
      ftext=PTN(i)%Seq(1:PTN(i)%Len)
      rtext=PTN(i)%SeqRC(1:PTN(i)%Len)
      deg: DO j=1,(DNAlen-PTN(i)%Len+1)
        IF (DegCmpr(ftext(1:PTN(i)%Len),CurrDNA%DNAseq(j:j+PTN(i)%Len-1))) THEN
          DO k=j,j+PTN(i)%Len-1
            CurrDNA%PScore(k)=CurrDNA%PScore(k)+1
          END DO
        END IF
        IF (.not.PTN(i)%SelfCompl) THEN
          IF (DegCmpr(rtext(1:PTN(i)%Len),CurrDNA%DNAseq(j:j+PTN(i)%Len-1))) THEN
            DO k=j,j+PTN(i)%Len-1
              CurrDNA%PScore(k)=CurrDNA%PScore(k)+1
            END DO
          END IF
        END IF
      END DO deg

    ELSE

! not degenerate

      curr=0
      start=1
      finis=DNAlen
! forward direction
      forward: DO n=1,DNAlen
        j=INDEX(CurrDNA%DNAseq(start:finis),PTN(i)%Seq(1:PTN(i)%Len))
        curr=curr+j
        IF (j.eq.0) THEN
          EXIT forward
        ELSE
          DO k=curr,(curr+PTN(i)%Len)
            CurrDNA%PScore(k)=CurrDNA%PScore(k)+1
          END DO
          start=curr+1
        END IF
      END DO forward
! reverse direction if needed
      IF (.not.PTN(i)%SelfCompl) THEN
        curr=0
        start=1
        reverse: DO n=1,DNAlen
          j=INDEX(CurrDNA%DNAseq(start:finis),PTN(i)%SeqRC(1:PTN(i)%Len))
          curr=curr+j
          IF (j.eq.0) THEN
            EXIT reverse
          ELSE
            DO k=curr,(curr+PTN(i)%Len)
              CurrDNA%PScore(k)=CurrDNA%PScore(k)+1
            END DO
            start=curr+1
          END IF
        END DO reverse
      END IF
    END IF
  END DO main

  DO i=1,DNAlen
    CurrDNA%TotalPScore=CurrDNA%TotalPScore+CurrDNA%PScore(i)
  END DO
  CurrDNA%TotalPScore=CurrDNA%TotalPScore*20/DNAlen

END SUBROUTINE Pattern_Score
SUBROUTINE Repeat_Score
!
! This subroutine calculates the score for finding tandem repeats, both
! direct and inverted (RC), in the trial DNA sequence.  The entire sequence
! is queried against itself in the first run.  If a codon has been mutated,
! only the small region around the mutation is queried against the sequence.
! The score is then applied to the nts it is found within.
!
! Only evaluating the mutation site speeds up the calculation more than 10-fold.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,j,pos

  IF (TEST1) PRINT *,"Repeat_Score" ! TEST1

  IF (MutProtPos.eq.0) THEN
    CALL Find_Repeats
  ELSE
    CALL Find_Mutated_Repeats
  END IF

  CurrDNA%TotalRScore = 0.0     ! Initialize the repeat scores
  DO i=1,DNAlen
    CurrDNA%TotalRScore=CurrDNA%TotalRScore+CurrDNA%RScore(i)
  END DO
  CurrDNA%TotalRScore=CurrDNA%TotalRScore*20/DNAlen

END SUBROUTINE Repeat_Score
SUBROUTINE Temp_Score
!
! This subroutine returns the Melting Temperature scores for each overlap as
! well as the total Tm score.
!
! The score is calculated as follow:
!   - for the temperatures within the range (MeltTemp-MeltTol ... MeltTemp +
!     MeltTol) score is calculated using a quadratic function of difference
!   - outside of this range also the 10 times square of the second difference
!     is added
!
  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: i,c
  REAL :: diff,maxT,minT
  REAL,EXTERNAL :: TmCalc

  IF (TEST1) PRINT *,"Temp_Score" !TEST1

  maxT=0
  minT=999

  CurrDNA%TotalTScore=0
  DO i = 1,CurrDNA%NumOlaps
    CurrDNA%TScore(i)=0
  END DO

  DO i=1,CurrDNA%NumOlaps
    CurrDNA%MeltT(i)=TmCalc(CurrDNA%OlapsPos(i,1),CurrDNA%OlapsPos(i,2))
  END DO

  DO i = 1,CurrDNA%NumOlaps
    diff=ABS(MeltTemp-CurrDNA%MeltT(i))
    IF(diff.gt.MeltTol)THEN
      diff=MAX(1.0,(diff-MeltTol))
      CurrDNA%TScore(i)=(diff**2)/10
    ELSE
      CurrDNA%TScore(i)=0
    END IF
  END DO

  DO i=1,CurrDNA%NumOlaps
    CurrDNA%TotalTScore=CurrDNA%TotalTScore+CurrDNA%TScore(i)
  END DO
  CurrDNA%TotalTScore=CurrDNA%TotalTScore*20/DNAlen

END SUBROUTINE Temp_Score
REAL FUNCTION TmCalc(start,finish)
!
! This function returns the melting temperature for an overlap of nucleotides
!  start to finish.
!
! The melting temperature is based on the paper by John SantaLucia
!  Jr., "A unified view of polymer, dumbbell, and oligonucleotide DNA
!  nearest-neighbor thermodynamics", Biochemistry Vol. 95, Issue 4,
!  1460-1465, 1998.  Tm in Celcius is calculated using the following formula:
!
!  Tm = [dH/(dS+(R*ln(OligoConc))+(0.368*ln(Na)*N)+(??*??MgConc*N))]-273.15
!
!      where
!              dH = sum of individual dH for each nucleotide pair, kcal/mol
!              dS = sum of individual dS for each nucleotide pair, cal/k*mol
!              RGasConstant = gas constant, 1.987 cal/K*mol
!              OligoConc = template concentration, mol/liter
!              Na = monovalent cation concentration (sodium), mol/liter
!              MgConc = magnesium concentration, mol/liter
!              N = number of backbone phosphates (number of nucleotides - 1)
!              Kelvin = convert Kelvin to Celsius, 273.15
!
! Additionally, the Tm is modified by the terminal nucleotides:
!
!              dH = dH + 2.2 for A/T
!              dS = dS + 6.935 for A/T
!
! The Tm is also modified by the presence of self-complementarity:
!
!              dS = dS - 1.4
!
! Numerical sequences are used instead of strings
!
! This subroutine is by far the most heavily used subroutine in the program.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: start,finish,i,i1,i2
  REAL :: dh,ds
  LOGICAL :: self_compl

! Lookup table: NUMseq values {-3,-1,1,3} -> indices {1,2,3,4} (C,A,T,G)
  INTEGER, SAVE :: nt_idx(-3:3)
  DATA nt_idx /1, 0, 2, 0, 3, 0, 4/

! Nearest-neighbor dH table (kcal/mol), indexed dh_nn(first_nt, second_nt)
!   Index order: C=1, A=2, T=3, G=4
!   DATA fills column-major: col 1 first, then col 2, etc.
  REAL, SAVE :: dh_nn(4,4)
  DATA dh_nn / &
    -8.0, -8.4, -8.2, -9.8, &   ! col 1 (2nd=C): CC,AC,TC,GC
    -8.5, -7.9, -7.2, -8.2, &   ! col 2 (2nd=A): CA,AA,TA,GA
    -7.8, -7.2, -7.9, -8.4, &   ! col 3 (2nd=T): CT,AT,TT,GT
   -10.6, -7.8, -8.5, -8.0  /   ! col 4 (2nd=G): CG,AG,TG,GG

! Nearest-neighbor dS table (cal/K*mol)
  REAL, SAVE :: ds_nn(4,4)
  DATA ds_nn / &
   -19.8612,  -22.44082, -22.24469, -24.37776, &  ! col 1 (2nd=C)
   -22.73082, -22.2473,  -21.34081, -22.24469, &  ! col 2 (2nd=A)
   -21.02469, -20.38082, -22.2473,  -22.44082, &  ! col 3 (2nd=T)
   -27.17776, -21.02469, -22.73082, -19.8612   /  ! col 4 (2nd=G)

  IF (TEST3) PRINT *,"TmCalc" !TEST3

! Initialize values

  self_compl=.FALSE.
  dh=0.2
  ds=-5.68

! Make sure the overlap is more than 7 nt long

   IF ((finish-start+1).le.7) THEN
     tmcalc=0
     RETURN
   END IF

! Sum the dH, dS values using lookup tables

  DO i=start,finish-1
    i1 = nt_idx(CurrDNA%NUMseq(i))
    i2 = nt_idx(CurrDNA%NUMseq(i+1))
    dh = dh + dh_nn(i1,i2)
    ds = ds + ds_nn(i1,i2)
  END DO

! Correct for A or T at the termini

  IF (ABS(CurrDNA%NUMseq(start)).eq.1) THEN
    dh=dh+2.2
    ds=ds+6.935
  END IF

  IF (ABS(CurrDNA%NUMseq(finish)).eq.1) THEN
    dh=dh+2.2
    ds=ds+6.935
  END IF

! Correct for self-complementarity

  inner1: DO i=start,finish
    IF ((CurrDNA%NUMseq(i)+CurrDNA%NUMseq(finish+start-i)).eq.0) THEN
      self_compl=.TRUE.
    ELSE
      self_compl=.FALSE.
      EXIT inner1
    END IF
  END DO inner1

  IF (self_compl) ds=ds-1.4

! Make corrections for oligo concentration

  IF (self_compl) THEN
    ds=ds+OligoCorrSC
  ELSE
    ds=ds+OligoCorr
  END IF

! Make corrections for cation concentrations

  ds=ds+(SaltCorr*(finish-start))

! Now find the actual Tm

  TmCalc=(1000*dh/ds)-Kelvin

END FUNCTION TmCalc
SUBROUTINE TmCorrect
!
! Create salt and oligo corrections for Tm
! Uses a 2D lookup table indexed by snapped sodium and magnesium concentrations.

  USE dnaworks_data
  USE dnaworks_test
  IMPLICIT NONE

  INTEGER :: si, mi, i
  REAL :: best_dist, dist

! Sodium levels (M): 21 values (fine step 0.005 from 0.01 to 0.075)
  REAL, SAVE :: na_lvl(21)
  DATA na_lvl /0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, &
               0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, &
               0.100, 0.150, 0.200, 0.250, 0.500, 0.750, 1.000/

! Mg levels (M): 16 values (fine step 0.0005 from 0.002 to 0.005)
  REAL, SAVE :: mg_lvl(16)
  DATA mg_lvl /0.0000, 0.0005, 0.0010, 0.0015, 0.0020, 0.0025, 0.0030, &
               0.0035, 0.0040, 0.0045, 0.0050, 0.0100, 0.0200, 0.0500, &
               0.1000, 0.2000/

! Salt correction table: salt_table(sodium_idx, mg_idx)
! Stored column-major (Fortran default).
! Each column is one Mg level; each row is one Na level.
! Intermediate values computed via bilinear interpolation from original data.
  REAL, SAVE :: salt_table(21,16)
  DATA salt_table / &
    -1.6960,-1.5831,-1.4703,-1.3574,-1.3065, &  ! Mg=0.0000
    -1.2555,-1.2046,-1.1536,-1.1027,-1.0730, &
    -1.0434,-1.0137,-0.9841,-0.9544,-0.8480, &
    -0.6964,-0.5933,-0.5094,-0.2547,-0.1064, &
     0.0000,                                  &
    -0.9125,-0.8921,-0.8716,-0.8512,-0.8351, &  ! Mg=0.0005
    -0.8190,-0.8028,-0.7867,-0.7706,-0.7564, &
    -0.7422,-0.7281,-0.7139,-0.6997,-0.6449, &
    -0.5514,-0.4772,-0.4159,-0.2031,-0.0709, &
     0.0258,                                  &
    -0.7996,-0.7846,-0.7695,-0.7545,-0.7410, &  ! Mg=0.0010
    -0.7274,-0.7139,-0.7003,-0.6868,-0.6758, &
    -0.6649,-0.6539,-0.6430,-0.6320,-0.5836, &
    -0.5030,-0.4385,-0.3805,-0.1838,-0.0580, &
     0.0355,                                  &
    -0.7287,-0.7158,-0.7029,-0.6900,-0.6790, &  ! Mg=0.0015
    -0.6681,-0.6571,-0.6462,-0.6352,-0.6255, &
    -0.6158,-0.6062,-0.5965,-0.5868,-0.5449, &
    -0.4707,-0.4095,-0.3579,-0.1709,-0.0484, &
     0.0451,                                  &
    -0.6803,-0.6696,-0.6588,-0.6481,-0.6378, &  ! Mg=0.0020
    -0.6275,-0.6171,-0.6068,-0.5965,-0.5881, &
    -0.5797,-0.5714,-0.5630,-0.5546,-0.5127, &
    -0.4449,-0.3901,-0.3386,-0.1612,-0.0387, &
     0.0516,                                  &
    -0.6449,-0.6352,-0.6255,-0.6159,-0.6062, &  ! Mg=0.0025
    -0.5965,-0.5868,-0.5772,-0.5675,-0.5598, &
    -0.5520,-0.5443,-0.5365,-0.5288,-0.4901, &
    -0.4256,-0.3724,-0.3241,-0.1516,-0.0323, &
     0.0565,                                  &
    -0.6094,-0.6008,-0.5922,-0.5836,-0.5746, &  ! Mg=0.0030
    -0.5656,-0.5565,-0.5475,-0.5385,-0.5314, &
    -0.5243,-0.5172,-0.5101,-0.5030,-0.4675, &
    -0.4063,-0.3547,-0.3095,-0.1419,-0.0258, &
     0.0613,                                  &
    -0.5836,-0.5755,-0.5675,-0.5594,-0.5510, &  ! Mg=0.0035
    -0.5426,-0.5343,-0.5259,-0.5175,-0.5107, &
    -0.5040,-0.4972,-0.4904,-0.4837,-0.4497, &
    -0.3917,-0.3418,-0.2983,-0.1338,-0.0209, &
     0.0661,                                  &
    -0.5578,-0.5503,-0.5427,-0.5352,-0.5275, &  ! Mg=0.0040
    -0.5197,-0.5120,-0.5042,-0.4965,-0.4901, &
    -0.4836,-0.4772,-0.4707,-0.4643,-0.4320, &
    -0.3772,-0.3289,-0.2870,-0.1258,-0.0161, &
     0.0709,                                  &
    -0.5384,-0.5309,-0.5234,-0.5159,-0.5088, &  ! Mg=0.0045
    -0.5017,-0.4946,-0.4875,-0.4804,-0.4739, &
    -0.4675,-0.4610,-0.4546,-0.4481,-0.4175, &
    -0.3643,-0.3192,-0.2773,-0.1194,-0.0113, &
     0.0741,                                  &
    -0.5191,-0.5116,-0.5040,-0.4965,-0.4901, &  ! Mg=0.0050
    -0.4836,-0.4772,-0.4707,-0.4643,-0.4578, &
    -0.4514,-0.4449,-0.4385,-0.4320,-0.4030, &
    -0.3514,-0.3095,-0.2676,-0.1129,-0.0064, &
     0.0774,                                  &
    -0.3966,-0.3912,-0.3859,-0.3805,-0.3753, &  ! Mg=0.0100
    -0.3702,-0.3650,-0.3599,-0.3547,-0.3502, &
    -0.3457,-0.3411,-0.3366,-0.3321,-0.3095, &
    -0.2708,-0.2354,-0.1999,-0.0677, 0.0290, &
     0.1064,                                  &
    -0.2741,-0.2698,-0.2655,-0.2612,-0.2573, &  ! Mg=0.0200
    -0.2534,-0.2496,-0.2457,-0.2418,-0.2386, &
    -0.2354,-0.2321,-0.2289,-0.2257,-0.2096, &
    -0.1773,-0.1483,-0.1225,-0.0129, 0.0709, &
     0.1419,                                  &
    -0.1064,-0.1043,-0.1021,-0.1000,-0.0974, &  ! Mg=0.0500
    -0.0948,-0.0923,-0.0897,-0.0871,-0.0852, &
    -0.0832,-0.0813,-0.0793,-0.0774,-0.0645, &
    -0.0451,-0.0226,-0.0032, 0.0774, 0.1451, &
     0.2031,                                  &
     0.0193, 0.0204, 0.0215, 0.0226, 0.0245, &  ! Mg=0.1000
     0.0264, 0.0284, 0.0303, 0.0322, 0.0341, &
     0.0361, 0.0380, 0.0400, 0.0419, 0.0484, &
     0.0645, 0.0806, 0.0935, 0.1612, 0.2160, &
     0.2644,                                  &
     0.1451, 0.1462, 0.1472, 0.1483, 0.1496, &  ! Mg=0.2000
     0.1509, 0.1522, 0.1535, 0.1548, 0.1561, &
     0.1574, 0.1586, 0.1599, 0.1612, 0.1677, &
     0.1805, 0.1902, 0.1999, 0.2515, 0.2934, &
     0.3353                                   /

  IF (TEST0) PRINT *,"TmCorrect" !TEST0

! Clamp concentrations to valid ranges

  IF (OligoConc.lt.1e-9) OligoConc=1e-9
  IF (OligoConc.gt.1e-4) OligoConc=1e-4
  IF (SodiumConc.gt.1.000) SodiumConc=1.000
  IF (SodiumConc.lt.0.010) SodiumConc=0.010
  IF (MgConc.gt.0.2000) MgConc=0.2000
  IF (MgConc.lt.0.0) MgConc=0.0

! Compute oligo corrections

  OligoCorr=RGasConstant*(LOG(((OligoConc/100)/2)))
  OligoCorrSC=RGasConstant*(LOG((OligoConc/100)))

! Compute salt correction

  IF (SLOW) THEN
!   Exact calculation: SantaLucia (1998) monovalent correction with
!   von Ahsen et al. (2001) Mg2+ equivalent ionic strength.
!   SaltCorr = 0.368 * ln([Na+] + 3.3 * sqrt([Mg2+]))
    SaltCorr = 0.368 * LOG(SodiumConc + 3.3 * SQRT(MgConc))
  ELSE
!   Snap to nearest table level and look up
    si = 1
    best_dist = ABS(SodiumConc - na_lvl(1))
    DO i = 2, 21
      dist = ABS(SodiumConc - na_lvl(i))
      IF (dist .lt. best_dist) THEN
        best_dist = dist
        si = i
      END IF
    END DO
    SodiumConc = na_lvl(si)

    mi = 1
    best_dist = ABS(MgConc - mg_lvl(1))
    DO i = 2, 16
      dist = ABS(MgConc - mg_lvl(i))
      IF (dist .lt. best_dist) THEN
        best_dist = dist
        mi = i
      END IF
    END DO
    MgConc = mg_lvl(mi)

    SaltCorr = salt_table(si, mi)
  END IF

END SUBROUTINE TmCorrect
