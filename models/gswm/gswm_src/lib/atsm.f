C        SUBROUTINE ATSM
C
C        PURPOSE
C           NDIM POINTS OF A GIVEN TABLE WITH MONOTONIC ARGUMENTS ARE
C           SELECTED AND ORDERED SUCH THAT
C           ABS(ARG(I)-X).GE.ABS(ARG(J)-X) IF I.GT.J.

C
C        USAGE
C           CALL ATSM (X,Z,F,IROW,ICOL,ARG,VAL,NDIM)
C
C        DESCRIPTION OF PARAMETERS
C           X      - THE SEARCH ARGUMENT.
C           Z      - THE VECTOR OF ARGUMENT VALUES (DIMENSION IROW).
C                    THE ARGUMENT VALUES MUST BE STORED IN INCREASING
C                    OR DECREASING SEQUENCE.
C           F      - IN CASE ICOL=1, F IS THE VECTOR OF FUNCTION VALUES
C                    (DIMENSION IROW).
C                    IN CASE ICOL=2, F IS AN IROW BY 2 MATRIX. THE FIRST
C                    COLUMN SPECIFIES THE VECTOR OF FUNCTION VALUES AND
C                    THE SECOND THE VECTOR OF DERIVATIVES.
C           IROW   - THE DIMENSION OF VECTOR Z AND OF EACH COLUMN
C                    IN MATRIX F.
C           ICOL   - THE NUMBER OF COLUMNS IN F (I.E. 1 OR 2).
C           ARG    - THE RESULTING VECTOR OF SELECTED AND ORDERED
C                    ARGUMENT VALUES (DIMENSION NDIM).
C           VAL    - THE RESULTING VECTOR OF SELECTED FUNCTION VALUES
C                    (DIMENSION NDIM) IN CASE ICOL=1. IN CASE ICOL=2,
C                    VAL IS THE VECTOR OF FUNCTION AND DERIVATIVE VALUES
C                    (DIMENSION 2*NDIM) WHICH ARE STORED IN PAIRS (I.E.
C                    EACH FUNCTION VALUE IS FOLLOWED BY ITS DERIVATIVE
C                    VALUE).
C           NDIM   - THE NUMBER OF POINTS WHICH MUST BE SELECTED OUT OF
C                    THE GIVEN TABLE (Z,F).
C
C        REMARKS
C           NO ACTION IN CASE IROW LESS THAN 1.
C           IF INPUT VALUE NDIM IS GREATER THAN IROW, THE PROGRAM
C           SELECTS ONLY A MAXIMUM TABLE OF IROW POINTS.  THEREFORE THE
C           USER OUGHT TO CHECK CORRESPONDENCE BETWEEN TABLE (ARG,VAL)
C           AND ITS DIMENSION BY COMPARISON OF NDIM AND IROW, IN ORDER
C           TO GET CORRECT RESULTS IN FURTHER WORK WITH TABLE (ARG,VAL).
C           THIS TEST MAY BE DONE BEFORE OR AFTER CALLING
C           SUBROUTINE ATSM.
C           SUBROUTINE ATSM ESPECIALLY CAN BE USED FOR GENERATING THE
C           TABLE (ARG,VAL) NEEDED IN SUBROUTINES ALI, AHI, AND ACFI.
C

C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           SELECTION IS DONE BY SEARCHING THE SUBSCRIPT J OF THAT
C           ARGUMENT, WHICH IS NEXT TO X (BINARY SEARCH).
C           AFTERWARDS NEIGHBOURING ARGUMENT VALUES ARE TESTED AND
C           SELECTED IN THE ABOVE SENSE.
C
C     ..................................................................
C
      SUBROUTINE ATSM(X,Z,F,IROW,ICOL,ARG,VAL,NDIM)
C
C
      DIMENSION Z(*),F(*),ARG(*),VAL(*)
C
C     CASE IROW=1 IS CHECKED OUT
      IF(IROW-1)23,21,1
    1 N=NDIM
C
C     IF N IS GREATER THAN IROW, N IS SET EQUAL TO IROW.
      IF(N-IROW)3,3,2
    2 N=IROW
C
C     CASE IROW.GE.2
C     SEARCHING FOR SUBSCRIPT J SUCH THAT Z(J) IS NEXT TO X.
    3 IF(Z(IROW)-Z(1))5,4,4
    4 J=IROW
      I=1
      GOTO 6
    5 I=IROW
      J=1
    6 K=(J+I)/2
      IF(X-Z(K))7,7,8
    7 J=K
      GOTO 9
    8 I=K
    9 IF(IABS(J-I)-1)10,10,6
   10 IF(ABS(Z(J)-X)-ABS(Z(I)-X))12,12,11
   11 J=I

C
C     TABLE SELECTION
   12 K=J
      JL=0
      JR=0
      DO 20 I=1,N
      ARG(I)=Z(K)
      IF(ICOL-1)14,14,13
   13 VAL(2*I-1)=F(K)
      KK=K+IROW
      VAL(2*I)=F(KK)
      GOTO 15
   14 VAL(I)=F(K)
   15 JJR=J+JR
      IF(JJR-IROW)16,18,18
   16 JJL=J-JL
      IF(JJL-1)19,19,17
   17 IF(ABS(Z(JJR+1)-X)-ABS(Z(JJL-1)-X))19,19,18
   18 JL=JL+1
      K=J-JL
      GOTO 20
   19 JR=JR+1
      K=J+JR
   20 CONTINUE
      RETURN
C
C     CASE IROW=1
   21 ARG(1)=Z(1)
      VAL(1)=F(1)
      IF(ICOL-2)23,22,23
   22 VAL(2)=F(2)
   23 RETURN
      END

