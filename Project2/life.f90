      SUBROUTINE DISPLAY(STATUS,N)
      DIMENSION STATUS(35,35)
      CHARACTER*70 STRING
      INTEGER neighborcount
      
      DO 7 I=1,N
      STRING=' '
      DO 14 J=1,N

      IF (STATUS(I,J).EQ.1.) THEN
      STRING=TRIM(STRING)//' @ '
      ELSE 
      STRING=TRIM(STRING)//' . '
      ENDIF   
      
   14 CONTINUE
      WRITE(*,*)STRING
    7 CONTINUE
! *** UNFORTUNATELY THIS IS QUANTIZED IN UNITS OF 1 SEC
!      CALL SLEEP(1)
! *** ANOTHER WAY TO DO THIS (WASTE SOME TIME)
      DO 17 J=1,9000000
      RN1=RAND()
   17 CONTINUE
      RETURN
      END
! **************************************************
      PROGRAM LIFE
      DIMENSION STATUS(35,35),NEXT(35,35)
      
      N=35
      
      WRITE(*,*)"ENTER 1 FOR MANUAL SETUP, 0 FOR RANDOM"
      READ(*,*)RAN
      
! *** SET UP THE BOARD
      DO 7 I=1,N
      DO 14 J=1,N
      STATUS(I,J)=0.
      NEXT(I,J)=0.
   14 CONTINUE
    7 CONTINUE

      IF(RAN.EQ.1.)THEN     
! *** INITIALIZE LIVE CELLS
   15 WRITE(*,*)'ENTER COORDINATES OF LIVE CELLS, 0,0 WHEN DONE'
      READ(*,*)I,J
      IF((I.EQ.0).AND.(J.EQ.0))GOTO 17
      STATUS(I,J)=1.
! *** DISPLAY MOVE
      CALL DISPLAY(STATUS,N)
      GOTO 15
      ELSE

      WRITE(*,*)"ENTER A FIVE DIGIT INTEGER TO CHANGE SEED"
      READ(*,*)NSEED
      CALL SRAND(NSEED)
         
      DO 57 I=1,N
      DO 54 J=1,N
      RN1=RAND()
      IF(RN1.LT.0.75)THEN
      STATUS(I,J)=0.
      ELSE
      STATUS(I,J)=1.
      ENDIF
   54 CONTINUE
   57 CONTINUE
    
      ENDIF
     
 17   CONTINUE
      
! *** MAKE IT ALIVE!!! ***
! *** SEE http://www.bitstorm.org/gameoflife/
! *** RULES:
! *** For a space that is 'populated':
! *** Each cell with one or no neighbors dies, as if by loneliness. 
! *** Each cell with four or more neighbors dies, as if by overpopulation. 
! *** Each cell with two or three neighbors survives. 
! *** For a space that is 'empty' or 'unpopulated':
! *** Each cell with three neighbors becomes populated. 
      
      WRITE(*,*)"ADJUST YOUR SCREEN TO BE 36 LINES LONG"



      DO WHILE (.true.)

      CALL SLEEP(1)
      

      DO 64 I=1,N
            DO 71 J=1,N

! Now make if statements for a given I,J pair
! First, we can make temp variable named neighbor and check number of neighbors
                  neighborcount = 0
!                 LHS Check
                  IF (J.GT.1) THEN
                        neighborcount = neighborcount + STATUS(I,J-1)
                        IF (I.GT.1) THEN
                              neighborcount = neighborcount + STATUS(I-1,J-1)
                              ENDIF
                        IF (I.LT.N) THEN
                              neighborcount = neighborcount + STATUS(I+1,J-1)
                              ENDIF
                  ENDIF
!                 Middle Check
                  IF (I.GT.1) THEN
                        neighborcount = neighborcount + STATUS(I-1,J)
                        ENDIF
                  IF (I.LT.N) THEN
                        neighborcount = neighborcount + STATUS(I+1,J)
                        ENDIF
!                 RHS Check
                  IF (J.LT.N) THEN
                        neighborcount = neighborcount + STATUS(I,J+1)
                        IF (I.GT.1) THEN
                              neighborcount = neighborcount + STATUS(I-1,J+1)
                              ENDIF
                        IF (I.LT.N) THEN
                              neighborcount = neighborcount + STATUS(I+1,J+1)
                              ENDIF
                  ENDIF


! Now we can directly implement the rules:

                  IF ( (STATUS(I,J)==0).AND.(neighborcount==3) ) THEN
                        NEXT(I,J)= 1
                  ENDIF

                  IF ( (STATUS(I,J) == 1).AND.(neighborcount.LE.1).OR.(STATUS(I,J) == 1).AND.(neighborcount.GE.4) ) THEN
                        NEXT(I,J) = 0
                  ELSE IF ( (STATUS(I,J) == 1).AND.(neighborcount == 2).OR.(STATUS(I,J) == 1).AND.(neighborcount == 3) ) THEN
                        NEXT(I,J) = 1
                  ENDIF


      71    CONTINUE
      64 CONTINUE
      STATUS = NEXT
      WRITE(*,*) !Empty line to separate each grid
      CALL DISPLAY(STATUS,N)

      END DO
      
      STOP 
      END

