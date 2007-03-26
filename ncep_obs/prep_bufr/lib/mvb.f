C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MVB(IB1,NB1,IB2,NB2,NBM)                               
                                                                        
      DIMENSION IB1(*),IB2(*),NVAL(24000)
                                                                        
C-----------------------------------------------------------------------
CFPP$ EXPAND(UPB,PKB)                                                   
C-----------------------------------------------------------------------
                                                                        
      IF(NBM.GT.24000) CALL BORT('MVB - NBM>24000')
      JB1 = 8*(NB1-1)                                                   
      JB2 = 8*(NB2-1)                                                   
                                                                        
      DO N=1,NBM                                                        
      CALL UPB(NVAL(N),8,IB1,JB1)
      ENDDO

      DO N=1,NBM
      CALL PKB(NVAL(N),8,IB2,JB2)
      ENDDO                                                             
                                                                        
      RETURN                                                            
      END                                                               
