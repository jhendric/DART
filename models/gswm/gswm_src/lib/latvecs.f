      SUBROUTINE LATVECS(IJ)

C changes to construct LATVECT arrays from the north to south pole
C	consistent with changes to ATMOS5	--M. Hagan (12/11/92)


      COMMON/LATVECT/DTHETA,CLATT(91),XLAT(91),SNCLAT(91),CSCLAT(91),
     1    TNCLAT(91),CTNCLAT(91),SIN2I(91),GMLAT(91),DIP(91)

      IEQ = 1 + ((IJ-1)/2)
      PI=ACOS(-1.)
      RADEG=180./PI
      COLAT = 0.0
      DTHETA = PI/FLOAT(IJ+1)

      DO 1 I=1,IJ

      COLAT = COLAT + DTHETA
      CLATT(I) = COLAT
      XLAT(I) = PI*.5-COLAT
      SNCLAT(I) = SIN(COLAT)
      CSCLAT(I) = COS(COLAT)
      DEGLAT = XLAT(I)*RADEG

C     GMLAT IS AN ARRAY OF GEOMAGNETIC LATITUDES FOR A GIVEN GEOGRAPHIC
C     LONGITUDE, HERE TAKEN TO BE -75 (75 DEG WEST) OR +285 DEG EAST
      CALL GEOGM(DEGLAT,285.,GMLT,GMLG)
      GMLAT(I) = GMLT/RADEG

      IF(I.EQ.IEQ) GO TO 210
      TNCLAT(I) = TAN(COLAT)
      TLAT=1./TNCLAT(I)
      GO TO 200
210   TLAT=0.0
      TNCLAT(I) = 0.0
200   CONTINUE
      CTNCLAT(I) = TLAT
      DIP(I) = ABS(ATAN(2.*TLAT))
      SIN2I(I) = SIN(DIP(I))*SIN(DIP(I))

1     CONTINUE

      RETURN
      END
