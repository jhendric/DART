	subroutine seth2onm(mois,nzonal)
C_________________________________________________________________________

C  This subroutine reads in arrays of H2O forcing provided by h2osub.f
C	which is based upon the Groves (JATP, 1982) parameterization.
C            35 latitudes: -85 to +85 in steps of 5 degrees
C    41 (irregularly spaced) altitudes: X=0--->10. in steps of dx=1.
C	Units are: J/m3/s; converted from J/kg/s in ~hagan/testlat/testh2o.f
C  		arrays used herein are subsequently created in
C			~hagan/testlat/plotout/h2onointerp.pro
C
C  The heating rates (J/m3/s) are stored in the array dh2o(37,41).
C  and the corresponding array of altitudes in zh2o(37,41)
C   N.B. Values at the poles are identical to those at +/-85
c
C changes to accommodate complex forcing and multiple wavenumbers
C					M. Hagan (3/9/95)
C changes to accommodate seasonal variations 
C					M. Hagan (5/16/95)
C_________________________________________________________________________

        dimension zx1(41),zxm(41),zy1(37),zyn(37),temp(119)
        dimension zx1a(41),zxma(41),zy1a(37),zyna(37),tempa(119)

	character*17 dchoice(4)
	character*17 dachoice(4)
	character*17 djchoice(4)
	character*17 dochoice(4)
	character*17 dfile

        real xh2o(41), xlih2o(35)
	complex dch2o(37,41),dch2oi(35,41)

        common/hth2o/zh2o(41),xlh2o(37),dh2o(37,41),zdh2o(37,41,3),
     +	sigma2,dah2o(37,41),zadh2o(37,41,3),sigma2a

        data dchoice/'h2o_jan.eastward1','h2o_jan.standing0',
     +             'h2o_jan.westward2','h2o_jan.westward3'/

        data dachoice/'h2o_apr.eastward1','h2o_apr.standing0',
     +             'h2o_apr.westward2','h2o_apr.westward3'/

        data djchoice/'h2o_jul.eastward1','h2o_jul.standing0',
     +             'h2o_jul.westward2','h2o_jul.westward3'/

        data dochoice/'h2o_oct.eastward1','h2o_oct.standing0',
     +             'h2o_oct.westward2','h2o_oct.westward3'/


       icount=(mois-1)/3+1
       print * ,"In SR SETH2O icount=",icount
       print * ," nzonal=",nzonal
       print * ," and mois=",mois

       if(mois.eq.1.and.nzonal.lt.1) dfile=dchoice(nzonal+2)
       if(mois.eq.1.and.nzonal.ge.1) dfile=dchoice(nzonal+1)
       if(mois.eq.4.and.nzonal.lt.1) dfile=dachoice(nzonal+2)
       if(mois.eq.4.and.nzonal.ge.1) dfile=dachoice(nzonal+1)
       if(mois.eq.7.and.nzonal.lt.1) dfile=djchoice(nzonal+2)
       if(mois.eq.7.and.nzonal.ge.1) dfile=djchoice(nzonal+1)
       if(mois.eq.10.and.nzonal.lt.1) dfile=dochoice(nzonal+2)
       if(mois.eq.10.and.nzonal.ge.1) dfile=dochoice(nzonal+1)

	pi=acos(-1.)
	pio2=pi/2.
        dtor=pi/180.

C  Read in array of latitudes defining H2O grid: -85 to +85 in 5 deg
C  Convert the values to colatitude
C  Write another array in reverse order for print (consistent with SETATMOS)

	OPEN(UNIT=63,file=dfile,form='formatted')

201	format(5(5X ,7(f10.5),/))
222	format(2(1x,f5.2))
202	format(5(5X ,3((f10.5,f10.5)),/))
203	format(/)

	read(63,203)
	read(63,201) (xlih2o(i),i=1,35)


C  Reverse the order of  the array (+80 - -80 - consistent setatmos)
C  *******need colat to be increasing for SURF1/SURFD********
 	do 2 j=1,35
 	il=37-j
 	xlh2o(il)=xlih2o(j)
    2  continue
	xlh2o(1)=90.
	xlh2o(37)=-90.
	print '(37(1x,f4.0))',(xlh2o(i),i=1,37)

C convert to colatitude:
	do 22 j=1,37
        xlh2o(j)=pio2-xlh2o(j)*dtor    
   22  continue

C  Read in x altitude and corresponding array of h2o forcing at 37 latitudes:

	do 3 j=1,41    
	read(63,222) xh2o(j),zh2o(j)
	read(63,202) (dch2oi(i,j),i=1,35)
    3  continue

C  Reverse the order of  the array (+90 - -90 - consistent setatmos)
C  *******need colat to be increasing for SURF1/SURFD********
 	do 34 j=1,41
 	do 33 i=1,35
 	ii=37-i
 	dh2o(ii,j)=real(dch2oi(i,j))
 	dah2o(ii,j)=aimag(dch2oi(i,j))
   33   continue
	dh2o(1,j)=dh2o(2,j)
	dh2o(37,j)=dh2o(36,j)
	dah2o(1,j)=dah2o(2,j)
	dah2o(37,j)=dah2o(36,j)
   34   continue

	do 36 j=1,41
	do 35 i=1,37
 	dch2o(i,j)=cmplx(dh2o(i,j),dah2o(i,j))
   35   continue
c	write(15,222) xh2o(j),zh2o(j)
c	write(15,*) (dch2o(i,j),i=1,37)
   36   continue

	CLOSE(UNIT=63)


C___________________________________________________________________________
C  Set up arrays for later interpolation in heat92:
C

        sigma2=1.

	call surf1(37,41,xlh2o,zh2o,dh2o,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zdh2o,temp(1),sigma2,ierr)

        sigma2a=1.

	call surf1(37,41,xlh2o,zh2o,dah2o,37,zx1a,zxma,zy1a,zyna,zxy11a,
     +  zxym1a,zxy1na,zxymna,255,zadh2o,tempa(1),sigma2a,ierr)

      return
      end

 
