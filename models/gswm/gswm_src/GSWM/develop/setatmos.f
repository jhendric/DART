	subroutine setatmos
c Modified to include HWM93 wind option		M. Hagan (10/20/98)
c
	real zx1(101),zxm(101),zy1(37),zyn(37),temp(239)
        real zxy11,zxym1,zxy1n,zxymn
c
	integer igm,ihd,iyp,ibr,nowind,ihm,latgrad,polebc
	integer wback
c
	dimension d(8),tt(2),w(2),sw(25),f(37),df(37),po(37)
	dimension ih(19)
	dimension ap(1),it(19),iu(19),ifd(19),il(19)
C	dimension zx1(101),zxm(101),zy1(37),zyn(37),temp(239)
	dimension ubr(19)

	character*10 gchoice(12)
	character*10 gfile
	character*13  hchoice(12)
	character*13  hfile
	character*10 rchoice(2)
	character*10  rfile

      common/interp1/zpt(37,101,3),zpu(37,101,3),zpr(37,101,3),
     +zpxh(37,101,3),do2(37,101,3)
      common/interp2/x(37),y(101),zt(37,101),zu(37,101),zr(37,101),
     +zxh(37,101),zph(37,101),o2(37,101),sigma
	COMMON/MODE/NZONAL,PERIOD,FREQ,MOIS,NSS,FLUX

c  The following common block will contain control flags for the
c  background atmosphere.  Put it in the relevant subroutines


	common /backatm/igm,ihd,iyp,ibr,nowind,ihm,latgrad,polebc
C
	common /write/wback
C
	data ap/4./
	data (sw(i),i=1,25)/25*1./
	data rs,bc,poo,go/8.31432,1.38054e-23,1.01325e+05,9.80665/
	data re,omega/6356766.,.72722e-04/
C
	data gchoice/'groves_jan','groves_jan','groves_mar',
     +             'groves_apr','groves_may','groves_jun',
     +             'groves_jul','groves_aug','groves_sep',
     +             'groves_oct','groves_nov','groves_dec'/
C
	data hchoice/'hrdi0-128_jan','hrdi0-128_feb','hrdi0-128_mar',
     +             'hrdi0-128_apr','hrdi0-128_may','hrdi0-128_jun',
     +             'hrdi0-128_jul','hrdi0-128_aug','hrdi0-128_sep',
     +             'hrdi0-128_oct','hrdi0-128_nov','hrdi0-128_dec'/
C
	data rchoice/'randel_jan','randel_apr'/
C*************************** OPTIONS ************************************     

C       default = parameters (below) = zero
C               = geostrophic winds calculated from MSISE90 up to
C                 108 km and from pole to pole (u=zero at poles)

C   igm = 1  read in Groves/MSIS winds and overwrite up to 108 km

C   igm = 2  read in Groves/MSIS winds and overwrite only troposphere/
C            lower stratosphere winds

C            Note: for igm.neq.0 , the correct Groves/MSIS data file
C                  must be included in OPEN statement and tarsrc file

C   ihd = 1  overwrite mean winds with interpolated/extrapolated HRDI data
C            52-96 km for appropriate month

C   ihd = 2  overwrite mean winds with interpolated/extrapolated HRDI data
C            0-128 km for VALIDATED FOR APRIL ONLY (M. Hagan 2/26/98)

C   ihm = 1  HWM93 mean winds at all altitudes
C
C   ihm = 2  HRDI/Groves winds below 128; HWM93 mean winds above (EUV standard)
C
C   iyp = 1  overwrite with Portnyagin winds above 80 km with smooth
C            merger between 72 and 80 km

C*************************** OPTIONS ************************************     

C  GSWM standard for tides:

C Filename with the Groves information

	gfile=gchoice(mois)

c Filename with the hrdi information

	hfile=hchoice(mois)

C Filename with the Randel (NMC) data:

	if(ibr.eq.1)then
	   if(mois.eq.1)then
	      rfile=rchoice(1)
	   elseif(mois.eq.4)then
	      rfile=rchoice(2)
	   else
	      print*,"Error finding NMC wind file. Stop."
	      stop
	   endif
	endif
c
	pi=acos(-1.)

	iyd=15+30*(mois-1)

C_________________________________________________________________________

C  This subroutine sets up arrays which ATMOS5 will later use for interpolation
C  and calculation of various partial derivatives.  The arrays set up here are:

C		     zt = zonal mean temp
C		     zu = mean zonal wind
C                    zr = log mass density
C 		     zxh = pressure height

C_________________________________________________________________________
c
c Do loop for all latitude: colat 0 to PI in radians

	do 2 i=1,37
	   colat=float(i-1)*5.*pi/180.
	   x(i)=colat
c
c Do loop for all altitude 0-400km

	   do 2 j=1,101
	      z=float(j-1)*4.
	      y(j)=z
	      
	      deglat=90.-float(i-1)*5.

	      sw(7)=0.
	      sw(8)=0.
	      sw(14)=0.
c
c New findings suggest using SW(10)=0. JKH 4/29/98

	      sw(10)=0.
	      call tselec(sw)
	      call wtselec(sw)

	      call gtd6(iyd,0.,z,deglat,285.,12.,flux,flux,ap,48,d,tt)

C Calculate # density and convert from cm-3 to m-3

	      xn=(d(1)+d(2)+d(3)+d(4)+d(5)+d(7))*1.0e+06
c
c Initialize pressure at ground for each latitude

	      if(j.gt.1) go to 20

	      po(i)=xn*bc*tt(2)
 20	      p=xn*bc*tt(2)
c
c Convert pressure from Pa to mbar for printout only

	      zph(i,j)=p/100.
	      zxh(i,j)=-alog(p/po(i))
c
c Calculate ln mass density and convert from g/cm3 to kg/m3

	      zr(i,j)=alog(d(6)*1.0e+03)
	      zt(i,j)=tt(2)
	      if(d(4).gt.0.) then
		 o2(i,j)=alog10(d(4))
	      else
		 o2(i,j)=0.0
	      endif
c
c End of both latitude and altitude loops

 2	   continue

CTEST 7/22/98 write out zph for use in interpolating some wind pressure alts
c	do j=1,18
c	   write(*,1111)(zph(i,j),i=1,37,2)
c	end do
c 1111	format(19(f6.1))
c	return
cchange cira_jan again
CTEST 7/22/98 end

C_________________________________________________________________________

C  Calculate geostrophic wind up to 108 km with exponential tail-off above
C  108 km; 108 km is chosen so that tables from the Portnyagin or other
C  observationally-based model can be inserted in place of geostrophic
C  winds above 80 km.  In this case some scheme is/should be introduced
C  to introduce weighted averages at 72 km and 76 km to effect a smooth
C  transition from msise90 at 68 km to the observationally-based model
C  at 80 km.

C ...Or else enter zeroes for background winds (nowind=1)

c
c jkh 9/4/98 initialize zonal mean wind to zero to 400km

	do j=1,101
	   do i=1,37
	      zu(i,j)=0.0
	   enddo
	enddo

C  Below 108 km:

      dtheta=5.*pi/180.
      omega=2.*pi/(24.*3600.)

C  Latitudes away from equator and poles: Altitude 0 to 108km

	do  j=1,28
	  z=float(j-1)*4.

	  do  i=1,37
	     f(i)=zxh(i,j)
 	  end do

          call det3(dtheta,f,df,37,ier)

C If the "NO WIND" flag is set to 0 (turned off)

	  if(nowind.eq.0)then

C North half

	     do  i=4,16
		colat=float(i-1)*5.*pi/180.
		coeff=-1./(2.*omega*re*cos(colat))
		zu(i,j)=coeff*df(i)*po(i)*exp(-zxh(i,j)-zr(i,j))
	     end do		!i loop of North colatitude (15-75 deg)

C South half

	     do  i=22,34
		colat=float(i-1)*5.*pi/180.
		coeff=-1./(2.*omega*re*cos(colat))
		zu(i,j)=coeff*df(i)*po(i)*exp(-zxh(i,j)-zr(i,j))
	     end do		!i loop of South colatitude (105-165 deg)

	  endif			!if NOWIND = 0

	end do			!j loop of altitude (0-108 km)

C Equatorial points; linear interpolation between +15,-15 degrees latitude

	do j=1,28

C If the "NO WIND" flag is set to 0 (turned off)

	   if(nowind.eq.0)then

	      do  i=17,21
		 zu(i,j)=zu(16,j)+float(i-16)*(zu(22,j)-zu(16,j))/6.
	      end do

	   elseif(nowind.eq.1)then !all background zonal winds =0
	      do i=1,37
		 zu(i,j)=0.0
	      end do

	   endif		!nowind=0

	end do			! altitude 0-108km


C_________________________________________________________________________

C  Read in GROVES-MSIS winds below 108 km.  
 	
c If user wants nonzero background wind

	if(nowind.eq.0)then

	   if(igm.ne.0) then

C Read in Groves/MSISE90 background winds to overwrite geostrophic winds

	      OPEN(UNIT=13,file=gfile,form='formatted',STATUS='old')

 202	      format(4x,19i4)
 203	      format(/)

	      read(13,203)
	      j=29
	      do 204 jj=1,28	!altitudes from 0 to 108km
		 j=j-1		!top to bottom
		 read(13,202) (iu(i),i=1,19)
		 i=1
		 do 206 k=1,37,2 !all latitudes, every 10 deg

C  Overwrite with Groves/MSIS winds at all altitudes:

		    if(igm.eq.1)then
		       zu(k,j)=float(iu(i))

C  Overwrite with Groves/MSIS winds to 20 km, smoothing

		    else	!igm=2
		       if(j.eq.6) zu(k,j)=.66*zu(k,j)+.33*float(iu(i))
		       if(j.eq.5) zu(k,j)=.33*zu(k,j)+.66*float(iu(i))
		       if(j.le.4) zu(k,j)=float(iu(i))

		    endif	!igm=1
		    i=i+1

 206		 continue	!all latitudes, k loop

C Fill in the 5 degree latitude grid

		 do 205 i=2,36,2
		    zu(i,j)=(zu(i+1,j)+zu(i-1,j))*.5
 205		 continue

 204	      continue
	      CLOSE(UNIT=13)
	   endif		!igm.ne.0
C_________________________________________________________________________

C  If the flag, iyp, is not zero, introduce Portnyagin zonal mean model
C  between 80 and 108 km, with
C  smoothed transition between 70 and 80 km:

	   if(iyp.ne.0)then

	      do j=19,28
		 z=float(j-1)*4.

		 do i=1,37
		    xlat=(90.-float(i-1)*5.)

		    call wind2(xlat,z,uu)
		    if(j.eq.19)then !smoothing at 72km
		       zu(i,19)=.66*zu(i,19)+.33*uu
		    elseif(j.eq.20)then	!smoothing at 76km
		       zu(i,20)=.33*zu(i,20)+.66*uu
		    else
		       zu(i,j)=uu
		    endif
		    
		 end do		!i=latitude loop

	      end do		!j=altitude loop

	   endif		!iyp.ne.0

C_________________________________________________________________________

C  Read in HRDI winds between 52 and 96 km (12 altitudes) if ihd=1,
C  or between 0 and 128 km (33 altitudes) if ihd=2.  New HRDI background
C  files contain groves to 12 km--thus, set igm=0 to run ihd=2.


	   if(ihd.ne.0)then

	      OPEN(UNIT=19,file=hfile,form='formatted',STATUS='old')
	      read(19,903)	!skip top two lines

	      if(ihd.eq.1)then	!traditional hrdi, 52-96km
		 jin=14
		 jend=25
		 do j=1,jin-1	!skip down to 52km
		    read(19,*)
		 end do
	      elseif(ihd.eq.2)then !new hrdi, whole background to 128km
		 jin=1
		 jend=33
	      endif
c
 902	      format(4x,19i4)
 903	      format(/)

C HRDI data on GSWM background grid 52km-14th height; 96km-25th height;
C 128km-33rd height:
C delta z=4km, delta lat=10 degrees (interpolated to std. 5 deg. later)
c
c Read in winds according to user's preference
c
	      do 906 j=jin,jend
		 read(19,902) (iu(i),i=1,19)
		 i=1
c
c Overwrite HRDI winds into background
c
		 do k=1,37,2	!latitude every 10 degrees
		    zu(k,j)=float(iu(i))
		    i=i+1
		 end do
 906	      continue
c
c Smooth above and below if ihd=1; else all background is hrdi (ihd=2)
c
	      if(ihd.eq.1)then
		 do 907 k=1,37,2 !smooth at top (100-104 km) and bottom(44-48)
		    zu(k,13)=.33*zu(k,13)+.66*zu(k,14)
		    zu(k,12)=.66*zu(k,12)+.33*zu(k,14)
		    zu(k,26)=.33*zu(k,26)+.66*zu(k,25)
		    zu(k,27)=.66*zu(k,27)+.33*zu(k,25)
		    jend=jend+2
 907		 continue
	      endif
c
c Interpolate to the GSWM grid:  to 128 if ihd=2, to 104 if ihd=1
c
	      do 904 j=jin,jend
		 do 905 i=2,36,2
		    zu(i,j)=(zu(i+1,j)+zu(i-1,j))*.5
 905		 continue
 904	      continue

	      CLOSE(UNIT=19)
	      
	   endif		!if ihd.ne.0

C_________________________________________________________________________

C  Read in NMC (Bill Randel) winds between 12(4th height) and
C       48 km(13th height)

	   if(ibr.ne.0) then

	      OPEN(UNIT=23,file=rfile,form='formatted',STATUS='old')
	      read(23,*)
	      read(23,*)
	      read(23,*)
	      read(23,*)
	      read(23,*)

 302	      format(5x,19(f4.0,1x))

	      do jj=1,10
		 j=jj+3
		 read(23,302) (ubr(i),i=1,19)

C  Overwrite with NMC winds in stratosphere; smooth merger between 4 - 12km
C                               and between 48 - 56 km
		 i=1
		 do k=1,37,2
		    zu(k,j)=(ubr(i))
		    i=i+1
		 end do		!i
	      end do		!jj

	      i=1
	      do 307 k=1,37,2
		 zu(k,14)=.33*zu(k,14)+.66*zu(k,13)
		 zu(k,15)=.66*zu(k,15)+.33*zu(k,13)
		 zu(k,3)=.33*zu(k,3)+.66*zu(k,4)
		 zu(k,2)=.66*zu(k,2)+.33*zu(k,4)
		 i=i+1
 307	      continue

	      do 304 j=2,15
		 do 305 i=2,36,2
		    zu(i,j)=(zu(i+1,j)+zu(i-1,j))*.5
 305		 continue
 304	      continue
	      CLOSE(UNIT=23)
	   endif		!if ibr.ne.0

	endif			!if nowind ne 0
C_________________________________________________________________________


C   Points above 70 degrees; smooth transition to zero at pole

	do 51 j=1,28

	   xlatn=70.*pi/180.
	   do 8 i=1,4
	      xlat=(90.-float(i-1)*5.)*pi/180.
	      zu(i,j)=zu(5,j)*cos(xlat)/cos(xlatn)
 8	      zu(38-i,j)=zu(33,j)*cos(xlat)/cos(xlatn)

 51	   continue


C_________________________________________________________________________

C  Above 108 km:

	   if(ihd.eq.1)then
	      do 9 j=jend+2,101	!29,101
		 z=float(j-1)*4.

		 do 9 i=1,37
		    zu(i,j)=zu(i,28)*exp((108.-z)/10.)
 9		 continue
	   endif
 
C_________________________________________________________________________

C  If the flag, ihm, is not zero, introduce HWM93 zonal mean model

	   if(ihm.ne.0)then

	      if(ihm.eq.1)then	!overwrite entire wind field with HWM93
		 hin=1
		 hend=101
	      elseif(ihd.eq.2)then !retain HRDI background to 120km
		 hin=32
		 hend=101
	      endif
	      do j=hin,hend
		 z=float(j-1)*4.

		 do i=1,37
	            deglat=90.-float(i-1)*5.

	      	    call gws5(iyd,0.,z,deglat,285.,12.,flux,flux,ap,w)
	            if (ihd.eq.2)then !transition from HRDI at 124-128 km
		       if(j.eq.32) zu(i,j)=.66*zu(i,j)+.33*w(2)
		       if(j.eq.33) zu(i,j)=.33*zu(i,j)+.66*w(2)
		       if(j.ge.34) zu(i,j)=w(2)
	            endif
		 end do		!i=latitude loop

	      end do		!j=altitude loop

	   endif		!ihm.ne.0

C_________________________________________________________________________

C  Pre-process arrays which are needed later for interpolation.  The following
C  SURF1 subroutine is from "FITPACK", obtained from Dick Valent 
C  (valent@bierstadt.ucar.edu) in NCAR/SCD consulting office.  The routine
C  that later (in ATMOS5) uses the arrays from SURF1 is called SURFD.

	sigma=1.

	call surf1(37,101,x,y,zt,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zpt,temp,sigma,ierr)

        call surf1(37,101,x,y,zu,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zpu,temp,sigma,ierr)

	call surf1(37,101,x,y,zxh,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zpxh,temp,sigma,ierr)

	call surf1(37,101,x,y,zr,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zpr,temp,sigma,ierr)

	call surf1(37,101,x,y,o2,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,do2,temp,sigma,ierr)

C_________________________________________________________________________

C  Print out array of zonal mean temperature

100	format(20i4)
111     format(4hlat=,19i4)
	ik=0
        do 1 k=1,37,2
	   ik=ik+1
	   il(ik)=(90.01-x(k)*180./pi)
1       continue

101	format(1x,'zonal mean temperature')
c
c Write to the background file if a new month

	if(wback.eq.1)then
 	       write(10, 101)
	       write(10, 111) (il(i),i=1,19)
	endif
c
c Modify loop to print out background up to 400km jkh 8/98

	ix=101
	do 30 i=1,101,2
	   iz=400-(i-1)*4
c	   ix=77
c	   do 30 i=1,77,2
c	      iz=304-(i-1)*4

	   kt=0
	   do 40 k=1,37,2
	      kt=kt+1
	      it(kt)=zt(k,ix)
 40	   continue

	   ix=ix-2
	   if(wback.eq.1)then
	      write(10, 100) iz,(it(k),k=1,19)
	   endif
 30	continue

C  Print out array of zonal mean winds

102	format(1x,'zonal mean winds')

	if(wback.eq.1)then
 	       write(10, 102)
	       write(10, 111) (il(i),i=1,19)
	endif
c
c Modify loop to print out background up to 400km jkh 8/98

	ix=101
	do 31 i=1,101,2
	   iz=400-(i-1)*4
c	   ix=77
c	   do 31 i=1,77,2
c	      iz=304-(i-1)*4

	   kt=0
	   do 41 k=1,37,2
	      kt=kt+1
	      iu(kt)=zu(k,ix)
 41	   continue

	   ix=ix-2
	   if(wback.eq.1)then
	      write(10, 100) iz,(iu(k),k=1,19)
	   endif
 31	continue

 133	format(1x,'pressure (mbar)')
c       write(10, 133)
c	write(10, 111) (il(i),i=1,19)

	ix=26
	do 34 i=1,26  
	   iz=100-(i-1)*4

	   kt=0
	   do 44 k=1,37,2
	      kt=kt+1
	      ih(kt)=zph(k,ix)
 44	   continue

	   ix=ix-1
c	write(10, 100) iz,(ih(k),k=1,19)

34	continue


C  Print out array of normalized intrinsic frequency

103	format(1x,'normalized intrinsic frequency')
c	if(wback.eq.1)then
c 	       write(10, 103)
c	       write(10, 111) (il(i),i=1,19)
c	endif
c
c	ix=77
c	do 32 i=1,77,2
c	iz=304-(i-1)*4
c
c	kt=0
c	do 42 k=1,37,2
c	kt=kt+1
c	sint=sin(x(k))
c	if(k.eq.1) sint=1.
c	if(k.eq.37) sint=1.
c	dfreq=freq+float(nzonal)*zu(k,ix)/(re*sint)
c	dfreq=dfreq/(2.*omega)
c	ifd(kt)=(dfreq*1000.)
c42 	continue
c
c	ix=ix-2
c	if(wback.eq.1)then
c 	write(10, 100) iz,(ifd(k),k=1,19)
c	endif
c32	continue


	return
	end

