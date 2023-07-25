double precision function difSDVe(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev,scale,pow,tpl)
!with accounting for the precise Bouger law
!strength_in is now real line strength multiplied to path length L
!scale is scaling factor due to the variable sensitivity of the RAD cell
!now we need to calc real power loss due to exp(-alpha*L) and corresponding difference
!and then multiply the result to the scaling factor
	
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow
double precision qRe1, qRe2, qIm, pow0 !, qRe0 !, qRe3, qRe4, sh, amp,g0s
integer tpl

call qSDV(sg0,GamD,Gam0,Gam2,Shift0,Shift2,sg+dev,qRe1,qIm)
call qSDV(sg0,GamD,Gam0,Gam2,Shift0,Shift2,sg-dev,qRe2,qIm)

!call HTP(sg0,GamD,Gam0,Gam2,Shift0,Shift2,NVC,0.D0,sg+dev,qRe1,qIm)
!call HTP(sg0,GamD,Gam0,Gam2,Shift0,Shift2,NVC,0.D0,sg-dev,qRe2,qIm)


!call qSDV(sg0,GamD,Gam0,Gam2,Shift0,Shift2,sg0,qRe0,qIm)
!qRe0 = qRe0*strength_in

pow0 = pow


!if (qRe0.gt.1.D-2) then !if absorption*length in the line center exceeds 0.001 use exponent
if (tpl.eq.1) then
	qRe1 = scale*(1.D0-dexp(-qRe1*strength_in))*(1+pow0*(sg+dev-sg0))
	qRe2 = scale*(1.D0-dexp(-qRe2*strength_in))*(1+pow0*(sg-dev-sg0))
else if (tpl.eq.2) then
    qRe1 = dexp(-qRe1*strength_in)*(scale+pow0*(sg+dev-sg0))
	qRe2 = dexp(-qRe2*strength_in)*(scale+pow0*(sg-dev-sg0))
end if

!else  !otherwise use linear approach
!	qRe1 = (qRe1+amp*qRe3)*strength_in*(1+pow0*(sg+dev-sg0))
!	qRe2 = (qRe2+amp*qRe4)*strength_in*(1+pow0*(sg-dev-sg0))
!endif

difSDVe = qRe1 - qRe2

end function difSDVe


double precision function difSDVe_add(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev,scale,pow,amp,NVC,tpl)
!with accounting for the precise Bouger law
!strength_in is now real line strength multiplied to path length L
!scale is scaling factor due to the variable sensitivity of the RAD cell
!now we need to calc real power loss due to exp(-alpha*L) and corresponding difference
!and then multiply the result to the scaling factor

!also accounts for 2-nd BWO harmonic interacting with the other CO line (R1 for R0, R3 for R1)
	
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev, scale,pow,NVC
double precision qRe1, qRe2, qIm, qRe3, qRe4, sh, amp,g0s, pow0, g2s
integer tpl

!sh = -2.0D0 !0.5*230538.D0-sg0	for 115 line
!sh = 230538.45D0
!sh = sg0-17.6D0 ! for 230 line
sh = sg0-2.0D0 !-17.6D0
!sh = 230520.384
!amp = 2.0D-3  ! CO-230 s1 ->2.5D-3 CO-230 s2 ->1.4D-3
g0s = 0.9D0  !0.951D0
g2s = 0.9

!call qSDV(sg0,GamD,Gam0,Gam2,Shift0,Shift2,sg+dev,qRe1,qIm)
!call qSDV(sg0,GamD,Gam0,Gam2,Shift0,Shift2,sg-dev,qRe2,qIm)

call HTP(sg0,GamD,Gam0,Gam2,Shift0,Shift2,NVC,0.D0,sg+dev,qRe1,qIm)
call HTP(sg0,GamD,Gam0,Gam2,Shift0,Shift2,NVC,0.D0,sg-dev,qRe2,qIm)


!call qSDV(sg0,GamD,Gam0,Gam2,Shift0,Shift2,sg0,qRe0,qIm)
!qRe0 = qRe0*strength_in

qRe3 = 0.D0
qRe4 = 0.D0

pow0 = pow

call qSDV(2.d0*sh,GamD*2.,Gam0*g0s,Gam2*g2s,0.d0,0.d0,2.d0*(sg+dev),qRe3,qIm)
call qSDV(2.d0*sh,GamD*2.,Gam0*g0s,Gam2*g2s,0.d0,0.d0,2.d0*(sg-dev),qRe4,qIm)


!if (qRe0.gt.1.D-2) then !if absorption*length in the line center exceeds 0.001 use exponent
if (tpl.eq.1) then
	qRe1 = scale*(1.D0-dexp(-(qRe1+amp*qRe3)*strength_in))*(1+pow0*(sg+dev-sg0))
	qRe2 = scale*(1.D0-dexp(-(qRe2+amp*qRe4)*strength_in))*(1+pow0*(sg-dev-sg0))
else if (tpl.eq.2) then
    qRe1 = dexp(-(qRe1+amp*qRe3)*strength_in)*(scale+pow0*(sg+dev-sg0))
	qRe2 = dexp(-(qRe2+amp*qRe4)*strength_in)*(scale+pow0*(sg-dev-sg0))
end if
!else  !otherwise use linear approach
!	qRe1 = (qRe1+amp*qRe3)*strength_in*(1+pow0*(sg+dev-sg0))
!	qRe2 = (qRe2+amp*qRe4)*strength_in*(1+pow0*(sg-dev-sg0))
!endif

difSDVe_add = qRe1 - qRe2

end function difSDVe_add


double precision function difSDV(sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev)
	
double precision sg,sg0,GamD,Gam0,Gam2,Shift0,Shift2,strength_in,dev
double precision qRe1, qRe2, qIm

call qSDV(sg0,GamD,Gam0,Gam2,Shift0,Shift2,sg+dev,qRe1,qIm)
call qSDV(sg0,GamD,Gam0,Gam2,Shift0,Shift2,sg-dev,qRe2,qIm)

difSDV = strength_in*(qRe1 - qRe2)

end function difSDV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine 	HTP(sg0,GamD,Gam0,Gam2,Shift0,Shift2,anuVC,eta,sg,LS_pCqSDHC_R,LS_pCqSDHC_I)
!-------------------------------------------------
!	"pCqSDHC": partially-Correlated quadratic-Speed-Dependent Hard-Collision
!	Subroutine to Compute the complex normalized spectral shape of an 
!	isolated line by the pCqSDHC model
!
!	Input/Output Parameters of Routine (Arguments or Common)
!	---------------------------------
!	T	    : Temperature in Kelvin (Input).
!	amM1	: Molar mass of the absorber in g/mol(Input).
!	sg0		: Unperturbed line position in cm-1 (Input).
!     GamD	: Doppler HWHM in cm-1 (Input)
!	Gam0	: Speed-averaged line-width in cm-1 (Input). 	
!	Gam2	: Speed dependence of the line-width in cm-1 (Input).
!	anuVC	: Velocity-changing frequency in cm-1 (Input).
!	eta		: Correlation parameter, No unit (Input).
!	Shift0	: Speed-averaged line-shift in cm-1 (Input).
!	Shift2	: Speed dependence of the line-shift in cm-1 (Input)	 
!	sg		: Current WaveNumber of the Computation in cm-1 (Input).
!
!	Output Quantities (through Common Statements)
!	-----------------
!	LS_pCqSDHC_R: Real part of the normalized spectral shape (cm)
!	LS_pCqSDHC_I: Imaginary part of the normalized spectral shape (cm)
!
!	Called Routines: 'CPF'	(Complex Probability Function)
!	---------------  'CPF3'	(Complex Probability Function for the region 3)
!
!	Called By: Main Program
!	---------
!
!     Double Precision Version
!
!-------------------------------------------------
	implicit none
	 double precision sg0,GamD
	 double precision Gam0,Gam2,anuVC,eta,Shift0,Shift2
	 double precision sg
	 double precision pi,rpi,cte
	 double precision xz1,xz2,yz1,yz2,xXb,yXb
	 double precision wr1,wi1,wr2,wi2,wrb,wib
	 double precision SZ1,SZ2,DSZ,SZmx,SZmn
	 double precision LS_pCqSDHC_R,LS_pCqSDHC_I
	double complex c0,c2,c0t,c2t
	double complex X,Y,iz,Z1,Z2
	double complex Aterm,Bterm,LS_pCqSDHC
!
!-------------------------------------------------
!
	cte=dsqrt(dlog(2.D0))/GamD
	pi=4.d0*datan(1.d0)
	rpi=dsqrt(pi)
	iz=dcmplx(0.d0,1.d0)
! Calculating the different parameters 
	c0=dcmplx(Gam0,-Shift0)
	c2=dcmplx(Gam2,-Shift2)
	c0t=(1.d0-eta)*(c0-1.5d0*c2)+anuVC
	c2t=(1.d0-eta)*c2
	Y=1.d0/((2.d0*cte*C2t))**2			
!
	
	X=(iz*(sg-sg0)+c0t)/c2t
!	
	if (cdabs(C2t).eq.0.d0) go to 110
	if (cdabs(X).le.3.d-8*cdabs(Y)) go to 120
	if (cdabs(Y).le.1.d-15*cdabs(X)) go to 140
! calculating Z1 and Z2
	Z1=cdsqrt(X+Y)-cdsqrt(Y)
	Z2=Z1+2.d0*cdsqrt(Y)
! calculating the real and imaginary parts of Z1 and Z2
	xZ1=-dimag(Z1)
	yZ1=dreal(Z1)
	xZ2=-dimag(Z2)
	yZ2=dreal(Z2)
! check if Z1 and Z2 are close to each other
	SZ1=dsqrt(xZ1*xZ1+yZ1*yZ1)
	SZ2=dsqrt(xZ2*xZ2+yZ2*yZ2)
	DSZ=dabs(SZ1-SZ2)
	SZmx=dmax1(SZ1,SZ2)
	SZmn=dmin1(SZ1,SZ2)
! when Z1 and Z2 are close to each other, ensure that they are in 
! the same interval of CPF 
	if (DSZ.le.1.d0.and.SZmx.gt.8.d0.and.SZmn.le.8.d0) then
	Call CPF3 ( xZ1, yZ1, WR1, WI1 ) 
	Call CPF3 ( xZ2, yZ2, WR2, WI2 ) 
	else	
	Call CPF ( xZ1, yZ1, WR1, WI1 ) 
	Call CPF ( xZ2, yZ2, WR2, WI2 ) 
	endif
! calculating the A and B terms of the profile
	Aterm=rpi*cte*(dcmplx(wr1,wi1)-dcmplx(wr2,wi2))
	Bterm=(-1.d0+rpi/(2.d0*cdsqrt(Y))*(1.d0-Z1**2)*dcmplx(wr1,wi1)-rpi/(2.d0*cdsqrt(Y))*(1.d0-Z2**2)*dcmplx(wr2,wi2))/C2t
	go to 10
! when C2t=0
110   continue
	Z1=(iz*(sg-sg0)+C0t)*cte
	xZ1=-dimag(Z1)
	yZ1=dreal(Z1)
	Call CPF ( xZ1, yZ1, WR1, WI1 )
	Aterm=rpi*cte*dcmplx(WR1,WI1)
	if (cdabs(Z1).le.4.d3) then
      Bterm=rpi*cte*((1.d0-Z1**2)*dcmplx(WR1,WI1)+Z1/rpi)
	else
	Bterm=cte*(rpi*dcmplx(WR1,WI1)+0.5d0/Z1-0.75d0/(Z1**3))
	endif
	go to 10
! when abs(Y) is much larger than abs(X)
120   continue
	Z1=(iz*(sg-sg0)+C0t)*cte
	Z2=cdsqrt(X+Y)+cdsqrt(Y)
	xZ1=-dimag(z1)
	yZ1=dreal(z1)
	xZ2=-dimag(z2)
	yZ2=dreal(z2)
	Call CPF ( xZ1, yZ1, WR1, WI1 )
	Call CPF ( xZ2, yZ2, WR2, WI2 ) 
	Aterm=rpi*cte*(dcmplx(WR1,WI1)-dcmplx(WR2,WI2))
	Bterm=(-1.d0+rpi/(2.d0*cdsqrt(Y))*(1.d0-Z1**2)*dcmplx(wr1,wi1)-rpi/(2.d0*cdsqrt(Y))*(1.d0-Z2**2)*dcmplx(wr2,wi2))/C2t
	go to 10
! when abs(X) is much larger than abs(Y)
140   continue
	xZ1=-dimag(cdsqrt(X+Y))
	yZ1=dreal(cdsqrt(X+Y))
	Call CPF ( xZ1, yZ1, WR1, WI1 ) 
	if (cdabs(cdsqrt(X)).le.4.d3) then
	  xXb=-dimag(cdsqrt(X))
	  yXb=dreal(cdsqrt(X))
	  Call CPF ( xXb, yXb, WRb, WIb ) 
	  Aterm=(2.d0*rpi/C2t)*(1.d0/rpi-cdsqrt(X)*dcmplx(WRb,WIb))
	  Bterm=(1.d0/C2t)*(-1.d0+2.d0*rpi*(1.d0-X-2.d0*Y)*(1.d0/rpi-cdsqrt(X)*dcmplx(wrb,wib))+2.d0*rpi*cdsqrt(X+Y)*dcmplx(wr1,wi1))
! and when abs(X) is much larger than 1
	else
	  Aterm=(1.d0/C2t)*(1.d0/X-1.5d0/(X**2))
	  Bterm=(1.d0/C2t)*(-1.d0+(1.d0-X-2.d0*Y)*(1.d0/X-1.5d0/(X**2))+2.d0*rpi*cdsqrt(X+Y)*dcmplx(wr1,wi1))
	endif
!
10    continue
!
	LS_pCqSDHC=(1.d0/pi)*(Aterm/(1.d0-(anuVC-eta*(C0-1.5d0*C2))*Aterm+eta*C2*Bterm))

	LS_pCqSDHC_R=dreal(LS_pCqSDHC)
	LS_pCqSDHC_I=dimag(LS_pCqSDHC)

   
Return
End Subroutine HTP


	
subroutine qSDV(sg0,GamD,Gam0,Gam2,Shift0,Shift2,sg,LS_qSDV_R,LS_qSDV_I)
!-------------------------------------------------
!	"qSDV": quadratic-Speed-Dependent Voigt
!	Subroutine to Compute the complex normalized spectral shape of an 
!	isolated line by the qSDV model following the two references:
!	[1] Ngo NH, Lisak D, Tran H, Hartmann J-M. An isolated line-shape model
!	to go beyond the Voigt profile in spectroscopic databases and radiative
!	transfer codes. J Quant Radiat Transfer 2013;129:89-100.	 	
!	[2] Tran H, Ngo NH, Hartmann J-M. Efficient computation of some speed-dependent 
!	isolated line profiles. J Quant Radiat Transfer 2013;129:199-203.
!
!	Input/Output Parameters of Routine (Arguments or Common)
!	---------------------------------
!	T	    : Temperature in Kelvin (Input).
!	amM1	: Molar mass of the absorber in g/mol(Input).
!	sg0		: Unperturbed line position in cm-1 (Input).
!	GamD	: Doppler HWHM in cm-1 (Input)
!	Gam0	: Speed-averaged line-width in cm-1 (Input). 	
!	Gam2	: Speed dependence of the line-width in cm-1 (Input).
!	Shift0	: Speed-averaged line-shift in cm-1 (Input).
!	Shift2	: Speed dependence of the line-shift in cm-1 (Input)	 
!	sg		: Current WaveNumber of the Computation in cm-1 (Input).
!
!	Output Quantities (through Common Statements)
!	-----------------
!	LS_qSDV_R: Real part of the normalized spectral shape (cm)
!	LS_qSDV_I: Imaginary part of the normalized spectral shape (cm)
!
!	Called Routines: 'CPF'	(Complex Probability Function)
!	---------------  'CPF3'	(Complex Probability Function for the region 3)
!
!	Called By: Main Program
!	---------
!
!	Double Precision Version
!
!-------------------------------------------------
	implicit none
	 double precision sg0,GamD
	 double precision Gam0,Gam2,Shift0,Shift2
	 double precision sg
	 double precision pi,rpi,cte
	 double precision xz1,xz2,yz1,yz2,xXb,yXb
	 double precision wr1,wi1,wr2,wi2,wrb,wib
	 double precision SZ1,SZ2,DSZ,SZmx,SZmn
	 double precision LS_qSDV_R,LS_qSDV_I
	double complex c0,c2,c0t,c2t
	double complex X,Y,iz,Z1,Z2,csqrtY
	double complex Aterm,LS_qSDV
!
!-------------------------------------------------
!
	cte=dsqrt(dlog(2.D0))/GamD
	pi=4.d0*datan(1.d0)
	rpi=dsqrt(pi)
	iz=dcmplx(0.d0,1.d0)
! Calculating the different parameters 
	c0=dcmplx(Gam0,Shift0)
	c2=dcmplx(Gam2,Shift2)
	c0t=(c0-1.5d0*c2)
	c2t=c2
!
	if (cdabs(C2t).eq.0.d0) go to 110	
	X=(iz*(sg0-sg)+c0t)/c2t
	Y=1.d0/((2.d0*cte*C2t))**2		
	csqrtY=(Gam2-iz*Shift2)/(2.d0*cte*(Gam2**2+Shift2**2))
	if (cdabs(X).le.3.d-8*cdabs(Y)) go to 120
	if (cdabs(Y).le.1.d-15*cdabs(X)) go to 140
! calculating Z1 and Z2
	Z1=cdsqrt(X+Y)-csqrtY
	Z2=Z1+2.d0*csqrtY
! calculating the real and imaginary parts of Z1 and Z2
	xZ1=-dimag(Z1)
	yZ1=dreal(Z1)
	xZ2=-dimag(Z2)
	yZ2=dreal(Z2)
! check if Z1 and Z2 are close to each other
	SZ1=dsqrt(xZ1*xZ1+yZ1*yZ1)
	SZ2=dsqrt(xZ2*xZ2+yZ2*yZ2)
	DSZ=dabs(SZ1-SZ2)
	SZmx=dmax1(SZ1,SZ2)
	SZmn=dmin1(SZ1,SZ2)
! when Z1 and Z2 are close to each other, ensure that they are in 
! the same interval of CPF 
	if (DSZ.le.1.d0.and.SZmx.gt.8.d0.and.SZmn.le.8.d0) then
	Call CPF3 ( xZ1, yZ1, WR1, WI1 ) 
	Call CPF3 ( xZ2, yZ2, WR2, WI2 ) 
	else	
	Call CPF ( xZ1, yZ1, WR1, WI1 ) 
	Call CPF ( xZ2, yZ2, WR2, WI2 ) 
	endif
! calculating the A term of the profile
	Aterm=rpi*cte*(dcmplx(wr1,wi1)-dcmplx(wr2,wi2))
!	write(15,*)sg,wr1,wr2
	go to 10
! when C2t=0
110   continue
	Z1=(iz*(sg0-sg)+C0t)*cte
	xZ1=-dimag(Z1)
	yZ1=dreal(Z1)
	Call CPF ( xZ1, yZ1, WR1, WI1 )
	Aterm=rpi*cte*dcmplx(WR1,WI1)
	go to 10
! when abs(Y) is much larger than abs(X)
120   continue
	Z1=(iz*(sg0-sg)+C0t)*cte
	Z2=cdsqrt(X+Y)+csqrtY
	xZ1=-dimag(z1)
	yZ1=dreal(z1)
	xZ2=-dimag(z2)
	yZ2=dreal(z2)
	Call CPF ( xZ1, yZ1, WR1, WI1 )
	Call CPF ( xZ2, yZ2, WR2, WI2 ) 
	Aterm=rpi*cte*(dcmplx(WR1,WI1)-dcmplx(WR2,WI2))
	go to 10
! when abs(X) is much larger than abs(Y)
140   continue
	if (cdabs(cdsqrt(X)).le.4.d3) then
	  xXb=-dimag(cdsqrt(X))
	  yXb=dreal(cdsqrt(X))
	  Call CPF ( xXb, yXb, WRb, WIb ) 
	  Aterm=(2.d0*rpi/C2t)*(1.d0/rpi-cdsqrt(X)*dcmplx(WRb,WIb))
!! and when abs(X) is much larger than 1
	else
	  Aterm=(1.d0/C2t)*(1.d0/X-1.5d0/(X**2))
	endif
!
10    continue
!
	LS_qSDV=(1.d0/pi)*Aterm

	LS_qSDV_R=dreal(LS_qSDV)
	LS_qSDV_I=dimag(LS_qSDV)

   
      Return
      End Subroutine qSDV


      Subroutine CPF(X,Y,WR,WI)
!-------------------------------------------------
! "CPF": Complex Probability Function
! .........................................................
!         .       Subroutine to Compute the Complex       .
!         .        Probability Function W(z=X+iY)         .
!         .     W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0      .
!         .    Which Appears when Convoluting a Complex   .
!         .     Lorentzian Profile by a Gaussian Shape    .
!         .................................................
!
!             WR : Real Part of W(z)
!             WI : Imaginary Part of W(z)
!
! This Routine was Taken from the Paper by J. Humlicek, which 
! is Available in Page 309 of Volume 21 of the 1979 Issue of
! the Journal of Quantitative Spectroscopy and Radiative Transfer
! Please Refer to this Paper for More Information
!
! Accessed Files:  None
! --------------
!
! Called Routines: None                               
! ---------------                                 
!
! Called By: 'CompAbs' (COMPute ABSorpton)
! ---------
!
! Double Precision Version
!
!-------------------------------------------------
!      
      Implicit None
        Integer I
	double complex zm1,zm2,zterm,zsum,zone,zi
      Double Precision X,Y,WR,WI
      Double Precision T,U,S,Y1,Y2,Y3,R,R2,D,D1,D2,D3,D4
      Double Precision TT(15),pipwoeronehalf
!      
      Dimension T(6),U(6),S(6)
      Data T/.314240376d0,.947788391d0,1.59768264d0,2.27950708d0,3.02063703d0,3.8897249d0/
      Data U/1.01172805d0,-.75197147d0,1.2557727d-2,1.00220082d-2,-2.42068135d-4,5.00848061d-7/
      Data S/1.393237d0,.231152406d0,-.155351466d0,6.21836624d-3,9.19082986d-5,-6.27525958d-7/
	Data zone,zi/(1.d0,0.D0),(0.d0,1.D0)/
	data tt/0.5d0,1.5d0,2.5d0,3.5d0,4.5d0,5.5d0,6.5d0,7.5d0,8.5d0,9.5d0,10.5d0,11.5d0,12.5d0,13.5d0,14.5d0/
	data pipwoeronehalf/0.564189583547756d0/

! new Region 3
	if(dsqrt(x*x+y*Y).gt.8.D0)then
	zm1=zone/dcmplx(x,y)
	zm2=zm1*zm1
	zsum=zone
	zterm=zone
	do i=1,15
	zterm=zterm*zm2*tt(i)
	zsum=zsum+zterm
	end do
	zsum=zsum*zi*zm1*pipwoeronehalf
	wr=dreal(zsum)
	wi=dimag(zsum)
	return
	end if
!
      WR=0.d0
      WI=0.d0
      Y1=Y+1.5d0
      Y2=Y1*Y1
      If( (Y.GT.0.85d0) .OR. (DABS(X).LT.(18.1d0*Y+1.65d0)) )GoTo 2
!
!       Region 2
!
      If( DABS(X).LT.12.d0 )WR=DEXP(-X*X)
      Y3=Y+3.d0
      Do 1 I=1,6
      R=X-T(I)
      R2=R*R
      D=1.d0/(R2+Y2)
      D1=Y1*D
      D2=R*D
      WR=WR+Y*(U(I)*(R*D2-1.5d0*D1)+S(I)*Y3*D2)/(R2+2.25d0)
      R=X+T(I)
      R2=R*R
      D=1.d0/(R2+Y2)
      D3=Y1*D
      D4=R*D
      WR=WR+Y*(U(I)*(R*D4-1.5d0*D3)-S(I)*Y3*D4)/(R2+2.25d0)
      WI=WI+U(I)*(D2+D4)+S(I)*(D1-D3)
 1    Continue  
      Return
!
!       Region 1
!
 2    Continue
      Do 3 I=1,6
      R=X-T(I)
      D=1.d0/(R*R+Y2)
      D1=Y1*D
      D2=R*D
      R=X+T(I)
      D=1.d0/(R*R+Y2)
      D3=Y1*D
      D4=R*D
      WR=WR+U(I)*(D1+D3)-S(I)*(D2-D4)
      WI=WI+U(I)*(D2+D4)+S(I)*(D1-D3)
 3    Continue  
      Return
      End


      Subroutine CPF3(X,Y,WR,WI)
!-------------------------------------------------
! "CPF": Complex Probability Function
! .........................................................
!         .       Subroutine to Compute the Complex       .
!         .        Probability Function W(z=X+iY)         .
!         .     W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0      .
!         .    Which Appears when Convoluting a Complex   .
!         .     Lorentzian Profile by a Gaussian Shape    .
!         .................................................
!
!             WR : Real Part of W(z)
!             WI : Imaginary Part of W(z)
!
! This Routine takes into account the region 3 only, i.e. when sqrt(x**2+y**2)>8. 
!
! Accessed Files:  None
! --------------
!
! Called Routines: None                               
! ---------------                                 
!
! Called By: 'pCqSDHC'
! ---------
!
! Double Precision Version
! 
!-------------------------------------------------
!      
      Implicit None
        Integer I
	double complex zm1,zm2,zterm,zsum,zone,zi
      Double Precision X,Y,WR,WI
      Double Precision TT(15),pipwoeronehalf
!      
	Data zone,zi/(1.d0,0.D0),(0.d0,1.D0)/
	data tt/0.5d0,1.5d0,2.5d0,3.5d0,4.5d0,5.5d0,6.5d0,7.5d0,8.5d0,9.5d0,10.5d0,11.5d0,12.5d0,13.5d0,14.5d0/
	data pipwoeronehalf/0.564189583547756d0/

! Region 3
	zm1=zone/dcmplx(x,y)
	zm2=zm1*zm1
	zsum=zone
	zterm=zone
	do i=1,15
	zterm=zterm*zm2*tt(i)
	zsum=zsum+zterm
	end do
	zsum=zsum*zi*zm1*pipwoeronehalf
	wr=dreal(zsum)
	wi=dimag(zsum)
	return
      End