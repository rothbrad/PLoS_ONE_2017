program S2AnodalorCathodal2msat200ms


implicit none
! Program to find the S2 (threshold) for Anodal or Cathodal stimulation
! Bidomain model with parsimonious ionic current model 
! Total time 300 ms
! Stimulius is set to zero initially and turned for only 2ms
! unequal space steps dx=0.005 cm (50um) dy=0.002 cm  (20um)
! time step dt=0.001ms (1us)


integer, parameter :: asX = 101 	! asX- array size
integer, parameter :: asY = 251 	! asY- array size
real(8), parameter :: dt = 0.001	! time step in ms (1us)
real(8), parameter :: dx=0.005   	! cm   (50um) tissue size will be 0.5cm by 0.5cm
real(8), parameter :: dy=0.002  	! cm   (20um)
integer, parameter :: nelecX=11		! electrode size 0.05 cm
integer, parameter :: nelecY=26		! electrode length in x and y direction 0.5mm by 0.5mm
real(8), parameter :: maxResiconstant = 0.001 ! mV (1 uV) maxresidue constant
integer, parameter :: maxResiarray=asX*asY
real, parameter :: s2starttime=200.0	! S2 start time (ms)
real, parameter :: s2stoptime=202.0		! S2 end time (ms)
integer, parameter :: jumpX=10
integer, parameter :: jumpY=25



real(8) :: dVmdt=0.0		! deravative of Vm
real(8) :: IionT=0.0		! total ion current		
real(8) :: timetaken=0.0

! Bidomain model parameters according to N. G. Sepulveda, B. J. Roth, and J. P. Wikswo Jr., “Current injection
! into a two-dimensional anisotropic bidomain,” Biophys. J., vol. 55, pp.987-999, 1989.


real(8), parameter :: Cm=1.0		! membrane capacitance uF/cm^2
real(8), parameter :: Beta=2000.0  	! Surface to volume ratio (per cm)
real(8), parameter :: gix=2.0   	! Conductivity in x direction in intracellular space (mS/cm)
real(8), parameter :: gex=8.0		! Conductivity in x direction in extracellular space (mS/cm)
real(8), parameter :: giy=0.2		! Conductivity in y direction in intracellular space (mS/cm)
real(8), parameter :: gey=2.0		! Conductivity in y direction in extracellular space (mS/cm)

real(8) :: Ve(asX,asY)  		! array to hold calculated Ve data (extracellular potential)
real(8) :: Vm(asX,asY) 			! array to hold calculated Vm data (membrane potential)
real(8) :: gixNd2Vmx2(asX,asY)	! second deravatives of Vm in x-direction
real(8) :: giyNd2Vmy2(asX,asY)	! second deravatives of Vm in y-direction
real(8)	:: Istm=0.0				! Stimulus current in uA/cm^2 
real(8) :: Veold(asX,asY)
real(8) :: Veoldold(asX,asY)

	
real(8), parameter :: w=1.8      		! the overrelaxation parameter. should be more than 1 and less than 2 for fast convergence. Ref. page 857. NR in Fortran 
real(8) :: D1,D2,D3,D5,D6,D7,D8,D9 		! constants used for calculation
real(8) :: dVe=0.0         				! defining the residual for Ve
real(8) :: absResiVe=0    				! variable to hold the absolute value of Residual Ve
real(8) :: tofindMaxVe(maxResiarray) 	! array to hold the max residual values for Ve
real(8) :: maxResiVe=0.0   				! varialbe to hold the maximum value of the residual Ve
 

real(8) :: dx2,dy2			! variables to hold square of dx and dy 	
real(8) :: gx,gy			! gx=gix+gex and gy=giy+gey defined below		
real(8) :: InvCm, InvBeta	
real(8) :: Invdx2, Invdy2			
real(8) :: d2Vex2, d2Vey2  	! second deravatives of Ve in x-direction and y-direction

integer :: k 				! time loop variable
integer :: icount=0			! counter write data to screen every 1ms
integer :: l=0,i=0,j=0,nout=0,mm=0,kk=101,ll=0,ww=0,validcounter=0

! Electrode data
real(8) :: currentinX=0.0,currentinY=0.0,totalCurrent=0.0	! currents around x and y sides of electrode
real(8) :: Velec=0.0		! voltage in electrode

! INa current variables and constants
real(8), parameter :: ENa = 65.0		! Reversal potential of Na in mV
real(8), parameter :: GNaMax = 11.0		! Max conductance of Na channel mS/cm2
real(8), parameter :: Em=-41.0			! mV
real(8), parameter :: km=4.0			! mV
real(8) :: m(asX,asY)					! Na activation gate
real(8) :: INa = 0.0					! Fast Na current uA/uF		
reaL(8) :: mInfinity=0.0				! steady state value for m 
real(8) :: tauM=0.12					! time constasnt for m in ms
real(8) :: dmdt=0.0

real(8), parameter :: Eh=-74.7	! mV
real(8), parameter :: kh=4.4	! mV
reaL(8) :: h(asX,asY)			! Na inactivation gate
real(8) :: hInfinity=0.0		! steady state value for h
real(8) :: tauh=0.0
real(8) :: tauHo=6.80738		! time constasnt for ho in ms
real(8) :: deltaH=0.799163		! constant dimensionless
real(8) :: dhdt=0.0

! Irep current variables and constants
real(8), parameter :: Vr=-83.0		! mV
real(8), parameter :: b1=0.3		! mS/cm2
real(8), parameter :: b2=0.047		! 	
real(8) :: Irep=0.0		! 

! Recording the time taken to run the program
real :: elapsed(2)
real :: total					!total time elapsed
integer :: starttime(3),endtime(3),today1(3),today2(3)

character(len=1024) :: filename(500),format_str1,format_str2,format_str3
integer :: UnitNo(500), ii=0

print *,'Program started'

! Storing the time program started and the date
call idate(today1)   		! today(1)=day, (2)=month, (3)=year
call itime(starttime)    	! starttime(1)=hour, (2)=minute, (3)=second

! formatting for out put files
format_str1="(a,i1,'.txt')"
format_str2="(a,i2,'.txt')"
format_str3="(a,i3,'.txt')"

! Opening files to write in every milli second 

do ii=1,9 
	UnitNo(ii) = ii+100
	write(filename(ii),format_str1)'Vm',ii
	open(UNIT=UnitNo(ii), FILE=trim(filename(ii)), FORM="FORMATTED", ACTION="WRITE", STATUS="REPLACE", POSITION="APPEND")
end do

do ii=10,99 
	UnitNo(ii) = ii+100
	write(filename(ii),format_str2)'Vm',ii
	open(UNIT=UnitNo(ii), FILE=trim(filename(ii)), FORM="FORMATTED", ACTION="WRITE", STATUS="REPLACE", POSITION="APPEND")
end do


do ii=100,300 
	UnitNo(ii) = ii+100
	write(filename(ii),format_str3)'Vm',ii
	open(UNIT=UnitNo(ii), FILE=trim(filename(ii)), FORM="FORMATTED", ACTION="WRITE", STATUS="REPLACE", POSITION="APPEND")
end do


dx2 = dx**2
dy2 = dy**2

Invdx2 = 1.0/dx2
Invdy2 = 1.0/dy2

gx = gix + gex
gy = giy + gey

InvCm = 1.0/Cm
InvBeta = 1.0/Beta

d2Vex2 = 0
d2Vey2 = 0
	

D1 = ((2.0*gx)/dx2) + ((2.0*gy)/dy2)	! units are mS cm^-3
D2 = gex*InvBeta*InvCm					! units 10^-3 cm^2 s^-1
D3 = gey*InvBeta*InvCm					! units 10^-3 cm^2 s^-1
D5 = gx*Invdx2							! units are mS cm^-3
D6 = gy*Invdy2							! units are mS cm^-3
D7 = InvBeta*InvCm						! units are 10^6 R cm^3 s^-1
D8 = gey*dx2							! units are mS cm
D9 = gex*dy2							! units are mS cm


! initializing Ve to zero and Vm resting potential
do j=1,asY                
	do i =1,asX
    Ve(i,j) = 0.0
    Vm(i,j) = -83.0		   	! setting the initial membrane potential to resting potential 
	gixNd2Vmx2(i,j) = 0.0
	giyNd2Vmy2(i,j) = 0.0
	Veold(i,j) = 0.0
	Veoldold(i,j) = 0.0
	m(i,j)=2.7540e-05	
	h(i,j)=0.8683
	end do
end do


do k=0,300000! time loop for 300 ms for testing purpose

	timetaken = (k*dt) - 5.0

! S1 stimulus
	if(timetaken>=0.0)then    	! cathodal Stimulus applied at the 0th ms, twice the threshold value
		Istm = -5996.094  		! mA/m  
		
	endif
	
! S1 stops
	if(timetaken>=2.0)then  	! Stimulus turned off after 2ms
		Istm = 0.0
	end if

! S2 stimulus
	if(timetaken>=s2starttime)then 	! Anodal Stimulus applied at the 200th ms
		Istm = 57923.828  			! mA/m
	endif
	
! S2 stops
	if(timetaken>=s2stoptime)then	! Stimulus is removed after 2ms
		Istm = 0.0
	end if
	
	

	do j=2,asY-1
	  	do i=2,asX-1

			! ignoring the region covered by the electrode
			if ((i>= 2 .AND. i<= nelecX) .AND. (j>=2 .AND. j<= nelecY)) then
				do ww=1,nelecY-1
					do ll=1,nelecX-1
						Vm(ll,ww)=-800.0
					end do
				enddo
			else

				d2Vex2 = (Ve(i+1,j) - (2*Ve(i,j)) + Ve(i-1,j))*Invdx2

				d2Vey2 = (Ve(i,j+1) - (2*Ve(i,j)) + Ve(i,j-1))*Invdy2

				! INa current calculation
			
				mInfinity = 1/(1+exp((-(Vm(i,j)-Em))/km))
			
				dmdt = (mInfinity - m(i,j))/tauM
				m(i,j) = m(i,j) + dt*dmdt	
			
				tauH = (2*tauHo*exp((deltaH*(Vm(i,j)-Eh))/kh))/(1+exp((Vm(i,j)-Eh)/kh))
				if (tauH<dt) then
					tauH=dt
				end if
			
				hInfinity = 1/(1+exp((Vm(i,j)-Eh)/kh))
				
				dhdt = (hInfinity - h(i,j))/tauH
				h(i,j) = h(i,j) + dt*dhdt	
			
				INa = GNamax*(m(i,j)**3)*h(i,j)*(Vm(i,j)-ENa)   !  Calculation of INa
					
				! Irep current calculation
				Irep = b1*(Vm(i,j)-Vr)*exp(-b2*(Vm(i,j)-Vr))    ! Calculation of Irep
				
				IionT = Irep+INa
						
				dVmdt = -((D2*d2Vex2) + (D3*d2Vey2) + (InvCm*IionT))
				
				Vm(i,j) = Vm(i,j) + dt*dVmdt
			
			end if
		end do
	end do

! Boundry conditions for Vm

	
	do j=1,nelecY
		Vm(nelecX,j) = Vm(nelecX+1,j) + Ve(nelecX+1,j) - Ve(nelecX,j)       ! Right side of the electrode
	end do
	
	do i=1,nelecX
		Vm(i,nelecY) = Vm(i,nelecY+1) + Ve(i,nelecY+1) - Ve(i,nelecY)		! Top of the electrode
	end do
			
	do j=1,asY
		Vm(asX,j) = Vm(asX-1,j)  
	end do
	
	do i=1,asX		
		Vm(i,asY) = Vm(i,asY-1)  
	end do
	
	do i=nelecX+1,asX
		Vm(i,1) = Vm(i,2) ! y=1 line starting from the end of the electrode
	end do
	
	do j=nelecY+1,asY
		Vm(1,j) = Vm(2,j) ! x=1 line starting from the end of the electrode
	end do
	
	do j=2,asY-1
		do i=2,asX-1
			gixNd2Vmx2(i,j) = gix*((Vm(i+1,j) - (2*Vm(i,j)) + Vm(i-1,j))*Invdx2)
			giyNd2Vmy2(i,j) = giy*((Vm(i,j+1) - (2*Vm(i,j)) + Vm(i,j-1))*Invdy2)
		end do
	end do
	
	do j=1,asY
		do i=1,asX
			Veoldold(i,j) = Veold(i,j)
			Veold(i,j) = Ve(i,j)
			Ve(i,j) = 2*Veold(i,j) - Veoldold(i,j)
		end do
	end do

do l=0,1000
	
	mm=1
	
	do j=2,asY-1
	do i=2,asX-1
		if ((i>= 2 .AND. i<= nelecX) .AND. (j>=2 .AND. j<= nelecY)) then
			tofindMaxVe(mm)=0.0
		else		
	dVe= ((((Ve(i+1,j)+Ve(i-1,j))*D5)+((Ve(i,j+1)+Ve(i,j-1))*D6)+gixNd2Vmx2(i,j)+giyNd2Vmy2(i,j))/D1)-Ve(i,j)				
	Ve(i,j) = w*dVe + Ve(i,j)
	absResiVe = ABS(dVe)	! finding the absolute value of the residual
	tofindMaxVe(mm) = absResiVe   ! assigning all absolute values of residual to an array 
		end if
		
	mm=mm+1
	end do
	end do	
	maxResiVe = MAXVAL(tofindMaxVe)    	! to find the MAXVALUE(absResi) note. maxResi depend on w value
	! Calculating current around electrode at the left edge of the tissue starting at (1,1)
		currentinX = 0.5*D8*(Ve(1,nelecY+1)+Ve(nelecX,nelecY+1))

		do ll=2,nelecX-1
			currentinX = currentinX + (D8*Ve(ll,nelecY+1))
		end do

		currentinY = 0.5*D9*(Ve(nelecX+1,1)+Ve(nelecX+1,nelecY))
			
		do ww=2,nelecY-1
			currentinY = currentinY + (D9*Ve(nelecX+1,ww)) 
		end do
		
		totalCurrent = currentinX+currentinY+((Istm/4)*dx*dy)
		Velec=totalCurrent/(((nelecX-1)*D8)+((nelecY-1)*D9))
		
		! Boundry conditions for the electrode

		do ww=1,nelecY
			Ve(nelecX,ww) = Velec 	! Right side
		end do
			
		do ll=nelecX+1,asX
			Ve(ll,1) = Ve(ll,2)		! bottom side of tissue from end of electrode
		end do
			
		do ww=nelecY+1,asY-1
			Ve(1,ww) = Ve(2,ww) 	! Left side of tissue from end of electrode
		end do
			
		do ll=1,nelecX
			Ve(ll,nelecY) = Velec	! top side
		end do
				
		do ww= 1,asY
			Ve(asX,ww) = 0 
		end do
				
		do ll= 1,asX
			Ve(ll,asY) = 0 
		end do
			
	if (maxResiVe < maxResiconstant ) then	  	! checking if residue is less than a certain tolerance value. (0.001 mV)
		goto 300						! if max residue is < 0.001 mV break out of iteration
	end if
	
end do

	
	300 continue 
	! Terminal out put of data
	if (nout==1000) then
		write(*,*)		
		write(6,*) timetaken,'ms',maxResiVe,l, Istm,k! write to the scrren every 1ms: timeloop no, max residue value, iteration no
		write(*,*)
		
		do j=asY,1,-jumpY								
			write(6,3000) (Vm(i,j),i=asX,1,-jumpX)
		end do

		nout = 0
	end if	
	nout=nout+1

	! Write to the files
	if(icount == 1000 ) then
		WRITE(UNIT=kk, FMT="(101('',e15.4))") Vm ! data to be written to the file with 100 columns
		close(kk)
		kk=kk+1
		icount=0
	endif
	
	icount = icount+1	
	
	3000 format (4096(f5.0,' '))	
		


end do


total = etime(elapsed)
print *,'Time elapsed: ',total,'s'
write (*,*)

write ( *, 1000 )  today1(2), today1(1), today1(3), starttime
1000 format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; Start time ', i2.2, ':', i2.2, ':', i2.2)

call idate(today2)
call itime(endtime)     ! endtimetime(1)=hour, (2)=minute, (3)=second
write(*,2000) today2(2), today2(1), today2(3), endtime
2000 format ( 'Date ', i2.2, '/', i2.2, '/', i4.4, '; End time ', i2.2, ':', i2.2, ':', i2.2)


stop


end program S2AnodalorCathodal2msat200ms

