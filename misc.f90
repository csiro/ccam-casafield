!
! THIS FILE CONTAINS MISC SUBROUTINES
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function converts between two sets of units
!

Subroutine calrescale(inunit,outunit,x)

Implicit None

Character(len=*), intent(in) :: inunit,outunit
Real, intent(inout) :: x
real, dimension(1) :: xd
Character*80 baseunit

baseunit=''

xd(1)=x
If (inunit.NE.outunit) Then
  Call baserescale(inunit,xd,1,baseunit,1)
  Call baserescale(outunit,xd,-1,baseunit,1)
End If
x=xd(1)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The function rescales to a base value for calrescale
!

Subroutine baserescale(inunit,x,inverse,baseunit,sz)

Implicit None

Integer, intent(in) :: inverse,sz
Character(len=*), intent(in) :: inunit
Character(len=*), intent(inout) :: baseunit
Real, dimension(sz), intent(inout) :: x
Character*80 actunit

actunit=''

Select Case(inunit)

  Case('minute')
    actunit='hour'
    If (inverse.EQ.1) Then
      x=x*60.
    Else
      x=x/60.
    End If

  Case('minutes')
    actunit='hour'
    If (inverse.EQ.1) Then
      x=x*60.
    Else
      x=x/60.
    End If
  
  Case('hour')
    actunit='hour'
    ! No change

  Case('hours')
    actunit='hour'
    ! No change

  Case('hrs')
    actunit='hour'  
    ! No change
    
  Case('Pa')
    actunit='Pa'  
    ! No change
    
  Case('hPa')
    actunit='Pa'  
    If (inverse.EQ.1) Then
      x=x*100.
    Else
      x=x/100.
    End If
  
  Case('meters')
    actunit='meters'  
    ! No change

  Case('none')
    actunit='none'
    ! No change

  Case('')
    actunit='none'  
    ! No change
    
  Case('g/g')
    actunit='ratio'  
    ! No change

  Case('kg/kg')
    actunit='ratio'  
    ! No change
    
  Case('W/m^2')
    actunit='W/m^2'
    ! No change
    
  Case('W/m2')
    actunit='W/m^2'
    ! No change
   
  Case('%')
    actunit='ratio'
    If (inverse.EQ.1) Then
      x=x*100.
    Else
      x=x/100.
    End If

  Case('ratio')
    actunit='ratio'
    ! No change

  Case('K')
    actunit='K'
    ! No change

  Case('Kelvin')
    actunit='K'
    ! No change

  Case('C')
    actunit='K'
    If (inverse.EQ.1) Then
      x=x+273.16
    Else
      x=x-273.16
    End If

  Case('gpm')
    actunit='gpm'
    ! No change
    
  Case('m2/s2')
    actunit='gpm'
    If (inverse.EQ.1) Then
      x=x/0.98
    Else
      x=x*0.98
    End If

  Case('kg/m^2/s')
    actunit='mm/hr'
    If (inverse.EQ.1) Then
      x=x*3600.
    Else
      x=x/3600.
    End If
  
  Case('mm/day')
    actunit='mm/hr'
    If (inverse.EQ.1) Then
      x=x/24.
    Else
      x=x*24.
    End If
 
  Case('mm/hr')
    actunit='mm/hr'
    ! No change

  Case DEFAULT
    Write(6,*) "ERROR: Unknown unit ",trim(inunit)," for conversion."
    Write(6,*) "       Please contact MJT and get him to fix this."
    Stop
  
End Select

If (inverse.EQ.1) Then
  baseunit=actunit
Else
  If ((baseunit.NE.actunit).AND.(actunit.NE.'none')) Then
    Write(6,*) "ERROR: Mismatched units ",trim(baseunit)," and ",trim(actunit)
    Stop
  End If
End If

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine cleans a string from 0's to spaces (32's)
!

Subroutine stringclean(instr)

Implicit None

Character(len=*), intent(inout) :: instr
Integer i,iend

iend=Len(instr)
Do i=1,iend
  If (Ichar(instr(i:i)).EQ.0) Then
    instr(i:i)=Char(32)
  End If
End Do

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine rescales a 1D field
!

Subroutine fieldrescale(inunit,outunit,f,asize)

Implicit None

Character(len=*), intent(in) :: inunit,outunit
Character*80 baseunit
Integer, intent(in) :: asize
Real, dimension(1:asize), intent(inout) :: f
Integer i


! Not a simple multiplication factor.  Need to rescale point-by-point.
If (inunit.NE.outunit) Then
  baseunit=''
  Call baserescale(inunit,f,1,baseunit,asize)
  Call baserescale(outunit,f,-1,baseunit,asize)
End If

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine rescales a 4D field
!

Subroutine arrfieldrescale(inunit,outunit,f,asize)

Implicit None

Character(len=*), intent(in) :: inunit,outunit
Integer, dimension(1:4), intent(in) :: asize
Real, dimension(1:asize(1),1:asize(2),1:asize(3),1:asize(4)), intent(inout) :: f
real, dimension(1:asize(1)*asize(2)*asize(3)*asize(4)) :: dum
Integer i

i=asize(1)*asize(2)*asize(3)*asize(4)
dum(:)=reshape(f,(/i/))
Call fieldrescale(inunit,outunit,dum(:),i)
f(:,:,:,:)=reshape(dum(:),asize(:))

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine advances the date by a specified number of minutes
!

Subroutine advdate(indate,outdate,ot)

Implicit None

Integer, intent(in) :: ot
Integer, dimension(1:6), intent(in) :: indate
Integer, dimension(1:6), intent(out) :: outdate
Integer i,maxday

outdate=indate
outdate(5)=outdate(5)+ot

i=Int(Real(outdate(5))/60.)
outdate(5)=mod(outdate(5),60)

outdate(4)=outdate(4)+i
i=Int(Real(outdate(4))/24.)
outdate(4)=mod(outdate(4),24)

maxday=0
outdate(3)=outdate(3)+i
Do While(outdate(3).GT.maxday)

  Select Case(outdate(2))
  
    Case(1)
      maxday=31
      
    Case(2)
      maxday=28
      If (Mod(Real(outdate(1)),4.).EQ.0.) Then
        maxday=29
      End If

    Case(3)
      maxday=31

    Case(4)
      maxday=30

    Case(5)
      maxday=31

    Case(6)
      maxday=30

    Case(7)
      maxday=31

    Case(8)
      maxday=31

    Case(9)
      maxday=30

    Case(10)
      maxday=31

    Case(11)
      maxday=30

    Case(12)
      maxday=31
  
    Case DEFAULT
      Write(6,*) "ERROR: Internal error in advdate"
      Stop
  
  End Select
  
  If (outdate(3).GT.maxday) Then
    outdate(3)=outdate(3)-maxday
    outdate(2)=outdate(2)+1
    maxday=0
    If (outdate(2).GT.12) Then
      outdate(1)=outdate(1)+1
      outdate(2)=1
    End If
  End If
  
End Do

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine interpolates field from sigma levels to meters
!

Subroutine convertlvl(arrdata,inlvl,outlvl,arrsize)

Implicit None

Integer, dimension(1:3), intent(in) :: arrsize
Real, dimension(1:arrsize(1),1:arrsize(2),1:arrsize(3)), intent(in) :: inlvl
Real, dimension(1:arrsize(3)), intent(in) :: outlvl
Real, dimension(1:arrsize(1),1:arrsize(2),1:arrsize(3)), intent(inout) :: arrdata
Real, dimension(1:arrsize(3)) :: tempdata
Integer i,j,k

Do i=1,arrsize(1)
  Do j=1,arrsize(2)
    ! Interpolate between levels
    Do k=1,arrsize(3)
      Call lineintp(arrdata(i,j,:),inlvl(i,j,:),outlvl(k),tempdata(k),arrsize(3))
    End Do
    arrdata(i,j,:)=tempdata
  End Do
End Do

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the integral of the specified array
!

Subroutine calintegrate(arrdata,xlvl,oval,num)

Implicit None

Integer, intent(in) :: num
Real, dimension(1:num), intent(in) :: xlvl
Real, dimension(1:num), intent(in) :: arrdata
Real, intent(out) :: oval
Integer i

! variable step length?  Use Trapezoidal

oval=0.
Do i=2,num
  oval=oval+(arrdata(i-1)+arrdata(i))*(xlvl(i)-xlvl(i-1))/2.
End Do

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine integrates temperature to determine height of each
! sigma level in meters.  Currently, we assume the hydrostatic approximation.
!

Subroutine calsigmalevel(slvl,olvl,arrsize,tempdata)

Implicit None

Integer, dimension(1:3), intent(in) :: arrsize
Real, dimension(1:arrsize(3)), intent(in) :: slvl
Real, dimension(1:arrsize(1),1:arrsize(2),1:arrsize(3)), intent(out) :: olvl
Real, dimension(1:arrsize(1),1:arrsize(2),1:arrsize(3)), intent(in) :: tempdata
Real, dimension(1:arrsize(3)) :: ilvl
Integer i,j,k

olvl=0.
ilvl=Log(slvl)
Do i=1,arrsize(1)
  Do j=1,arrsize(2)
    Do k=2,arrsize(3)
      Call calintegrate(tempdata(i,j,1:k),ilvl(1:k),olvl(i,j,k),k)
    End Do
  End Do
End Do
olvl=(-287./9.81)*olvl

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine linearly interpolates along an array
!

Subroutine lineintp(dataval,posin,posout,oval,num)

Implicit None

Integer, intent(in) :: num
Real, intent(in) :: posout
Real, intent(out) :: oval
Real, dimension(1:num), intent(in) :: dataval,posin
Integer apos(1),bpos(1)
Real m

If (posout.LT.Minval(posin)) Then
  Write(6,*) "ERROR: Must extrapolate below lowest value"
  Stop
Else If (posout.GT.Maxval(posin)) Then
  Write(6,*) "ERROR: Must extrapolate above highest value"
  Stop
Else
  ! interpolate
  bpos=Minloc(posin,posin.GE.posout)
  apos=Maxloc(posin,posin.LE.posout)
  If (apos(1).EQ.bpos(1)) Then
    oval=dataval(apos(1))
  Else 
    m=(dataval(bpos(1))-dataval(apos(1)))/(posin(bpos(1))-posin(apos(1)))
    oval=m*(posout-posin(apos(1)))+dataval(apos(1))
  End If

End If

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is a bilinear interpolator
!

Real function ipol(dat,serlon,serlat)

Implicit None

Real, intent(in), dimension(2,2) :: dat
Real, intent(in) :: serlon,serlat
Real b,c,d

! Very basic for now.  Later replace with a bicubic.
! Current : z = a + b x + c y + d x y

If ((serlon.EQ.0.).AND.(serlat.EQ.0.)) then

  ipol=dat(1,1)

Else

! After a lot of matrix inverse solving, it turns out that the answer is
! trival (assuming a regular grid).

  d=dat(2,2)-dat(1,2)-dat(2,1)+dat(1,1)
  b=dat(2,1)-dat(1,1)
  c=dat(1,2)-dat(1,1)

  ipol=b*serlon+c*serlat+d*serlon*serlat+dat(1,1)

Endif

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function converts a string to a real (with error checking)
!

Real function sr(instr)

Implicit None

Character(len=*), intent(in) :: instr
Integer ierr

Read(instr,*,iostat=ierr) sr

if (ierr.ne.0) then
  Write(6,*) "ERROR: String "//trim(instr)//" is not a number."
  Stop
end if

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the saturation vapor pressure
!

subroutine getqsat(qsat,temp,ps)

implicit none

real, intent(in) :: temp,ps
real, intent(out) :: qsat
real, dimension(0:220), save :: table
real esatf,tdiff
logical, save :: first=.true.

if (first) then
  table(0:4)=    (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9 /)                                !-146C
  table(5:9)=    (/ 6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9 /)                             !-141C
  table(10:14)=  (/ 36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9 /)                          !-136C
  table(15:19)=  (/ 0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648 /)  !-131C
  table(20:24)=  (/ 0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774 /)  !-126C
  table(25:29)=  (/ 0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081 /)   !-121C
  table(30:34)=  (/ 0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866 /)       !-116C
  table(35:39)=  (/ 0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280 /)         !-111C
  table(40:44)=  (/ 0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951 /)            !-106C
  table(45:49)=  (/ 0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143 /)             !-101C
  table(50:55)=  (/ .001403, .001719, .002101, .002561, .003117, .003784 /)             !-95C
  table(56:63)=  (/ .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658 /) !-87C
  table(64:72)=  (/ .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577 /) !-78C
  table(73:81)=  (/ .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032 /)   !-69C
  table(82:90)=  (/ .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080 /)    !-60C
  table(91:99)=  (/ 1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476 /)    !-51C
  table(100:107)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098 /)            !-43C
  table(108:116)=(/ 10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88 /)    !-34C
  table(117:126)=(/ 27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85 /) !-24C 
  table(127:134)=(/ 77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67 /)      !-16C
  table(135:142)=(/ 171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78 /)  !-8C
  table(143:150)=(/ 353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78 /)   !0C
  table(151:158)=(/ 656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2 /)   !8C
  table(159:166)=(/ 1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3 /)   !16C
  table(167:174)=(/ 1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1 /)   !24C
  table(175:182)=(/ 3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1 /)   !32C
  table(183:190)=(/ 5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7 /)   !40C
  table(191:197)=(/ 7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0 /)         !47C
  table(198:204)=(/ 11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0 /)    !54C
  table(205:211)=(/ 15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0 /)    !61C
  table(212:218)=(/ 21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0 /)    !68C
  table(219:220)=(/ 29845.0, 31169.0 /)
  first=.false.
end if

tdiff=min(max( temp-123.16, 0.), 219.)
esatf=(1.-(tdiff-aint(tdiff)))*table(int(tdiff))+ (tdiff-aint(tdiff))*table(int(tdiff)+1)
qsat=.622*esatf/(ps-esatf)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the nearest valid value (i.e., for filling)
!

subroutine findnear(pxy,i,j,sermsk,rlld,ecodim)

implicit none

integer, dimension(2), intent(out) :: pxy
integer, dimension(2), intent(in) :: ecodim
integer, intent(in) ::i,j
integer nn,nxa,nxb,nya,nyb
integer nxc,nxd
real, dimension(ecodim(1),ecodim(2),2), intent(in) :: rlld
Real, dimension(ecodim(1),ecodim(2)) :: dismsk
logical, dimension(ecodim(1),ecodim(2)), intent(in) :: sermsk
logical, dimension(ecodim(1),ecodim(2)) :: xmsk

nn=1
nxa=i
nxb=i
nxc=i
nxd=i
nya=j
nyb=j
do while (.not.(any(sermsk(nxa:nxb,nya:nyb)).or.any(sermsk(nxc:nxd,nya:nyb))))
  nn=nn*2
  nxa=max(i-nn,1)
  nxb=min(i+nn,ecodim(1))
  nxc=min(i-nn+ecodim(1),ecodim(1))
  nxd=max(i+nn-ecodim(1),1)
  nya=max(j-nn,1)
  nyb=min(j+nn,ecodim(2))
end do
xmsk=.false.
xmsk(nxa:nxb,nya:nyb)=sermsk(nxa:nxb,nya:nyb)
if (nxc<nxd) then
  xmsk(nxc:nxd,nya:nyb)=xmsk(nxc:nxd,nya:nyb).or.sermsk(nxc:nxd,nya:nyb)
  nxa=min(nxa,nxc)
  nxb=max(nxb,nxd)
end if
where(xmsk(nxa:nxb,nya:nyb))
  dismsk(nxa:nxb,nya:nyb)=abs(rlld(i,j,1)-rlld(nxa:nxb,nya:nyb,1))
end where
where ((dismsk(nxa:nxb,nya:nyb)>180.).and.xmsk(nxa:nxb,nya:nyb))
  dismsk(nxa:nxb,nya:nyb)=abs(360.-dismsk(nxa:nxb,nya:nyb))
end where
where (xmsk(nxa:nxb,nya:nyb))
  dismsk(nxa:nxb,nya:nyb)=dismsk(nxa:nxb,nya:nyb)**2+(rlld(i,j,2)-rlld(nxa:nxb,nya:nyb,2))**2
end where
pxy=Minloc(dismsk(nxa:nxb,nya:nyb),xmsk(nxa:nxb,nya:nyb))
pxy(1)=pxy(1)+nxa-1
pxy(2)=pxy(2)+nya-1

return
end
