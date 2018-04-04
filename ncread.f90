! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
!
! THESE SUBROUTINES READ DATA IN NETCDF FORMAT
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determines the size of the array contained in a nc
! file
!

Subroutine getncdims(ncid,ncdim)

use netcdf_m

Implicit None

Integer, intent(in) :: ncid
Integer, dimension(4), intent(out) :: ncdim
Character*80 outname
Character*4, dimension(4) :: varnamelist
Integer ncstatus,i

varnamelist(1)="lon"
varnamelist(2)="lat"
varnamelist(3)="lev"
varnamelist(4)="time"
Do i=1,4
  Call ncfinddimlen(ncid,varnamelist(i),outname,ncdim(i))
Enddo


Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine extracts attribute data from the nc file
!

Subroutine getncdata(ncid,varname,datalab,dataval)

use netcdf_m

Implicit None

Integer, intent(in) :: ncid
Character(len=*), intent(in) :: varname,datalab
Character(len=*), intent(out) :: dataval
Character*80 outname
Integer ncstatus,varid

dataval=''

Call ncfindvarid(ncid,varname,outname,varid)
If (outname.EQ.'') Then
  !Write(6,*) "WARN: Cannot determine variable id ",trim(varname)," (",ncstatus,")"
  return
End If

ncstatus = nf_get_att_text(ncid,varid,datalab,dataval)
If (ncstatus.NE.nf_noerr) Then
  !Write(6,*) "WARN: Cannot read attribute ",trim(datalab)," (",ncstatus,")"
  return
End If

Call stringclean(dataval)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine extracts numerical attribute data from the nc file
!

Subroutine ncgetnumval(ncid,varname,datalab,dataval,ncstatus)

use netcdf_m

Implicit None

Integer, intent(in) :: ncid
Integer, intent(out) :: ncstatus
Character(len=*), intent(in) :: varname,datalab
Real, intent(out) :: dataval
real, dimension(1) :: rvals
Character*80 outname
Integer varid

dataval=0.
Call ncfindvarid(ncid,varname,outname,varid)
ncstatus = nf_get_att_real(ncid,varid,datalab,rvals(1))
dataval=rvals(1)

Return
End
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores (nc) dates in an array
!

Subroutine ncdateconvert(timedate,datearray)

Implicit None

Character(len=*), intent(in) :: timedate
Integer, dimension(1:6), intent(out) :: datearray
integer ierr
Character*3 mthlab
Character*1 tmp1,tmp2,tmp3,tmp4,tmp5
Character*5 tmp6  

mthlab=timedate(4:6)

Select Case(mthlab)

  Case('jan')
    datearray(2)=1

  Case('feb')
    datearray(2)=2

  Case('mar')
    datearray(2)=3

  Case('apr')
    datearray(2)=4

  Case('may')
    datearray(2)=5

  Case('jun')
    datearray(2)=6

  Case('jul')
    datearray(2)=7

  Case('aug')
    datearray(2)=8

  Case('sep')
    datearray(2)=9

  Case('oct')
    datearray(2)=10

  Case('nov')
    datearray(2)=11

  Case('dec')
    datearray(2)=12

  Case DEFAULT
    datearray(2)=-1

End Select

If (datearray(2).EQ.-1) then
  Read(timedate,'(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)',iostat=ierr) datearray(1),tmp1,datearray(2),tmp2,datearray(3), &
      tmp3,datearray(4),tmp4,datearray(5),tmp5,datearray(6)
  if (ierr.ne.0) then
    Read(timedate,'(I4.4,A1,I1.1,A1,I1.1)',iostat=ierr) datearray(1),tmp1,datearray(2),tmp2,datearray(3)
    datearray(4:6)=0
  end if
Else  
  Read(timedate,'(I2.2,A5,I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)',iostat=ierr) datearray(3),tmp6,datearray(1),tmp3,datearray(4), &
      tmp4,datearray(5),tmp5,datearray(6)  
End If

if (ierr.ne.0) then
  Write(6,*) "ERROR: Cannot read nc date."
  Write(6,*) "       Please contact MJT and get him to fix this."
  call finishbanner
  Stop -1
end if

Return
End



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads a field from a nc file
!

Subroutine getncarray(ncid,varname,arrsize,arrdata)

use netcdf_m

Implicit None

Integer, intent(in) :: ncid
Integer, dimension(1:4,1:2), intent(in) :: arrsize
Integer, parameter :: maxdim=4
Integer ndims,varid,ncstatus,i,li,lj,lk,ll,tinx(1),dumnum
Integer, dimension(:), allocatable :: startpos,npos,cid,inx
Integer, dimension(1:maxdim) :: aid,dumcount,duminx
Integer, dimension(0:maxdim) :: ii,ij
Real, dimension(1:arrsize(1,2),1:arrsize(2,2),1:arrsize(3,2),1:arrsize(4,2)), intent(out) :: arrdata
Real, dimension(:,:,:,:), allocatable :: dumarr
real, dimension(1) :: rvals
Real scale,offset
Character(len=*), intent(in) :: varname

Character*80 outname

arrdata = 0.

Call ncfindvarid(ncid,varname,outname,varid)
If (outname.EQ.'') Then
  Write(6,*) "ERROR: Cannot determine var id ",trim(varname)
  call finishbanner
  Stop -1
End If

ncstatus = nf_inq_varndims(ncid,varid,ndims)
If (ncstatus.NE.nf_noerr) Then
  Write(6,*) "ERROR: Cannot determine var dims ",trim(varname)," (",ncstatus,")"
  call finishbanner
  Stop -1
End If

If (ndims.GT.maxdim) Then
  Write(6,*) "WARN: ",trim(varname)," has more than ",maxdim," dimensions."
  Write(6,*) "       This code is untested."
End If

! Need to determine which dimension is which.
! Determine var dim index
Allocate(cid(1:ndims))
ncstatus=nf_inq_vardimid(ncid,varid,cid)
If (ncstatus.NE.nf_noerr) Then
  Write(6,*) "ERROR: Cannot determine var dimid ",trim(varname)," (",ncstatus,")"
  call finishbanner
  Stop -1
End If

! Determine lon, lat, level and time index
Call ncfinddimid(ncid,"lon",outname,aid(1))
Call ncfinddimid(ncid,"lat",outname,aid(2))
Call ncfinddimid(ncid,"level",outname,aid(3))
Call ncfinddimid(ncid,"time",outname,aid(4))

! Cross-reference indices and define slab size
dumcount=1
dumnum=1
duminx=0 ! ii(0) and ij(0) are set to 1.
Allocate(startpos(1:ndims),npos(1:ndims),inx(1:ndims))
Do i=1,ndims
  tinx=Maxloc(aid,aid.EQ.cid(i))
  inx(i)=tinx(1)
  If ((inx(i).LT.1).OR.(inx(i).GT.maxdim)) Then
    Write(6,*) "Non-standard dimension found in ",trim(varname)
    startpos(i)=1
    npos(i)=1
    inx(i)=0 ! ii(0) and ij(0) are set to 1.
  Else
    startpos(i)=arrsize(inx(i),1)
    npos(i)=arrsize(inx(i),2)
    ! Shuffle indices to handle ndims > maxdim
    If (dumnum.GT.maxdim) Then
      Write(6,*) "ERROR: Internal error in getncarray.  dumnum > maxdim."
      Write(6,*) "       Please contact MJT and get him to fix this."
      call finishbanner
      Stop -1
    End If
    dumcount(dumnum)=npos(i)
    duminx(dumnum)=inx(i)
    dumnum=dumnum+1
  End If
End Do

! Create array to accept the slab
Allocate(dumarr(1:dumcount(1),1:dumcount(2),1:dumcount(3),1:dumcount(4)))

! Get the slab
ncstatus = nf_get_vara_real(ncid,varid,startpos,npos,dumarr)
If (ncstatus.NE.nf_noerr) Then
  Write(6,*) "ERROR: Cannot read var ",trim(varname)," data (",ncstatus,")"
  call finishbanner
  Stop -1
End If


! Map slab onto output array
ii=1
ij=1
Do li=1,arrsize(1,2)
  ii(1)=li
  Do lj=1,arrsize(2,2)
    ii(2)=lj  
    Do lk=1,arrsize(3,2)
      ii(3)=lk    
      Do ll=1,arrsize(4,2)
        ii(4)=ll
      
        ij(1:maxdim)=ii(duminx(1:maxdim))
      
        arrdata(ii(1),ii(2),ii(3),ii(4))=dumarr(ij(1),ij(2),ij(3),ij(4))
      End Do
    End Do
  End Do
End Do

Deallocate(startpos,npos,cid,dumarr,inx)

ncstatus = nf_get_att_real(ncid,varid,'scale_factor',rvals(1))
If (ncstatus.NE.nf_noerr) Then
  scale=1.
else
  scale=rvals(1)
End If

ncstatus = nf_get_att_real(ncid,varid,'add_offset',rvals(1))
If (ncstatus.NE.nf_noerr) Then
  offset=0.
else
  offset=rvals(1)
End If

arrdata=scale*arrdata+offset

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determines the date of a netCDF file
!

Subroutine getnctime(ncid,nctmunit,nctmdate)

Implicit None

Integer, intent(in) :: ncid
integer tst
Character(len=*), intent(out) :: nctmunit,nctmdate
character*80 cdum

nctmunit=''
nctmdate=''

Call getncdata(ncid,"time","units",nctmunit)
if (nctmunit.eq.'') return

tst=index(nctmunit,'since')
if (tst.ne.0) then
  tst=scan(nctmunit,' ')-1
  nctmdate=''
  nctmdate(1:19)=nctmunit(tst+8:tst+26)
  cdum=''
  cdum(1:tst)=nctmunit(1:tst)
  nctmunit=cdum
else
  Call getncdata(ncid,"time","time_origin",nctmdate)
end if

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determines the index of a specified time in the nc
! file
!

Subroutine getncindex(ncid,valname,rval,valindex)

Implicit None

Integer, intent(in) :: ncid
Real, intent(in) :: rval
Character(len=*), intent(in) :: valname
Integer, intent(out) :: valindex
Character*80 outname
Real, dimension(:), allocatable :: vallab
Integer, dimension(1) :: vi
Integer valnum,valident


! Find id
Call ncfinddimlen(ncid,valname,outname,valnum)
Allocate(vallab(1:valnum))
Call getncval(ncid,trim(outname),vallab,valnum)
vi=maxloc(vallab,vallab.LE.rval)
If ((vi(1).LT.1).OR.(vi(1).GT.valnum)) then
  vi=minloc(vallab,vallab.GE.rval)
End If
If (vallab(vi(1)).NE.rval) Then
  Write(6,*) "Replace ",trim(valname)," ",rval," with ",vallab(vi(1))
End If

valindex=vi(1)

Deallocate(vallab)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns the dimid for lon, lat, lvl and tim
!

Subroutine ncfinddimlen(ncid,valname,outname,valnum)

use netcdf_m

Implicit none

Integer, intent(in) :: ncid
Character(len=*), intent(in) :: valname
Character(len=*), intent(out) :: outname
Integer, intent(out) :: valnum
Integer ncstatus,valident

Call ncfinddimid(ncid,valname,outname,valident)
If (outname=='') Then
  !Write(6,*) "WARN: Cannot find dimension ",trim(valname)
  valnum=1
  Return
End If

! Find number of elements
ncstatus = nf_inq_dimlen(ncid,valident,valnum)
If (ncstatus/=nf_noerr) Then
  Write(6,*) "ERROR: Cannot determine length of ",trim(valname)," (",ncstatus,")"
  call finishbanner
  Stop -1
End If

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determines the min and max of the lon and lat values
!

Subroutine getnclonlat(ncid,lonlat)

Implicit None

Integer, intent(in) :: ncid
Real, dimension(1:2,1:2), intent(out) :: lonlat
Real, dimension(:), allocatable :: vallab
Integer valnum,valident
Character*80 outname

Call ncfinddimlen(ncid,"lon",outname,valnum)
Allocate(vallab(1:valnum))
Call getncval(ncid,outname,vallab,valnum)
lonlat(1,1)=vallab(1)
lonlat(1,2)=vallab(valnum)
Deallocate(vallab)

Call ncfinddimlen(ncid,"lat",outname,valnum)
Allocate(vallab(1:valnum))
Call getncval(ncid,outname,vallab,valnum)
lonlat(2,1)=vallab(1)
lonlat(2,2)=vallab(valnum)
Deallocate(vallab)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine generates a list of values for a dimension
!

Subroutine getncval(ncid,outname,vallab,valnum)

use netcdf_m

Implicit None

Integer, intent(in) :: ncid,valnum
Character(len=*), intent(in) :: outname
Real, dimension(1:valnum), intent(out) :: vallab
Character*80 actname
Integer ncstatus,valident
integer, dimension(1) :: nstart,ncount

Call ncfindvarid(ncid,outname,actname,valident)

nstart=1
ncount=valnum
ncstatus = nf_get_vara_real(ncid,valident,nstart,ncount,vallab)
If (ncstatus/=nf_noerr) Then
  Write(6,*) "ERROR: Cannot read ",trim(outname)," data (",ncstatus,")"
  call finishbanner
  Stop -1
End If

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determines the min time, max time and time step
! from the nc file.
!

Subroutine getncminmaxtime(ncid,mintime,maxtime,timestep)

Implicit None

Integer, intent(in) :: ncid
Real, intent(out) :: mintime,maxtime,timestep
Real, dimension(:), allocatable :: vallab
Character*80 outname
Integer valnum,valident


Call ncfinddimlen(ncid,"time",outname,valnum)
Allocate(vallab(1:valnum))
Call getncval(ncid,outname,vallab,valnum)
mintime=vallab(1)
maxtime=vallab(valnum)
timestep=(maxtime-mintime)/Real(valnum-1)
Deallocate(vallab)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns the varid for lon, lat, lvl and tim
!

Subroutine ncfindvarid(ncid,valname,outname,valident)

use netcdf_m

Implicit none

Integer, intent(in) :: ncid
Character(len=*), intent(in) :: valname
Character(len=*), intent(out) :: outname
Integer, intent(out) :: valident
Integer, parameter :: maxname=4
Integer, parameter :: maxlist=14
Character*80, dimension(1:maxlist,1:maxname) :: varnamelist
Integer ncstatus,ierr,i,j

varnamelist(:,:)=""
varnamelist(1,1)='u'
varnamelist(1,2)='U'
varnamelist(1,3)='zonal_wnd'
varnamelist(2,1)='v'
varnamelist(2,2)='V'
varnamelist(2,3)='merid_wnd'
varnamelist(3,1)='omega'
varnamelist(3,2)='W'
varnamelist(4,1)='temp'
varnamelist(4,2)='TA'
varnamelist(4,3)='air_temp'
varnamelist(5,1)='mixr'
varnamelist(5,2)='H'
varnamelist(5,3)='rh'
varnamelist(5,4)='mix_rto'
varnamelist(6,1)='zs'
varnamelist(6,2)='TOPO'
varnamelist(6,3)='topo'
varnamelist(7,1)='pblh'
varnamelist(7,2)='ZI'
varnamelist(8,1)='fg'
varnamelist(8,2)='HFLX'
varnamelist(9,1)='zolnd'
varnamelist(9,2)='ZRUF'
varnamelist(10,1)='alb'
varnamelist(10,2)='ALBEDO'
varnamelist(11,1)='pmsl'
varnamelist(11,2)='mslp'
varnamelist(12,1)='tss'
varnamelist(12,2)='sfc_temp'
varnamelist(12,3)='tsu'
varnamelist(13,1)='ps'
varnamelist(13,2)='sfc_pres'
varnamelist(14,1)='hgt'
varnamelist(14,2)='zht'

outname=""

ncstatus = nf_inq_varid(ncid,trim(valname),valident)
If (ncstatus.EQ.nf_noerr) Then
  outname=valname
Else 

  Call ncfinddimid(ncid,valname,outname,valident)
  If (outname.NE.'') Then
    ncstatus = nf_inq_varid(ncid,outname,valident)
    If (ncstatus.NE.nf_noerr) Then
      Write(6,*) "ERROR: Cannot read id for ",trim(outname)
      call finishbanner
      Stop -1
    End If
  Else

    i=1
    j=1
    ierr=1
    Do While (ierr.NE.0)
      If (valname.EQ.varnamelist(i,j)) Then
        ierr=0
      Else
        i=i+1
      End If
  
      If (i.GT.maxlist) Then
        i=1
        j=j+1
      End If
  
      If (j.GT.maxname) Then
        outname=''
        valident=-1
        Return
      End If

    End Do


    j=1
    ierr=1
    Do While (ierr.NE.0)
      ncstatus = nf_inq_varid(ncid,trim(varnamelist(i,j)),valident)
      If (ncstatus.EQ.nf_noerr) Then
        ierr=0
        outname=varnamelist(i,j)
      Else
        j=j+1
      End If
  
      If (j.GT.maxname) Then
        outname=''
        valident=-1
        ierr=0
      End If
  
    End Do  
  End If
End If

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns the varid for lon, lat, lvl and tim
!

Subroutine ncfinddimid(ncid,valname,outname,valident)

use netcdf_m

Implicit none

Integer, intent(in) :: ncid
Character(len=*), intent(in) :: valname
Character(len=*), intent(out) :: outname
Integer, intent(out) :: valident
Integer, parameter :: maxname=5
Character*80, dimension(1:4,1:maxname) :: varnamelist
Integer ncstatus,ierr,i,j

varnamelist(:,:)=""
varnamelist(1,1)="lon"
varnamelist(2,1)="lat"
varnamelist(3,1)="lev"
varnamelist(4,1)="time"
varnamelist(1,2)="longitude"
varnamelist(2,2)="latitude"
varnamelist(3,2)="level"
varnamelist(3,3)="pres"
varnamelist(3,4)="lvl"
varnamelist(3,5)="sigma_level"

ncstatus = nf_inq_dimid(ncid,valname,valident)
If (ncstatus==nf_noerr) Then
  outname=valname
  Return
End if

outname=""

i=1
j=1
ierr=1
Do While (ierr.NE.0)
  If (valname==varnamelist(i,j)) Then
    ierr=0
  Else
    i=i+1
  End If
  
  If (i==5) Then
    i=1
    j=j+1
  End If
  
  If (j>maxname) Then
    outname=''
    valident=-1
    Return
  End If

End Do

j=1
ierr=1
Do While (ierr/=0)

  ncstatus = nf_inq_dimid(ncid,trim(varnamelist(i,j)),valident)
  If (ncstatus==nf_noerr) Then
    ierr=0
    outname=varnamelist(i,j)
  Else
    j=j+1
  End If
  
  If (j>maxname) Then
    ierr=0
    outname=''
    valident=-1
  End If
  
End Do  

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine rescales the levels from sigma units to meters above
! the terrain height
!

Subroutine nccallvlheight(ncid,outlvl,lvlnum,otype)

Implicit None

Integer, intent(in) :: ncid,lvlnum
Real, dimension(1:lvlnum), intent(out) :: outlvl
Character(len=*), intent(in) :: otype
Real, dimension(:,:,:), allocatable :: arrdata
Real, dimension(1:lvlnum) :: tlvl
Integer, dimension(1:4) :: ncdim
Integer it,i,j,itop(1),ibot(1)
Real x,holdlvl

Call getncdims(ncid,ncdim)
If (lvlnum.NE.ncdim(3)) Then
  Write(6,*) 'ERROR: Internal error in nccallvlheight'
  call finishbanner
  Stop -1
End If

! Define height array
Allocate(arrdata(1:ncdim(1),1:ncdim(2),1:ncdim(3)))

! Step through all time steps (skip first time step)
outlvl=0.
Do it=2,ncdim(4)
  Call ncgetlvlheight(ncid,it,arrdata,ncdim,otype)
  If (it.EQ.2) Then
    itop=Maxloc(arrdata(1,1,:))
    ibot=Minloc(arrdata(1,1,:))
    outlvl(ibot)=0.
    outlvl(itop)=9999999
  End If
   
  ! Calculate average levels
  Do i=1,ncdim(1)
    Do j=1,ncdim(2)
      tlvl=outlvl+arrdata(i,j,:)/Real(ncdim(1)*ncdim(2)*(ncdim(4)-1))
      tlvl(ibot(1))=Max(outlvl(ibot(1)),arrdata(i,j,ibot(1)))
      tlvl(itop(1))=Min(outlvl(itop(1)),arrdata(i,j,itop(1)))
      outlvl=tlvl
    End Do
  End Do
        
End Do
Deallocate(arrdata)

! Round levels
Do i=1,ncdim(3)
  holdlvl=outlvl(i)
  x=Int(Log10(outlvl(i))-2.)
  If (x.LT.0.) x=0.
  x=10**x
  outlvl(i)=Nint(outlvl(i)/x)*x
  If ((i.EQ.ibot(1)).AND.(holdlvl.GT.outlvl(i))) outlvl(i)=outlvl(i)+x
  If ((i.EQ.itop(1)).AND.(holdlvl.LT.outlvl(i))) outlvl(i)=outlvl(i)-x
End Do

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine converts from sigma levels to meters for a specified
! time step.
!

Subroutine nccalsigmaheight(ncid,it,nclvl,outdata,arrsize)

Implicit None

Integer, intent(in) :: ncid,it
Integer, dimension(1:3), intent(in) :: arrsize
Real, dimension(1:arrsize(3)), intent(in) :: nclvl
Real, dimension(1:arrsize(1),1:arrsize(2),1:arrsize(3)), intent(out) :: outdata
Integer, dimension(1:4,1:2) :: tempsize
Real, dimension(1:arrsize(1),1:arrsize(2),0:arrsize(3),1) :: tempdata
Real, dimension(1:arrsize(1),1:arrsize(2),0:arrsize(3),1) :: mixrdata
Real, dimension(1:arrsize(1),1:arrsize(2),0:arrsize(3)) :: tempout
Real, dimension(0:arrsize(3)) :: lvldata
Character*80 outname,inunit
Integer valid,ncstatus
Real offset

! Define hyperslab
tempsize=1.
tempsize(1:3,2)=arrsize
tempsize(4,1)=it

tempdata = 0. ! for Cray compiler
mixrdata = 0. ! for Cray compiler

! get temp data
Call getncarray(ncid,'temp',tempsize,tempdata(:,:,1:tempsize(3,2),:))
Call getncdata(ncid,'temp','units',inunit)
Call arrfieldrescale(inunit,'K',tempdata(:,:,1:tempsize(3,2),:),tempsize(:,2))

! Get mixr data
Call getncarray(ncid,'mixr',tempsize,mixrdata(:,:,1:tempsize(3,2),:))
Call getncdata(ncid,'mixr','units',inunit)
Call arrfieldrescale(inunit,'kg/kg',mixrdata(:,:,1:tempsize(3,2),:),tempsize(:,2))

Call ncfindvarid(ncid,'tscrn',outname,valid)
If (valid.NE.-1) Then
  tempsize(3,2)=1
  Call getncarray(ncid,'tscrn',tempsize,tempdata(:,:,0,:))
  Call getncdata(ncid,'tscrn','units',inunit)
  Call arrfieldrescale(inunit,'K',tempdata(:,:,0,:),tempsize(:,2))
  
  Call ncgetnumval(ncid,'tscrn','add_offset',offset,ncstatus)
  If (tempdata(1,1,0,1).EQ.offset) Then
    Write(6,*) "tscrn is not defined at time step ",it   
    tempdata(:,:,0,:)=tempdata(:,:,1,:)
  End If
Else
  Write(6,*) "tscrn is not defined in nc file"
  tempdata(:,:,0,:)=tempdata(:,:,1,:)
End If

Call ncfindvarid(ncid,'qgscrn',outname,valid)
If (valid.NE.-1) Then
  tempsize(3,2)=1
  Call getncarray(ncid,'qgscrn',tempsize,mixrdata(:,:,0,:))
  Call getncdata(ncid,'qgscrn','units',inunit)
  Call arrfieldrescale(inunit,'kg/kg',mixrdata(:,:,0,:),tempsize(:,2))
  
  Call ncgetnumval(ncid,'qgscrn','add_offset',offset,ncstatus)
  If (tempdata(1,1,0,1).EQ.offset) Then
    Write(6,*) "qgscrn is not defined at time step ",it   
    mixrdata(:,:,0,:)=mixrdata(:,:,1,:)
  End If
Else
  Write(6,*) "qgscrn is not defined in nc file"
  mixrdata(:,:,0,:)=mixrdata(:,:,1,:)
End If


! Convert temperture to virtual temperature
tempdata=tempdata*(1.+(461.5/287.-1.)*mixrdata)

! include sigma=1
tempsize(3,2)=arrsize(3)+1
lvldata(1:arrsize(3))=nclvl
lvldata(0)=1.
Call calsigmalevel(lvldata,tempout,tempsize(1:3,2),tempdata(:,:,:,1))
outdata=tempout(:,:,1:arrsize(3))

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the levels from a nc file
!

Subroutine ncgetlvlheight(ncid,it,outlvl,ncdim,otype)

Implicit None

Integer, intent(in) :: ncid,it
Integer, dimension(1:4), intent(in) :: ncdim
Real, dimension(1:ncdim(1),1:ncdim(2),1:ncdim(3)), intent(out) :: outlvl
Character(len=*), intent(in) :: otype
Real, dimension(:), allocatable :: nclvl
Integer i,j
Character*80 utype
   
Call getncdata(ncid,'lev','units',utype)
Allocate(nclvl(1:ncdim(3)))
Call getncval(ncid,'lev',nclvl,ncdim(3))

If (utype.NE.otype) Then

  outlvl=0.

  Select Case(utype)

    Case('sigma_level')
      ! Determine current sigma levels
      Call nccalsigmaheight(ncid,it,nclvl,outlvl,ncdim(1:3))
  
      Case DEFAULT
      Write(6,*) "ERROR: Cannot convert from unit ",trim(utype)
      Write(6,*) "       Please contact MJT and get him to fix this."
      call finishbanner
      Stop -1
      
  End Select

  Select Case(otype)
  
    Case('meters')
      ! Do nothing.
    
    Case DEFAULT
      Write(6,*) "ERROR: Cannot convert to unit ",trim(otype)
      Write(6,*) "       Please contact MJT and get him to fix this."
      call finishbanner
      Stop -1
  
  End Select

Else
  Do i=1,ncdim(1)
    Do j=1,ncdim(2)
      outlvl(i,j,:)=nclvl
    End Do
  End Do
End If

Deallocate(nclvl)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine converts velocity units (e.g., Pa/s to m/s).
!

Subroutine ncvelconvert(ncid,arrdata,arrsize,utype,otype)

Implicit None

Integer, intent(in) :: ncid
Integer, dimension(1:4,1:2), intent(in) :: arrsize
Real, dimension(1:arrsize(1,2),1:arrsize(2,2),1:arrsize(3,2),1:arrsize(4,2)), intent(inout) :: arrdata
Character(len=*), intent(in) :: utype,otype
Real, dimension(1:arrsize(1,2),1:arrsize(2,2),1:arrsize(3,2),1:arrsize(4,2)) :: tempdata
!Real, dimension(1:arrsize(1,2),1:arrsize(2,2),1:arrsize(3,2),1:arrsize(4,2)) :: mixrdata
Real, dimension(1:arrsize(1,2),1:arrsize(2,2),1,1:arrsize(4,2)) :: presdata
Real, dimension(1:arrsize(1,2),1:arrsize(2,2),1,1:arrsize(4,2)) :: dotpresdata
Real, dimension(1:arrsize(3,2)) :: sigmalvl
Integer, dimension(1:4,1:2) :: tsize
Integer i,j,k,l
Character*80 inunit

If (utype.NE.otype) Then
  Select Case(utype)

    Case('Pa/s')
      Call getncarray(ncid,'temp',arrsize,tempdata)
      Call getncdata(ncid,'temp','units',inunit)
      Call arrfieldrescale(inunit,'K',tempdata,arrsize(:,2))      
      !Call getncarray(ncid,'mixr',arrsize,mixrdata)
      !Call getncdata(ncid,'mixr','units',inunit)
      !Call arrfieldrescale(inunit,'kg/kg',mixrdata,arrsize(:,2))      

!     !Calculate virutal temperature ...
      !tempdata=tempdata*(1.+(461.5/287.-1.)*mixrdata)

      tsize=arrsize
      tsize(3,:)=1
      Call getncarray(ncid,'ps',tsize,presdata)
      Call getncdata(ncid,'ps','units',inunit)
      Call arrfieldrescale(inunit,'Pa',presdata,tsize(:,2))
      Call getncarray(ncid,'dpsdt',tsize,dotpresdata)
      Call getncdata(ncid,'dpsdt','units',inunit)
      Call arrfieldrescale(inunit,'hPa/day',dotpresdata,tsize(:,2))
      dotpresdata=dotpresdata/864. ! Convert to Pa/s
      Call getncval(ncid,'lev',sigmalvl,arrsize(3,2))
      Do k=1,arrsize(3,2)
	    arrdata(:,:,k,:)=(-287./9.81)*tempdata(:,:,k,:)*arrdata(:,:,k,:)/(sigmalvl(k)*presdata(:,:,1,:)) &
                        +sigmalvl(k)*dotpresdata(:,:,1,:)/presdata(:,:,1,:)
      End Do
  
    Case DEFAULT
      Write(6,*) "ERROR: Cannot convert from velocity unit ",trim(utype)
      Write(6,*) "       Please contact MJT and get him to fix this"
      call finishbanner
      Stop -1
  
  End Select
  
  Select Case(otype)
  
    Case('m/s')
      ! Do nothing
      
    Case DEFAULT
      Write(6,*) "ERROR: Cannot convert to velocity unit ",trim(otype)
      Write(6,*) "       Please contact MJT and get him to fix this"
      call finishbanner
      Stop -1
  
  End Select
End If

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns data from nc file as specified by the HPAC
! metadata
!

Subroutine getmeta(ncid,varname,arrdata,arrsize)

Implicit None

Integer, intent(in) :: ncid
Integer, dimension(1:4,1:2), intent(in) :: arrsize
Character(len=*), dimension(1:2), intent(in) :: varname
Real, dimension(1:arrsize(1,2),1:arrsize(2,2),1:arrsize(3,2),1:arrsize(4,2)), intent(out) :: arrdata
Real, dimension(:,:,:,:), allocatable :: tdata
Real, dimension(:,:,:), allocatable :: temp
Real, dimension(:,:), allocatable :: zs
Real, dimension(:), allocatable :: sigma,sigmaout
Character*80 outname,inunit,outunit,utype
Character*80, dimension(1:2) :: chartemp
Integer, dimension(1:4) :: ncsize,tempsize
Integer valid,i,j,k,indx

Select Case(varname(1))

  Case('BOWEN')
    Allocate(tdata(arrsize(1,2),arrsize(2,2),arrsize(3,2),arrsize(4,2)))
    Call getncarray(ncid,'fg',arrsize,arrdata)
    Call getncarray(ncid,'eg',arrsize,tdata)
    Call getncdata(ncid,'fg','units',outunit)
    Call getncdata(ncid,'eg','units',inunit)
    Call arrfieldrescale(inunit,outunit,tdata,arrsize(:,2))
    where (tdata.ne.0.)
      arrdata=arrdata/tdata
    elsewhere
      arrdata=arrdata/0.001
    end where
    Deallocate(tdata)
    Return
    
  Case('zht')
    Call ncfindvarid(ncid,varname(1),outname,valid)
    If (outname.NE.'') then
      Call getncarray(ncid,outname,arrsize,arrdata)
      Call getncdata(ncid,outname,'units',inunit)
    Else
      Call ncfindvarid(ncid,'zs',outname,valid)
      If (outname.NE.'') then
        If (arrsize(3,2).EQ.1) Then
	  Write(6,*) "Using topography to calculate surface geopotential height."
	  chartemp(1)=outname
      chartemp(2)='m'
	  call getncarray(ncid,outname,arrsize,arrdata)
	  arrdata=arrdata*9.80616
	  inunit='m2/s2'
	  Return   
        Else
          Write(6,*) "ERROR: Surface geopotential height requested"
          Write(6,*) "       at multiple levels."
          call finishbanner
          Stop -1
	End if
      Else
        Write(6,*) "ERROR: Cannot find surface geopotential height"
        Write(6,*) "       or topography in input file."
        call finishbanner
        Stop -1
      End if
    End if
    
  Case DEFAULT
    Call ncfindvarid(ncid,varname(1),outname,valid)
    If (outname.EQ.'') Then
      Write(6,*) "ERROR: Cannot locate "//trim(varname(1))//" in nc file."
      call finishbanner
      Stop -1
    End if
    If (outname.NE.varname(1)) Write(6,*) "Located "//trim(varname(1))//" as "//trim(outname)
    Call getncarray(ncid,outname,arrsize,arrdata)
    Call getncdata(ncid,outname,'units',inunit)
    if (inunit.eq.'') inunit=varname(2) ! bug fix

End Select

!*******************************************************************
! Special cases (to handle bugs....)
If ((outname.EQ.'topo').AND.(inunit.EQ.'degrees_east')) then
  Write(6,*) "Replace topo unit degrees_east with m."
  inunit='m'
End if
    
If ((outname.EQ.'hgt').AND.(inunit.EQ.'m')) then
  Write(6,*) "Replace hgt unit m with gpm."    
  inunit='gpm'
End if
!*******************************************************************


If (inunit.NE.varname(2)) Then
  Write(6,*) "Converting ",trim(inunit)," to ",trim(varname(2))
  If (inunit.EQ.'Pa/s') Then
    Call ncvelconvert(ncid,arrdata,arrsize,inunit,'m/s')
  Else
    ! Convert units
    Call arrfieldrescale(inunit,varname(2),arrdata,arrsize(:,2))
  End If
End If
   

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine matches metadata in the MEDOC file to metadata in
! the nc file
!

Subroutine matchmeta(ncid,varname,varnum)

Implicit None

Integer, intent(in) :: ncid
Integer, intent(inout) :: varnum
Character(len=*), dimension(1:varnum,1:2), intent(inout) :: varname
Character*80 outname
Integer i,oldvarnum,valid

Do i=1,varnum

  Select Case(varname(i,1))
  
    Case('BOWEN')
      Call ncfindvarid(ncid,'fg',outname,valid)
      If (outname.NE.'') Call ncfindvarid(ncid,'eg',outname,valid)
      
    Case DEFAULT 
      Call ncfindvarid(ncid,varname(i,1),outname,valid)

  End Select

  If (outname.EQ.'') varname(i,:)=''
End Do

! Clean-up
oldvarnum=varnum
varnum=0
Do i=1,oldvarnum
  If (varname(i,1).NE.'') then
    varnum=varnum+1
    If (varnum.NE.i) Then
      varname(varnum,:)=varname(i,:)
      varname(i,:)=''
    End If
  End If
End Do

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads topography data
!

Subroutine readtopography(topounit,toponame,ecodim,lonlat,schmidt,dsx,header)

use netcdf_m

Implicit None

Integer, intent(in) :: topounit
Integer, dimension(2), intent(out) :: ecodim
Integer ierr,ncid,varid
Character(len=*), intent(in) :: toponame
Character*47, intent(out) :: header
Real, dimension(1:2), intent(out) :: lonlat
Real, intent(out) :: schmidt,dsx
real, dimension(1) :: rvals

ierr=nf_open(toponame,nf_nowrite,ncid)
if (ierr==nf_noerr) then
  ierr=nf_get_att_real(ncid,nf_global,'lon0',lonlat(1))
  if (ierr/=nf_noerr) then
    write(6,*) "ERROR reading lon0"
    call finishbanner
    stop -1
  end if
  ierr=nf_get_att_real(ncid,nf_global,'lat0',lonlat(2))
  if (ierr/=nf_noerr) then
    write(6,*) "ERROR reading lat0"
    call finishbanner
    stop -1
  end if
  ierr=nf_get_att_real(ncid,nf_global,'schmidt',rvals(1))
  if (ierr/=nf_noerr) then
    write(6,*) "ERROR reading schmidt"
    call finishbanner
    stop -1
  end if
  schmidt=rvals(1)
  ierr=nf_inq_dimid(ncid,'longitude',varid)
  ierr=nf_inq_dimlen(ncid,varid,ecodim(1))
  ierr=nf_inq_dimid(ncid,'latitude',varid)
  ierr=nf_inq_dimlen(ncid,varid,ecodim(2))
  ierr=nf_close(ncid)
else
  Open(topounit,FILE=toponame,FORM='formatted',STATUS='old',IOSTAT=ierr)
  Read(topounit,*,IOSTAT=ierr) ecodim(1),ecodim(2),lonlat(1),lonlat(2),schmidt,dsx,header
  Close(topounit)

  If (ierr.NE.0) then
    Write(6,*) "ERROR: Cannot read file ",trim(toponame)
    call finishbanner
    Stop -1
  End if
end if

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads the land/sea mask from topography data
!

Subroutine gettopols(topounit,toponame,lsmsk,ecodim)

use netcdf_m

Implicit None

Integer, intent(in) :: topounit
Integer, dimension(2), intent(in) :: ecodim
integer, dimension(3) :: spos,npos
Integer ierr,ilout,varid,ncid
Character(len=*), intent(in) :: toponame
Character*47 :: dc
Character*9 formout
Real, dimension(6) :: dr
Real, dimension(ecodim(1),ecodim(2)), intent(out) :: lsmsk

ierr=nf_open(toponame,nf_nowrite,ncid)
if (ierr==0) then
  spos=1
  npos(1)=ecodim(1)
  npos(2)=ecodim(2)
  npos(3)=1
  ierr=nf_inq_varid(ncid,'lsm',varid)
  ierr=nf_get_vara_real(ncid,varid,spos,npos,lsmsk)
  ierr=nf_close(ncid)
else
  ilout=Min(ecodim(1),30) ! To be compatiable with terread

  Open(topounit,FILE=toponame,FORM='formatted',STATUS='old',IOSTAT=ierr)
  Read(topounit,*,IOSTAT=ierr) dr(1:6),dc
  Write(formout,'("(",i3,"f7.0)")',IOSTAT=ierr) ilout
  Read(topounit,formout,IOSTAT=ierr) lsmsk ! Dummy topo data
  Write(formout,'("(",i3,"f4.1)")',IOSTAT=ierr) ilout
  Read(topounit,formout,IOSTAT=ierr) lsmsk ! Read ls mask
  Close(topounit)
end if

If (ierr/=0) then
  Write(6,*) "ERROR: Cannot read file ",trim(toponame)
  call finishbanner
  Stop -1
End if

lsmsk=1.-lsmsk

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads the land/sea mask from topography data
!

Subroutine gettopohgt(topounit,toponame,hgt,lsmsk,ecodim)

use netcdf_m

Implicit None

Integer, intent(in) :: topounit
Integer, dimension(2), intent(in) :: ecodim
integer, dimension(3) :: spos,npos
Integer ierr,ilout,varid,ncid
Character(len=*), intent(in) :: toponame
Character*47 :: dc
Character*9 formout
Real, dimension(6) :: dr
Real, dimension(ecodim(1),ecodim(2)), intent(out) :: hgt,lsmsk

ierr=nf_open(toponame,nf_nowrite,ncid)
if (ierr==0) then
  spos=1
  npos(1)=ecodim(1)
  npos(2)=ecodim(2)
  npos(3)=1
  ierr=nf_inq_varid(ncid,'zs',varid)
  if (ierr/=0) then
    write(6,*) "ERROR reading zs ",ierr
    call finishbanner
    stop -1
  end if
  ierr=nf_get_vara_real(ncid,varid,spos,npos,hgt)
  if (ierr/=0) then
    write(6,*) "ERROR reading zs ",ierr
    call finishbanner
    stop -1
  end if
  hgt=9.80616*hgt
  ierr=nf_inq_varid(ncid,'lsm',varid)
  ierr=nf_get_vara_real(ncid,varid,spos,npos,lsmsk)
  ierr=nf_close(ncid)
else
  ilout=Min(ecodim(1),30) ! To be compatiable with terread

  Open(topounit,FILE=toponame,FORM='formatted',STATUS='old',IOSTAT=ierr)
  Read(topounit,*,IOSTAT=ierr) dr(1:6),dc
  Write(formout,'("(",i3,"f7.0)")',IOSTAT=ierr) ilout
  Read(topounit,formout,IOSTAT=ierr) hgt
  Write(formout,'("(",i3,"f4.1)")',IOSTAT=ierr) ilout
  Read(topounit,formout,IOSTAT=ierr) lsmsk ! Read ls mask
  Close(topounit)
end if

If (ierr/=0) then
  Write(6,*) "ERROR: Cannot read file ",trim(toponame)
  call finishbanner
  Stop -1
End if

lsmsk=1.-lsmsk

Return
End