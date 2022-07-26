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
! THESE SUBROUTINES WRITE DATA IN NETCDF FORMAT
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises the output NetCDF file for 3D
!

Subroutine ncinit(ncidarr,outfile,dimnum,dimvar,outputunit,adate)

Integer, dimension(1:3), intent(in) :: dimnum
Integer, dimension(0:4), intent(out) :: ncidarr
Integer, dimension(1:4), intent(out) :: dimvar
Integer, dimension(1:6), intent(in) :: adate
Character(len=*), dimension(1:2), intent(in) :: outputunit
Character(len=*), intent(in) :: outfile

call ncinitgen(ncidarr,outfile,dimnum,dimvar,outputunit,adate,'even')

return
end

Subroutine ncinitgen(ncidarr,outfile,dimnum,dimvar,outputunit,adate,latsp)

use netcdf_m

Implicit None

Integer, dimension(1:3), intent(in) :: dimnum
Integer, dimension(0:4), intent(out) :: ncidarr
Integer, dimension(1:4), intent(out) :: dimvar
Integer, dimension(1:6), intent(in) :: adate
Character(len=*), intent(in) :: outfile,latsp
Character(len=*), dimension(1:2), intent(in) :: outputunit
Integer status,i,j,strlen
Integer, dimension(1:4) :: dims
Character*80, dimension(1:4,1:4) :: desc
Character*80, dimension(2:4) :: vardesc
Character*80 timedesc
Character*3 mthnam
Character*3 findmonth

ncidarr=0
dims(1:3)=dimnum(1:3)
dims(4)=nf_unlimited

vardesc(2)="long_name"
vardesc(3)="units"
vardesc(4)="point_spacing"

desc(1,1)="lon"
desc(1,2)="longitude"
desc(1,3)="degrees_east"
desc(1,4)="even"
desc(2,1)="lat"
desc(2,2)="latitude"
desc(2,3)="degrees_north"
desc(2,4)="even"
desc(3,1)="lev"
desc(3,2)="level"
desc(3,3)="layer"
desc(3,4)="uneven"
desc(4,1)="time"
desc(4,2)="time"
desc(4,3)="hours"
desc(4,4)="even"

desc(3,1)=outputunit(1)
desc(3,3)=outputunit(2)
desc(2,4)=latsp

If (adate(1)/=0) then
  Write(timedesc,'("hours since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') adate(:)
  desc(4,3)=timedesc
End if

! Create NetCDF file
#ifdef usenc3
status=nf_create(outfile,nf_64bit_offset,ncidarr(0))
#else
status=nf_create(outfile,nf_netcdf4,ncidarr(0))
#endif
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error opening NetCDF file (",status,")"
  call finishbanner
  Stop -1
End If

! Define dimensions

Do i=1,4
  If (dims(i)/=1) Then
    status=nf_def_dim(ncidarr(0),desc(i,1),dims(i),ncidarr(i))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      call finishbanner
      Stop -1
    End If
    status=nf_def_var(ncidarr(0),trim(desc(i,1)),nf_float,1,ncidarr(i),dimvar(i))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      call finishbanner
      Stop -1
    End If
    Do j=2,4
      strlen=Len_trim(desc(i,j))
      status=nf_put_att_text(ncidarr(0),dimvar(i),trim(vardesc(j)),strlen,trim(desc(i,j)))
      If (status /= nf_noerr) Then
        Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,j))
        call finishbanner
        Stop -1
      End If
    End Do
  End If
End Do

If (adate(1)/=0) Then
  mthnam=findmonth(adate(2))
  Write(timedesc,'(I2.2,"-",A3,"-",I4.4,":",I2.2,":",I2.2,":",I2.2)') adate(3),mthnam,adate(1),adate(4),adate(5),adate(6)
  strlen=Len_trim(timedesc)
  status=nf_put_att_text(ncidarr(0),dimvar(4),"time_origin",strlen,trim(timedesc))
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining time in NetCDF file (",status,")"
    call finishbanner
    Stop -1
  End If
End If


Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises the output NetCDF file for CC
!

Subroutine ncinitcc(ncidarr,outfile,dimnum,dimvar,adate)

Implicit None

Integer, dimension(1:3), intent(in) :: dimnum
Integer, dimension(0:4), intent(out) :: ncidarr
Integer, dimension(1:4), intent(out) :: dimvar
Integer, dimension(1:6), intent(in) :: adate
Character(len=*), intent(in) :: outfile

call ncinitccgen(ncidarr,outfile,dimnum,dimvar,adate,'real')

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises the output NetCDF file for CC
!

Subroutine ncinitccgen(ncidarr,outfile,dimnum,dimvar,adate,mode)

use netcdf_m

Implicit None

Integer, dimension(1:3), intent(in) :: dimnum
Integer, dimension(0:4), intent(out) :: ncidarr
Integer, dimension(1:4), intent(out) :: dimvar
Integer, dimension(1:6), intent(in) :: adate
Character(len=*), intent(in) :: outfile,mode
Integer status,i,strlen,vtype
Integer, dimension(1:4) :: dims
Character*80, dimension(1:4,1:3) :: desc
Character*80 timedesc
Character*3 mthnam
Character*3 findmonth

ncidarr=0
dims(1:3)=dimnum(1:3)
dims(4)=nf_unlimited

desc(1,1)="longitude"
desc(1,2)="degrees_east"
desc(1,3)="even"
desc(2,1)="latitude"
desc(2,2)="degrees_north"
desc(2,3)="even"
desc(3,1)="lev"
desc(3,2)="sigma_level"
desc(3,3)="uneven"
desc(4,1)="time"
desc(4,2)="---"
desc(4,3)="even"

If (adate(1).NE.0) then
  Write(timedesc,'("minutes since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') adate(:)
  desc(4,2)=timedesc
Else
  If (adate(2)==1) desc(4,2)='months'
  If (adate(3)==1) desc(4,2)='days'
  If (adate(4)==1) desc(4,2)='hours'
  If (adate(5)==1) desc(4,2)='minutes'
  If (adate(6)==1) desc(4,2)='seconds'  
End if

! Create NetCDF file
#ifdef usenc3
status=nf_create(outfile,nf_64bit_offset,ncidarr(0))
#else
status=nf_create(outfile,nf_netcdf4,ncidarr(0))
#endif
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error opening NetCDF file (",status,")"
  call finishbanner
  Stop -1
End If

! Define dimensions
Do i=1,4
  If (dims(i)/=1) Then
    status=nf_def_dim(ncidarr(0),desc(i,1),dims(i),ncidarr(i))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      call finishbanner
      Stop -1
    End If
    if ((i==4).and.(mode=='int')) then
      vtype=nf_int
    else
      vtype=nf_float
    end if
    status=nf_def_var(ncidarr(0),desc(i,1),vtype,1,ncidarr(i),dimvar(i))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      call finishbanner
      Stop -1
    End If
    If (i/=3) Then
      strlen=Len_trim(desc(i,2))
      status=nf_put_att_text(ncidarr(0),dimvar(i),"long_name",strlen,trim(desc(i,2)))
      If (status /= nf_noerr) Then
        Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
        call finishbanner
        Stop -1
      End If
      !status=nf_put_att_text(ncidarr(0),dimvar(i),"positive",4,"down")
      !If (status /= nf_noerr) Then
      !  Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      !  call finishbanner
      !  Stop -1
      !End If
    End if
    
    strlen=Len_trim(desc(i,2))
    status=nf_put_att_text(ncidarr(0),dimvar(i),"units",strlen,trim(desc(i,2)))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      call finishbanner
      Stop -1
    End If
    strlen=Len_trim(desc(i,3))
    status=nf_put_att_text(ncidarr(0),dimvar(i),"point_spacing",strlen,trim(desc(i,3)))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,3))
      call finishbanner
      Stop -1
    End If
    
  End If
End Do

If (adate(1)/=0) Then
  mthnam=findmonth(adate(2))
  Write(timedesc,'(I2.2,"-",A3,"-",I4.4,":",I2.2,":",I2.2,":",I2.2)') adate(3),mthnam,adate(1),adate(4),adate(5),adate(6)
  strlen=Len_trim(timedesc)
  status=nf_put_att_text(ncidarr(0),dimvar(4),"time_origin",strlen,trim(timedesc))
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining time in NetCDF file (",status,")"
    call finishbanner
    Stop -1
  End If
End If


Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine defines fields in the NetCDF file
!

Subroutine ncaddvargen(ncidarr,elemdesc,numtype,numdim,varid,sc,of)

use netcdf_m

Implicit None

Integer, intent(in) :: numtype,numdim
Integer, dimension(0:4), intent(in) :: ncidarr
Integer, intent(out) :: varid
Character(len=*), dimension(1:3), intent(in) :: elemdesc
Integer, dimension(1:numdim) :: dimtype
Integer strlen,ierr
Integer status
integer*2, dimension(1) :: ivals2
Real, intent(in) :: sc,of
real, dimension(1) :: rvals

Select Case (numdim)
  Case(2)
    dimtype(1:2)=ncidarr(1:2)
  Case(3)
    dimtype(1:2)=ncidarr(1:2)
    dimtype(3)=ncidarr(4)
  Case(4)
    dimtype(1:4)=ncidarr(1:4)
  Case DEFAULT
    Write(6,*) "ERROR: Unsupported number of field dimensions :",numdim
    Write(6,*) "       Please contact MJT and get him to fix it"
    call finishbanner
    Stop -1
End Select

status=nf_def_var(ncidarr(0),elemdesc(1),numtype,numdim,dimtype,varid)
If (status /= nf_noerr) Then
  If (status == nf_enameinuse) then
    Write(6,*) "WARN:  Variable ",trim(elemdesc(1))," already exists."
    status = nf_inq_varid(ncidarr(0),elemdesc(1),varid)
    If (status == nf_noerr) Return
  End if
  Write(6,*) "ERROR: Error defining variable in NetCDF file (",status,") : ",trim(elemdesc(1))
  call finishbanner
  Stop -1
End If

strlen=Len_trim(elemdesc(2))
status=nf_put_att_text(ncidarr(0),varid,"long_name",strlen,trim(elemdesc(2)))
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
  call finishbanner
  Stop -1
End If

strlen=Len_trim(elemdesc(3))
status=nf_put_att_text(ncidarr(0),varid,"units",strlen,trim(elemdesc(3)))
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
  call finishbanner
  Stop -1
End If

!status=nf_put_att_int2(ncidarr(0),varid,"missing_value",nf_int2,1,-32500)
!If (status /= nf_noerr) Then
!  Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
!  call finishbanner
!  Stop -1
!End If
!

if (numtype==nf_short) then
  ivals2=-32500
  status=nf_put_att_int2(ncidarr(0),varid,"valid_min",nf_int2,1,ivals2)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
    call finishbanner
    Stop -1
  End If
  ivals2=32500
  status=nf_put_att_int2(ncidarr(0),varid,"valid_max",nf_int2,1,ivals2)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
    call finishbanner
    Stop -1
  End If
  rvals=of
  status=nf_put_att_real(ncidarr(0),varid,"add_offset",nf_float,1,rvals)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
    call finishbanner
    Stop -1
  End If
  rvals=sc
  status=nf_put_att_real(ncidarr(0),varid,"scale_factor",nf_float,1,rvals)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
    call finishbanner
    Stop -1
  End If
end if

status=nf_put_att_text(ncidarr(0),varid,"FORTRAN_format",5,"G11.4")
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
  call finishbanner
  Stop -1
End If


Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine defines a dimension
!

subroutine ncadd_dimension(ncidarr,elemdesc,dimlen,dimid)

use netcdf_m

implicit None

integer, dimension(0:4), intent(in) :: ncidarr
integer, intent(in) :: dimlen
integer, intent(out) :: dimid
character(len=*), intent(in) :: elemdesc
integer status

status = nf_def_dim(ncidarr(0),elemdesc,dimlen,dimid)
if (status /= nf_noerr) Then
  write(6,*) "ERROR: Error defining dimension in NetCDF file (",status,"): ",trim(elemdesc)
  call finishbanner
  stop -1
end If

return
end
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine defines a 1D variable
!
    
subroutine ncadd_1dvar(ncidarr,elemdesc,numtype,dimid)

use netcdf_m

implicit none

integer, dimension(0:4), intent(in) :: ncidarr
integer, intent(in) :: numtype, dimid
character(len=*), dimension(1:3), intent(in) :: elemdesc
integer status, strlen, varid
integer ch_dimid
integer, dimension(2) :: dimid_a
integer, parameter :: ch_len = 50

if ( numtype == nf_char ) then
    
  status = nf_inq_dimid(ncidarr(0),'chid',ch_dimid)
  if ( status /= nf_noerr ) then
    status = nf_def_dim(ncidarr(0),'chid',ch_len,ch_dimid)
    if ( status /= nf_noerr ) then
      write(6,*) "ERROR: Error defining dimension chid in NetCDF file (",status,")"
      call finishbanner
      stop -1
    end if
  end if
  dimid_a(1) = ch_dimid
  dimid_a(2) = dimid
  status = nf_def_var(ncidarr(0),elemdesc(1),nf_char,2,dimid_a(1:2),varid)
  if (status /= nf_noerr) Then
    if (status == nf_enameinuse) then
      write(6,*) "WARN:  Variable ",trim(elemdesc(1))," already exists."
      status = nf_inq_varid(ncidarr(0),elemdesc(1),varid)
      if (status == nf_noerr) Return
    end if
    write(6,*) "ERROR: Error defining variable in NetCDF file (",status,") : ",trim(elemdesc(1))
    call finishbanner
    stop -1
  end If
  
else

  status = nf_def_var(ncidarr(0),elemdesc(1),numtype,1,dimid,varid)
  if (status /= nf_noerr) Then
    if (status == nf_enameinuse) then
      write(6,*) "WARN:  Variable ",trim(elemdesc(1))," already exists."
      status = nf_inq_varid(ncidarr(0),elemdesc(1),varid)
      if (status == nf_noerr) Return
    end if
    write(6,*) "ERROR: Error defining variable in NetCDF file (",status,") : ",trim(elemdesc(1))
    call finishbanner
    stop -1
  end If
  
end if

strlen = Len_trim(elemdesc(2))
status = nf_put_att_text(ncidarr(0),varid,"long_name",strlen,trim(elemdesc(2)))
if (status /= nf_noerr) Then
  write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
  call finishbanner
  stop -1
end If

strlen = Len_trim(elemdesc(3))
status = nf_put_att_text(ncidarr(0),varid,"units",strlen,trim(elemdesc(3)))
if (status /= nf_noerr) Then
  write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
  call finishbanner
  stop -1
end If

return
end
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine changes from define mode to data mode
!

Subroutine ncenddef(ncidarr)

use netcdf_m

Implicit None

Integer, dimension(0:4), intent(in) :: ncidarr
Integer status

status = nf_enddef(ncidarr(0))
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error ending define mode in NetCDF file (",status,")"
  call finishbanner
  Stop -1
End If


Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores lon, lat, lvl data
!

Subroutine nclonlat(ncidarr,dimid,alonlat,alvl,dimnum)

use netcdf_m

Implicit None

Integer, dimension(0:4), intent(in) :: ncidarr
Integer, dimension(1:4), intent(in) :: dimid
Integer, dimension(1:3), intent(in) :: dimnum
Integer, dimension(1:4) :: dimnumout
Real, dimension(1:3,1:2), intent(in) :: alonlat
Real, dimension(1:dimnum(3)), intent(in) :: alvl
Real, dimension(1) :: atime

atime=0.
dimnumout(1:3)=dimnum(1:3)
dimnumout(4)=1

Call nclonlatgen(ncidarr,dimid,alonlat,alvl,atime,dimnumout)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores lon, lat, lvl data (general version)
!

Subroutine nclonlatgen(ncidarr,dimid,alonlat,alvl,atime,dimnum)

use netcdf_m

Implicit None

Integer, dimension(0:4), intent(in) :: ncidarr
Integer, dimension(1:4), intent(in) :: dimid
Integer, dimension(1:4), intent(in) :: dimnum
Real, dimension(1:3,1:2), intent(in) :: alonlat
Real, dimension(1:dimnum(3)), intent(in) :: alvl
Real, dimension(1:dimnum(4)), intent(in) :: atime
Real, dimension(:), allocatable :: ldata
Real sgn
Integer i,j,status,vtype
integer, dimension(1) :: nstart

Do i=1,2
  sgn=Abs(alonlat(3,i))
  If (alonlat(1,i)>alonlat(2,i)) sgn=-sgn
  Allocate(ldata(1:dimnum(i)))
  Do j=1,dimnum(i)
    ldata(j)=alonlat(1,i)+sgn*Real(j-1)
  End Do
  
  nstart=1
  status = nf_put_vara_real(ncidarr(0),dimid(i),nstart,dimnum(i:i),ldata)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing lon and lat data (",status,")"
    call finishbanner
    Stop -1
  End If
  Deallocate(ldata)
End Do

If (dimnum(3)/=1) Then
  nstart=1
  status = nf_put_vara_real(ncidarr(0),dimid(3),nstart,dimnum(3:3),alvl)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing lvl data (",status,")"
    call finishbanner
    Stop -1
  End If
End If

status=nf_inq_vartype(ncidarr(0),dimid(4),vtype)
select case(vtype)
  case(nf_float)
    nstart=1
    status = nf_put_vara_real(ncidarr(0),dimid(4),nstart,dimnum(4:4),atime)
  case(nf_int)
    nstart=1
    status = nf_put_vara_int(ncidarr(0),dimid(4),nstart,dimnum(4:4),nint(atime))
  case DEFAULT
    write(6,*) 'ERROR: Unsupported time vartype ',vtype
    call finishbanner
    stop -1
end select
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error writing time data (",status,")"
  call finishbanner
  Stop -1
End If

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores lon, lat, lvl data (array version)
!

Subroutine nclonlatarr(ncidarr,dimid,alonlat,latarr,alvl,atime,dimnum)

use netcdf_m

Implicit None

Integer, dimension(0:4), intent(in) :: ncidarr
Integer, dimension(1:4), intent(in) :: dimid
Integer, dimension(1:4), intent(in) :: dimnum
Real, dimension(1:3), intent(in) :: alonlat
real, dimension(1:dimnum(2)), intent(in) :: latarr
Real, dimension(1:dimnum(3)), intent(in) :: alvl
Real, dimension(1:dimnum(4)), intent(in) :: atime
Real, dimension(:), allocatable :: ldata
Real sgn
real, dimension(1) :: rvals
Integer i,j,status,vtype
integer, dimension(1) :: nstart,ncount

sgn=Abs(alonlat(3))
If (alonlat(1)>alonlat(2)) sgn=-sgn
Allocate(ldata(1:dimnum(1)))
Do j=1,dimnum(1)
  ldata(j)=alonlat(1)+sgn*Real(j-1)
End Do

nstart=1
status = nf_put_vara_real(ncidarr(0),dimid(1),nstart,dimnum(1:1),ldata)
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error writing lon data (",status,")"
  call finishbanner
  Stop -1
End If
Deallocate(ldata)

nstart=1
status = nf_put_vara_real(ncidarr(0),dimid(2),nstart,dimnum(2:2),latarr)
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error writing lat data (",status,")"
  call finishbanner
  Stop -1
End If

If (dimnum(3)/=1) Then
  nstart=1
  status = nf_put_vara_real(ncidarr(0),dimid(3),nstart,dimnum(3:3),alvl)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing lvl data (",status,")"
    call finishbanner
    Stop -1
  End If
End If

If (dimnum(4)/=1) Then
  status=nf_inq_vartype(ncidarr(0),dimid(4),vtype)
  select case(vtype)
    case(nf_float)
      nstart=1
      status = nf_put_vara_real(ncidarr(0),dimid(4),nstart,dimnum(4:4),atime)
    case(nf_int)
      nstart=1
      status = nf_put_vara_int(ncidarr(0),dimid(4),nstart,dimnum(4:4),nint(atime))
    case DEFAULT
      write(6,*) 'ERROR: Unsupported time vartype ',vtype
      call finishbanner
      stop -1
  end select
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing time data (",status,")"
    call finishbanner
    Stop -1
  End If
Else
  nstart=1
  ncount=1
  rvals=0.
  status = nf_put_vara_real(ncidarr(0),dimid(4),nstart,ncount,rvals)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing time data (",status,")"
    call finishbanner
    Stop -1
  End If
End If

Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine writes data to a NetCDF file (generalised ncwritedat)
! without position
!

Subroutine ncwritedatgen(ncidarr,dataout,dimnum,varid)

Implicit None

Integer, dimension(0:4), intent(in) :: ncidarr
Integer, intent(in) :: varid
Integer, dimension(1:4), intent(in) :: dimnum
Integer, dimension(1:4) :: startpos
Real, dimension(1:dimnum(1),1:dimnum(2),1:dimnum(3),1:dimnum(4)), intent(in) :: dataout

startpos=1
call ncwritedatgen2(ncidarr,dataout,startpos,dimnum,varid)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine writes data to a NetCDF file (generalised ncwritedat)
! and position
!

Subroutine ncwritedat4(ncidarr,dataout,dimnum,varid)

Implicit None

Integer, dimension(0:4), intent(in) :: ncidarr
Integer, intent(in) :: varid
Integer, dimension(1:4,1:2), intent(in) :: dimnum
Integer, dimension(1:4) :: startpos,size
Real, dimension(1:dimnum(1,2),1:dimnum(2,2),1:dimnum(3,2),1:dimnum(4,2)), intent(in) :: dataout

startpos=dimnum(:,1)
size=dimnum(:,2)
call ncwritedatgen2(ncidarr,dataout,startpos,size,varid)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine writes data to a NetCDF file (generalised ncwritedat)
! and position
!

Subroutine ncwritedatgen2(ncidarr,dataout,startpos,dimnum,varid)

use netcdf_m

Implicit None

Integer, dimension(0:4), intent(in) :: ncidarr
Integer, intent(in) :: varid
Integer, dimension(4), intent(in) :: startpos,dimnum
Real, dimension(dimnum(1),dimnum(2),dimnum(3),dimnum(4)), intent(in) :: dataout
real, dimension(:,:,:,:), allocatable :: dum
real offset,scale
real, dimension(1) :: rvals
Integer, dimension(4) :: start,ncount
integer, dimension(:,:,:,:), allocatable :: idum
Integer status,xtype,numofdim

status = nf_inq_varndims(ncidarr(0),varid,numofdim)
status = nf_inq_vartype(ncidarr(0),varid,xtype)
status = nf_get_att_real(ncidarr(0),varid,'scale_factor',rvals(1))
if (status.ne.nf_noerr) then
  scale=1.
else
  scale=rvals(1)
end if
status = nf_get_att_real(ncidarr(0),varid,'add_offset',rvals(1))
if (status.ne.nf_noerr) then
  offset=0.
else
  offset=rvals(1)
end if

If (numofdim.GT.4) Then
  Write(6,*) "ERROR Max number of dimensions reached."
  Write(6,*) "      Please contact MJT and get him to"
  Write(6,*) "      fix this."
  call finishbanner
  Stop -1
End if

start(1:numofdim)=startpos(1:numofdim)
ncount(1:numofdim)=dimnum(1:numofdim)

Select Case(xtype)

  Case(nf_float)
    allocate( dum(dimnum(1),dimnum(2),dimnum(3),dimnum(4)) )
    dum=(dataout-offset)/scale
    status = nf_put_vara_real(ncidarr(0),varid,start(1:numofdim),ncount(1:numofdim),dum)
    deallocate( dum )

  Case(nf_short)
    allocate( idum(dimnum(1),dimnum(2),dimnum(3),dimnum(4)) )  
    idum=nint((dataout-offset)/scale)
    status = nf_put_vara_int(ncidarr(0),varid,start(1:numofdim),ncount(1:numofdim),idum)
    deallocate( idum )

  Case (nf_int)
    allocate( idum(dimnum(1),dimnum(2),dimnum(3),dimnum(4)) )    
    idum=nint((dataout-offset)/scale)
    status = nf_put_vara_int(ncidarr(0),varid,start(1:numofdim),ncount(1:numofdim),idum)
    deallocate( idum )

  Case DEFAULT
    Write(6,*) "ERROR: Internal error in ncwritedat.  Unknown vartype ",xtype
    call finishbanner
    Stop -1
  
End select

If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error writing data to NetCDF file (",status,")"
  call finishbanner
  Stop -1
End If

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine defines a 1D variable (real)
!
    
subroutine ncput_1dvar_real(ncidarr,elemdesc,nsize,ndata)

use netcdf_m

implicit none

integer, dimension(0:4), intent(in) :: ncidarr
integer, intent(in) :: nsize
real, dimension(nsize) :: ndata
character(len=*), intent(in) :: elemdesc
integer status, varid, vtype
integer, dimension(1) :: nstart, ncount

status = nf_inq_varid(ncidarr(0),elemdesc,varid)
if ( status /= nf_noerr ) then
  write(6,*) "ERROR: Cannot find variable in NetCDF file (",status,") : ",trim(elemdesc)
  call finishbanner
  stop -1
end if

status=nf_inq_vartype(ncidarr(0),varid,vtype)
select case(vtype)
  case(nf_float)
    nstart(1) = 1
    ncount(1) = nsize
    status = nf_put_vara_real(ncidarr(0),varid,nstart,ncount,ndata)
  case(nf_int)
    nstart(1) = 1
    ncount(1) = nsize
    status = nf_put_vara_int(ncidarr(0),varid,nstart,ncount,nint(ndata))
  case DEFAULT
    write(6,*) 'ERROR: Unsupported vartype for ncput_1dvar_real',vtype
    call finishbanner
    stop -1
end select
if (status /= nf_noerr) Then
  write(6,*) "ERROR: Error writing 1dvar data (",status,")"
  call finishbanner
  stop -1
end If

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine defines a 1D variable (char)
!
    
subroutine ncput_1dvar_text(ncidarr,elemdesc,nsize,ndata)

use netcdf_m

implicit none

integer, dimension(0:4), intent(in) :: ncidarr
integer, intent(in) :: nsize
character(len=*), dimension(nsize) :: ndata
character(len=*), intent(in) :: elemdesc
integer status, varid, vtype, n
integer, dimension(2) :: nstart, ncount, nstride

status = nf_inq_varid(ncidarr(0),elemdesc,varid)
if ( status /= nf_noerr ) then
  write(6,*) "ERROR: Cannot find variable in NetCDF file (",status,") : ",trim(elemdesc)
  call finishbanner
  stop -1
end if

status = nf_inq_vartype(ncidarr(0),varid,vtype)
select case(vtype)
  case(nf_char)
    do n = 1,nsize
      nstart(1) = 1
      nstart(2) = n
      ncount(1) = len_trim(ndata(n))
      ncount(2) = 1
      nstride(1) = 1
      nstride(2) = 1
      status = nf_put_vars_text(ncidarr(0),varid,nstart,ncount,nstride,ndata(n))
      if (status /= nf_noerr) Then
        write(6,*) "ERROR: Error writing 1dvar data (",status,")"
        call finishbanner
        stop -1
      end If
    end do
  case DEFAULT
    write(6,*) 'ERROR: Unsupported vartype for ncput_1dvar_char ',vtype
    call finishbanner
    stop -1
end select

return
end

    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine closes the NetCDF file
!

Subroutine ncclose(ncidarr)

use netcdf_m

Implicit None

Integer, dimension(0:4), intent(in) :: ncidarr
Integer status

status=nf_close(ncidarr(0))
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error closing NetCDF file (",status,")"
  call finishbanner
  Stop -1
End If


Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function returns the month name for a given number
!

Character*3 Function findmonth(mthnum)

Implicit None

Integer, intent(in) :: mthnum

Select Case(mthnum)

  Case(1)
    findmonth="jan"

  Case(2)
    findmonth="feb"

  Case(3)
    findmonth="mar"

  Case(4)
    findmonth="apr"

  Case(5)
    findmonth="may"

  Case(6)
    findmonth="jun"

  Case(7)
    findmonth="jul"

  Case(8)
    findmonth="aug"

  Case(9)
    findmonth="sep"

  Case(10)
    findmonth="oct"

  Case(11)
    findmonth="nov"

  Case(12)
    findmonth="dec"
  
  Case DEFAULT
    Write(6,*) "ERROR: Invaid month"
    call finishbanner
    Stop -1
    
End Select

Return
End

subroutine ncatt(ncidarr,desc,rval)

use netcdf_m

implicit none

integer, dimension(0:4), intent(in) :: ncidarr
integer ncstatus
real, intent(in) :: rval
real, dimension(1) :: rvals
character(len=*), intent(in) :: desc

rvals=rval
ncstatus=nf_put_att_real(ncidarr(0),nf_global,desc,nf_real,1,rvals)

return
end
