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

Implicit None

Include "netcdf.inc"

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

desc(1,:)=(/ "lon",  "longitude", "degrees_east",  "even" /)
desc(2,:)=(/ "lat",  "latitude",  "degrees_north", "even" /)
desc(3,:)=(/ "lev",  "level",     "layer",         "uneven" /)
desc(4,:)=(/ "time", "time",      "hours",         "even" /)

desc(3,1)=outputunit(1)
desc(3,3)=outputunit(2)
desc(2,4)=latsp

If (adate(1).NE.0) then
  Write(timedesc,'("hours since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') adate(:)
  desc(4,3)=timedesc
End if

! Create NetCDF file
status=nf_create(outfile,nf_clobber,ncidarr(0))
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error opening NetCDF file (",status,")"
  Stop
End If

! Define dimensions

Do i=1,4
  If (dims(i).NE.1) Then
    status=nf_def_dim(ncidarr(0),desc(i,1),dims(i),ncidarr(i))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      Stop
    End If
    status=nf_def_var(ncidarr(0),trim(desc(i,1)),nf_float,1,ncidarr(i),dimvar(i))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      Stop
    End If
    Do j=2,4
      strlen=Len_trim(desc(i,j))
      status=nf_put_att_text(ncidarr(0),dimvar(i),trim(vardesc(j)),strlen,trim(desc(i,j)))
      If (status /= nf_noerr) Then
        Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,j))
        Stop
      End If
    End Do
  End If
End Do

If (adate(1).NE.0) Then
  mthnam=findmonth(adate(2))
  Write(timedesc,'(I2.2,"-",A3,"-",I4.4,":",I2.2,":",I2.2,":",I2.2)') adate(3),mthnam,adate(1),adate(4),adate(5),adate(6)
  strlen=Len_trim(timedesc)
  status=nf_put_att_text(ncidarr(0),dimvar(4),"time_origin",strlen,trim(timedesc))
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining time in NetCDF file (",status,")"
    Stop
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

Implicit None

Include "netcdf.inc"

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

desc(1,:)=(/ "longitude", "degrees_east",  "even" /)
desc(2,:)=(/ "latitude",  "degrees_north", "even" /)
desc(3,:)=(/ "lev",       "sigma_level",   "uneven" /)
desc(4,:)=(/ "time",      "---",           "even" /)

If (adate(1).NE.0) then
  Write(timedesc,'("minutes since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') adate(:)
  desc(4,2)=timedesc
Else
  If (adate(2).EQ.1) desc(4,2)='months'
  If (adate(3).EQ.1) desc(4,2)='days'
  If (adate(4).EQ.1) desc(4,2)='hours'
  If (adate(5).EQ.1) desc(4,2)='minutes'
  If (adate(6).EQ.1) desc(4,2)='seconds'  
End if

! Create NetCDF file
status=nf_create(outfile,nf_clobber,ncidarr(0))
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error opening NetCDF file (",status,")"
  Stop
End If

! Define dimensions
Do i=1,4
  If (dims(i).NE.1) Then
    status=nf_def_dim(ncidarr(0),desc(i,1),dims(i),ncidarr(i))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      Stop
    End If
    if ((i.eq.4).and.(mode.eq.'int')) then
      vtype=nf_int
    else
      vtype=nf_float
    end if
    status=nf_def_var(ncidarr(0),desc(i,1),vtype,1,ncidarr(i),dimvar(i))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      Stop
    End If
    If (i.EQ.3) Then
      strlen=Len_trim(desc(i,2))
      status=nf_put_att_text(ncidarr(0),dimvar(i),"long_name",strlen,trim(desc(i,2)))
      If (status /= nf_noerr) Then
        Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
        Stop
      End If
      status=nf_put_att_text(ncidarr(0),dimvar(i),"positive",4,"down")
      If (status /= nf_noerr) Then
        Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
        Stop
      End If
    End if
    
    strlen=Len_trim(desc(i,2))
    status=nf_put_att_text(ncidarr(0),dimvar(i),"units",strlen,trim(desc(i,2)))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,2))
      Stop
    End If
    strlen=Len_trim(desc(i,3))
    status=nf_put_att_text(ncidarr(0),dimvar(i),"point_spacing",strlen,trim(desc(i,3)))
    If (status /= nf_noerr) Then
      Write(6,*) "ERROR: Error defining dim in NetCDF file (",status,"): ",trim(desc(i,3))
      Stop
    End If
    
  End If
End Do

If (adate(1).NE.0) Then
  mthnam=findmonth(adate(2))
  Write(timedesc,'(I2.2,"-",A3,"-",I4.4,":",I2.2,":",I2.2,":",I2.2)') adate(3),mthnam,adate(1),adate(4),adate(5),adate(6)
  strlen=Len_trim(timedesc)
  status=nf_put_att_text(ncidarr(0),dimvar(4),"time_origin",strlen,trim(timedesc))
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining time in NetCDF file (",status,")"
    Stop
  End If
End If


Return
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine defines fields in the NetCDF file
!

Subroutine ncaddvargen(ncidarr,elemdesc,numtype,numdim,varid,sc,of)

Implicit None

Include "netcdf.inc"

Integer, intent(in) :: numtype,numdim
Integer, dimension(0:4), intent(in) :: ncidarr
Integer, intent(out) :: varid
Character(len=*), dimension(1:3), intent(in) :: elemdesc
Integer, dimension(1:numdim) :: dimtype
Integer strlen,ierr
Integer status
Real, intent(in) :: sc,of

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
    Stop
End Select

status=nf_def_var(ncidarr(0),elemdesc(1),numtype,numdim,dimtype,varid)
If (status /= nf_noerr) Then
  If (status == nf_enameinuse) then
    Write(6,*) "WARN:  Variable ",trim(elemdesc(1))," already exists."
    status = nf_inq_varid(ncidarr(0),elemdesc(1),varid)
    If (status == nf_noerr) Return
  End if
  Write(6,*) "ERROR: Error defining variable in NetCDF file (",status,") : ",trim(elemdesc(1))
  Stop
End If

strlen=Len_trim(elemdesc(2))
status=nf_put_att_text(ncidarr(0),varid,"long_name",strlen,trim(elemdesc(2)))
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
Stop
End If

strlen=Len_trim(elemdesc(3))
status=nf_put_att_text(ncidarr(0),varid,"units",strlen,trim(elemdesc(3)))
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
Stop
End If

!status=nf_put_att_int2(ncidarr(0),varid,"missing_value",nf_int2,1,-32500)
!If (status /= nf_noerr) Then
!  Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
!Stop
!End If
!

if (numtype.eq.nf_short) then
  status=nf_put_att_int2(ncidarr(0),varid,"valid_min",nf_int2,1,-32500)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
    Stop
  End If
  status=nf_put_att_int2(ncidarr(0),varid,"valid_max",nf_int2,1,32500)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
    Stop
  End If
  status=nf_put_att_real(ncidarr(0),varid,"add_offset",nf_float,1,of)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
  Stop
  End If

  status=nf_put_att_real(ncidarr(0),varid,"scale_factor",nf_float,1,sc)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
  Stop
  End If
end if

status=nf_put_att_text(ncidarr(0),varid,"FORTRAN_format",5,"G11.4")
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
Stop
End If


Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine defines misc fields in the NetCDF file
!
!
!Subroutine ncaddmisc(ncidarr,elemdesc,numtype,dimtype,varid)
!
!Implicit None
!
!Include "netcdf.inc"
!
!Integer, intent(in) :: numtype,dimtype
!Integer, dimension(0:4), intent(in) :: ncidarr
!Integer, intent(out) :: varid
!Character(len=*), dimension(1:3), intent(in) :: elemdesc
!Integer strlen,ierr,numdim
!Integer status
!
!If ((dimtype.LT.1).OR.(dimtype.GT.4)) Then
!  numdim=0
!Else
!  numdim=1
!End If
!
!
!status=nf_def_var(ncidarr(0),elemdesc(1),numtype,numdim,ncidarr(dimtype),varid)
!If (status /= nf_noerr) Then
!  Write(6,*) "ERROR: Error defining variable in NetCDF file (",status,") : ",trim(elemdesc(1))
!  Stop
!End If
!
!strlen=Len_trim(elemdesc(2))
!If (strlen.GT.0) Then
!  status=nf_put_att_text(ncidarr(0),varid,"long_name",strlen,trim(elemdesc(2)))
!  If (status /= nf_noerr) Then
!    Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
!  Stop
!  End If
!End If
!
!
!strlen=Len_trim(elemdesc(3))
!If (strlen.GT.0) Then
!  status=nf_put_att_text(ncidarr(0),varid,"positive",strlen,trim(elemdesc(3)))
!  If (status /= nf_noerr) Then
!    Write(6,*) "ERROR: Error defining var in NetCDF file (",status,"): ",trim(elemdesc(1))
!  Stop
!  End If
!End If
!
!
!Return
!End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine changes from define mode to data mode
!

Subroutine ncenddef(ncidarr)

Implicit None

Include "netcdf.inc"

Integer, dimension(0:4), intent(in) :: ncidarr
Integer status

status = nf_enddef(ncidarr(0))
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error ending define mode in NetCDF file (",status,")"
  Stop
End If


Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores lon, lat, lvl data
!

Subroutine nclonlat(ncidarr,dimid,alonlat,alvl,dimnum)

Implicit None

Include "netcdf.inc"

Integer, dimension(0:4), intent(in) :: ncidarr
Integer, dimension(1:4), intent(in) :: dimid
Integer, dimension(1:3), intent(in) :: dimnum
Integer, dimension(1:4) :: dimnumout
Real, dimension(1:3,1:2), intent(in) :: alonlat
Real, dimension(1:dimnum(3)), intent(in) :: alvl
Real, dimension(1) :: atime

atime=0.
dimnumout=1
dimnumout(1:3)=dimnum(1:3)

Call nclonlatgen(ncidarr,dimid,alonlat,alvl,atime,dimnumout)

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores lon, lat, lvl data (general version)
!

Subroutine nclonlatgen(ncidarr,dimid,alonlat,alvl,atime,dimnum)

Implicit None

Include "netcdf.inc"

Integer, dimension(0:4), intent(in) :: ncidarr
Integer, dimension(1:4), intent(in) :: dimid
Integer, dimension(1:4), intent(in) :: dimnum
Real, dimension(1:3,1:2), intent(in) :: alonlat
Real, dimension(1:dimnum(3)), intent(in) :: alvl
Real, dimension(1:dimnum(4)), intent(in) :: atime
Real, dimension(:), allocatable :: ldata
Real sgn
Integer i,j,status,vtype

Do i=1,2
  sgn=Abs(alonlat(3,i))
  If (alonlat(1,i).GT.alonlat(2,i)) sgn=-sgn
  Allocate(ldata(1:dimnum(i)))
  Do j=1,dimnum(i)
    ldata(j)=alonlat(1,i)+sgn*Real(j-1)
  End Do
  
  status = nf_put_vara_real(ncidarr(0),dimid(i),1,dimnum(i),ldata)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing lon and lat data (",status,")"
    Stop
  End If
  Deallocate(ldata)
End Do

If (dimnum(3).NE.1) Then
  status = nf_put_vara_real(ncidarr(0),dimid(3),1,dimnum(3),alvl)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing lvl data (",status,")"
    Stop
  End If
End If

If (dimnum(4).NE.1) Then
  status=nf_inq_vartype(ncidarr(0),dimid(4),vtype)
  select case(vtype)
    case(nf_float)
      status = nf_put_vara_real(ncidarr(0),dimid(4),1,dimnum(4),atime)
    case(nf_int)
      status = nf_put_vara_int(ncidarr(0),dimid(4),1,dimnum(4),nint(atime))
    case DEFAULT
      write(6,*) 'ERROR: Unsupported time vartype ',vtype
      stop
  end select
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing time data (",status,")"
    Stop
  End If
Else
  status = nf_put_vara_real(ncidarr(0),dimid(4),1,1,0.)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing time data (",status,")"
    Stop
  End If
End If

Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores lon, lat, lvl data (array version)
!

Subroutine nclonlatarr(ncidarr,dimid,alonlat,latarr,alvl,atime,dimnum)

Implicit None

Include "netcdf.inc"

Integer, dimension(0:4), intent(in) :: ncidarr
Integer, dimension(1:4), intent(in) :: dimid
Integer, dimension(1:4), intent(in) :: dimnum
Real, dimension(1:3), intent(in) :: alonlat
real, dimension(1:dimnum(2)), intent(in) :: latarr
Real, dimension(1:dimnum(3)), intent(in) :: alvl
Real, dimension(1:dimnum(4)), intent(in) :: atime
Real, dimension(:), allocatable :: ldata
Real sgn
Integer i,j,status,vtype

sgn=Abs(alonlat(3))
If (alonlat(1).GT.alonlat(2)) sgn=-sgn
Allocate(ldata(1:dimnum(1)))
Do j=1,dimnum(1)
  ldata(j)=alonlat(1)+sgn*Real(j-1)
End Do
  
status = nf_put_vara_real(ncidarr(0),dimid(1),1,dimnum(1),ldata)
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error writing lon data (",status,")"
  Stop
End If
Deallocate(ldata)
  
status = nf_put_vara_real(ncidarr(0),dimid(2),1,dimnum(2),latarr)
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error writing lat data (",status,")"
  Stop
End If

If (dimnum(3).NE.1) Then
  status = nf_put_vara_real(ncidarr(0),dimid(3),1,dimnum(3),alvl)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing lvl data (",status,")"
    Stop
  End If
End If

If (dimnum(4).NE.1) Then
  status=nf_inq_vartype(ncidarr(0),dimid(4),vtype)
  select case(vtype)
    case(nf_float)
      status = nf_put_vara_real(ncidarr(0),dimid(4),1,dimnum(4),atime)
    case(nf_int)
      status = nf_put_vara_int(ncidarr(0),dimid(4),1,dimnum(4),nint(atime))
    case DEFAULT
      write(6,*) 'ERROR: Unsupported time vartype ',vtype
      stop
  end select
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing time data (",status,")"
    Stop
  End If
Else
  status = nf_put_vara_real(ncidarr(0),dimid(4),1,1,0.)
  If (status /= nf_noerr) Then
    Write(6,*) "ERROR: Error writing time data (",status,")"
    Stop
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

Implicit None

Include "netcdf.inc"

Integer, dimension(0:4), intent(in) :: ncidarr
Integer, intent(in) :: varid
Integer, dimension(1:4), intent(in) :: startpos,dimnum
Real, dimension(1:dimnum(1),1:dimnum(2),1:dimnum(3),1:dimnum(4)), intent(in) :: dataout
real offset,scale
Integer start(1:4),ncount(1:4)
Integer status,xtype,numofdim

status = nf_inq_varndims(ncidarr(0),varid,numofdim)
status = nf_inq_vartype(ncidarr(0),varid,xtype)
status = nf_get_att_real(ncidarr(0),varid,'scale_factor',scale)
if (status.ne.nf_noerr) scale=1.
status = nf_get_att_real(ncidarr(0),varid,'add_offset',offset)
if (status.ne.nf_noerr) offset=0.

If (numofdim.GT.4) Then
  Write(6,*) "ERROR Max number of dimensions reached."
  Write(6,*) "      Please contact MJT and get him to"
  Write(6,*) "      fix this."
  Stop
End if

start(1:numofdim)=startpos(1:numofdim)
ncount(1:numofdim)=dimnum(1:numofdim)

Select Case(xtype)

  Case(nf_float)
    status = nf_put_vara_real(ncidarr(0),varid,start(1:numofdim),ncount(1:numofdim),(dataout-offset)/scale)

  Case(nf_short)
    status = nf_put_vara_int(ncidarr(0),varid,start(1:numofdim),ncount(1:numofdim),nint((dataout-offset)/scale))

  Case (nf_int)
    status = nf_put_vara_int(ncidarr(0),varid,start(1:numofdim),ncount(1:numofdim),nint((dataout-offset)/scale))

  Case DEFAULT
    Write(6,*) "ERROR: Internal error in ncwritedat.  Unknown vartype ",xtype
    Stop
  
End select

If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error writing data to NetCDF file (",status,")"
  Stop
End If

Return
End



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine closes the NetCDF file
!

Subroutine ncclose(ncidarr)

Implicit None

Include "netcdf.inc"

Integer, dimension(0:4), intent(in) :: ncidarr
Integer status

status=nf_close(ncidarr(0))
If (status /= nf_noerr) Then
  Write(6,*) "ERROR: Error closing NetCDF file (",status,")"
  Stop
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
    Stop
    
End Select

Return
End

subroutine ncatt(ncidarr,desc,rval)

implicit none

include 'netcdf.inc'

integer, dimension(0:4), intent(in) :: ncidarr
integer ncstatus
real, intent(in) :: rval
character(len=*), intent(in) :: desc

ncstatus=nf_put_att_real(ncidarr(0),nf_global,desc,nf_real,1,rval)

return
end