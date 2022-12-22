! Conformal Cubic Atmospheric Model
    
! Copyright 2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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


program casafield

! This code maps CASA field data to conformal cubic coordinates

implicit none

include 'version.h'

character*1024, dimension(:,:), allocatable :: options
integer nopts

! Start banner
write(6,*) "=============================================================================="
write(6,*) "CCAM: Starting casafield"
write(6,*) "=============================================================================="


#ifndef stacklimit
! For linux only - removes stacklimit on all processors
call setstacklimit(-1)
#endif 


write(6,*) 'CASAFIELD - CASA fields to CC grid'
write(6,*) version

! Read switches
nopts=3
allocate (options(nopts,2))
options(:,1) = (/ '-t', '-i', '-o' /)
options(:,2) = ''

call readswitch(options,nopts)
call defaults(options,nopts)

call createcasa(options,nopts)

deallocate(options)

! Complete
write(6,*) "CCAM: casafield completed successfully"
call finishbanner

stop
end

subroutine finishbanner

implicit none

! End banner
write(6,*) "=============================================================================="
write(6,*) "CCAM: Finished casafield"
write(6,*) "=============================================================================="

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine displays the help message
!

subroutine help()

implicit None

write(6,*)
write(6,*) "Usage:"
write(6,*) "  casafield -t topoin -i casain -o casaout"
write(6,*)
write(6,*) "Options:"
write(6,*) "  -t topoin    CCAM topography file"
write(6,*) "  -i casain    CASA field input dataset"
write(6,*) "  -o casaout   CASA field output filename"
write(6,*)

call finishbanner
stop

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determins the default values for the switches
!

subroutine defaults(options,nopts)

implicit none

integer nopts
character(len=*), dimension(nopts,2), intent(inout) :: options
integer out,in,topo
integer locate

topo=locate('-t',options(:,1),nopts)
if (options(topo,2).EQ.'') then
  write(6,*) "ERROR: Must specify topography filename"
  stop
end if

in=locate('-i',options(:,1),nopts)
if (options(in,2).EQ.'') then
  write(6,*) "ERROR: Must specify input filename"
  stop
end if

out=locate('-o',options(:,1),nopts)
if (options(out,2).EQ.'') then
  write(6,*) "ERROR: Must specify output filename"
  stop
end if

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes the CASA data
!

subroutine createcasa(options,nopts)

use ccinterp

implicit none

integer, parameter :: nfield = 5
integer, intent(in) :: nopts
character(len=*), dimension(nopts,2), intent(in) :: options
character*80, dimension(3) :: outputdesc
character*1024 returnoption,outfile,infile,topofile
character*47 header
real, dimension(:,:,:), allocatable :: rlld,cfield
real, dimension(:,:), allocatable :: gridout,lsdata,topdata
real, dimension(2) :: lonlat
real, dimension(3,2) :: alonlat
real, dimension(1) :: alvl,atime
real schmidt,dsx,ds
integer, dimension(2) :: sibdim
integer, dimension(4) :: dimnum,dimid,dimcount
integer, dimension(0:4) :: ncidarr
integer, dimension(6) :: adate
integer, dimension(nfield) :: varid
integer sibsize,tunit,i,j,k,ierr

outfile=returnoption('-o',options,nopts)
infile=returnoption('-i',options,nopts)
topofile=returnoption('-t',options,nopts)

! Read topography file
tunit=1
call readtopography(tunit,topofile,sibdim,lonlat,schmidt,dsx,header)

write(6,*) "Dimension : ",sibdim
write(6,*) "lon0,lat0 : ",lonlat
write(6,*) "Schmidt   : ",schmidt
allocate(gridout(sibdim(1),sibdim(2)),rlld(sibdim(1),sibdim(2),2))
allocate(topdata(sibdim(1),sibdim(2)))
allocate(lsdata(sibdim(1),sibdim(2)),cfield(sibdim(1),sibdim(2),nfield))

call gettopols(tunit,topofile,lsdata,sibdim)
lsdata=1.-lsdata

! Determine lat/lon to CC mapping
call ccgetgrid(rlld,gridout,sibdim,lonlat,schmidt,ds)

! Read CASA field data
call getdata(cfield,gridout,lsdata,rlld,sibdim,infile,nfield)

! Prep nc output
dimnum(1:2)=sibdim(1:2) ! CC grid dimensions
dimnum(3)=1 ! Turn off level
dimnum(4)=1 ! single month
adate=0 ! Turn off date
adate(2)=1 ! time units=months
call ncinitcc(ncidarr,outfile,dimnum(1:3),dimid,adate)
call ncatt(ncidarr,'lat0',lonlat(2))
call ncatt(ncidarr,'lon0',lonlat(1))
call ncatt(ncidarr,'schmidt0',schmidt)
outputdesc(1)='sorder'
outputdesc(2)='Soil order'
outputdesc(3)='none'
call ncaddvargen(ncidarr,outputdesc,5,2,varid(1),1.,0.)
outputdesc(1)='ndep'
outputdesc(2)='N deposition rate'
outputdesc(3)='g m-2 year-1'
call ncaddvargen(ncidarr,outputdesc,5,2,varid(2),1.,0.)
outputdesc(1)='nfix'
outputdesc(2)='N fixation rate'
outputdesc(3)='g m-2 year-1'
call ncaddvargen(ncidarr,outputdesc,5,2,varid(3),1.,0.)
outputdesc(1)='pdust'
outputdesc(2)='P dust deposition rate'
outputdesc(3)='g m-2 year-1'
call ncaddvargen(ncidarr,outputdesc,5,2,varid(4),1.,0.)
outputdesc(1)='pweather'
outputdesc(2)='P weathering rate'
outputdesc(3)='g m-2 year-1'
call ncaddvargen(ncidarr,outputdesc,5,2,varid(5),1.,0.)
call ncenddef(ncidarr)
alonlat(:,1)=(/ 1., real(sibdim(1)), 1. /)
alonlat(:,2)=(/ 1., real(sibdim(2)), 1. /)
alvl=1.
atime(1)=0
call nclonlatgen(ncidarr,dimid,alonlat,alvl,atime,dimnum)

write(6,*) 'Write CASA field data'
dimcount=(/ sibdim(1), sibdim(2), 1, 1 /)
do i=1,nfield
  call ncwritedatgen(ncidarr,cfield(:,:,i),dimcount,varid(i))
end do
call ncclose(ncidarr)

deallocate(gridout,rlld,topdata,lsdata,cfield)

return
end

