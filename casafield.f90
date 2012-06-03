program casafield

! This code maps CASA field data to conformal cubic coordinates

implicit none

character*80, dimension(:,:), allocatable :: options
integer nopts

write(6,*) 'CASAFIELD - CASA fields to CC grid (APR-12)'

! Read switches
nopts=3
allocate (options(nopts,2))
options(:,1) = (/ '-t', '-i', '-o' /)
options(:,2) = ''

call readswitch(options,nopts)
call defaults(options,nopts)

call createcasa(options,nopts)

deallocate(options)

stop
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
character*80 returnoption,outfile,infile,topofile
character*45 header
character*9 formout
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
integer sibsize,tunit,i,j,k,ierr,ilout

outfile=returnoption('-o',options,nopts)
infile=returnoption('-i',options,nopts)
topofile=returnoption('-t',options,nopts)

! Read topography file
tunit=1
open(tunit,file=topofile,form='formatted',status='old',iostat=ierr)
read(tunit,*,iostat=ierr) sibdim(1),sibdim(2),lonlat(1),lonlat(2),schmidt,dsx,header

write(6,*) "Dimension : ",sibdim
write(6,*) "lon0,lat0 : ",lonlat
write(6,*) "Schmidt   : ",schmidt
allocate(gridout(sibdim(1),sibdim(2)),rlld(sibdim(1),sibdim(2),2))
allocate(topdata(sibdim(1),sibdim(2)))
allocate(lsdata(sibdim(1),sibdim(2)),cfield(sibdim(1),sibdim(2),nfield))

ilout=min(sibdim(1),30) ! To be compatiable with terread
write(formout,'("(",i3,"f7.0)")',IOSTAT=ierr) ilout
read(tunit,formout,iostat=ierr) topdata ! read topography
write(formout,'("(",i3,"f4.1)")',IOSTAT=ierr) ilout
read(tunit,formout,iostat=ierr) lsdata ! Read ls mask
close(tunit)

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
outputdesc=(/ 'sorder', 'Soil order', 'none' /)
call ncaddvargen(ncidarr,outputdesc,5,2,varid(1),1.,0.)
outputdesc=(/ 'ndep', 'N deposition rate', 'g m-2 year-1' /)
call ncaddvargen(ncidarr,outputdesc,5,2,varid(2),1.,0.)
outputdesc=(/ 'nfix', 'N fixation rate', 'g m-2 year-1' /)
call ncaddvargen(ncidarr,outputdesc,5,2,varid(3),1.,0.)
outputdesc=(/ 'pdust', 'P dust deposition rate', 'g m-2 year-1' /)
call ncaddvargen(ncidarr,outputdesc,5,2,varid(4),1.,0.)
outputdesc=(/ 'pweather', 'P weathering rate', 'g m-2 year-1' /)
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

