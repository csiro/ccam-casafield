! This subroutine is to extract (in memory) data from the CASA dataset.
!

Subroutine getdata(dataout,grid,lsdata,rlld,sibdim,casafile,nfield)

Use ccinterp

Implicit None

include 'netcdf.inc'

integer, intent(in) :: nfield
integer, dimension(2), intent(in) :: sibdim
integer, dimension(sibdim(1),sibdim(2)) :: countt
integer, dimension(4,2) :: arrsize
integer, dimension(4) :: ncsize
integer, dimension(2) :: pxy
integer, dimension(1) :: pos
integer ncstatus,ncid
integer j,n,ii,jj,lci,lcj,nface
real, dimension(sibdim(1),sibdim(2),nfield), intent(out) :: dataout
real, dimension(sibdim(1),sibdim(2),0:12) :: datatmp
real, dimension(sibdim(1),sibdim(2)), intent(in) :: grid,lsdata
real, dimension(sibdim(1),sibdim(2),2), intent(in) :: rlld
real, dimension(:,:), allocatable :: coverout
real, dimension(2,2) :: emlonlat
real, dimension(:), allocatable :: dis
real aglon,aglat,alci,alcj,ssum
character(len=*), intent(in) :: casafile
character*160, dimension(2) :: varname
logical, dimension(sibdim(1),sibdim(2)) :: sermask

dataout=0.

write(6,*) 'Process CASA field datasets'

! read size and coordinates
ncstatus=nf_open(casafile,nf_nowrite,ncid)
if (ncstatus.NE.nf_noerr) Then
  write(6,*) "ERROR: Cannot open NetCDF file ",trim(casafile)," (",ncstatus,")"
  stop
end If
call getncdims(ncid,ncsize)
call getnclonlat(ncid,emlonlat)
arrsize=1
arrsize(1:2,2)=ncsize(1:2)


!--------------------------------------------------------------------      
! allocate arrays    
allocate(coverout(arrsize(1,2),arrsize(2,2)))

do n=1,5
  coverout=0.
  countt=0

  select case(n)
    case(1)
      write(6,*) "Processing sorder"
      varname=(/ 'sorder', 'integer 1 to 12' /)
    case(2)
      write(6,*) "Processing ndep"
      varname=(/ 'ndep', 'N g/m2 per year' /)
    case(3)
      write(6,*) "Processing nfix"
      varname=(/ 'nfix', 'N g/m2 per year' /)
    case(4)
      write(6,*) "Processing pdust"
      varname=(/ 'pdust', 'P g/m2 per year' /)
    case(5)
      write(6,*) "Processing pwea"
      varname=(/ 'pwea', 'P g/m2 per year' /)
  end select
  call getmeta(ncid,varname,coverout,arrsize)

  select case(n)
    case(1)
      datatmp=0.
      ! bin land points
      do jj=1,arrsize(2,2)
        aglat=(emlonlat(2,2)-emlonlat(2,1))*real(jj-1)/real(arrsize(2,2)-1)+emlonlat(2,1)
        do ii=1,arrsize(1,2)
          aglon=(emlonlat(1,2)-emlonlat(1,1))*real(ii-1)/real(arrsize(1,2)-1)+emlonlat(1,1)
          ! find cc grid point
          call lltoijmod(aglon,aglat,alci,alcj,nface)
          lci = nint(alci)
          lcj = nint(alcj)
          lcj = lcj+nface*sibdim(1)
          ! bin soil order
          j=nint(coverout(ii,jj))
          datatmp(lci,lcj,j)=datatmp(lci,lcj,j)+1.
          countt(lci,lcj)=countt(lci,lcj)+1
        end do
      end do
      
      ! fill land values
      do lcj=1,sibdim(2)
        do lci=1,sibdim(1)
          if (countt(lci,lcj).eq.0) then
            aglon=rlld(lci,lcj,1)
            aglat=rlld(lci,lcj,2)
            if (aglon.lt.emlonlat(1,1)) aglon=aglon+360.
            if (aglon.gt.emlonlat(1,1)+360.) aglon=aglon-360.
            ii=nint((aglon-emlonlat(1,1))*real(arrsize(1,2)-1)/(emlonlat(1,2)-emlonlat(1,1)))+1
            if (ii>arrsize(1,2)) ii=ii-arrsize(1,2)
            jj=nint((aglat-emlonlat(2,1))*real(arrsize(2,2)-1)/(emlonlat(2,2)-emlonlat(2,1)))+1
            jj=min(max(jj,1),arrsize(2,2))
            j=nint(coverout(ii,jj))
            datatmp(lci,lcj,j)=1.
            countt(lci,lcj)=1
          end if      
        end do
      end do
      
      ! assign soil order
      do lcj=1,sibdim(2)
        do lci=1,sibdim(1)
          if (nint(lsdata(lci,lcj)).eq.1) then
            ssum=sum(datatmp(lci,lcj,1:12))
            if (ssum.gt.0.001) then
              pos=maxloc(datatmp(lci,lcj,1:12))
              dataout(lci,lcj,n)=pos(1)
              countt(lci,lcj)=1
            else
              countt(lci,lcj)=0
            end if
          else
            dataout(lci,lcj,n)=0.
            countt(lci,lcj)=1
          end if
        end do
      end do 

      ! replace land values
      sermask=countt.gt.0.and.nint(lsdata).eq.1
      if (any(sermask)) then
        do lci=1,sibdim(1)
          do lcj=1,sibdim(2)
            if (countt(lci,lcj).eq.0) then
              call findnear(pxy,lci,lcj,sermask,rlld,sibdim)
              dataout(lci,lcj,n)=dataout(pxy(1),pxy(2),n)
     	      countt(lci,lcj)=countt(pxy(1),pxy(2))
            end if
          end do
        end do
      else
        write(6,*) 'WARN: Cannot find any non-trivial points'
        write(6,*) '      Assume data is trivial'
        dataout(:,:,n)=0.
      end if	
      
    case(5)
      ! water points
      where(lsdata.eq.0)
        dataout(:,:,n)=0.
        countt(:,:)=1
      end where
      
      ! bin land points
      do jj=1,arrsize(2,2)
        aglat=(emlonlat(2,2)-emlonlat(2,1))*real(jj-1)/real(arrsize(2,2)-1)+emlonlat(2,1)
        do ii=1,arrsize(1,2)
          aglon=(emlonlat(1,2)-emlonlat(1,1))*real(ii-1)/real(arrsize(1,2)-1)+emlonlat(1,1)
          ! find cc grid point
          call lltoijmod(aglon,aglat,alci,alcj,nface)
          lci = nint(alci)
          lcj = nint(alcj)
          lcj = lcj+nface*sibdim(1)
          if (nint(lsdata(lci,lcj)).eq.1) then
            if (coverout(ii,jj).gt.0.001) then
              dataout(lci,lcj,n)=dataout(lci,lcj,n)+coverout(ii,jj)
              countt(lci,lcj)=countt(lci,lcj)+1
            end if
          end if
        end do
      end do

      ! fill land values
      do lcj=1,sibdim(2)
        do lci=1,sibdim(1)
          if (countt(lci,lcj).eq.0) then
            aglon=rlld(lci,lcj,1)
            aglat=rlld(lci,lcj,2)
            if (aglon.lt.emlonlat(1,1)) aglon=aglon+360.
            if (aglon.gt.emlonlat(1,1)+360.) aglon=aglon-360.
            ii=nint((aglon-emlonlat(1,1))*real(arrsize(1,2)-1)/(emlonlat(1,2)-emlonlat(1,1)))+1
            if (ii>arrsize(1,2)) ii=ii-arrsize(1,2)
            jj=nint((aglat-emlonlat(2,1))*real(arrsize(2,2)-1)/(emlonlat(2,2)-emlonlat(2,1)))+1
            jj=min(max(jj,1),arrsize(2,2))
            if (coverout(ii,jj).gt.0.001) then
              dataout(lci,lcj,n)=coverout(ii,jj)
              countt(lci,lcj)=1
            end if
          end if      
        end do
      end do
      
      ! replace land values
      sermask=countt.gt.0
      if (any(sermask)) then
        do lci=1,sibdim(1)
          do lcj=1,sibdim(2)
            if (countt(lci,lcj).eq.0) then
              call findnear(pxy,lci,lcj,sermask,rlld,sibdim)
              dataout(lci,lcj,n)=dataout(pxy(1),pxy(2),n)
     	      countt(lci,lcj)=countt(pxy(1),pxy(2))
            end if
          end do
        end do
      else
        write(6,*) 'WARN: Cannot find any non-trivial points'
        write(6,*) '      Assume data is trivial'
        dataout(:,:,n)=0.
        countt=1
      end if	

      dataout(:,:,n)=dataout(:,:,n)/real(countt)

    case default
      ! bin land points
      do jj=1,arrsize(2,2)
        aglat=(emlonlat(2,2)-emlonlat(2,1))*real(jj-1)/real(arrsize(2,2)-1)+emlonlat(2,1)
        do ii=1,arrsize(1,2)
          aglon=(emlonlat(1,2)-emlonlat(1,1))*real(ii-1)/real(arrsize(1,2)-1)+emlonlat(1,1)
          ! find cc grid point
          call lltoijmod(aglon,aglat,alci,alcj,nface)
          lci = nint(alci)
          lcj = nint(alcj)
          lcj = lcj+nface*sibdim(1)
          ! bin emission
          if (nint(lsdata(lci,lcj)).eq.1) then
            dataout(lci,lcj,n)=dataout(lci,lcj,n)+coverout(ii,jj)
          end if
          countt(lci,lcj)=countt(lci,lcj)+1
        end do
      end do

      ! fill missing values
      do lcj=1,sibdim(2)
        do lci=1,sibdim(1)
          if (countt(lci,lcj).eq.0) then
            if (nint(lsdata(lci,lcj)).eq.1) then
              aglon=rlld(lci,lcj,1)
              aglat=rlld(lci,lcj,2)
              if (aglon.lt.emlonlat(1,1)) aglon=aglon+360.
              if (aglon.gt.emlonlat(1,1)+360.) aglon=aglon-360.
              ii=nint((aglon-emlonlat(1,1))*real(arrsize(1,2)-1)/(emlonlat(1,2)-emlonlat(1,1)))+1
              if (ii>arrsize(1,2)) ii=ii-arrsize(1,2)
              jj=nint((aglat-emlonlat(2,1))*real(arrsize(2,2)-1)/(emlonlat(2,2)-emlonlat(2,1)))+1
              jj=min(max(jj,1),arrsize(2,2))
              dataout(lci,lcj,n)=coverout(ii,jj)
            end if      
            countt(lci,lcj)=1
          end if
        end do
      end do

      dataout(:,:,n)=dataout(:,:,n)/real(countt)

  end select
  
end do

ncstatus=nf_close(ncid)

deallocate(coverout)

write(6,*) "Task complete"

return
end
