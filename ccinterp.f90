!
! THIS MODULE INTERFACES WITH SETXYZ.F (to isolate the common blocks)
!

Module ccinterp

Use parm_m, only : rlong0, rlat0, schmidt
use precis_m, only : rx

Private

Real(kind=rx), parameter :: schm13 = 0.1

Public ccgetgrid, lltoijmod, getcc, cgg2

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine determines the lon, lat and grid spacing on a CC
! grid (using setxyz_m.f90 and xyzinfo_m.f90).
!
Subroutine ccgetgrid(rlld,gridout,ecodim,lonlat,schmidtin,dst)

Use newmpar_m
Use setxyz_m
Use xyzinfo_m, only : ds,rlong,rlat,em

Implicit None

Integer, dimension(1:2), intent(in) :: ecodim
Real, dimension(1:ecodim(1),1:ecodim(2)), intent(out) :: gridout
Real, dimension(1:ecodim(1),1:ecodim(2),1:2), intent(out) :: rlld
Real, dimension(1:2), intent(in) :: lonlat
Real, intent(in) :: schmidtin
Real, intent(out) :: dst
Real, parameter :: pi = 3.1415926536
Real(kind=rx), parameter :: erad = 6.37122e6
Integer i,j,n
Integer, parameter :: diag = 0
Integer, parameter :: id = 13
Integer, parameter :: jd = 63
Integer, parameter :: ntang = 2

schmidt=schmidtin      
rlong0=lonlat(1)
rlat0=lonlat(2)
il=ecodim(1)
jl=ecodim(2)
ifull=il*jl
npanels=jl/il-1
iquad=1+il*((8*npanels)/(npanels+4))

Write(6,*) 'Start setxyz'
Call setxyz(il,jl,kl,npanels,ifull,iquad,diag,id,jd,rlong0,rlat0,schmidt,schm13,ntang,erad)
Write(6,*) 'End setxyz'
      
Do j=1,ecodim(2)
  Do i=1,ecodim(1)
    n=i+(j-1)*ecodim(1)
    gridout(i,j)=(ds/em(n))/1000. ! km
    rlld(i,j,1)=rlong(n)*180./pi
    If (rlld(i,j,1).GT.180.) rlld(i,j,1)=rlld(i,j,1)-360.
    rlld(i,j,2)=rlat(n)*180./pi
  End Do
End Do

dst=ds

Return
End subroutine ccgetgrid

Subroutine getcc(rlld,gridout,xyz,axyz,bxyz,ecodim,lonlat,schmidtin,dst)

Use newmpar_m
Use setxyz_m
Use xyzinfo_m, only : ds,rlong,rlat,em,x,y,z,ax,ay,az,bx,by,bz

Implicit None

Integer, dimension(1:2), intent(in) :: ecodim
Real, dimension(1:ecodim(1),1:ecodim(2)), intent(out) :: gridout
Real, dimension(1:ecodim(1),1:ecodim(2),1:2), intent(out) :: rlld
real, dimension(ecodim(1),ecodim(2),3), intent(out) :: xyz,axyz,bxyz
Real, dimension(1:2), intent(in) :: lonlat
Real, intent(in) :: schmidtin
Real, intent(out) :: dst
Real, parameter :: pi = 3.1415926536
Real(kind=rx), parameter :: erad = 6.37122e6
Integer i,j,n
Integer, parameter :: diag = 0
Integer, parameter :: id = 13
Integer, parameter :: jd = 63
Integer, parameter :: ntang = 2

schmidt=schmidtin      
rlong0=lonlat(1)
rlat0=lonlat(2)
il=ecodim(1)
jl=ecodim(2)
ifull=il*jl
npanels=jl/il-1
iquad=1+il*((8*npanels)/(npanels+4))

Write(6,*) 'Start setxyz'
Call setxyz(il,jl,kl,npanels,ifull,iquad,diag,id,jd,rlong0,rlat0,schmidt,schm13,ntang,erad)
Write(6,*) 'End setxyz'
      
Do j=1,ecodim(2)
  Do i=1,ecodim(1)
    n=i+(j-1)*ecodim(1)
    gridout(i,j)=(ds/em(n))/1000. ! km
    rlld(i,j,1)=rlong(n)*180./pi
    If (rlld(i,j,1).GT.180.) rlld(i,j,1)=rlld(i,j,1)-360.
    rlld(i,j,2)=rlat(n)*180./pi
    xyz(i,j,1)=x(n)
    xyz(i,j,2)=y(n)
    xyz(i,j,3)=z(n)
    axyz(i,j,1)=ax(n)
    axyz(i,j,2)=ay(n)
    axyz(i,j,3)=az(n)
    bxyz(i,j,1)=bx(n)
    bxyz(i,j,2)=by(n)
    bxyz(i,j,3)=bz(n)
  End Do
End Do

dst=ds

Return
End subroutine getcc

Subroutine cgg2(rlld,gridout,ecodim,lonlat,schmidtin,dst,in,ie,is,iw)

Use newmpar_m
Use setxyz_m
use indices_m
Use xyzinfo_m, only : ds,rlong,rlat,em

Implicit None

Integer, dimension(1:2), intent(in) :: ecodim
Real, dimension(1:ecodim(1),1:ecodim(2)), intent(out) :: gridout
Real, dimension(1:ecodim(1),1:ecodim(2),1:2), intent(out) :: rlld
Real, dimension(1:2), intent(in) :: lonlat
Real, intent(in) :: schmidtin
Real, intent(out) :: dst
Real, parameter :: pi = 3.1415926536
Real(kind=rx), parameter :: erad = 6.37122e6
integer, dimension(1:ecodim(1),1:ecodim(2)), intent(out) :: in,ie,is,iw
Integer i,j,n
Integer, parameter :: diag = 0
Integer, parameter :: id = 13
Integer, parameter :: jd = 63
Integer, parameter :: ntang = 2

schmidt=schmidtin      
rlong0=lonlat(1)
rlat0=lonlat(2)
il=ecodim(1)
jl=ecodim(2)
ifull=il*jl
npanels=jl/il-1
iquad=1+il*((8*npanels)/(npanels+4))

Write(6,*) 'Start setxyz'
Call setxyz(il,jl,kl,npanels,ifull,iquad,diag,id,jd,rlong0,rlat0,schmidt,schm13,ntang,erad)
Write(6,*) 'End setxyz'
      
Do j=1,ecodim(2)
  Do i=1,ecodim(1)
    n=i+(j-1)*ecodim(1)
    gridout(i,j)=(ds/em(n))/1000. ! km
    rlld(i,j,1)=rlong(n)*180./pi
    If (rlld(i,j,1).GT.180.) rlld(i,j,1)=rlld(i,j,1)-360.
    rlld(i,j,2)=rlat(n)*180./pi
    in(i,j)=i_n(n)
    ie(i,j)=i_e(n)
    is(i,j)=i_s(n)
    iw(i,j)=i_w(n)
  End Do
End Do

dst=ds

Return
End subroutine cgg2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calls latltoij through the module (i.e., so the
! global variables are avaliable).
!

Subroutine lltoijmod(aglon,aglat,alci,alcj,nface)

use latltoij_m

Implicit None

Real, intent(in) :: aglon,aglat
Real, intent(out) :: alci,alcj
Integer, intent(out) :: nface

Call latltoij(aglon,aglat,alci,alcj,nface,rlong0,rlat0,schmidt,schm13)

Return
End subroutine lltoijmod

End module ccinterp

