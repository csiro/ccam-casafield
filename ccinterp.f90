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

use indices_m
Use newmpar_m
Use setxyz_m
Use xyzinfo_m, only : ds,rlong,rlat,em,f,fu,fv,dmdx,dmdy,dmdxv,dmdyu,wts,x,y,z

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
!npanels=jl/il-1
iquad=1+il*((8*npanels)/(npanels+4))

Write(6,*) 'Start setxyz'
Call setxyz(il,jl,kl,npanels,ifull,iquad,diag,id,jd,rlong0,rlat0,schmidt,schm13,ntang,erad)
Write(6,*) 'End setxyz'

deallocate(i_nn,i_ss,i_ww,i_ee,i_ne)
deallocate(i_se,i_en,i_wn,i_wu,i_sv,i_wu2,i_sv2,i_eu2)
deallocate(i_nv2,i_ev2,i_nu2,i_eu,i_nv)
deallocate(lwws,lws,lwss,les,lees,less,lwwn,lwnn,leen)
deallocate(lenn,lsww,lsw,lssw,lsee,lsse,lnww,lnw,lnnw)
deallocate(lnee,lnne)
deallocate(f,fu,fv,dmdx,dmdy,dmdxv,dmdyu)
deallocate(wts)

Do j=1,ecodim(2)
  Do i=1,ecodim(1)
    n=i+(j-1)*ecodim(1)
    gridout(i,j)=(ds/em(n))/1000. ! km
    rlld(i,j,1)=rlong(n)*180./pi
    If (rlld(i,j,1)>180.) rlld(i,j,1)=rlld(i,j,1)-360.
    rlld(i,j,2)=rlat(n)*180./pi
  End Do
End Do

dst=ds

deallocate(em,rlat,rlong,x,y,z)
deallocate(i_n,i_s,i_w,i_e)

Return
End subroutine ccgetgrid

Subroutine getcc(rlld,gridout,xyz,axyz,bxyz,ecodim,lonlat,schmidtin,dst)

use indices_m
Use newmpar_m
Use setxyz_m
Use xyzinfo_m, only : ds,rlong,rlat,em,x,y,z,ax,ay,az,bx,by,bz,f,fu,fv,dmdx,dmdy,dmdxv,dmdyu,wts

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
!npanels=jl/il-1
iquad=1+il*((8*npanels)/(npanels+4))

Write(6,*) 'Start setxyz'
Call setxyz(il,jl,kl,npanels,ifull,iquad,diag,id,jd,rlong0,rlat0,schmidt,schm13,ntang,erad)
Write(6,*) 'End setxyz'

deallocate(i_nn,i_ss,i_ww,i_ee,i_ne)
deallocate(i_se,i_en,i_wn,i_wu,i_sv,i_wu2,i_sv2,i_eu2)
deallocate(i_nv2,i_ev2,i_nu2,i_eu,i_nv)
deallocate(lwws,lws,lwss,les,lees,less,lwwn,lwnn,leen)
deallocate(lenn,lsww,lsw,lssw,lsee,lsse,lnww,lnw,lnnw)
deallocate(lnee,lnne)
deallocate(f,fu,fv,dmdx,dmdy,dmdxv,dmdyu)
deallocate(wts)

Do j=1,ecodim(2)
  Do i=1,ecodim(1)
    n=i+(j-1)*ecodim(1)
    gridout(i,j)=(ds/em(n))/1000. ! km
    rlld(i,j,1)=rlong(n)*180./pi
    If (rlld(i,j,1)>180.) rlld(i,j,1)=rlld(i,j,1)-360.
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

deallocate(em,rlat,rlong,x,y,z)
deallocate(i_n,i_s,i_w,i_e)

Return
End subroutine getcc

Subroutine cgg2(rlld,gridout,ecodim,lonlat,schmidtin,dst,in,ie,is,iw)

Use newmpar_m
Use setxyz_m
use indices_m
Use xyzinfo_m, only : ds,rlong,rlat,em,f,fu,fv,dmdx,dmdy,dmdxv,dmdyu,wts,x,y,z

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
!npanels=jl/il-1
iquad=1+il*((8*npanels)/(npanels+4))

Write(6,*) 'Start setxyz'
Call setxyz(il,jl,kl,npanels,ifull,iquad,diag,id,jd,rlong0,rlat0,schmidt,schm13,ntang,erad)
Write(6,*) 'End setxyz'

deallocate(i_nn,i_ss,i_ww,i_ee,i_ne)
deallocate(i_se,i_en,i_wn,i_wu,i_sv,i_wu2,i_sv2,i_eu2)
deallocate(i_nv2,i_ev2,i_nu2,i_eu,i_nv)
deallocate(lwws,lws,lwss,les,lees,less,lwwn,lwnn,leen)
deallocate(lenn,lsww,lsw,lssw,lsee,lsse,lnww,lnw,lnnw)
deallocate(lnee,lnne)
deallocate(f,fu,fv,dmdx,dmdy,dmdxv,dmdyu)
deallocate(wts)
      
Do j=1,ecodim(2)
  Do i=1,ecodim(1)
    n=i+(j-1)*ecodim(1)
    gridout(i,j)=(ds/em(n))/1000. ! km
    rlld(i,j,1)=rlong(n)*180./pi
    If (rlld(i,j,1)>180.) rlld(i,j,1)=rlld(i,j,1)-360.
    rlld(i,j,2)=rlat(n)*180./pi
    in(i,j)=i_n(n)
    ie(i,j)=i_e(n)
    is(i,j)=i_s(n)
    iw(i,j)=i_w(n)
  End Do
End Do

dst=ds

deallocate(em,rlat,rlong,x,y,z)
deallocate(i_n,i_s,i_w,i_e)

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

