FF = ifort
XFLAGS = -O3 -Qipo -QaxSSSE3,SSE3,SSE2 /arch:IA32 /Gm /static
LIBS = netcdf.lib
LDFLAGS = /link netcdf.lib /stack:9999999999

OBJT = casafield.obj casaread.obj setxyz_m.obj ccinterp.obj readswitch.obj jimcc_m.obj \
       latltoij_m.obj xyzinfo_m.obj newmpar_m.obj indices_m.obj \
       parm_m.obj precis_m.obj ind_m.obj jimco_m.obj jim_utils.obj nfft_m.obj \
       ncread.obj ncwrite.obj misc.obj

casafield :$(OBJT)
	$(FF) $(XFLAGS) $(OBJT) $(LIBS) $(LDFLAGS) -out:casafield.exe

clean:
	del *.obj core *.mod *.exe
# This section gives the rules for building object modules.

.SUFFIXES:.f90
.f90.obj:
	$(FF) -c $(XFLAGS) $<
.f.obj:
	$(FF) -c $(XFLAGS) $<

casafield.obj : ccinterp.obj
casaread.obj : ccinterp.obj
ccinterp.obj : ccinterp.f90 setxyz_m.obj xyzinfo_m.obj latltoij_m.obj newmpar_m.obj
latltoij_m.obj : latltoij_m.f90 xyzinfo_m.obj newmpar_m.obj
setxyz_m.obj : setxyz_m.f90 newmpar_m.obj indices_m.obj parm_m.obj precis_m.obj ind_m.obj xyzinfo_m.obj jimco_m.obj jimcc_m.obj 
xyzinfo_m.obj : xyzinfo_m.f90 precis_m.obj
newmpar_m.obj : newmpar_m.f90 
precis_m.obj : precis_m.f90
indices_m.obj : indices_m.f90
parm_m.obj : parm_m.f90 precis_m.obj 
ind_m.obj : ind_m.f90 newmpar_m.obj 
jimcc_m.obj : jimcc_m.f90 parm_m.obj precis_m.obj 
jimco_m.obj : jimco_m.f90 precis_m.obj jim_utils.obj nfft_m.obj 
jim_utils.obj : jim_utils.f90 precis_m.obj 
nfft_m.obj : nfft_m.f90 precis_m.obj 