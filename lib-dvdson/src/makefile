#
FC=xlf
FFLAGS= -O
#
OBJECTS=      addabs.o   dvdson.o   newvec.o   ovflow.o   \
              setup.o    dvdrvr.o   multbc.o   orthnrm.o  \
              tstsel.o

#
dvdson :  $(OBJECTS)
	ar ur dvdson.a addabs.o   dvdson.o   newvec.o   ovflow.o   \
              setup.o    dvdrvr.o   multbc.o   orthnrm.o  \
              tstsel.o
	ranlib dvdson.a
	mv dvdson.a ../dvdson.a
#

addabs.o :  addabs.f
	$(FC)  -c  $(FFLAGS)  addabs.f


dvdrvr.o :  dvdrvr.f
	$(FC)  -c  $(FFLAGS)  dvdrvr.f


dvdson.o :  dvdson.f
	$(FC)  -c  $(FFLAGS)  dvdson.f


multbc.o :  multbc.f
	$(FC)  -c  $(FFLAGS)  multbc.f


newvec.o :  newvec.f
	$(FC)  -c  $(FFLAGS)  newvec.f


orthnrm.o :  orthnrm.f
	$(FC)  -c  $(FFLAGS)  orthnrm.f


ovflow.o :  ovflow.f
	$(FC)  -c  $(FFLAGS)  ovflow.f


setup.o :  setup.f
	$(FC)  -c  $(FFLAGS)  setup.f


tstsel.o :  tstsel.f
	$(FC)  -c  $(FFLAGS)  tstsel.f

