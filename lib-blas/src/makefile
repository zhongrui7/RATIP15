#
FC=xlf
FFLAGS= -O
#
OBJECTS=      dasum.o    dgemv.o    drotg.o    dswap.o    \
              dsyrk.o    dtrmv.o    daxpy.o    dger.o     \
              dsbmv.o    dsymm.o    dtbmv.o    dtrsm.o    \
              dcopy.o    dinit.o    dscal.o    dsymv.o    \
              dtbsv.o    dtrsv.o    ddot.o     dmach.o    \
              dspmv.o    dsyr.o     dtpmv.o    gather.o   \
              dgbmv.o    dnrm2.o    dspr.o     dsyr2.o    \
              dtpsv.o    idamax.o   dgemm.o    drot.o     \
              dspr2.o    dsyr2k.o   dtrmm.o    scatter.o  \
              iinit.o    daxpyi.o   dgthrz.o   ddoti.o    \
              dgthr.o droti.o dsctr.o

#
blas :  $(OBJECTS)
	ar r blas.a    dasum.o    dgemv.o    drotg.o    dswap.o    \
              dsyrk.o    dtrmv.o    daxpy.o    dger.o     \
              dsbmv.o    dsymm.o    dtbmv.o    dtrsm.o    \
              dcopy.o    dinit.o    dscal.o    dsymv.o    \
              dtbsv.o    dtrsv.o    ddot.o     dmach.o    \
              dspmv.o    dsyr.o     dtpmv.o    gather.o   \
              dgbmv.o    dnrm2.o    dspr.o     dsyr2.o    \
              dtpsv.o    idamax.o   dgemm.o    drot.o     \
              dspr2.o    dsyr2k.o   dtrmm.o    scatter.o  \
              iinit.o    daxpyi.o   dgthrz.o   ddoti.o    \
              dgthr.o droti.o dsctr.o
	ranlib blas.a
	mv blas.a ../blas.a
#

dasum.o :  dasum.f
	$(FC)  -c  $(FFLAGS)  dasum.f


daxpy.o :  daxpy.f
	$(FC)  -c  $(FFLAGS)  daxpy.f


daxpyi.o :  daxpyi.f
	$(FC)  -c  $(FFLAGS)  daxpyi.f


dcopy.o :  dcopy.f
	$(FC)  -c  $(FFLAGS)  dcopy.f


ddot.o :  ddot.f
	$(FC)  -c  $(FFLAGS)  ddot.f


ddoti.o :  ddoti.f
	$(FC)  -c  $(FFLAGS)  ddoti.f


dgbmv.o :  dgbmv.f
	$(FC)  -c  $(FFLAGS)  dgbmv.f


dgemm.o :  dgemm.f
	$(FC)  -c  $(FFLAGS)  dgemm.f


dgemv.o :  dgemv.f
	$(FC)  -c  $(FFLAGS)  dgemv.f


dger.o :  dger.f
	$(FC)  -c  $(FFLAGS)  dger.f


dgthr.o :  dgthr.f
	$(FC)  -c  $(FFLAGS)  dgthr.f


dgthrz.o :  dgthrz.f
	$(FC)  -c  $(FFLAGS)  dgthrz.f


dinit.o :  dinit.f
	$(FC)  -c  $(FFLAGS)  dinit.f


dmach.o :  dmach.f
	$(FC)  -c  $(FFLAGS)  dmach.f


dnrm2.o :  dnrm2.f
	$(FC)  -c  $(FFLAGS)  dnrm2.f


drot.o :  drot.f
	$(FC)  -c  $(FFLAGS)  drot.f


drotg.o :  drotg.f
	$(FC)  -c  $(FFLAGS)  drotg.f


droti.o :  droti.f
	$(FC)  -c  $(FFLAGS)  droti.f


dsbmv.o :  dsbmv.f
	$(FC)  -c  $(FFLAGS)  dsbmv.f


dscal.o :  dscal.f
	$(FC)  -c  $(FFLAGS)  dscal.f


dsctr.o :  dsctr.f
	$(FC)  -c  $(FFLAGS)  dsctr.f


dspmv.o :  dspmv.f
	$(FC)  -c  $(FFLAGS)  dspmv.f


dspr.o :  dspr.f
	$(FC)  -c  $(FFLAGS)  dspr.f


dspr2.o :  dspr2.f
	$(FC)  -c  $(FFLAGS)  dspr2.f


dswap.o :  dswap.f
	$(FC)  -c  $(FFLAGS)  dswap.f


dsymm.o :  dsymm.f
	$(FC)  -c  $(FFLAGS)  dsymm.f


dsymv.o :  dsymv.f
	$(FC)  -c  $(FFLAGS)  dsymv.f


dsyr.o :  dsyr.f
	$(FC)  -c  $(FFLAGS)  dsyr.f


dsyr2.o :  dsyr2.f
	$(FC)  -c  $(FFLAGS)  dsyr2.f


dsyr2k.o :  dsyr2k.f
	$(FC)  -c  $(FFLAGS)  dsyr2k.f


dsyrk.o :  dsyrk.f
	$(FC)  -c  $(FFLAGS)  dsyrk.f


dtbmv.o :  dtbmv.f
	$(FC)  -c  $(FFLAGS)  dtbmv.f


dtbsv.o :  dtbsv.f
	$(FC)  -c  $(FFLAGS)  dtbsv.f


dtpmv.o :  dtpmv.f
	$(FC)  -c  $(FFLAGS)  dtpmv.f


dtpsv.o :  dtpsv.f
	$(FC)  -c  $(FFLAGS)  dtpsv.f


dtrmm.o :  dtrmm.f
	$(FC)  -c  $(FFLAGS)  dtrmm.f


dtrmv.o :  dtrmv.f
	$(FC)  -c  $(FFLAGS)  dtrmv.f


dtrsm.o :  dtrsm.f
	$(FC)  -c  $(FFLAGS)  dtrsm.f


dtrsv.o :  dtrsv.f
	$(FC)  -c  $(FFLAGS)  dtrsv.f


gather.o :  gather.f
	$(FC)  -c  $(FFLAGS)  gather.f


iinit.o :  iinit.f
	$(FC)  -c  $(FFLAGS)  iinit.f


idamax.o :  idamax.f
	$(FC)  -c  $(FFLAGS)  idamax.f


scatter.o :  scatter.f
	$(FC)  -c  $(FFLAGS)  scatter.f

