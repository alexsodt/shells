
CC = g++
F77 = gfortran 
CC2 = gcc
INCLUDE = -I /opt/local/include

# uncomment these to use your system lapack/blas
#LAPACK = -llapack -lblas
#LAPACK_BUILD = 

# comment these out to use your system lapack/blas
LAPACK = -llapack -lblas
LAPACK_BUILD = $(LAPACK)

GSLLIB = -L/opt/local/lib #-lgsl
#-L /usr/lib -L /usr/lib64 -lgsl 
LDFLAGS = -L $(LD_RUN_PATH)  $(LAPACK) -lgfortran $(GSLLIB)
#LDFLAGS =  -framework Accelerate

QHULL = -I qhull-2011.1/src/libqhull

FLAGS = -g 

.C.o: 
	$(CC) $(INCLUDE) $(GSLINC) -o $*.o -c $(FLAGS) $*.C

.c.o:
	$(CC2)	$(INCLUDE) $(QHULL) $(COPT) $(CDEFN)   -g	-c	$*.c

.f.o:
	$(F77) -O3 -c $*.f -o $*.o

EXECS = hopBorders processHop rehop histo lifetime outputIndex avhist

all: $(EXECS)

clean: 
	rm -f lapack.a blas.a
	rm  *.o $(EXECS)
	rm -f qhull-2011.1/lib/libqhullstatic.a
	cd qhull-2011.1 ; make clean
	cd qhull-2011.1/build ; make clean


BLASOBJS = \
BLAS/caxpy.o  \
BLAS/ccopy.o  \
BLAS/cdotc.o  \
BLAS/cdotu.o  \
BLAS/cgbmv.o  \
BLAS/cgemm.o  \
BLAS/cgemv.o  \
BLAS/cgerc.o  \
BLAS/cgeru.o  \
BLAS/chbmv.o  \
BLAS/chemm.o  \
BLAS/chemv.o  \
BLAS/cher.o  \
BLAS/cher2.o  \
BLAS/cher2k.o  \
BLAS/cherk.o  \
BLAS/chpmv.o  \
BLAS/chpr.o  \
BLAS/chpr2.o  \
BLAS/crotg.o  \
BLAS/cscal.o  \
BLAS/csrot.o  \
BLAS/csscal.o  \
BLAS/cswap.o  \
BLAS/csymm.o  \
BLAS/csyr2k.o  \
BLAS/csyrk.o  \
BLAS/ctbmv.o  \
BLAS/ctbsv.o  \
BLAS/ctpmv.o  \
BLAS/ctpsv.o  \
BLAS/ctrmm.o  \
BLAS/ctrmv.o  \
BLAS/ctrsm.o  \
BLAS/ctrsv.o  \
BLAS/dasum.o  \
BLAS/daxpy.o  \
BLAS/dcabs1.o  \
BLAS/dcopy.o  \
BLAS/ddot.o  \
BLAS/dgbmv.o  \
BLAS/dgemm.o  \
BLAS/dgemv.o  \
BLAS/dger.o  \
BLAS/dnrm2.o  \
BLAS/drot.o  \
BLAS/drotg.o  \
BLAS/drotm.o  \
BLAS/drotmg.o  \
BLAS/dsbmv.o  \
BLAS/dscal.o  \
BLAS/dsdot.o  \
BLAS/dspmv.o  \
BLAS/dspr.o  \
BLAS/dspr2.o  \
BLAS/dswap.o  \
BLAS/dsymm.o  \
BLAS/dsymv.o  \
BLAS/dsyr.o  \
BLAS/dsyr2.o  \
BLAS/dsyr2k.o  \
BLAS/dsyrk.o  \
BLAS/dtbmv.o  \
BLAS/dtbsv.o  \
BLAS/dtpmv.o  \
BLAS/dtpsv.o  \
BLAS/dtrmm.o  \
BLAS/dtrmv.o  \
BLAS/dtrsm.o  \
BLAS/dtrsv.o  \
BLAS/dzasum.o  \
BLAS/dznrm2.o  \
BLAS/icamax.o  \
BLAS/idamax.o  \
BLAS/isamax.o  \
BLAS/izamax.o  \
BLAS/lsame.o  \
BLAS/sasum.o  \
BLAS/saxpy.o  \
BLAS/scabs1.o  \
BLAS/scasum.o  \
BLAS/scnrm2.o  \
BLAS/scopy.o  \
BLAS/sdot.o  \
BLAS/sdsdot.o  \
BLAS/sgbmv.o  \
BLAS/sgemm.o  \
BLAS/sgemv.o  \
BLAS/sger.o  \
BLAS/snrm2.o  \
BLAS/srot.o  \
BLAS/srotg.o  \
BLAS/srotm.o  \
BLAS/srotmg.o  \
BLAS/ssbmv.o  \
BLAS/sscal.o  \
BLAS/sspmv.o  \
BLAS/sspr.o  \
BLAS/sspr2.o  \
BLAS/sswap.o  \
BLAS/ssymm.o  \
BLAS/ssymv.o  \
BLAS/ssyr.o  \
BLAS/ssyr2.o  \
BLAS/ssyr2k.o  \
BLAS/ssyrk.o  \
BLAS/stbmv.o  \
BLAS/stbsv.o  \
BLAS/stpmv.o  \
BLAS/stpsv.o  \
BLAS/strmm.o  \
BLAS/strmv.o  \
BLAS/strsm.o  \
BLAS/strsv.o  \
BLAS/xerbla.o  \
BLAS/zaxpy.o  \
BLAS/zcopy.o  \
BLAS/zdotc.o  \
BLAS/zdotu.o  \
BLAS/zdrot.o  \
BLAS/zdscal.o  \
BLAS/zgbmv.o  \
BLAS/zgemm.o  \
BLAS/zgemv.o  \
BLAS/zgerc.o  \
BLAS/zgeru.o  \
BLAS/zhbmv.o  \
BLAS/zhemm.o  \
BLAS/zhemv.o  \
BLAS/zher.o  \
BLAS/zher2.o  \
BLAS/zher2k.o  \
BLAS/zherk.o  \
BLAS/zhpmv.o  \
BLAS/zhpr.o  \
BLAS/zhpr2.o  \
BLAS/zrotg.o  \
BLAS/zscal.o  \
BLAS/zswap.o  \
BLAS/zsymm.o  \
BLAS/zsyr2k.o  \
BLAS/zsyrk.o  \
BLAS/ztbmv.o  \
BLAS/ztbsv.o  \
BLAS/ztpmv.o  \
BLAS/ztpsv.o  \
BLAS/ztrmm.o  \
BLAS/ztrmv.o  \
BLAS/ztrsm.o  \
BLAS/ztrsv.o  


LAPACKOBJS = \
lapack/util/xerbla.o \
lapack/util/dlamch.o \
lapack/util/ieeeck.o \
lapack/util/iladlc.o \
lapack/util/iladlr.o \
lapack/util/ilaenv.o \
lapack/util/iparmq.o \
lapack/util/lsame.o \
lapack/double/dgebak.o \
lapack/double/dgebal.o \
lapack/double/dgeev.o \
lapack/double/dgehd2.o \
lapack/double/dgehrd.o \
lapack/double/dgesv.o \
lapack/double/dgetf2.o \
lapack/double/dgetrf.o \
lapack/double/dgetri.o \
lapack/double/dgetrs.o \
lapack/double/dhseqr.o \
lapack/double/disnan.o \
lapack/double/dlabad.o \
lapack/double/dlacpy.o \
lapack/double/dladiv.o \
lapack/double/dlae2.o \
lapack/double/dlaev2.o \
lapack/double/dlaexc.o \
lapack/double/dlahqr.o \
lapack/double/dlahr2.o \
lapack/double/dlaisnan.o \
lapack/double/dlaln2.o \
lapack/double/dlange.o \
lapack/double/dlanst.o \
lapack/double/dlansy.o \
lapack/double/dlanv2.o \
lapack/double/dlapy2.o \
lapack/double/dlaqr0.o \
lapack/double/dlaqr1.o \
lapack/double/dlaqr2.o \
lapack/double/dlaqr3.o \
lapack/double/dlaqr4.o \
lapack/double/dlaqr5.o \
lapack/double/dlarfb.o \
lapack/double/dlarfg.o \
lapack/double/dlarf.o \
lapack/double/dlarft.o \
lapack/double/dlarfx.o \
lapack/double/dlartg.o \
lapack/double/dlascl.o \
lapack/double/dlaset.o \
lapack/double/dlasr.o \
lapack/double/dlasrt.o \
lapack/double/dlassq.o \
lapack/double/dlaswp.o \
lapack/double/dlasy2.o \
lapack/double/dlasyf.o \
lapack/double/dlatrd.o \
lapack/double/dorg2l.o \
lapack/double/dorg2r.o \
lapack/double/dorghr.o \
lapack/double/dorgql.o \
lapack/double/dorgqr.o \
lapack/double/dorgtr.o \
lapack/double/dorm2r.o \
lapack/double/dormhr.o \
lapack/double/dormqr.o \
lapack/double/dsteqr.o \
lapack/double/dsterf.o \
lapack/double/dsyev.o \
lapack/double/dsytd2.o \
lapack/double/dsytf2.o \
lapack/double/dsytrd.o \
lapack/double/dsytrf.o \
lapack/double/dsytri.o \
lapack/double/dtrevc.o \
lapack/double/dtrexc.o \
lapack/double/dtrti2.o \
lapack/double/dtrtri.o


blas.a: $(BLASOBJS)
	ar -r blas.a $(BLASOBJS)

lapack.a: $(LAPACKOBJS)
	ar -r lapack.a $(LAPACKOBJS)

qhull-2011.1/lib/libqhullstatic.a:
	cd qhull-2011.1/build ; cmake .. ; cd .. ; make lib/libqhullstatic.a

amoeba.a: amoeba.o amotry.o nrutil.o mnbrak.o frprmn.o brent.o dbrent.o linmin.o dlinmin.o f1dim.o df1dim.o ft.o
	ar -r amoeba.a amoeba.o amotry.o nrutil.o frprmn.o brent.o dbrent.o linmin.o dlinmin.o f1dim.o mnbrak.o df1dim.o ft.o


hopBorders: hopBorders.o util.o dcd.o pdb.o alignSet.o surfacesPS.o dovoronoi.o process_lipids.o qhull-2011.1/lib/libqhullstatic.a
	$(CC) -o hopBorders hopBorders.o util.o pdb.o dcd.o alignSet.o surfacesPS.o dovoronoi.o process_lipids.o $(LDFLAGS) $(GSLLIB) qhull-2011.1/lib/libqhullstatic.a

processHop: processHop.o util.o $(LAPACK_BUILD)
	$(CC) -o processHop processHop.o util.o $(LDFLAGS) $(GSLLIB)

rehop: rehop.o util.o dcd.o pdb.o alignSet.o surfacesPS.o dovoronoi.o orderp.o qhull-2011.1/lib/libqhullstatic.a
	$(CC) -o rehop rehop.o util.o pdb.o dcd.o alignSet.o surfacesPS.o orderp.o dovoronoi.o $(LDFLAGS) $(GSLLIB) qhull-2011.1/lib/libqhullstatic.a

histo: histo.o util.o
	$(CC) -o histo histo.o util.o $(LDFLAGS)

lifetime: lifetime.o util.o
	$(CC) -o lifetime lifetime.o util.o $(LDFLAGS)

outputIndex: outputIndex.o process_lipids.o util.o dcd.o pdb.o alignSet.o 
	$(CC) -o outputIndex outputIndex.o process_lipids.o util.o dcd.o pdb.o alignSet.o $(LDFLAGS)

avhist: avhist.o util.o
	$(CC) -o avhist avhist.o util.o
