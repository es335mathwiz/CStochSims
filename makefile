
#identify operating system
UNAME= $(shell uname)
NUWEBFLAGS = -t
SPAMADIR=../sparseAMA

ifeq ($(UNAME),Linux)
#compilers
CC = gcc
FCFLAGS = -c -O2 -I $(SPAMADIR)/src/main/include   -I /msu/res5/software/myUsr/include/ -I/msu/res5/software/myUsr/include/ -I /msu/res1/Software/matio-1.5.1/src
FCFLAGS = -c -g -Wall -I $(SPAMADIR)/src/main/include   -I /msu/res5/software/myUsr/include -I/msu/res5/software/myUsr/include/ -I /msu/res1/Software/matio-1.5.1/src
#lapack
LAPACKLIBS=   -L /msu/res5/software/ARPACK96forCluster -larpack_linux -L/msu/res5/software/lapackGithubForCluster -llapack -lrefblas
CUNITLIBS= -L/msu/res5/software/myUsr/lib/ -l cunit
MATIOLIBS= -L/msu/res1/Software/matio-1.5.1/src/.libs/ -lmatio  -lhdf5

endif

ifeq ($(UNAME),Darwin)
#compilers
CC = gcc-6
FCFLAGS = -c -O2 -I$(SPAMADIR)/src/main/include   -I /Users/garyanderson/myUsr/include/ -I/Users/garyanderson/myUsr/include/\
-I /usr/local/Cellar/libmatio/1.5.10/include
FCFLAGS = -c -Wall -g -I $(SPAMADIR)/src/main/include   -I /Users/garyanderson/myUsr/include/ -I/Users/garyanderson/myUsr/include/\
-I /usr/local/Cellar/libmatio/1.5.10/include
#lapack
LAPACKLIBS=  -L /Users/garyanderson/ARPACK96/  -larpack_MACOS -L /Users/garyanderson/lapack-release/ -llapack -lrefblas
CUNITLIBS= -L /Users/garyanderson/myUsr/lib -l cunit
MATIOLIBS= -L/usr/local/Cellar/libmatio/1.5.10/lib -lmatio 

endif


#compilers
FC = gfortran
SPARSEAMALIB= -L../sparseAMA -lsparseAMA
STOCHSIMSLIB= -L./ -lstochSims
.w.c:
	nuweb $(NUWEBFLAGS) $*

.w.o:
	make $*.c
	make $*.o

.c.o:
	$(CC) $(FCFLAGS) -c $*.c

.f.o:
	$(FC) $(FCFLAGS) -c $*.f

.PHONY: Build

Build: stochRun runrbcTryC
#	$(FC) -o stochRun -g  stochRun.o juillard.o $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS) $(CUNITLIBS) $(MATIOLIBS)
#	$(FC) -o runrbcTryC -g  runrbcTryC.o $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS) $(CUNITLIBS) $(MATIOLIBS)
#	$(FC) -o stochSimsUnitTests -g  stochSimsUnitTests.o juillard.o $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS) $(CUNITLIBS) $(MATIOLIBS)

myNewt.o:			 stackC.w
		nuweb $(NUWEBFLAGS)  stackC.w
	$(CC) $(FCFLAGS) -c myNewt.c

juillard.o:	juillard.c
	$(CC) $(FCFLAGS) -c juillard.c


stochRun:	stochRun.o  juillard.o libstochSims.a
	$(FC) -o stochRun -g  stochRun.o juillard.o $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS)  $(CUNITLIBS) $(MATIOLIBS)

stochSimsUnitTests:	stochSimsUnitTests.o  rbcTryC.o rbcTryCDrv.o rbcTryCData.o rbcTryCShocks.o rbcTryCSupport.o libstochSims.a
	$(FC) -o stochSimsUnitTests -g  stochSimsUnitTests.o  rbcTryC.o rbcTryCDrv.o rbcTryCData.o rbcTryCShocks.o rbcTryCSupport.o  $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS)  $(CUNITLIBS) $(MATIOLIBS)

runrbcTryC:	runrbcTryC.o  rbcTryC.o rbcTryCDrv.o rbcTryCData.o rbcTryCShocks.o rbcTryCSupport.o libstochSims.a
	$(FC) -o runrbcTryC -g  runrbcTryC.o  rbcTryC.o rbcTryCDrv.o rbcTryCData.o rbcTryCShocks.o rbcTryCSupport.o  $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS)  $(CUNITLIBS) $(MATIOLIBS)

libstochSims.a:	myNewt.o \
		stackC.o ranlib.o compXEtm1.o  \
	generateNextXTMinusOne.o generatePathX.o  \
	generateNextXT.o stochSims.o generateDraws.o ranlib.o
	ar -cvq libstochSims.a myNewt.o \
		stackC.o  ranlib.o compXEtm1.o \
	generateNextXTMinusOne.o generatePathX.o  \
	generateNextXT.o  stochSims.o generateDraws.o ranlib.o


clean: 
	rm -f *.o stochRun stochSimsUnitTests libstochSims.a

