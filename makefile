#identify operating system
UNAME= $(shell uname)
NUWEBFLAGS = -t
SPAMADIR=../sparseAMA

ifeq ($(UNAME),Linux)
#compilers
CC = gcc
FCFLAGS = -c -O2 -I $(SPAMADIR)/src/main/include   -I /msu/res5/software/myUsr/include/ -I../stackStochSims/
FCFLAGS = -c -g -Wall -I $(SPAMADIR)/src/main/include   -I /msu/res5/software/myUsr/include -I../stackStochSims/
#lapack
LAPACKLIBS=   -L /msu/res5/software/ARPACK96forCluster -larpack_linux -L/msu/res5/software/lapackGithubForCluster -llapack -lrefblas
CUNITLIBS= -L /msu/res5/software/myUsr/lib/ -l cunit
endif

ifeq ($(UNAME),Darwin)
#compilers
CC = gcc-6
FCFLAGS = -c -O2 -I$(SPAMADIR)/src/main/include   -I /Users/garyanderson/myUsr/include/ -I../stackStochSims/
FCFLAGS = -c -Wall -g -I $(SPAMADIR)/src/main/include   -I /Users/garyanderson/myUsr/include/ -I../stackStochSims/
#lapack
LAPACKLIBS=  -L /Users/garyanderson/ARPACK96/  -larpack_MACOS -L /Users/garyanderson/lapack-release/ -llapack -lrefblas
CUNITLIBS= -L /Users/garyanderson/myUsr/lib -l cunit
endif

echo:
	@echo $(FCFLAGS)

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
	nuweb $(NUWEBFLAGS)  stochProto.w
	nuweb $(NUWEBFLAGS)  stochRun.w
	nuweb $(NUWEBFLAGS)  stackC.w
	$(CC) $(FCFLAGS) -c $*.c

.f.o:
	nuweb $(NUWEBFLAGS)  stochProto.w
	nuweb $(NUWEBFLAGS)  stochRun.w
	nuweb $(NUWEBFLAGS)  stackC.w
	$(FC) $(FCFLAGS) -c $*.f



myNewt.o:			 stackC.w
		nuweb $(NUWEBFLAGS)  stackC.w
	$(CC) $(FCFLAGS) -c myNewt.c

juillard.o:	juillard.c
	$(CC) $(FCFLAGS) -c juillard.c


stochRun:	stochRun.o libstochSims.a juillard.o
	$(FC) -o stochRun -g  stochRun.o juillard.o $(STOCHSIMSLIB) $(SPARSEAMALIB) $(LAPACKLIBS)

libstochSims.a:	myNewt.o dtime.o\
		stackC.o stochProto.o ranlib.o
	ar -cvq libstochSims.a myNewt.o dtime.o\
		stackC.o stochProto.o ranlib.o


clean: 
	rm -f *.o 
