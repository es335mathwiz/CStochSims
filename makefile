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


stochRun:	stochRun.o myNewt.o dtime.o\
		stackC.o juillard.o stochProto.o ranlib.o
	$(FC) -o stochRun -g  stochRun.o myNewt.o \
		stackC.o juillard.o stochProto.o ranlib.o dtime.o $(SPARSEAMALIB) $(LAPACKLIBS)



clean: 
	rm -f *.o 
