SHELL=/bin/sh
BENCHMARK=lu-mz
BENCHMARKU=LU-MZ

include ../config/make.def

OBJS = lu.o lu_data.o read_input.o setup_omp.o ompnpb.o \
       setcoeff.o setbv.o exact.o setiv.o \
       erhs.o ssor.o zone_setup.o exch_qbc.o \
       rhs.o l2norm.o jacld.o blts.o jacu.o buts.o \
       syncs.o error.o pintgr.o verify.o ${COMMON}/print_results.o \
       ${COMMON}/timers.o ${COMMON}/wtime.o

include ../sys/make.common


# npbparams.h is included by lu_data module (via lu_data.o)

${PROGRAM}: config
	@if [ x$(VERSION) = xvec -o x$(VERSION) = xVEC ] ; then	\
		${MAKE} VEC=_vec exec;		\
	else					\
		${MAKE} exec;			\
	fi

exec: $(OBJS)
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o :
	${FCOMPILE} $<

lu.o:		lu.f lu_data.o ompnpb.o
blts.o:		blts.f
buts.o:		buts.f
erhs.o:		erhs.f lu_data.o
zone_setup.o:   zone_setup.f lu_data.o ompnpb.o
exch_qbc.o:     exch_qbc.f lu_data.o ompnpb.o
error.o:	error.f
exact.o:	exact.f lu_data.o
jacld.o:	jacld.f lu_data.o
jacu.o:		jacu.f lu_data.o
l2norm.o:	l2norm.f
pintgr.o:	pintgr.f lu_data.o
read_input.o:	read_input.f lu_data.o ompnpb.o
rhs.o:		rhs.f lu_data.o
setbv.o:	setbv.f lu_data.o
setiv.o:	setiv.f lu_data.o
setcoeff.o:	setcoeff.f lu_data.o
ssor.o:		ssor.f lu_data.o
verify.o:	verify.f lu_data.o
setup_omp.o:	setup_omp.f lu_data.o ompnpb.o
syncs.o:	syncs.f
ompnpb.o:	ompnpb.f  lu_data.o
lu_data.o:	lu_data.f  npbparams.h

clean:
	- rm -f npbparams.h
	- rm -f *.o *~ *.mod
