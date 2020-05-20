SHELL=/bin/sh
BENCHMARK=lu-mz
BENCHMARKU=LU-MZ

include ../config/make.def

OBJS = lu.o lu_data.o read_input.o setup_mpi.o \
       setcoeff.o setbv.o exact.o setiv.o \
       erhs.o ssor.o zone_setup.o exch_qbc.o \
       rhs.o l2norm.o jacld.o blts.o jacu.o buts.o \
       syncs.o error.o pintgr.o verify.o proc_grid.o neighbors.o \
       subdomain.o exchange_1.o exchange_3.o  exchange_4.o \
       error_cond.o mpinpb.o \
       ${COMMON}/print_results.o ${COMMON}/timers.o

include ../sys/make.common


# npbparams.h is included by lu_data module (via lu_data.o)

${PROGRAM}: config $(OBJS)
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o :
	${FCOMPILE} $<

lu.o:		lu.f lu_data.o mpinpb.o
blts.o:		blts.f
buts.o:		buts.f
erhs.o:		erhs.f lu_data.o
zone_setup.o:   zone_setup.f lu_data.o mpinpb.o
exch_qbc.o:     exch_qbc.f lu_data.o mpinpb.o
error.o:	error.f lu_data.o mpinpb.o
exact.o:	exact.f lu_data.o
jacld.o:	jacld.f lu_data.o
jacu.o:		jacu.f lu_data.o
l2norm.o:	l2norm.f lu_data.o mpinpb.o
pintgr.o:	pintgr.f lu_data.o mpinpb.o
read_input.o:	read_input.f lu_data.o mpinpb.o
rhs.o:		rhs.f lu_data.o
setbv.o:	setbv.f lu_data.o
setiv.o:	setiv.f lu_data.o
setcoeff.o:	setcoeff.f lu_data.o
ssor.o:		ssor.f lu_data.o mpinpb.o
verify.o:	verify.f lu_data.o
setup_mpi.o:	setup_mpi.f lu_data.o mpinpb.o
proc_grid.o:	proc_grid.f lu_data.o mpinpb.o
neighbors.o:	neighbors.f lu_data.o mpinpb.o
subdomain.o:	subdomain.f lu_data.o mpinpb.o
exchange_1.o:	exchange_1.f lu_data.o mpinpb.o
exchange_3.o:	exchange_3.f lu_data.o mpinpb.o
exchange_4.o:	exchange_4.f lu_data.o mpinpb.o
error_cond.o:   error_cond.f  mpinpb.o
syncs.o:	syncs.f lu_data.o
mpinpb.o:       mpinpb.f lu_data.o
lu_data.o:      lu_data.f npbparams.h

clean:
	- rm -f *.o *~ *.mod mputil*
	- rm -f npbparams.h core
