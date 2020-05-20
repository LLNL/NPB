SHELL=/bin/sh
BENCHMARK=lu
BENCHMARKU=LU

include ../config/make.def

OBJS = lu.o lu_data.o init_comm.o read_input.o bcast_inputs.o \
       proc_grid.o neighbors.o nodedim.o subdomain.o setcoeff.o \
       setbv.o exact.o setiv.o erhs.o ssor.o exchange_1.o exchange_3.o \
       exchange_4.o exchange_5.o exchange_6.o rhs.o l2norm.o \
       jacld.o blts.o jacu.o buts.o mpinpb.o \
       error.o pintgr.o verify.o ${COMMON}/print_results.o ${COMMON}/timers.o

include ../sys/make.common


# npbparams.h is included by lu_data module (via lu_data.o)

${PROGRAM}: config $(OBJS)
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${FMPI_LIB}

.f.o :
	${FCOMPILE} $<

lu.o:		lu.f lu_data.o mpinpb.o
bcast_inputs.o:	bcast_inputs.f lu_data.o mpinpb.o
blts.o:	        blts.f lu_data.o
buts.o:	        buts.f lu_data.o
erhs.o:		erhs.f lu_data.o
error.o:	error.f lu_data.o mpinpb.o
exact.o:	exact.f lu_data.o
exchange_1.o:	exchange_1.f lu_data.o mpinpb.o
exchange_3.o:	exchange_3.f lu_data.o mpinpb.o
exchange_4.o:	exchange_4.f lu_data.o mpinpb.o
exchange_5.o:	exchange_5.f lu_data.o mpinpb.o
exchange_6.o:	exchange_6.f lu_data.o mpinpb.o
init_comm.o:	init_comm.f lu_data.o mpinpb.o 
jacld.o:	jacld.f lu_data.o
jacu.o:		jacu.f lu_data.o
l2norm.o:	l2norm.f lu_data.o mpinpb.o
neighbors.o:	neighbors.f lu_data.o
nodedim.o:	nodedim.f
pintgr.o:	pintgr.f lu_data.o mpinpb.o
proc_grid.o:	proc_grid.f lu_data.o mpinpb.o
read_input.o:	read_input.f lu_data.o mpinpb.o
rhs.o:		rhs.f lu_data.o
setbv.o:	setbv.f lu_data.o
setiv.o:	setiv.f lu_data.o
setcoeff.o:	setcoeff.f lu_data.o
ssor.o:		ssor.f lu_data.o mpinpb.o
subdomain.o:	subdomain.f lu_data.o mpinpb.o
verify.o:	verify.f lu_data.o
lu_data.o:      lu_data.f mpinpb.o npbparams.h 
mpinpb.o:       mpinpb.f

clean:
	- rm -f npbparams.h
	- rm -f *.o *.mod *~
