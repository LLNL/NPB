SHELL=/bin/sh
BENCHMARK=cg
BENCHMARKU=CG

include ../config/make.def

OBJS = cg.o cg_data.o mpinpb.o ${COMMON}/print_results.o  \
       ${COMMON}/${RAND}.o ${COMMON}/timers.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${FMPI_LIB}

.f.o:
	${FCOMPILE} $<

cg.o:		cg.f  cg_data.o mpinpb.o
cg_data.o:	cg_data.f mpinpb.o npbparams.h
mpinpb.o:	mpinpb.f

clean:
	- rm -f *.o *.mod *~ 
	- rm -f npbparams.h core



