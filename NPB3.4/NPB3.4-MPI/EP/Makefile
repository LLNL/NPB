SHELL=/bin/sh
BENCHMARK=ep
BENCHMARKU=EP

include ../config/make.def

OBJS = ep.o ep_data.o mpinpb.o \
	${COMMON}/print_results.o ${COMMON}/${RAND}.o ${COMMON}/timers.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${FMPI_LIB}

.f.o:
	${FCOMPILE} $<

ep.o:		ep.f ep_data.o mpinpb.o
ep_data.o:	ep_data.f npbparams.h
mpinpb.o:	mpinpb.f

clean:
	- rm -f *.o *~ *.mod
	- rm -f npbparams.h core



