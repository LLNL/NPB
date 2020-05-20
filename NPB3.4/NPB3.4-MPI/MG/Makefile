SHELL=/bin/sh
BENCHMARK=mg
BENCHMARKU=MG

include ../config/make.def

OBJS = mg.o mg_data.o mpinpb.o ${COMMON}/print_results.o  \
       ${COMMON}/${RAND}.o ${COMMON}/timers.o

include ../sys/make.common

${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${FMPI_LIB}

.f.o:
	${FCOMPILE} $<

mg.o:		mg.f  mg_data.o mpinpb.o
mg_data.o:	mg_data.f mpinpb.o npbparams.h
mpinpb.o:	mpinpb.f

clean:
	- rm -f *.o *.mod *~ 
	- rm -f npbparams.h core



