SHELL=/bin/sh
BENCHMARK=ft
BENCHMARKU=FT

include ../config/make.def

include ../sys/make.common

OBJS = ft.o ft_data.o ${COMMON}/${RAND}.o ${COMMON}/print_results.o \
       ${COMMON}/timers.o ${COMMON}/wtime.o

${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}



.f.o:
	${FCOMPILE} $<

ft.o:		ft.f  ft_data.o
ft_data.o:	ft_data.f  npbparams.h

clean:
	- rm -f *.o *~ mputil* *.mod
	- rm -f ft npbparams.h core
	- if [ -d rii_files ]; then rm -r rii_files; fi
