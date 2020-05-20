SHELL=/bin/sh
BENCHMARK=sp-mz
BENCHMARKU=SP-MZ

include ../config/make.def


OBJS = sp.o sp_data.o initialize.o exact_solution.o exact_rhs.o \
       set_constants.o adi.o rhs.o zone_setup.o \
       x_solve.o ninvr.o y_solve.o pinvr.o exch_qbc.o   \
       z_solve.o tzetar.o add.o txinvr.o error.o verify.o \
       setup_omp.o ompnpb.o \
       ${COMMON}/print_results.o ${COMMON}/timers.o ${COMMON}/wtime.o

include ../sys/make.common

# npbparams.h is included by sp_data module (via sp_data.o)

${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o:
	${FCOMPILE} $<

sp.o:             sp.f  sp_data.o ompnpb.o
initialize.o:     initialize.f  sp_data.o
exact_solution.o: exact_solution.f  sp_data.o
exact_rhs.o:      exact_rhs.f  sp_data.o
set_constants.o:  set_constants.f  sp_data.o
adi.o:            adi.f  sp_data.o
rhs.o:            rhs.f  sp_data.o
zone_setup.o:     zone_setup.f sp_data.o ompnpb.o
x_solve.o:        x_solve.f  sp_data.o
ninvr.o:          ninvr.f  sp_data.o
y_solve.o:        y_solve.f  sp_data.o
pinvr.o:          pinvr.f  sp_data.o
z_solve.o:        z_solve.f  sp_data.o
tzetar.o:         tzetar.f  sp_data.o
add.o:            add.f  sp_data.o
txinvr.o:         txinvr.f  sp_data.o
error.o:          error.f  sp_data.o
verify.o:         verify.f  sp_data.o ompnpb.o
exch_qbc.o:       exch_qbc.f sp_data.o ompnpb.o
setup_omp.o:      setup_omp.f  sp_data.o ompnpb.o
ompnpb.o:         ompnpb.f  sp_data.o
sp_data.o:        sp_data.f  npbparams.h

clean:
	- rm -f *.o *~ mputil* *.mod
	- rm -f npbparams.h core
