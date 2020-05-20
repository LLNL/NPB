SHELL=/bin/sh
BENCHMARK=bt-mz
BENCHMARKU=BT-MZ

include ../config/make.def


OBJS = bt.o bt_data.o initialize.o exact_solution.o exact_rhs.o \
       set_constants.o adi.o  rhs.o zone_setup.o  \
       x_solve.o y_solve.o  exch_qbc.o solve_subs.o \
       z_solve.o add.o error.o verify.o setup_mpi.o \
       mpinpb.o error_cond.o \
       ${COMMON}/print_results.o ${COMMON}/timers.o

include ../sys/make.common

# npbparams.h is included by bt_data module (via bt_data.o)

${PROGRAM}: config ${OBJS}
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o:
	${FCOMPILE} $<

bt.o:             bt.f  bt_data.o mpinpb.o
initialize.o:     initialize.f  bt_data.o
exact_solution.o: exact_solution.f  bt_data.o
exact_rhs.o:      exact_rhs.f  bt_data.o
set_constants.o:  set_constants.f  bt_data.o
adi.o:            adi.f  bt_data.o
rhs.o:            rhs.f  bt_data.o
zone_setup.o:     zone_setup.f  bt_data.o mpinpb.o
x_solve.o:        x_solve.f  bt_data.o
y_solve.o:        y_solve.f  bt_data.o
z_solve.o:        z_solve.f  bt_data.o
solve_subs.o:     solve_subs.f
add.o:            add.f  bt_data.o
error.o:          error.f  bt_data.o
verify.o:         verify.f  bt_data.o mpinpb.o
exch_qbc.o:       exch_qbc.f  bt_data.o mpinpb.o
setup_mpi.o:      setup_mpi.f  bt_data.o mpinpb.o
error_cond.o:     error_cond.f  mpinpb.o
mpinpb.o:         mpinpb.f  bt_data.o
bt_data.o:        bt_data.f  npbparams.h

clean:
	- rm -f *.o *~ *.mod mputil*
	- rm -f npbparams.h core
