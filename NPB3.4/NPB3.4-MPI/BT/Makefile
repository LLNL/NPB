SHELL=/bin/sh
BENCHMARK=bt
BENCHMARKU=BT

include ../config/make.def


OBJS = bt.o bt_data.o make_set.o initialize.o exact_solution.o \
       exact_rhs.o set_constants.o adi.o define.o copy_faces.o \
       rhs.o solve_subs.o x_solve.o y_solve.o z_solve.o add.o \
       error.o verify.o setup_mpi.o mpinpb.o \
       ${COMMON}/print_results.o ${COMMON}/timers.o

include ../sys/make.common

# npbparams.h is included by bt_data module (via bt_data.o)

${PROGRAM}: config $(OBJS)
	@if [ x$(SUBTYPE) = xfull -o x$(SUBTYPE) = xFULL ] ; then	\
		${MAKE} bt-full;		\
	elif [ x$(SUBTYPE) = xsimple -o x$(SUBTYPE) = xSIMPLE ] ; then	\
		${MAKE} bt-simple;		\
	elif [ x$(SUBTYPE) = xfortran -o x$(SUBTYPE) = xFORTRAN ] ; then \
		${MAKE} bt-fortran;		\
	elif [ x$(SUBTYPE) = xepio -o x$(SUBTYPE) = xEPIO ] ; then	\
		${MAKE} bt-epio;		\
	else					\
		${MAKE} bt-bt;			\
	fi

bt-bt: ${OBJS} btio.o
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} btio.o ${FMPI_LIB}

bt-full: ${OBJS} full_mpiio.o btio_common.o
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM}.mpi_io_full ${OBJS} btio_common.o full_mpiio.o ${FMPI_LIB}

bt-simple: ${OBJS} simple_mpiio.o btio_common.o
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM}.mpi_io_simple ${OBJS} btio_common.o simple_mpiio.o ${FMPI_LIB}

bt-fortran: ${OBJS} fortran_io.o btio_common.o
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM}.fortran_io ${OBJS} btio_common.o fortran_io.o ${FMPI_LIB}

bt-epio: ${OBJS} epio.o btio_common.o
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM}.ep_io ${OBJS} btio_common.o epio.o ${FMPI_LIB}

.f.o:
	${FCOMPILE} $<

.c.o:
	${CCOMPILE} $<


bt.o:             bt.f  bt_data.o mpinpb.o
make_set.o:       make_set.f  bt_data.o mpinpb.o
initialize.o:     initialize.f  bt_data.o
exact_solution.o: exact_solution.f  bt_data.o
exact_rhs.o:      exact_rhs.f  bt_data.o
set_constants.o:  set_constants.f  bt_data.o
adi.o:            adi.f  bt_data.o
define.o:         define.f  bt_data.o
copy_faces.o:     copy_faces.f  bt_data.o mpinpb.o
rhs.o:            rhs.f  bt_data.o
x_solve.o:        x_solve.f  bt_data.o mpinpb.o
y_solve.o:        y_solve.f  bt_data.o mpinpb.o
z_solve.o:        z_solve.f  bt_data.o mpinpb.o
solve_subs.o:     solve_subs.f
add.o:            add.f  bt_data.o
error.o:          error.f  bt_data.o mpinpb.o
verify.o:         verify.f  bt_data.o mpinpb.o
setup_mpi.o:      setup_mpi.f  bt_data.o mpinpb.o
btio.o:           btio.f  bt_data.o
btio_common.o:    btio_common.f  bt_data.o mpinpb.o
fortran_io.o:     fortran_io.f  bt_data.o mpinpb.o
simple_mpiio.o:   simple_mpiio.f  bt_data.o mpinpb.o
full_mpiio.o:     full_mpiio.f  bt_data.o mpinpb.o
epio.o:           epio.f  bt_data.o mpinpb.o
bt_data.o:        bt_data.f mpinpb.o npbparams.h
mpinpb.o:         mpinpb.f

clean:
	- rm -f *.o *.mod *~ mputil*
	- rm -f  npbparams.h core
