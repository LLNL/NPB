SHELL=/bin/sh
BENCHMARK=sp
BENCHMARKU=SP
BLK=
BLKFAC=0

include ../config/make.def


OBJS = sp.o sp_data.o initialize.o exact_solution.o exact_rhs.o \
       set_constants.o adi.o rhs.o work_lhs$(BLK).o     \
       x_solve$(BLK).o ninvr.o y_solve$(BLK).o pinvr.o    \
       z_solve$(BLK).o tzetar.o add.o txinvr.o error.o verify.o  \
       ${COMMON}/print_results.o ${COMMON}/timers.o ${COMMON}/wtime.o

include ../sys/make.common

# npbparams.h is included by sp_data module (via sp_data.o)

${PROGRAM}: config
	@ver=$(VERSION); bfac=$${ver:3}; \
	if [ x$${ver:0:3} = xblk -o x$${ver:0:3} = xBLK ] ; then	\
		${MAKE} BLK=_blk BLKFAC=$${bfac:-8} exec;		\
	else					\
		${MAKE} exec;			\
	fi

exec: $(OBJS)
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o:
	${FCOMPILE} $<

blk_par.h: FORCE
	sed -e 's/= 0/= $(BLKFAC)/' blk_par0.h > blk_par.h_wk
	@ if ! `diff blk_par.h_wk blk_par.h > /dev/null 2>&1`; then \
	mv -f blk_par.h_wk blk_par.h; else rm -f blk_par.h_wk; fi
FORCE:

sp.o:             sp.f  sp_data.o blk_par.h
initialize.o:     initialize.f  sp_data.o
exact_solution.o: exact_solution.f  sp_data.o
exact_rhs.o:      exact_rhs.f  sp_data.o
set_constants.o:  set_constants.f  sp_data.o
adi.o:            adi.f  sp_data.o
rhs.o:            rhs.f  sp_data.o
x_solve$(BLK).o:  x_solve$(BLK).f  sp_data.o work_lhs$(BLK).o
ninvr.o:          ninvr.f  sp_data.o
y_solve$(BLK).o:  y_solve$(BLK).f  sp_data.o work_lhs$(BLK).o
pinvr.o:          pinvr.f  sp_data.o
z_solve$(BLK).o:  z_solve$(BLK).f  sp_data.o work_lhs$(BLK).o
tzetar.o:         tzetar.f  sp_data.o
add.o:            add.f  sp_data.o
txinvr.o:         txinvr.f  sp_data.o
error.o:          error.f  sp_data.o
verify.o:         verify.f  sp_data.o
work_lhs$(BLK).o: work_lhs$(BLK).f  sp_data.o blk_par.h
sp_data.o:        sp_data.f  npbparams.h

clean:
	- rm -f *.o *~ *.mod mputil*
	- rm -f npbparams.h core blk_par.h
	- if [ -d rii_files ]; then rm -r rii_files; fi
