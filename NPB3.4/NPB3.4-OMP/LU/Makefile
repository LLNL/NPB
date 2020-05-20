SHELL=/bin/sh
BENCHMARK=lu
BENCHMARKU=LU
SVER=

include ../config/make.def

OBJS = lu.o lu_data.o read_input.o \
       domain.o setcoeff.o setbv.o exact.o setiv.o \
       erhs.o ssor$(SVER).o rhs.o l2norm.o \
       jacld.o blts.o jacu.o buts.o error.o \
       pintgr.o verify.o ${COMMON}/print_results.o \
       ${COMMON}/timers.o ${COMMON}/wtime.o

include ../sys/make.common

# npbparams.h is included by lu_data module (via lu_data.o)

${PROGRAM}: config
	@if [ x$(VERSION) = xdoac -o x$(VERSION) = xDOAC ] ; then	\
		${MAKE} SVER=_doac exec;			\
	elif [ x$(VERSION) = xhp -o x$(VERSION) = xHP ] ; then	\
		${MAKE} SVER=_hp exec;				\
	else							\
		${MAKE} exec-sync;				\
	fi

exec-sync: $(OBJS) syncs.o
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} syncs.o ${F_LIB}

exec: $(OBJS)
	${FLINK} ${FLINKFLAGS} -o ${PROGRAM} ${OBJS} ${F_LIB}

.f.o :
	${FCOMPILE} $<

lu.o:		lu.f lu_data.o
blts.o:		blts.f
buts.o:		buts.f
erhs.o:		erhs.f lu_data.o
error.o:	error.f lu_data.o
exact.o:	exact.f lu_data.o
jacld.o:	jacld.f lu_data.o
jacu.o:		jacu.f lu_data.o
l2norm.o:	l2norm.f
pintgr.o:	pintgr.f lu_data.o
read_input.o:	read_input.f lu_data.o
rhs.o:		rhs.f lu_data.o
setbv.o:	setbv.f lu_data.o
setiv.o:	setiv.f lu_data.o
setcoeff.o:	setcoeff.f lu_data.o
ssor$(SVER).o:	ssor$(SVER).f lu_data.o
domain.o:	domain.f lu_data.o
verify.o:	verify.f lu_data.o
syncs.o:	syncs.f lu_data.o
lu_data.o:	lu_data.f npbparams.h

clean:
	- rm -f *.o *~ *.mod npbparams.h
	- if [ -d rii_files ]; then rm -r rii_files; fi
