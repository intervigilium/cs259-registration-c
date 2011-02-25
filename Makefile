# Makefile

CC=gcc
SRC=reg_vec_mi.c
OUT=registration

PAPIDIR=/mnt/jc5/CS259/papi
INCLUDEDIR=-I${PAPIDIR}

CFLAGS=-g -pg
LDFLAGS=-L${PAPIDIR} -lm -lpapi -lutil_papi

all:
	${CC} ${CFLAGS} ${INCLUDEDIR} ${SRC} -o ${OUT} ${LDFLAGS}

clean:
	rm ${OUT}
