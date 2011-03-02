# Makefile

CC=gcc
SRC=convolution2d.c fft.c registration.c main.c
OBJ=${${SRC}=%.o}
OUT=registration

PAPIDIR=/mnt/jc5/CS259/papi

INCFLAGS=-I${PAPIDIR}
CFLAGS=-g -pg
LDFLAGS=-L${PAPIDIR} -lm -lpapi -lutil_papi

default: bin

bin: ${OBJ}
	${CC} ${CFLAGS} ${INCFLAGS} ${OBJ} -o ${OUT} ${LDFLAGS}

${%.o}: ${SRC}
	${CC} ${CFLAGS} ${INCFLAGS} -c $< -o $@

clean:
	rm ${OUT} ${OBJ}
