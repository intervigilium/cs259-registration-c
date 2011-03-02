# Makefile

PAPIDIR=/mnt/jc5/CS259/papi

CC=gcc
SRC=convolution2d.c fft.c registration.c main.c
OBJ=$(SRC:.c=.o)
OUT=registration

INCFLAGS=-I$(CURDIR) -I$(PAPIDIR)
CFLAGS=-g -pg
LDFLAGS=-L$(PAPIDIR) -lm -lpapi -lutil_papi

default: bin

bin: $(OBJ)
	$(CC) $(CFLAGS) $(INCFLAGS) $(OBJ) -o $(OUT) $(LDFLAGS)

%.o: $(SRC)
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(OUT)
