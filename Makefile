# Makefile

PAPIDIR=/mnt/jc5/CS259/papi

CC=gcc
SRC=fft.c convolution2d.c registration.c main.c
OBJ=$(SRC:.c=.o)
OUT=registration

INCFLAGS=-I$(CURDIR) -I$(PAPIDIR)
CFLAGS=-g -pg
LDFLAGS=-L$(PAPIDIR) -lm -lpapi -lutil_papi

default: $(OUT)

$(OUT): $(OBJ)
	$(CC) $(CFLAGS) $(INCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: $(SRC)
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(OUT)
