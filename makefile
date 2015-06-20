# CC = gcc-4.4
CC = gcc
CFLAGS = -Wall -g -O3 
LDFLAGS = -lm -g  -O3
all: genso

clean:
	rm genso.o genso

genso: genso.o
	$(CC) $(CFLAGS) genso.o -o genso $(LDFLAGS)

genso.o: genso.c soumod.c readidem.c readidem2.c readidem3.c velocity.c rupfront.c disazi.c risetime.c
	$(CC) $(CFLAGS) -c genso.c 
