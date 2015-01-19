
CFLAGS = -D_GNU_SOURCE=1 -W -Wall -O3 -std=c99
LDLIBS = -lm

fir: fir.o window.o

clean:
	rm -f fir *.o

