    CC     = gcc -g -std=c11
    CFLAGS = -O3 -mavx2 -march=native -DLIKWID_PERFMON -I/home/soft/likwid/include
    LFLAGS = -L/home/soft/likwid/lib -lm -llikwid

      PROG = InterpolationLinear 
      OBJS = utils.o \
             fatorLU.o

.PHONY: limpa faxina clean purge all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(LFLAGS)

limpa clean:
	@rm -f *~ *.bak

faxina purge:   limpa
	@rm -f *.o core a.out
	@rm -f $(PROG)

