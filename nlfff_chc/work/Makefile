OBJS = relax1.o init.o loop.o bfield.o optimization.o

CC = gcc
LDFLAGS = -lm 

relax4: $(OBJS)
	 $(CC)  $(LDFLAGS) $(OBJS) -o $@
relax1.o: relax1.c
init.o: init.c globals.h
loop.o: loop.c globals.h
bfield.o: bfield.c globals.h
optimization.o: optimization.c globals.h
