include ../Makefile.in

INCLUDES = -I./ 

CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)
LD = $(CC) -L./ 


OBJS = coarsen.o fm.o initpart.o match.o ccgraph.o memory.o \
       pmetis.o pqueue.o refine.o util.o timing.o debug.o \
       bucketsort.o graph.o stat.o kmetis.o kwayrefine.o \
       kwayfm.o balance.o ometis.o srefine.o sfm.o separator.o \
       mincover.o mmd.o mesh.o meshpart.o frename.o fortran.o \
       myqsort.o compress.o ometispar.o estmem.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c

../libmetis.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

checkin:
	@for file in *.[c,h]; \
	do \
	ci -u -m'Maintance' $$file;\
	done 

checkin2:
	@for file in *.[c,h]; \
	do \
	ci $$file;\
	rcs -U $$file;\
	co $$file;\
	done

clean:
	rm -f *.o

realclean:
	rm -f *.o ; rm -f ../libmetis.a
