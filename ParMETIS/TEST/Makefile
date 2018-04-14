include ../Makefile.in

BINDIR = ../Graphs

INCLUDES = -I../ParMETISLib $(INCDIR) 
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)


LIBSDIR = -L.. $(LIBDIR) 
LIBS = -lparmetis -lmetis $(XTRALIBS) -lm
PARMETISLIB = ../libparmetis.a ../libmetis.a


PTESTOBJS = main.o pio.o adaptgraph.o ptest.o
MESHTESTOBJS = mmain.o meshio.o pio.o

          
.c.o:
	$(CC) $(CFLAGS) -c $*.c


default: $(BINDIR)/ptest $(BINDIR)/meshtest 

$(BINDIR)/meshtest: $(MESHTESTOBJS) $(PARMETISLIB)
	$(LD)  -o $@ $(MESHTESTOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@

$(BINDIR)/ptest: $(PTESTOBJS) $(PARMETISLIB)
	$(LD)  -o $@ $(PTESTOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@

clean:
	rm -f *.o ;\
	rm -f $(BINDIR)/ptest$(VERNUM) ;\
	rm -f $(BINDIR)/meshtest$(VERNUM) 
