include ../Makefile.in

BINDIR = ./

INCLUDES = -I../Lib 
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)

LIBS = -lmetis -lm
LD = $(CC) $(LDOPTIONS) -L../ 
METISLIB = ../libmetis.a


NDCOMPAREOBJ 	= ndcompare.o ../Programs/io.o ../Programs/smbfactor.o
PWTESTOBJ 	= pwtest.o ../Programs/io.o 
KWTESTOBJ 	= kwtest.o ../Programs/io.o 
MESHTESTOBJ 	= meshtest.o 

.c.o:
	$(CC) $(CFLAGS) -c $*.c

default: $(BINDIR)/ndcompare$(VERNUM) $(BINDIR)/pwtest$(VERNUM) $(BINDIR)/kwtest$(VERNUM)


$(BINDIR)/pwtest$(VERNUM): $(PWTESTOBJ) $(METISLIB)
	$(LD) -o $@ $(PWTESTOBJ) $(LIBS)
	chmod 755 $@

$(BINDIR)/kwtest$(VERNUM): $(KWTESTOBJ) $(METISLIB)
	$(LD) -o $@ $(KWTESTOBJ) $(LIBS)
	chmod 755 $@

$(BINDIR)/ndcompare$(VERNUM): $(NDCOMPAREOBJ) $(METISLIB)
	$(LD) -o $@ $(NDCOMPAREOBJ) $(LIBS)
	chmod 755 $@

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
	rm -f *.o ;\
	rm -f $(BINDIR)/ndcompare$(VERNUM) ;\
	rm -f $(BINDIR)/pwtest$(VERNUM) ;\
	rm -f $(BINDIR)/kwtest$(VERNUM)

