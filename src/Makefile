

BINDIR = $(PWD)
CC = cc 
LIBS = -lm 

INCLUDES = -I. 
OPTIMIZATION = -O2

DEBUGFLAGS = 
EXTRADEFINES =  $(DEBUGFLAGS) 
CFLAGS = $(OPTIMIZATION) $(EXTRADEFINES) $(INCLUDES)
LD = $(CC) -L. $(LIBDIR)


MLEVEL_OBJS = main.o io.o util.o match.o coarsen.o memory.o fm.o  \
              GKlib.o htable.o list.o refine.o bucketlist.o stat.o \
              mlevelpart.o separator.o mlevelorder.o smbfactor.o \
              smbfct.o mincover.o initpart.o lanczos.o mmd.o balpart.o \
              kwaypart.o kwayfm.o kwaygreedy.o 

default:	
	touch *.c
	make $(BINDIR)/metis

$(BINDIR)/metis: $(MLEVEL_OBJS)
	$(LD) -o $@ $(MLEVEL_OBJS) $(LIBS)
	rm -f $(BINDIR)/pmetis $(BINDIR)/ometis $(BINDIR)/kmetis
	ln -s $@ $(BINDIR)/pmetis
	ln -s $@ $(BINDIR)/ometis
	ln -s $@ $(BINDIR)/kmetis

lib:	
	touch *.c
	make -f Makefile.lib

all: 
	make 
	make lib

clean:
	rm *.o

realclean:
	rm -f *.o $(BINDIR)/metis $(BINDIR)/pmetis $(BINDIR)/ometis $(BINDIR)/kmetis 
	make -f Makefile.lib realclean

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
