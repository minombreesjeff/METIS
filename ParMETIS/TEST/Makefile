INCLUDES = -I. -I..
LIBDIR = -L../

# The following is for mpich
LIBS = -lm -lparmetis
CC = mpicc 

# On T3D you may want to use the following
# LIBS = -lm -lmpi -lparmetis
# CC = mppcc

# On SP2 you may want to use the following
# LIBS = -lm -lmpi -lparmetis
# CC = mpcc


CFLAGS = $(OPTIMIZATION) $(INCLUDES)
LD = $(CC) $(LIBDIR)


P2P_OBJS = p2ptest.o 
KMETIS_OBJS = main.o io.o
          
.c.o:
	$(CC) $(CFLAGS) -c $*.c

all: p2ptest parmetis.test

p2ptest: $(P2P_OBJS)
	$(LD)  -o $@ $(P2P_OBJS) $(LIBS)
	chmod 744 $@

parmetis.test: $(KMETIS_OBJS)
	$(LD)  -o $@ $(KMETIS_OBJS) $(LIBS) 
	chmod 744 $@

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

