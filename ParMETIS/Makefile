OPTIMIZATION = -O3  
INCLUDES = -I. 

# The following is for mpich
CC = mpicc 

# On T3D/T3E you may want to use the following
# CC = mppcc 

# On SP2 you may want to use the following
# CC = mpcc

CFLAGS = $(OPTIMIZATION) $(INCLUDES)
AR = ar rv
RANLIB = ranlib


LIB_OBJS = comm.o util.o myqsort.o memory.o debug.o kmetis.o setup.o \
       coarsen.o edge_refine.o initpart.o serkmetis.o drivers.o fold.o initdiff.o \
       move.o node_refine.o initmsection.o order.o pqueue.o serometis.o mmd.o \
       xyzpart.o
          
.c.o:
	$(CC) $(CFLAGS) -c $*.c

libparmetis.a: $(LIB_OBJS)
	$(AR) $@ $(LIB_OBJS)
	$(RANLIB) $@


