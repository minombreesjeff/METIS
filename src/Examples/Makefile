

CC = cc 
LIBS = -lmetis -lm 

INCLUDES = -I. 
OPTIMIZATION = -O2

DEBUGFLAGS = 
EXTRADEFINES =  $(DEBUGFLAGS) 
CFLAGS = $(OPTIMIZATION) $(EXTRADEFINES) $(INCLUDES)
LD = $(CC) -L../ $(LIBDIR)


all: pmetis kmetis ometis

pmetis: pmetis.o
	$(LD) -o $@ pmetis.o $(LIBS)

kmetis: kmetis.o
	$(LD) -o $@ kmetis.o $(LIBS)

ometis: ometis.o
	$(LD) -o $@ ometis.o $(LIBS)



