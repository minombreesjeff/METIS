CC = /usr/bin/gcc
RM = rm -rvf

DEBUGOPTS = -g3 -O0 -DDEBUG

ifdef ASSERT
    COPTS = $(DEBUGOPTS)
else
  ifdef DEBUG
      COPTS = -rdynamic $(DEBUGOPTS) -DNDEBUG 
  else
      COPTS = -O3 -DNDEBUG
  endif
endif

METISROOT="metis"
METISINC="$(METISROOT)/include"
METISSOURCE="$(METISROOT)/libmetis"
METISBIN="$(METISROOT)/build/Linux-x86_64"

GKLIBSRC="$(METISROOT)/GKlib"

JUNK = -L"$(METISBIN)/libmetis/" -lmetis -I"$(METISSOURCE)" -I"$(METISINC)" \
    -I"$(GKLIBSRC)" -lm -fopenmp

HEADERS = includes.h proto.h defs.h macros.h struct.h
SOURCE = io.c kpart.c coarsen.c kwayrefine.c kwayfm.c initpart.c mem.c \
  graph.c check.c util.c

EXE = mt-metis

all: $(EXE)

$(EXE): animikii.c $(SOURCE) $(HEADERS) 
	$(CC) $(COPTS) animikii.c $(SOURCE) $(JUNK) -o $@

clean: .PHONY
	$(RM) animikii openmp_test *.o

.PHONY:
