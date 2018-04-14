CC = /usr/bin/gcc
RM = rm -rvf
ARCHIVE = ar rscv

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

METISROOT="../metis"
METISINC="$(METISROOT)/include"
METISSOURCE="$(METISROOT)/libmetis"

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)
METISBUILDDIR = build/$(systype)-$(cputype)

METISBIN="$(METISROOT)/$(METISBUILDDIR)"
MTMETISINCLUDE="../include"

GKLIBSRC="$(METISROOT)/GKlib"

INCS = -I"$(METISSOURCE)" -I"$(METISINC)" -I"$(GKLIBSRC)" -I"$(MTMETISINCLUDE)"
LIBS = -L"$(METISBIN)/libmetis/" -lmetis -lm 

HEADERS = includes.h \
          proto.h \
          defs.h \
          macros.h \
          struct.h

SOURCE = io.c \
         kpart.c \
         coarsen.c \
         kwayrefine.c \
         kwayfm.c \
         initpart.c \
         mem.c \
         graph.c \
         check.c \
         util.c

OFILES = io.o \
         kpart.o \
         coarsen.o \
         kwayrefine.o \
         kwayfm.o \
         initpart.o \
         mem.o \
         graph.o \
         check.o \
         mtmetis.o \
         util.o

EXE = ../mt-metis
LIB = ../libmtmetis.a

all: $(EXE) $(LIB)

$(EXE): animikii.c $(OFILES) $(HEADERS) 
	$(CC) $(COPTS) $(INCS) animikii.c $(OFILES) $(LIBS) -fopenmp -o $@

$(LIB): $(OFILES)
	$(ARCHIVE) $@ $^
 
%.o: %.c $(HEADERS)
	$(CC) $(COPTS) $(INCS) -c $< $(LIBS) -fopenmp -o $@

clean: .PHONY
	$(RM) *.o

distclean: clean
	$(RM) $(EXE) $(LIB)


.PHONY:
