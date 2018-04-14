CC = gcc
RM = rm -rfv
ARCHIVE = ar rscv

ifdef DEBUG
  OPTOPTS = -g -DDEBUG
else
  OPTOPTS = -mtune=native -march=native -O3
endif

ifdef USE_GKLIB
  GKOPT = -DUSE_GKLIB -I./GKlib
else
  GKOPT =
endif

ifdef NOMATH 
  MOPT = -DNOMATH 
else
  MOPT =
endif

LIBS = 
CWARN = -Wall -Wmissing-prototypes -Wmissing-declarations \
        -Wstrict-aliasing -Wextra -Werror -Wno-unused-parameter \
        -Wno-unused-function
CSTD = -std=gnu99
CSPEC = -fstrict-aliasing -static

COPTS = $(OPTOPTS) $(CSTD) $(CWARN) $(GKOPT) $(CSPEC) $(MOPT)

HEADERS = domlib.h \
          dlfile.h \
          dlutil.h \
          dlsort.h \
          dlprint.h \
          dlomp.h \
          dldefs.h \
          dlenv.h \
          dlstring.h \
          dlterm.h \
          dlcmdline.h \
          dlbuffer_headers.h \
          dlbuffer_funcs.h \
          dlmem_headers.h \
          dlmem_funcs.h \
          dldpq_headers.h \
          dldpq_funcs.h \
          dlht_headers.h \
          dlht_funcs.h \
          dlcb_headers.h \
          dlcb_funcs.h \
          dliset_headers.h \
          dliset_funcs.h \
          dldjset_headers.h \
          dldjset_funcs.h \
          dlmath_headers.h \
          dlmath_funcs.h \
          dlsort_headers.h \
          dlsort_funcs.h \
          dllist_headers.h \
          dllist_funcs.h \
          dlheap_headers.h \
          dlheap_funcs.h \
          dltree_headers.h \
          dltree_funcs.h

BASEOBJS = dlfile.o dlutil.o dlenv.o dlstring.o dlcmdline.o dlterm.o

all: libdom.a

test: libdom.a
	make -C test test

# Libraries ###################################################################

libdom.a: domlib.o $(BASEOBJS) 
	$(RM) $@
	$(ARCHIVE) $@ $^

# BASE OBJECTS ################################################################

%.o: %.c $(HEADERS)
	$(CC) $(COPTS) -c $< $(LIBS) -o $@

clean:
	$(RM) *.o && \
	$(RM) *.swp 

distclean: clean
	$(RM) libdom.a

rebuild: distclean all 

.PHONY: test
