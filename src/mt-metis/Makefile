RM = rm -rvf

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)
METISBUILDDIR = build/$(systype)-$(cputype)

all: metis/$(METISBUILDDIR)/libmetis.a
	$(MAKE) $(MFLAGS) -C src

metis/$(METISBUILDDIR)/libmetis.a:
	$(MAKE) $(MFLAGS) -C metis config
	$(MAKE) $(MFLAGS) -C metis

pkg: distclean
	./mkpkg.sh

clean:
	$(MAKE) $(MFLAGS) -C src clean
	$(MAKE) $(MFLAGS) -C metis clean

distclean:
	$(MAKE) $(MFLAGS) -C src distclean
	$(RM) mt-metis-*.tar.gz
