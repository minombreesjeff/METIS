
default:
	more Makefile.help

all:
	$(MAKE) -C libmetis
	$(MAKE) -C programs
	$(MAKE) -C test

clean:
	$(MAKE) -C libmetis clean
	$(MAKE) -C programs clean
	$(MAKE) -C test clean

realclean:
	$(MAKE) -C libmetis realclean
	$(MAKE) -C programs realclean
	$(MAKE) -C test realclean

dist:
	$(MAKE) -C libmetis dist
