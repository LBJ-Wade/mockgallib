
.PHONY: all libs check clean

libs:
	cd libs && $(MAKE)

pypackage:
	cd python && $(MAKE)

check:
	cd tests && sh test.sh

all:
	make libs
	make python
	make check


clean:
	cd libs && $(MAKE) clean
	cd python && $(MAKE) clean
