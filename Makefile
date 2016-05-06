.PHONY: all libs python check clean

all:
	make libs
	make python
	make check

libs:
	cd libs && $(MAKE)

python:
	cd python && $(MAKE)

check:
	cd tests && $(MAKE)

clean:
	cd libs && $(MAKE) clean
	cd python && $(MAKE) clean
