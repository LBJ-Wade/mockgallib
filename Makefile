.PHONY: all lib py test check clean count

WITHMPI := 1
export WITHMPI

ifdef WITHMPI
CC      := mpic++
CXX     := mpic++
endif

all:
	make lib
	make py
	make check

lib:
	cd lib && $(MAKE) lib

py:
	cd py && $(MAKE) py

test:
	cd test && $(MAKE) test

check:
	cd py && $(MAKE) check

clean:
	cd lib && $(MAKE) clean
	cd py && $(MAKE) clean

count:
	wc lib/*.cpp py/*.cpp > log/`date +%Y-%m-%d`.txt
