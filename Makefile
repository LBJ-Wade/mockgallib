.PHONY: all libs python check clean countm

all:
	make lib
	make py
	make check

libs:
	cd lib && $(MAKE)

python:
	cd py && $(MAKE)

check:
	cd test && $(MAKE)

clean:
	cd lib && $(MAKE) clean
	cd py && $(MAKE) clean

count:
	wc lib/*.cpp py/*.cpp > log/`date +%Y-%m-%d`.txt
