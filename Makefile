.PHONY: all lib py test check clean count

all:
	make lib
	make py
	make check

lib:
	cd lib && $(MAKE)

py:
	cd py && $(MAKE)

test:
	cd test && $(MAKE)

check:
	cd py && $(MAKE) check

clean:
	cd lib && $(MAKE) clean
	cd py && $(MAKE) clean

count:
	wc lib/*.cpp py/*.cpp > log/`date +%Y-%m-%d`.txt
