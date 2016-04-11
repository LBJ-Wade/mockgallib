
.PHONY: all

all:
	cd libs && $(MAKE)
	cd python && $(MAKE)
	cd tests && sh test.sh

clean:
	cd libs && $(MAKE) clean
	cd python && $(MAKE) clean
