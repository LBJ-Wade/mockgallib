all: check

TESTS := $(wildcard test_*.py)

test_%:

.PHONY: check

check:
	for t in $(TESTS); do python3 $$t || exit 1; done


