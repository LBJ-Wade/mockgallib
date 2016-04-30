SRC  := $(wildcard *.md)
HTML := $(addsuffix .html, $(basename $(SRC)))

all: $(HTML)

$(HTML): Makefile doc.css syntax.css

%.html: %.md
	@echo $<
	pandoc --mathjax -c doc.css -c syntax.css $< -o $@

.PHONY: clean

clean:
	rm -f $(HTML)
