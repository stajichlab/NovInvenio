REPORTS := $(wildcard */report.html)

.PHONY: all clean

all: index.html

index.html: generate_index.py $(REPORTS)
	./generate_index.py --output $@

clean:
	rm -f index.html
