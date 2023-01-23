CC=gcc
CFLAGS=-I.

build: lll_algorithm.c
	cc lll_algorithm.c -o lll_algorithm -lm

run: build
	./lll_algorithm

.PHONY: clean

clean:
	rm -f lll_algorithm