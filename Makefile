default:
	make bison && make prep && scons
	./_build/newick_parser data/four_taxon.tre | paste data/four_taxon.tre -
	./_build/doctest
	./_build/test/test

bison: src/parser.yy src/scanner.ll
	bison -o src/parser.cpp --defines=src/parser.hpp src/parser.yy
	flex -o src/scanner.cpp src/scanner.ll

prep:
	python test/prep/doctest.py

format:
	clang-format -i -style=file src/* test/test.c

clean:
	rm -rf _build

edit:
	vim -O2 src/sbn.hpp src/libsbn.cpp src/sbn.h src/driver.cpp src/driver.hpp src/parser.yy src/scanner.ll src/newick_parser.cpp test/prep/doctest.py

lint:
	cpplint src/sbn.hpp src/libsbn.cpp src/sbn.h src/driver.cpp src/driver.hpp src/newick_parser.cpp

.PHONY: bison prep format clean edit lint
