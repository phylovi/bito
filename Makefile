default:
	make bison && scons
	./_build/newick_parser data/four_taxon.tre | paste data/four_taxon.tre -
	./_build/doctest
	./_build/test/test

bison: src/parser.yy src/scanner.ll
	bison -o src/parser.cpp --defines=src/parser.hpp src/parser.yy
	flex -o src/scanner.cpp src/scanner.ll

format:
	clang-format -i -style=file src/* test/test.c


prep:
	python test/prep/doctest.py

clean:
	rm -rf _build

edit:
	vim src/doctest.cpp src/driver.cpp src/driver.hpp src/libsbn.cpp src/newick_parser.cpp src/parser.yy src/sbn.h src/sbn.hpp src/scanner.ll
