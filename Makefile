default:
	scons
	./_build/doctest
	pytest -s

bison: src/parser.yy src/scanner.ll
	bison -o src/parser.cpp --defines=src/parser.hpp src/parser.yy
	flex -o src/scanner.cpp src/scanner.ll

prep:
	python test/prep/doctest.py

format:
	clang-format -i -style=file src/*

clean:
	rm -rf _build

edit:
	vim -O2 src/tree.hpp src/libsbn.cpp src/sbn.h src/driver.cpp src/driver.hpp src/parser.yy src/scanner.ll test/prep/doctest.py

lint:
	cpplint src/tree.hpp src/libsbn.cpp src/sbn.h src/driver.cpp src/driver.hpp src/bitset.cpp src/bitset.hpp

.PHONY: bison prep format clean edit lint
