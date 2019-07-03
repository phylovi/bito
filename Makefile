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
	cpplint src/bitset.cpp src/bitset.hpp src/build.hpp src/default_dict.hpp src/doctest.cpp src/driver.cpp src/driver.hpp src/intpack.hpp src/libsbn.cpp src/libsbn.hpp src/tree.hpp

.PHONY: bison prep format clean edit lint
