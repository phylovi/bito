our_files = src/beagle.cpp src/beagle.hpp src/bitset.cpp src/bitset.hpp src/sbn_maps.hpp src/sbn_maps.cpp src/default_dict.hpp src/doctest.cpp src/driver.cpp src/driver.hpp src/intpack.hpp src/libsbn.cpp src/libsbn.hpp src/tree.cpp src/tree.hpp src/node.hpp src/node.cpp src/tree_collection.cpp src/tree_collection.hpp src/alignment.cpp src/alignment.hpp src/task_processor.hpp src/beagle.hpp src/beagle.cpp src/site_pattern.hpp src/site_pattern.cpp src/pylibsbn.cpp src/psp_indexer.cpp src/psp_indexer.hpp

default:
	scons
	./_build/doctest
	pytest -s

bison: src/parser.yy src/scanner.ll
	bison -o src/parser.cpp --defines=src/parser.hpp src/parser.yy
	flex -o src/scanner.cpp src/scanner.ll

prep:
	python test/prep/doctest.py
	clang-format -i -style=file src/doctest.cpp

format:
	clang-format -i -style=file $(our_files)
	yapf -i test_instance.py

clean:
	rm -rf _build

# We follow C++ core guidelines by allowing passing by non-const reference.
lint:
	cpplint --filter=-runtime/references,-build/c++11 $(our_files) \
		&& echo "LINTING PASS"

.PHONY: bison prep format clean edit lint
