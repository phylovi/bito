our_files = src/bitset.cpp src/bitset.hpp src/sbn_maps.hpp src/sbn_maps.cpp src/default_dict.hpp src/doctest.cpp src/driver.cpp src/driver.hpp src/intpack.hpp src/libsbn.cpp src/libsbn.hpp src/tree.cpp src/tree.hpp src/node.hpp src/node.cpp src/tree_collection.cpp src/tree_collection.hpp src/alignment.cpp src/alignment.hpp src/task_processor.hpp src/site_pattern.hpp src/site_pattern.cpp src/pylibsbn.cpp src/psp_indexer.cpp src/psp_indexer.hpp src/substitution_model.hpp src/substitution_model.cpp src/site_model.hpp src/site_model.cpp src/clock_model.hpp src/clock_model.cpp src/engine.hpp src/engine.cpp src/fat_beagle.cpp src/fat_beagle.hpp src/phylo_model.cpp src/phylo_model.hpp src/block_model.hpp src/block_model.cpp src/block_specification.cpp src/block_specification.hpp src/beagle_accessories.hpp

default:
	scons
	pip install -U dist/libsbn-*.whl
	./_build/doctest
	pytest
	./_build/noodle

bison: src/parser.yy src/scanner.ll
	bison -o src/parser.cpp --defines=src/parser.hpp src/parser.yy
	flex -o src/scanner.cpp src/scanner.ll

prep:
	python test/prep/doctest.py
	clang-format -i -style=file src/doctest.cpp

format:
	black vip/*py test/*py SConstruct
	docformatter --in-place vip/*py test/*py
	clang-format -i -style=file $(our_files)

clean:
	rm -rf _build

# We follow C++ core guidelines by allowing passing by non-const reference.
lint:
	cpplint --filter=-runtime/references,-build/c++11 $(our_files) \
		&& echo "LINTING PASS"

.PHONY: bison prep format clean edit lint
