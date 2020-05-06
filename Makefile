our_files := $(wildcard src/*.cpp) $(wildcard src/*.hpp)
our_files := $(filter-out src/doctest.h src/noodle.cpp src/parser.cpp src/parser.hpp src/prettyprint.hpp src/scanner.cpp, $(our_files))

default:
	scons
	pip install -U dist/libsbn-*.whl

test:
	make
	./_build/doctest
	./_build/gp_doctest
	pytest
	./_build/noodle

bison: src/parser.yy src/scanner.ll
	bison -o src/parser.cpp --defines=src/parser.hpp src/parser.yy
	flex -o src/scanner.cpp src/scanner.ll

prep:
	python test/prep/doctest.py
	clang-format -i -style=file src/doctest.cpp

docs:
	make -C doc clean
	PYTHONPATH=. sphinx-autogen doc/index.rst
	make -C doc html

deploy:
	make
	make docs
	git checkout gh-pages
	cp -a doc/_build/html/* .
	git add .
	git commit --amend -av -m "update docs"
	git push -f

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

.PHONY: bison prep format clean edit lint deploy docs test
