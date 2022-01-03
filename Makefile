our_files := $(wildcard src/*.cpp) $(wildcard src/*.hpp)
our_files := $(filter-out src/csv.h src/doctest.h src/noodle.cpp src/parser.cpp src/parser.hpp src/prettyprint.hpp src/scanner.cpp, $(our_files))
j_flags = $(shell echo "${MAKEFLAGS}" | grep -o -- "-j[0-9]\+" || true)

default:
	mkdir -p _ignore  # Needed for test outputs.
	scons ${j_flags}
	pip install -U dist/bito-*.whl

cmakebuild:
	@mkdir -p build
	@cd build && \
		cmake -DCMAKE_BUILD_TYPE=Release .. && \
		cmake --build . ${j_flags} && \
		ln -sf ../data . && \
		ln -sf libbito.so bito.so && \
		mkdir -p _ignore

cmakebuildtest:
	@mkdir -p build_test
	@cd build_test && \
		cmake -DCMAKE_BUILD_TYPE=Debug .. && \
		cmake --build . ${j_flags} && \
		ln -sf ../data . && \
		ln -sf libbito.so bito.so && \
		mkdir -p _ignore
	pip install ./build_test

cmakefasttest: cmakebuildtest
	@cd build_test && \
		./doctest --test-case-exclude="* tree sampling" && \
		./gp_doctest --test-case-exclude="UnrootedSBNInstance*" && \
		./noodle
	pytest

rungptest:
	./_build/gp_doctest --test-case-exclude="UnrootedSBNInstance*"

gptest:
	make ${j_flags}
	make rungptest

runfasttest:
	_build/doctest --test-case-exclude="* tree sampling"
	make rungptest
	pytest
	./_build/noodle

fasttest:
	make ${j_flags}
	make runfasttest

test:
	make ${j_flags}
	make runfasttest
	./_build/doctest

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
	rm -rf _build build dist bito.*.so $(find . -name __pycache)

# We follow C++ core guidelines by allowing passing by non-const reference.
lint:
	cpplint --filter=-runtime/references,-build/c++11 $(our_files) \
		&& echo "LINTING PASS"

.PHONY: bison prep format clean edit lint deploy docs fasttest runfasttest runtest test
