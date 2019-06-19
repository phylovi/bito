default:
	scons && ./_build/newick_parser _ignore/ex.nwk
	#scons && ./_build/test/doctest

format:
	clang-format -i -style=file src/* test/test.c

