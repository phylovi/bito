default:
	scons && ./_build/lexer
	#scons && ./_build/test/doctest

format:
	clang-format -i -style=file src/* test/test.c

