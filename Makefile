default:
	scons && ./_build/test/test

format:
	clang-format -i -style=file src/* test/test.c

