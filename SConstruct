import os

env = Environment(ENV=os.environ, CPPPATH=['include', 'src'], CPPFLAGS=['-Wall', '-Wextra', '-Wconversion'])

env.VariantDir('_build', 'src')
env.VariantDir('_build/test', 'test')

env.Program('_build/newick_parser', ['_build/newick_parser.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp'], LIBS=['fl','sbn'], LIBPATH='.')
env.Library('sbn', ['_build/libsbn.cpp', '_build/driver.o', '_build/parser.o', '_build/scanner.o'])

doctest = env.Program('_build/doctest.cpp')
test = env.Program(['_build/test/test.c', '_build/test/munit.c'], LIBS=['sbn'], LIBPATH='.', LINK='g++')
