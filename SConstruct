import os

env = Environment(ENV=os.environ, CPPPATH=['include', 'src']) #, CPPFLAGS=['-Wall', '-Wextra', '-Werror', '-Wconversion'])

env.VariantDir('_build', 'src')
env.Library('sbn', ['_build/libsbn.cpp', '_build/sbn.cpp'])

env.Command(['_build/parser.cpp','_build/parser.hpp', '_build/location.hh'], '_build/parser.yy', 'bison -o _build/parser.cpp --defines=_build/parser.hpp $SOURCE')
env.Command(['_build/scanner.cpp'], '_build/scanner.ll', 'flex -o _build/scanner.cpp $SOURCE')
env.Program('_build/newick_parser', ['_build/newick_parser.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp'], LIBS=['fl'])

# env.VariantDir('_build/test', 'test')
# doctest = env.Program(['_build/test/doctest.cpp', '_build/sbn.cpp'])
# test = env.Program(['_build/test/test.c', '_build/test/munit.c'], LIBS=['sbn'], LIBPATH='.', LINK='g++')

