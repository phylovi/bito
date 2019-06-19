import os

env = Environment(ENV=os.environ, CPPPATH=['include', 'src']) #, CPPFLAGS=['-Wall', '-Wextra', '-Werror', '-Wconversion'])

env.VariantDir('_build', 'src')
env.Library('sbn', ['_build/libsbn.cpp', '_build/sbn.cpp'])

env.Command(['_build/parser.cpp','_build/parser.h'], '_build/parser.yy', 'which bison; bison -o _build/parser.cpp --defines=_build/parser.h $SOURCE')
env.Command(['_build/lexer.cpp', '_build/stack.hh'], '_build/lexer.ll', 'flex -o _build/lexer.cpp --header-file=_build/stack.hh $SOURCE')

env.Program('_build/lexer', ['_build/lexer.cpp', '_build/parser.cpp'], LIBS=['fl'])

# env.VariantDir('_build/test', 'test')
# doctest = env.Program(['_build/test/doctest.cpp', '_build/sbn.cpp'])
# test = env.Program(['_build/test/test.c', '_build/test/munit.c'], LIBS=['sbn'], LIBPATH='.', LINK='g++')

