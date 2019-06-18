env = Environment(CPPPATH=['include', 'src']) #, CPPFLAGS=['-Wall', '-Wextra', '-Werror', '-Wconversion'])

env.VariantDir('_build', 'src')
env.Library('sbn', ['_build/libsbn.cpp', '_build/sbn.cpp'])

env.Command(['_build/parser.c','_build/parser.h'], '_build/parser.y', 'bison -o _build/parser.c --defines=_build/parser.h $SOURCE')
env.CFile('_build/lexer.c', '_build/lexer.l')
env.Program('_build/lexer', ['_build/lexer.c', '_build/parser.c'], LIBS=['fl'])

env.VariantDir('_build/test', 'test')
doctest = env.Program(['_build/test/doctest.cpp', '_build/sbn.cpp'])
# test = env.Program(['_build/test/test.c', '_build/test/munit.c'], LIBS=['sbn'], LIBPATH='.', LINK='g++')

