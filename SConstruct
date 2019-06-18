env = Environment(CPPPATH=['include', 'src'], CPPFLAGS=['-Wall', '-Wextra', '-Werror', '-Wconversion'])

env.VariantDir('_build', 'src')
env.Library('sbn', ['_build/libsbn.cpp', '_build/sbn.cpp'])

env.VariantDir('_build/test', 'test')
doctest = env.Program(['_build/test/doctest.cpp', '_build/sbn.cpp'])
# test = env.Program(['_build/test/test.c', '_build/test/munit.c'], LIBS=['sbn'], LIBPATH='.', LINK='g++')

