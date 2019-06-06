env = Environment(CPPPATH="include")

env.VariantDir('_build', 'src')
env.Library('sbn', '_build/libsbn.c')

env.VariantDir('_build/test', 'test')
env.Program('_build/test/test.cpp', LIBS=['sbn'], LIBPATH='.')

env.Program('_build/main.c', LIBS=['sbn'], LIBPATH='.')
