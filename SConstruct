env = Environment(CPPPATH=['include', 'src'])

env.VariantDir('_build', 'src')
env.Library('sbn', ['_build/libsbn.cpp', '_build/sbn.cpp'])

env.VariantDir('_build/test', 'test')
test = env.Program(['_build/test/test.c', '_build/test/munit.c'], LIBS=['sbn'], LIBPATH='.', LINK='g++')
