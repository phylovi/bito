env = Environment(CPPPATH="include")

env.VariantDir('_build', 'src')
env.Library('sbn', '_build/libsbn.c')

env.VariantDir('_build/test', 'test')
test = env.Program(['_build/test/test.c', '_build/test/munit.c'], LIBS=['sbn'], LIBPATH='.')

env.Program('_build/main.c', LIBS=['sbn'], LIBPATH='.')
