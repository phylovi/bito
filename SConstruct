import pybind11
import os

env = Environment(
    ENV=os.environ,
    CPPPATH=['include', 'src', pybind11.get_include()], # , '/home/ematsen/miniconda3/envs/libsbn/include/libhmsbeagle-1'],
    CCFLAGS=['-O3', '-Wall', '-Wextra', '-Wconversion'],
    CXXFLAGS=['-std=c++14', '-fPIC', '-shared'])

env.VariantDir('_build', 'src')
env.VariantDir('_build/test', 'test')

env.Program('_build/newick_parser', ['_build/newick_parser.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp'], LIBS=['fl'], LIBPATH='_build')
env.SharedLibrary("sbn"+os.popen("python3-config --extension-suffix").read().rstrip(), ['_build/libsbn.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp'], LINK='g++', SHLIBPREFIX='')
doctest = env.Program(['_build/doctest.cpp', '_build/driver.o', '_build/parser.o', '_build/scanner.o', '_build/bitset.cpp'])
