import pybind11
import os

env = Environment(
    ENV=os.environ,
    CPPPATH=['include', 'src', pybind11.get_include()], # , '/home/ematsen/miniconda3/envs/libsbn/include/libhmsbeagle-1'],
    CCFLAGS=['-O3', '-Wall', '-Wextra', '-Wconversion'],
    CXXFLAGS=['-std=c++14', '-fPIC', '-shared'])

env.VariantDir('_build', 'src')
env.VariantDir('_build/test', 'test')

#env.Program('_build/newick_parser', ['_build/newick_parser.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp'], LIBS=['fl','sbn'], LIBPATH='_build')
#env.SharedLibrary('_build/sbn', ['_build/libsbn.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp'])
#env.Library('_build/sbn', ['_build/libsbn.cpp', '_build/driver.o', '_build/parser.o', '_build/scanner.o'])

#doctest = env.Program(['_build/doctest.cpp', '_build/driver.o', '_build/parser.o', '_build/scanner.o', '_build/bitset.cpp'])
# test = env.Program(['_build/test/test.c', '_build/test/munit.c'], LIBS=['sbn','hmsbeagle'], LIBPATH=['_build','/home/ematsen/miniconda3/envs/libsbn/include/libhmsbeagle-1'], LINK='g++')

env.SharedLibrary("example"+os.popen("python3-config --extension-suffix").read().rstrip(), ['_build/example.cpp'], LIBS=[], LIBPATH=['_build'], LINK='g++', SHLIBPREFIX='')

# c++  `python3 -m pybind11 --includes` example.cpp -o example`python3-config --extension-suffix` -I../src -lsbn -L ../_build && python -c "import example; print(example.add(1,2))"

