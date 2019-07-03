import os
import pybind11

if 'CC' not in os.environ or 'CXX' not in os.environ:
    sys.exit("Need to have compilers defined by $CC and $CXX shell variables. This is done by conda.")

env = Environment(
    ENV=os.environ,
    CPPPATH=['include', 'src', pybind11.get_include()], # , '/home/ematsen/miniconda3/envs/libsbn/include/libhmsbeagle-1'],
    CCFLAGS=['-O3', '-Wall', '-Wextra', '-Wconversion'],
    CXXFLAGS=['-std=c++14', '-fPIC', '-shared'],
    CC = os.environ['CC'],
    CXX = os.environ['CXX']
    )

sources = [
    '_build/libsbn.cpp',
    '_build/bitset.cpp'
    '_build/driver.cpp',
    '_build/parser.cpp',
    '_build/scanner.cpp',
]

env.VariantDir('_build', 'src')

env.SharedLibrary("sbn"+os.popen("python3-config --extension-suffix").read().rstrip(), sources, SHLIBPREFIX='')
doctest = env.Program(['_build/doctest.cpp'+sources)
