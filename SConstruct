import os
import pybind11

env = Environment(
    ENV=os.environ,
    CPPPATH=['include', 'src', pybind11.get_include()], # , '/home/ematsen/miniconda3/envs/libsbn/include/libhmsbeagle-1'],
    CCFLAGS=['-O3', '-Wall', '-Wextra', '-Wconversion'],
    CXXFLAGS=['-std=c++14', '-fPIC', '-shared'],
    # CC = os.environ['CC'],
    # CXX = os.environ['CXX']
    CC = 'x86_64-conda_cos6-linux-gnu-gcc',
    CXX = 'x86_64-conda_cos6-linux-gnu-g++'
    )

env.VariantDir('_build', 'src')

env.SharedLibrary("sbn"+os.popen("python3-config --extension-suffix").read().rstrip(), ['_build/libsbn.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp', '_build/bitset.cpp'], SHLIBPREFIX='')
doctest = env.Program(['_build/doctest.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp', '_build/bitset.cpp'])
