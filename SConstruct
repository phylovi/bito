import glob
import os
import platform
import pybind11
import re
import sys

if 'CC' not in os.environ or 'CXX' not in os.environ:
    sys.exit("Need to have compilers defined by $CC and $CXX shell variables. This is done by conda.")

env = Environment(
    ENV=os.environ,
    CPPPATH=['include', 'src', pybind11.get_include()],
    CCFLAGS=['-O3', '-Wall', '-Wextra', '-Wconversion'],
    CXXFLAGS=['-std=c++14', '-fPIC', '-shared']
    )

conda_env_dir = env['ENV']['CONDA_PREFIX']
if not conda_env_dir:
    sys.exit("$CONDA_PREFIX is not set. This SConstruct is meant to be run in the libsbn conda environment.")
conda_base_dir = re.search('(.*)/envs/.*', conda_env_dir).group(1)

def find_conda_pkg_dir_containing(glob_str):
   for possible_location in [conda_base_dir, conda_env_dir]:
       g = glob.glob(possible_location+glob_str)
       if g:
           return g[0]
   sys.exit("I can't find the package directory containing "+glob_str)

beagle_pkg = find_conda_pkg_dir_containing("/pkgs/beagle-lib*/")
beagle_lib = beagle_pkg + 'lib'
beagle_include = beagle_pkg + 'include/libhmsbeagle-1'

for path in [beagle_lib, beagle_include]:
    if not os.path.isdir(path):
        sys.exit("Couldn't find:\n"+path)

env.Append(LIBPATH = beagle_lib)
env.Append(CPPPATH = beagle_include)

if platform.system() == 'Darwin':
    if 'DYLD_LIBRARY_PATH' in os.environ:
        os.environ['DYLD_LIBRARY_PATH'] += ':'+beagle_lib
    else:
        os.environ['DYLD_LIBRARY_PATH'] = beagle_lib
    env.Append(LINKFLAGS = ['-undefined', 'dynamic_lookup'])
elif platform.system() == 'Linux':
    env['CC'] = os.environ['CC']
    env['CXX'] = os.environ['CXX']
else:
    sys.exit("Sorry, we don't support "+platform.system()+".")

env.VariantDir('_build', 'src')
sources = [
    '_build/bitset.cpp',
    '_build/driver.cpp',
    '_build/parser.cpp',
    '_build/scanner.cpp'
]
env.SharedLibrary(
    "sbn"+os.popen("python3-config --extension-suffix").read().rstrip(),
    ['_build/libsbn.cpp'] + sources,
    SHLIBPREFIX='',
    LIBS=['hmsbeagle'])
doctest = env.Program(
    ['_build/doctest.cpp'] + sources,
    LIBS=['hmsbeagle'])
