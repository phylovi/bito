import glob
import os
import pybind11
import re
import sys
from conda.cli.python_api import Commands, run_command

if 'CC' not in os.environ or 'CXX' not in os.environ:
    sys.exit("Need to have compilers defined by $CC and $CXX shell variables. This is done by conda.")

env = Environment(
    ENV=os.environ,
    CPPPATH=['include', 'src', pybind11.get_include()],
    CCFLAGS=['-O3', '-Wall', '-Wextra', '-Wconversion'],
    CXXFLAGS=['-std=c++17', '-fPIC', '-shared'],
    CC = os.environ['CC'],
    CXX = os.environ['CXX']
    )

for s in run_command(Commands.INFO)[0].split('\n'):
    match = re.search('.*active env location : (.*)', s)
    if match:
        break
if not match:
    sys.exit("Could not find active env location.")
conda_env_dir = match.group(1)
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

sources = [
    '_build/bitset.cpp',
    '_build/driver.cpp',
    '_build/parser.cpp',
    '_build/scanner.cpp'
]

env.VariantDir('_build', 'src')

env.SharedLibrary(
    "sbn"+os.popen("python3-config --extension-suffix").read().rstrip(),
    ['_build/libsbn.cpp'] + sources,
    SHLIBPREFIX='',
    LIBS=['hmsbeagle'])
doctest = env.Program(
    ['_build/doctest.cpp'] + sources,
    LIBS=['hmsbeagle'])
