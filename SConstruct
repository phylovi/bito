import glob
import os
import pybind11
import sys
import re
from conda.cli.python_api import run_command, Commands

env = Environment(
    ENV=os.environ,
    CPPPATH=['include', 'src', pybind11.get_include()], # , '/home/ematsen/miniconda3/envs/libsbn/include/libhmsbeagle-1'],
    CCFLAGS=['-O3', '-Wall', '-Wextra', '-Wconversion'],
    CXXFLAGS=['-std=c++14', '-fPIC', '-shared'],
    CC = 'x86_64-conda_cos6-linux-gnu-gcc',
    CXX = 'x86_64-conda_cos6-linux-gnu-g++'
    )
    # LIBPATH='/home/ematsen/miniconda3/pkgs/flex-2.6.4-ha10e3a4_1/lib:/home/ematsen/re/libsbn/_build:/usr/local/lib:/usr/lib:/lib:/usr/local/lib:/usr/lib:/lib')

# Find the active env directory and the conda base directory.
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


# env.Append(LIBPATH = " "+find_conda_pkg_dir_containing("/pkgs/flex*/lib"))

env.VariantDir('_build', 'src')
env.VariantDir('_build/test', 'test')

env.SharedLibrary("sbn"+os.popen("python3-config --extension-suffix").read().rstrip(), ['_build/libsbn.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp', '_build/bitset.cpp'], LINK='g++', SHLIBPREFIX='')
doctest = env.Program(['_build/doctest.cpp', '_build/driver.cpp', '_build/parser.cpp', '_build/scanner.cpp', '_build/bitset.cpp'])
