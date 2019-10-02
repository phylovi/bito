import enscons, enscons.cpyext
import glob
import os
import platform
import pybind11
import pytoml as toml
import re
import sys

if 'CONDA_PREFIX' not in os.environ:
    sys.exit("\nThis SConstruct is meant to be run in the libsbn conda environment; see README for installation process.")
if 'CC' not in os.environ:
    sys.exit("\nDid you install compilers using conda? See README for installation process.")

metadata = dict(toml.load(open("pyproject.toml")))["tool"]["enscons"]
full_tag = enscons.get_abi3_tag()

env = Environment(
    tools=["default", "packaging", enscons.generate, enscons.cpyext.generate],
    PACKAGE_METADATA=metadata,
    WHEEL_TAG=full_tag,
    ENV=os.environ,
    CPPPATH=['include', 'src', pybind11.get_include()],
    # CCFLAGS=['-g', '-Wall', '-Wextra', '-Wconversion', '-pthread'],
    CCFLAGS=['-O3', '-pthread'],
    CXXFLAGS=['-std=c++14'],
    CC = os.environ['CC'],
    CXX = os.environ['CXX']
    )

# Sometimes conda installs the pybind11 headers inside a pythonXXX directory, so we get
# them here.
for d in glob.glob(pybind11.get_include() + "/python*/"):
    env.Append(CPPPATH = d)

conda_env_dir = env['ENV']['CONDA_PREFIX']
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

def set_library_path(ld_variable_name):
    conda_vars_dir =conda_env_dir+'/etc/conda'
    os.makedirs(conda_vars_dir+'/activate.d', exist_ok=True)
    with open(conda_vars_dir+'/activate.d/vars.sh', 'w') as fp:
        fp.write(f'export {ld_variable_name}='+beagle_lib+'\n')
    os.makedirs(conda_vars_dir+'/deactivate.d', exist_ok=True)
    with open(conda_vars_dir+'/deactivate.d/vars.sh', 'w') as fp:
        fp.write(f'unset {ld_variable_name}\n')

if platform.system() == 'Darwin':
    set_library_path('DYLD_LIBRARY_PATH')
    env.Append(LINKFLAGS = ['-undefined', 'dynamic_lookup'])
elif platform.system() == 'Linux':
    set_library_path('LD_LIBRARY_PATH')
else:
    sys.exit("Sorry, we don't support "+platform.system()+".")

env.VariantDir('_build', 'src')
sources = [
    '_build/alignment.cpp',
    '_build/beagle.cpp',
    '_build/bitset.cpp',
    '_build/driver.cpp',
    '_build/libsbn.cpp',
    '_build/node.cpp',
    '_build/parser.cpp',
    '_build/psp_indexer.cpp',
    '_build/sbn_maps.cpp',
    '_build/scanner.cpp',
    '_build/site_pattern.cpp',
    '_build/tree.cpp',
    '_build/tree_collection.cpp'
]
extension = env.SharedLibrary(
    "libsbn"+os.popen("python3-config --extension-suffix").read().rstrip(),
    ['_build/pylibsbn.cpp'] + sources,
    SHLIBPREFIX='',
    LIBS=['hmsbeagle'])
doctest = env.Program(
    ['_build/doctest.cpp'] + sources,
    LIBS=['hmsbeagle', 'pthread'])
noodle = env.Program(
    ['_build/noodle.cpp'] + sources,
    LIBS=['hmsbeagle', 'pthread'])

py_source = Glob("vip/*.py")

# "platlib" is Python-packaging-speak for a platform-specific library (not pure Python).
platlib = env.Whl("platlib", py_source + extension, root="")
whl = env.WhlFile(source=platlib)
print("\nTo install python wheel, execute:")
print(f"pip install -U {whl[0]}\n")

env.Default(doctest, whl, noodle)
