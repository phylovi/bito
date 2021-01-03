import enscons
import enscons.cpyext
import glob
import os
import platform
import pybind11
import pytoml as toml
import re
import sys


def check_file_exists(path):
    if not os.path.isfile(path):
        sys.exit("Couldn't find " + path)


if "CONDA_PREFIX" not in os.environ:
    sys.exit(
        "\nThis SConstruct is meant to be run in the libsbn conda environment; "
        "see README for installation process."
    )


conda_env_dir = os.environ["CONDA_PREFIX"]
conda_base_dir = re.search("(.*)/envs/.*", conda_env_dir).group(1)


def record_beagle_prefix(beagle_prefix, ld_variable_name):
    """
    Record the value of the BEAGLE_PREFIX environment variable so that it
    gets set upon `conda activate`.
    """
    conda_vars_dir = os.path.join(conda_env_dir, "etc/conda")
    conda_activate_dir = os.path.join(conda_vars_dir, "activate.d")
    conda_deactivate_dir = os.path.join(conda_vars_dir, "deactivate.d")
    os.makedirs(conda_activate_dir, exist_ok=True)
    os.makedirs(conda_deactivate_dir, exist_ok=True)
    with open(os.path.join(conda_activate_dir, "vars.sh"), "w") as fp:
        fp.write("export BEAGLE_PREFIX=" + beagle_prefix + "\n")
        fp.write(
            f"export {ld_variable_name}=" + os.path.join(beagle_prefix, "lib") + "\n"
        )
    with open(os.path.join(conda_deactivate_dir, "vars.sh"), "w") as fp:
        fp.write(f"unset BEAGLE_PREFIX\n")
        fp.write(f"unset {ld_variable_name}\n")


if "BEAGLE_PREFIX" in os.environ:
    beagle_prefix = os.environ["BEAGLE_PREFIX"]
else:
    print(
        """
It looks like we need to configure your BEAGLE install.

Please enter the prefix directory path where BEAGLE (specifically,
the `hmc-clock` branch) has been installed.
For example, if you supply `/usr/local`, then we should find

    /usr/local/lib/libhmsbeagle.[so|dylib]
    /usr/local/include/libhmsbeagle-1/libhmsbeagle/beagle.h

If your compilation works after this configuration step, run
    conda activate libsbn
which will mean that you don't have to do this configuration again.
Note that you also probably need to do this to run the tests.

If you want to change your prefix directory, do
    rm $CONDA_PREFIX/etc/conda/activate.d/vars.sh
    unset BEAGLE_PREFIX
and you'll get this prompt again when you `scons` or `make`.)
"""
    )

    beagle_prefix = input(">>> ").rstrip()

    check_file_exists(
        os.path.join(beagle_prefix, "include/libhmsbeagle-1/libhmsbeagle/beagle.h")
    )
    if platform.system() == "Darwin":
        record_beagle_prefix(beagle_prefix, "DYLD_LIBRARY_PATH")
    elif platform.system() == "Linux":
        record_beagle_prefix(beagle_prefix, "LD_LIBRARY_PATH")


metadata = dict(toml.load(open("pyproject.toml")))["tool"]["enscons"]
full_tag = enscons.get_abi3_tag()


def python3_get_includes():
    py_includes = set(os.popen("python3-config --includes").read().rstrip().split())
    return list(map(lambda s: s.replace("-I", ""), py_includes))


env = Environment(
    tools=["default", "packaging", enscons.generate, enscons.cpyext.generate],
    PACKAGE_METADATA=metadata,
    WHEEL_TAG=full_tag,
    ENV=os.environ,
    CPPPATH=["src", "lib/eigen", pybind11.get_include()] + python3_get_includes(),
    # CCFLAGS=["-g", "-pthread", "-fno-omit-frame-pointer"],
    CCFLAGS=["-O3", "-pthread", "-fno-omit-frame-pointer"],
    CXXFLAGS=["-std=c++17"],
)

# Sometimes conda installs the pybind11 headers inside a pythonXXX directory, so we get
# them here.
for d in glob.glob(os.path.join(pybind11.get_include(), "python*/")):
    env.Append(CPPPATH=d)

beagle_lib = os.path.join(beagle_prefix, "lib")
beagle_include = os.path.join(beagle_prefix, "include/libhmsbeagle-1")

for path in [beagle_lib, beagle_include]:
    if not os.path.isdir(path):
        sys.exit("Couldn't find:\n" + path)

env.Append(LIBPATH=beagle_lib)
env.Append(CPPPATH=beagle_include)


def perhaps_set_env(key, value):
    if key not in env:
        env[key] = value


if platform.system() == "Darwin":
    perhaps_set_env("CC", "clang")
    perhaps_set_env("CXX", "clang")
    env.Append(DYLD_LIBRARY_PATH=beagle_lib)
    env.Append(LINKFLAGS=["-undefined", "dynamic_lookup"])
    env.Append(CXXFLAGS=["-D_LIBCPP_DISABLE_AVAILABILITY"])
elif platform.system() == "Linux":
    perhaps_set_env("CC", "gcc")
    perhaps_set_env("CXX", "g++")
    env.Append(LD_LIBRARY_PATH=beagle_lib)
else:
    sys.exit("Sorry, we don't support " + platform.system() + ".")


env.VariantDir("_build", "src")
sources = [
    "_build/alignment.cpp",
    "_build/bitset.cpp",
    "_build/block_model.cpp",
    "_build/block_specification.cpp",
    "_build/clock_model.cpp",
    "_build/combinatorics.cpp",
    "_build/csv.cpp",
    "_build/driver.cpp",
    "_build/engine.cpp",
    "_build/fat_beagle.cpp",
    "_build/mersenne_twister.cpp",
    "_build/node.cpp",
    "_build/numerical_utils.cpp",
    "_build/parser.cpp",
    "_build/phylo_model.cpp",
    "_build/psp_indexer.cpp",
    "_build/rooted_gradient_transforms.cpp",
    "_build/rooted_tree.cpp",
    "_build/rooted_tree_collection.cpp",
    "_build/rooted_sbn_instance.cpp",
    "_build/sbn_maps.cpp",
    "_build/sbn_probability.cpp",
    "_build/sbn_support.cpp",
    "_build/scanner.cpp",
    "_build/site_model.cpp",
    "_build/site_pattern.cpp",
    "_build/subsplit_dag.cpp",
    "_build/subsplit_dag_node.cpp",
    "_build/substitution_model.cpp",
    "_build/taxon_name_munging.cpp",
    "_build/tree.cpp",
    "_build/tree_collection.cpp",
    "_build/unrooted_sbn_instance.cpp",
    "_build/unrooted_tree.cpp",
    "_build/unrooted_tree_collection.cpp",
]
gp_sources = [
    "_build/gp_dag.cpp",
    "_build/gp_engine.cpp",
    "_build/gp_instance.cpp",
    "_build/gp_operation.cpp",
]
extension = env.SharedLibrary(
    "libsbn" + os.popen("python3-config --extension-suffix").read().rstrip(),
    ["_build/pylibsbn.cpp"] + sources + gp_sources,
    SHLIBPREFIX="",
    LIBS=["hmsbeagle"],
)
doctest = env.Program(["_build/doctest.cpp"] + sources, LIBS=["hmsbeagle", "pthread"])
noodle = env.Program(["_build/noodle.cpp"] + sources, LIBS=["hmsbeagle", "pthread"])
gp_doctest = env.Program(
    ["_build/gp_doctest.cpp"] + sources + gp_sources, LIBS=["hmsbeagle", "pthread"]
)

py_source = Glob("vip/*.py")

# "platlib" is Python-packaging-speak for a platform-specific library (not pure Python).
platlib = env.Whl("platlib", py_source + extension, root="")
whl = env.WhlFile(source=platlib)
print("\nTo install python wheel, execute:")
print(f"pip install -U {whl[0]}\n")

env.Default(doctest, gp_doctest, whl, noodle)
