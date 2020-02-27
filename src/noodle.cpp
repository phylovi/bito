#include "libsbn.hpp"
#include "rooted_sbn_instance.hpp"

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

int main() {
  RootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/hello.nwk");
  inst.UseCurrentTreesAsRooted();
  inst.ReadFastaFile("data/hello.fasta");
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  inst.PrepareForPhyloLikelihood(simple_specification, 2);
  std::cout << inst.LogLikelihoods() << std::endl;
}
