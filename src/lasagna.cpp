#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include "eigen_sugar.hpp"

// Most complete example: https://gist.github.com/marcetcheverry/991042
// https://www.youtube.com/watch?v=m7E9piHcfr4
// https://jameshfisher.com/2017/01/28/mmap-file-write/
// https://stackoverflow.com/a/43301909/467327

int main(void) {
  std::string fname("_ignore/test_file");
  int fd = open(fname.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
  if (fd == -1) {
    std::cout << "could not create file\n";
    return 1;
  }
  size_t desired_vector_length = 5;
  auto v_size = desired_vector_length * sizeof(double);
  ftruncate(fd, v_size);
  double *mmapped_file = (double *)mmap(  //
      NULL,    // This address is ignored because we are using MAP_SHARED.
      v_size,  // Size of map.
      PROT_READ | PROT_WRITE,  // We want to read and write
      MAP_SHARED,              // We need MAP_SHARED to actually write to this.
      fd,                      // File descriptor.
      0                        // Offset.
  );
  Eigen::Map<EigenVectorXd> v(mmapped_file, desired_vector_length);
  std::cout << v << std::endl;
  v(1) += 1.;
  msync(mmapped_file, v_size, MS_SYNC);
  munmap(mmapped_file, v_size);
  close(fd);
  return 0;
}
