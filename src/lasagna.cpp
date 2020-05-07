#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include "eigen_sugar.hpp"

// Most complete example: https://gist.github.com/marcetcheverry/991042
// https://www.youtube.com/watch?v=m7E9piHcfr4
// https://jameshfisher.com/2017/01/28/mmap-file-write/
// https://stackoverflow.com/a/43301909/467327

int main(void) {
  std::string fname("_ignore/test_file");
  int fd = open(
      fname.c_str(),
      O_RDWR | O_CREAT,  // Open for reading and writing; create if it doesn't exit.
      S_IRUSR | S_IWUSR  // Make the file readable and writable by the user.
  );
  if (fd == -1) {
    std::cout << "could not create file\n";
    return 1;
  }
  size_t desired_vector_length = 5;
  auto v_size = desired_vector_length * sizeof(double);
  // Resizes fd so it's just right for our vector.
  ftruncate(fd, v_size);
  double *mmapped_memory = (double *)mmap(  //
      NULL,                    // This address is ignored as we are using MAP_SHARED.
      v_size,                  // Size of map.
      PROT_READ | PROT_WRITE,  // We want to read and write.
      MAP_SHARED,              // We need MAP_SHARED to actually write to memory.
      fd,                      // File descriptor.
      0                        // Offset.
  );
  Eigen::Map<EigenVectorXd> v(mmapped_memory, desired_vector_length);
  std::cout << v << std::endl;
  v(1) += 1.;
  // Synchronize memory with physical storage.
  msync(mmapped_memory, v_size, MS_SYNC);
  // Unmap memory mapped with mmap.
  munmap(mmapped_memory, v_size);
  close(fd);
  return 0;
}
