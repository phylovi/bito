#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include "eigen_sugar.hpp"

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
  lseek(fd, desired_vector_length - 1, SEEK_SET);
  write(fd, "", 1);
  double *mmapped_file = (double *)mmap(  //
      NULL,    // This address is ignored because we are using MAP_SHARED.
      v_size,  // Size of map.
      PROT_READ | PROT_WRITE,  // We want to read and write
      MAP_SHARED,              // We need MAP_SHARED to actually write to this.
      fd,                      // File descriptor.
      0                        // Offset.
  );
  Eigen::Map<EigenVectorXd> v(mmapped_file, desired_vector_length);
  v(0) = 5.;

  msync(mmapped_file, v_size, MS_SYNC);
  munmap(mmapped_file, v_size);
  close(fd);

  // reading file back

  fd = open(fname.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
  struct stat sb;

  if (fstat(fd, &sb) == -1) {
    perror("couldn't get file size\n");
  }

  printf("file size is %ld\n", sb.st_size);

  double *file_in_memory = (double *)mmap(NULL, v_size, PROT_READ, MAP_SHARED, fd, 0);

  Eigen::Map<EigenVectorXd> reconstituted_v(file_in_memory, desired_vector_length);
  std::cout << reconstituted_v << std::endl;
  // printf("%g", file_in_memory[0]);
  printf("\n");
  munmap(file_in_memory, v_size);
  close(fd);

  return 0;
}

int string_version(void) {
  int fd = open("_ignore/test_file", O_RDWR | O_CREAT, (mode_t)0600);
  const char *text = "hello";
  size_t textsize = strlen(text) + 1;
  // This writes a string termination character.
  lseek(fd, textsize - 1, SEEK_SET);
  write(fd, "", 1);
  char *map = (char *)mmap(  //
      NULL,                  // This address is ignored because we are using MAP_SHARED.
      textsize,              // Size of map.
      PROT_READ | PROT_WRITE,  // We want to read and write
      MAP_SHARED,              // We need MAP_SHARED to actually write to this.
      fd,                      // File descriptor.
      0                        // Offset.
  );
  memcpy(map, text, strlen(text));
  msync(map, textsize, MS_SYNC);
  for (int i = 0; i < strlen(text); i++) {
    printf("%c", map[i]);
  }
  munmap(map, textsize);
  close(fd);
  return 0;
}

