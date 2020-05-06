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
  int fd = open(fname.c_str(), O_RDWR | O_CREAT, (mode_t)0600);
  EigenVectorXd v(5);
  v << 100, 2, 3, 4, 5;
  std::cout << v << std::endl;
  auto v_size = v.size() * sizeof(double);
  double *map = (double *)mmap(  //
      NULL,    // This address is ignored because we are using MAP_SHARED.
      v_size,  // Size of map.
      PROT_READ | PROT_WRITE,  // We want to read and write
      MAP_SHARED,              // We need MAP_SHARED to actually write to this.
      fd,                      // File descriptor.
      0                        // Offset.
  );
  memcpy(map, v.data(), v_size);
  map[0] = 3.14;
  msync(map, v_size, MS_SYNC);
  munmap(map, v_size);
  close(fd);

  // reading file back

  fd = open(fname.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
  struct stat sb;

  if (fstat(fd, &sb) == -1) {
    perror("couldn't get file size\n");
  }

  printf("file size is %ld\n", sb.st_size);

  double *file_in_memory = (double *)mmap(NULL, v_size, PROT_READ, MAP_SHARED, fd, 0);

  printf("%g", file_in_memory[0]);
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

