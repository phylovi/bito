#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include "eigen_sugar.hpp"

// https://www.youtube.com/watch?v=m7E9piHcfr4
// https://jameshfisher.com/2017/01/28/mmap-file-write/
// https://stackoverflow.com/a/43301909/467327

int main(void) {
  std::string fname("_ignore/test_file");
  int fd = open(fname.c_str(), O_RDWR | O_CREAT, (mode_t)0600);
  EigenVectorXd v(4);
  v << 1, 2, 3, 4;
  auto v_size = v.size() * sizeof(double);
  char *map = (char *)mmap(  //
      NULL,                  // This address is ignored because we are using MAP_SHARED.
      v_size,                // Size of map.
      PROT_READ | PROT_WRITE,  // We want to read and write
      MAP_SHARED,              // We need MAP_SHARED to actually write to this.
      fd,                      // File descriptor.
      0                        // Offset.
  );
  memcpy(map, v.data(), v_size);
  msync(map, v_size, MS_SYNC);
  //  for (int i = 0; i < strlen(text); i++) {
  //    printf("%c", map[i]);
  //  }
  munmap(map, v_size);
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

