#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

// https://www.youtube.com/watch?v=m7E9piHcfr4
// https://jameshfisher.com/2017/01/28/mmap-file-write/

int main(void) {
  int fd = open("_ignore/test_file", O_RDWR | O_CREAT, (mode_t)0600);
  const char *text = "hello";
  size_t textsize = strlen(text) + 1;
  lseek(fd, textsize-1, SEEK_SET);
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

