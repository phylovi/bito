#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

int main(void)
{
  int fd = open("./README.md", O_RDONLY, S_IRUSR | S_IWUSR);
  struct stat sb;

  if (fstat(fd, &sb) == -1) {
    perror("couldn't get file size\n");
  }

  printf("file size is %ld\n", sb.st_size);

  char *file_in_memory = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);

  for (int i = 0; i < sb.st_size; i++) {
    printf("%c", file_in_memory[i]);
  }
  printf("\n");
  munmap(file_in_memory, sb.st_size);
  close(fd);

  return 0;
}
