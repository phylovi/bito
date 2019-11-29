#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>


const std::string FILE_NAME = "number.txt";
const std::string RESULT_FILE_NAME = "result.txt";

int main(void)
{
    int fd = ::open(FILE_NAME.c_str(), O_RDONLY);
    if (fd < 0) {
        std::cerr << "Error opeing " << FILE_NAME << std::endl;
        exit(1);
    }
    std::ofstream ofs(RESULT_FILE_NAME, std::ios_base::binary);
    if (!ofs) {
        std::cerr << "Error opeing " << RESULT_FILE_NAME << std::endl;
        exit(1);
    }

    auto file_size = lseek(fd, 0, SEEK_END);

    std::vector<std::chrono::milliseconds> duration_vec(5);
    for (std::vector<std::chrono::milliseconds>::size_type i = 0; i < duration_vec.size(); i++) {
        lseek(fd, 0, SEEK_SET);
        unsigned long long res = 0;
        auto begin = std::chrono::system_clock::now();

        char *chunk = reinterpret_cast<char*>(mmap(NULL, file_size, PROT_READ, MAP_FILE | MAP_SHARED, fd, 0));
        char *addr = chunk;

        for (size_t j = 0; j < file_size; j++) {
            res += *chunk++;
        }
        ofs << res;

        munmap(addr, file_size);

        duration_vec[i] = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - begin);
        std::cout<< duration_vec[i].count() << std::endl;
    }

    std::chrono::milliseconds total_time{0};
    for (auto const& v : duration_vec) {
        total_time += v;
    }
    std::cout << "Average exec time: " << total_time.count() / duration_vec.size() << std::endl;

    ::close(fd);
    return 0;
}
