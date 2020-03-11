#include "libsbn.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

// TagStringMap TransformStringValues(std::function<std::string(std::string)>,
std::string TransformString(std::function<std::string(const std::string &)> f,
                            std::string s) {
  return f(s);
}

// TagStringMap DequoteTagStringMap(TagStringMap tag_string_map);

// Parsing:
// Quotation marks are always taken off, and an error is thrown if quotation marks are
// not matching. If preserve_underscores is True, then underscores are kept. If
// preserve_underscores is False, then underscores are replaced with spaces. If a string
// is quoted,

std::string QuoteString(const std::string &in_str) {
  std::stringstream ss;
  ss << std::quoted(in_str);
  return ss.str();
}

std::string DequoteString(const std::string &in_str) {
  if (in_str.empty()) {
    return in_str;
  }
  // If the first character is an open single quote, then we assume that this is the
  // delimiter. Otherwise we assume that it's a double quote.
  char delim = (in_str.at(0) == '\'') ? '\'' : '"';
  std::stringstream ss(in_str);
  std::string out_str;
  ss >> std::quoted(out_str, delim);
  return out_str;
}

int main() {
  std::string unquoted_test(R"raw(hello 'there" friend)raw");
  std::string double_quoted_test(R"raw("this is a \" test")raw");
  std::string double_quoted_dequoted(R"raw(this is a " test)raw");
  std::string single_quoted_test(R"raw('this is a \' test')raw");
  std::string single_quoted_dequoted(R"raw(this is a ' test)raw");

  std::cout << QuoteString(unquoted_test) << std::endl;
  std::cout << (DequoteString(double_quoted_test) == double_quoted_dequoted)
            << std::endl;
  std::cout << (DequoteString(single_quoted_test) == single_quoted_dequoted)
            << std::endl;

  // R"raw("hi 'there\" ")raw"
  // R"raw(this is a " test)raw"

  return 0;

  uint32_t leaf_count = 10000;

  Node::NodePtr topology = Node::Ladder(leaf_count);

  std::vector<size_t> ids;
  ids.reserve(1 + 2 * leaf_count);

  auto t_start = now();
  for (int i = 0; i < 100; i++) {
    ids.clear();
    topology->PreOrder([&ids](const Node* node) { ids.push_back(node->Id()); });
  }

  std::chrono::duration<double> duration = now() - t_start;
  std::cout << "time: " << duration.count() << " seconds\n";
}
