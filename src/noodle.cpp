#include "libsbn.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

std::string QuoteString(const std::string &in_str) {
  std::stringstream ss;
  ss << std::quoted(in_str);
  return ss.str();
}

std::string DequoteString(const std::string &in_str) {
  if (in_str.empty()) {
    return std::string();
  }
  char delimiter = in_str.at(0);
  if (delimiter != '\'' && delimiter != '"') {
    return std::string(in_str);
  }
  std::stringstream ss(in_str);
  std::string out_str;
  ss >> std::quoted(out_str, delimiter);
  return out_str;
}

TagStringMap TransformStringValues(std::function<std::string(const std::string &)> f,
                                   const TagStringMap &in_map) {
  TagStringMap out_map;
  for (const auto &[tag, value] : in_map) {
    out_map.insert({tag, f(value)});
  }
  return out_map;
}

TagStringMap DequoteTagStringMap(const TagStringMap &tag_string_map) {
  return TransformStringValues(DequoteString, tag_string_map);
}

int main() {
  std::string unquoted_test(R"raw(hello 'there" friend)raw");
  std::string double_quoted_test(R"raw("this is a \" test")raw");
  std::string double_quoted_dequoted(R"raw(this is a " test)raw");
  std::string single_quoted_test(R"raw('this is a \' test')raw");
  std::string single_quoted_dequoted(R"raw(this is a ' test)raw");

  assert(QuoteString(unquoted_test) == R"raw("hello 'there\" friend")raw");
  assert(DequoteString(double_quoted_test) == double_quoted_dequoted);
  assert(DequoteString(single_quoted_test) == single_quoted_dequoted);
  assert(DequoteString(QuoteString(unquoted_test)) == unquoted_test);

  TagStringMap test_map(
      {{2, unquoted_test}, {3, double_quoted_test}, {5, single_quoted_test}});
  TagStringMap expected_test_map(
      {{2, unquoted_test}, {3, double_quoted_dequoted}, {5, single_quoted_dequoted}});

  assert(expected_test_map == DequoteTagStringMap(test_map));

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
