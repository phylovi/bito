#pragma once

#include "../src/task_processor.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TaskProcessor") {
  std::queue<int> executor_queue;
  std::queue<size_t> work_queue;
  std::vector<float> results(8);
  // Say we have 4 executors.
  for (auto i = 0; i < 4; i++) {
    executor_queue.push(i);
  }
  // Our Work in this example is just size_t's.
  for (size_t i = 0; i < results.size(); i++) {
    work_queue.push(i);
  }
  // And our task is just to cast this size_t to a float and store it in the
  // corresponding location of the results array.
  auto task = [&results](int /*executor*/, size_t work) {
    // std::cout << "work " << work << " on " << executor << std::endl;
    results[work] = static_cast<float>(work);
  };
  TaskProcessor<int, size_t> processor(executor_queue, work_queue, task);
  processor.Wait();
  std::vector<float> correct_results({0, 1, 2, 3, 4, 5, 6, 7});
  CHECK_EQ(results, correct_results);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
