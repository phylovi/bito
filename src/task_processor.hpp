// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

// This code started life as the example from the excellent blog post at
// https://embeddedartistry.com/blog/2017/2/1/c11-implementing-a-dispatch-queue-using-stdfunction
//
// The setup is that we have a Task that we want to run on a bunch of units of
// Work. The Task needs an Executor to run, and we have a pool of those
// available. We want to run the Tasks in parallel on as many threads as we have
// Executors.
//
// For this, build queues of Executors and Work, then make a TaskProcessor out
// of them that does the work as specified in the queues. Work begins right away
// in the constructor. To get the work to complete, you can just let the
// TaskProcessor go out of scope, or explicitly call the Wait method for it to
// complete.
//
// Note that we don't make any effort to ensure safety. For example, there is
// nothing keeping you from having abundant data races if your Executors have
// non-independent state.
//
// The fact that there are fewer Executors than Tasks is what requires some
// design like this-- we can't just use something like a C++17 parallel for
// loop.
//
// I realize that it's not recommended to write your own thread-handling
// library, and for good reason: see https://www.youtube.com/watch?v=QIHy8pXbneI
// https://sean-parent.stlab.cc/presentations/2016-08-08-concurrency/2016-08-08-concurrency.pdf
// However, here the tasks are few and big, the time required during locking is
// small (put an integer in a dequeue). The overhead of including a true
// threading library wouldn't be worth it for this example.

#ifndef SRC_TASK_PROCESSOR_HPP_
#define SRC_TASK_PROCESSOR_HPP_

#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

template <class Executor, class Work>
class TaskProcessor {
 public:
  using Task = std::function<void(Executor, Work)>;
  using ExecutorQueue = std::queue<Executor>;
  using WorkQueue = std::queue<Work>;

  TaskProcessor(ExecutorQueue executor_queue, WorkQueue work_queue, Task task)
      : executor_queue_(executor_queue),
        work_queue_(work_queue),
        task_(task),
        threads_(executor_queue.size()) {
    // Make as many threads as there are executors.
    for (size_t i = 0; i < threads_.size(); i++) {
      threads_[i] = std::thread(&TaskProcessor::thread_handler, this);
    }

    Wait();
  }

  // Delete (copy + move) x (constructor + assignment)
  TaskProcessor(const TaskProcessor &) = delete;
  TaskProcessor(const TaskProcessor &&) = delete;
  TaskProcessor &operator=(const TaskProcessor &) = delete;
  TaskProcessor &operator=(const TaskProcessor &&) = delete;

  ~TaskProcessor() {}

  void Wait() {
    condition_variable_.notify_all();
    // Wait for threads to finish before we exit
    for (size_t i = 0; i < threads_.size(); i++) {
      if (threads_[i].joinable()) {
        threads_[i].join();
      }
    }
    threads_rejoined_ = true;
    if (exception_ptr_ != nullptr) {
      std::rethrow_exception(exception_ptr_);
    }
  }

 private:
  ExecutorQueue executor_queue_;
  WorkQueue work_queue_;
  Task task_;
  std::vector<std::thread> threads_;
  std::mutex lock_;
  std::condition_variable condition_variable_;
  bool threads_rejoined_ = false;

  std::exception_ptr exception_ptr_ = nullptr;
  std::mutex exception_lock_;

  void thread_handler() {
    std::unique_lock<std::mutex> lock(lock_);
    bool exception_occurred = false;
    // Continue allocating free executors for work until no more work.
    while (work_queue_.size()) {
      // Check if any thread has recorded an exception. If so, end work.
      exception_lock_.lock();
      exception_occurred = (exception_ptr_ != nullptr);
      exception_lock_.unlock();
      if (exception_occurred) {
        break;
      }
      // Wait until we have an executor available. This right here is the key of
      // the whole implementation, giving a nice way to wait until the resources
      // are available to run the next thing in the queue.
      condition_variable_.wait(lock, [this] { return executor_queue_.size(); });
      // After wait, we own the lock.
      if (work_queue_.size()) {
        auto work = work_queue_.front();
        work_queue_.pop();
        auto executor = executor_queue_.front();
        executor_queue_.pop();
        // Unlock now that we're done messing with the queues.
        lock.unlock();
        // Run task.
        try {
          task_(executor, work);
        } catch (...) {
          // If the task throws an exception, make record it for the master thread.
          exception_lock_.lock();
          if (exception_ptr_ == nullptr) {
            exception_ptr_ = std::current_exception();
          }
          exception_lock_.unlock();
        }

        // Lock again so that we can mess with the queues.
        lock.lock();
        // We're done with the executor so we can put it back on the queue.
        executor_queue_.push(executor);
      }
    }
  }
};

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
#endif  // SRC_TASK_PROCESSOR_HPP_
