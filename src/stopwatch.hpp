// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Simple stopwatch class is for performing runtime speed tests.
// Can specify upon construction the time units (in seconds, milliseconds, or
// nanoseconds). Stopwatch uses milliseconds by default. Stopwatch has two states:
// either running or not. User can only capture times through GetElapsed...(), Lap() or
// Stop() while Stopwatch is running. Laps capture the elapsed time between each Lap()
// call (as well as Start()/Lap() and Lap()/Stop() intervals). GetTotal() returns a
// culumative total of all Start()/Stop() time intervals. User can't call Stop() on a
// stopped watch or Start() on a running watch. Clear() removes all memory and stops
// watch by default, but can optionally restart watch immediatedly.

#pragma once

// ** Doctest include must go first for all header tests to run.
#include "doctest.h"
// **

#include <stdio.h>

#include <chrono>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

#include "sugar.hpp"

class Stopwatch {
 public:
  using clock = std::chrono::high_resolution_clock;
  using time_point = std::chrono::time_point<clock>;
  using time_duration = std::chrono::duration<uint16_t, std::nano>;

  enum class TimeScale { SecondScale, MillisecondScale, NanosecondScale };
  using second_double = std::chrono::duration<double, std::chrono::seconds::period>;
  using nanosec_double =
      std::chrono::duration<double, std::chrono::nanoseconds::period>;
  using millisec_double =
      std::chrono::duration<double, std::chrono::milliseconds::period>;

  // Constructor:
  // Option to start stopwatch running upon initialization, defaults to not starting.
  // Option to specify time scale, defaults to milliseconds.
  Stopwatch(bool start_on_init = false, TimeScale scale = TimeScale::MillisecondScale)
      : scale_(scale), lap_seconds_(), interval_starts_(), is_running_(false) {
    interval_starts_.push_back(0);
    if (start_on_init) {
      Start();
    }
  };
  // Return the state of the watch.
  bool IsRunning() { return is_running_; }

  // Returns time elapsed from ith start/stop interval.
  double GetInterval(size_t which_interval) {
    Assert(which_interval < GetIntervalCount(),
           "Stopwatch.GetInterval(): ith start/stop interval does not exist.");
    if (is_running_ && (which_interval == GetIntervalCount() - 1)) {
      return GetElapsedOfCurrentInterval();
    }
    size_t start_point = interval_starts_[which_interval];
    size_t end_point = interval_starts_[which_interval + 1];
    double interval_total = 0.0;
    for (size_t i = start_point; i < end_point; i++) {
      interval_total += lap_seconds_[i];
    }
    return interval_total;
  };
  // Returns time elapsed during the latest completed start/stop interval.
  double GetLatestInterval() {
    Assert(GetCompleteIntervalCount() > 0,
           "Stopwatch.GetLatestInterval() cannot be called when no complete intervals "
           "exist.");
    size_t last_pos = GetCompleteIntervalCount() - 1;
    return GetInterval(last_pos);
  };
  // Returns time elapsed since last start. Retains no memory of this in Stopwatch.
  double GetElapsedOfCurrentInterval() {
    Assert(is_running_,
           "Stopwatch.GetElapsedOfCurrentInterval() cannot be called while Stopwatch "
           "is not running.");
    time_point now = GetCurrentTime();
    double elapsed_seconds = 0.0;
    size_t start_lap = interval_starts_[interval_starts_.size() - 1];
    for (size_t i = start_lap; i < lap_seconds_.size(); i++) {
      elapsed_seconds += lap_seconds_[i];
    }
    elapsed_seconds += GetElapsedOfCurrentLap();
    return elapsed_seconds;
  };
  // Returns number of time intervals that have occurred, including currently running.
  size_t GetIntervalCount() {
    return GetCompleteIntervalCount() + static_cast<size_t>(is_running_);
  }
  // Returns number of time intervals that have completed.
  size_t GetCompleteIntervalCount() {
    Assert(interval_starts_.size() >= 1,
           "Stopwatch.interval_starts_ should not be empty.");
    return interval_starts_.size() - 1;
  }

  // Returns total time elapsed over all start/stop intervals.
  // If watch is currently running, returns total over all previous interval plus time
  // elapsed in current interval.
  double GetTotal() {
    float total = 0.0;
    for (size_t i = 0; i < lap_seconds_.size(); i++) {
      total += lap_seconds_[i];
    }
    return (is_running_ ? total + GetElapsedOfCurrentLap() : total);
  };

  // Appends the time elapsed since latest lap (or start if on first lap) to laps and
  // returns time.
  double Lap() {
    Assert(is_running_,
           "Stopwatch.Lap() cannot be called while Stopwatch is not running.");
    time_point now = GetCurrentTime();
    double lap_seconds = DurationToSeconds(latest_lap_, now);
    lap_seconds_.push_back(lap_seconds);
    latest_lap_ = now;
    return lap_seconds;
  };
  // Get time elapsed over the ith lap.
  double GetLap(size_t which_lap) {
    Assert(which_lap < GetLapCount(), "Stopwatch.GetLap(): ith lap does not exist.");
    return lap_seconds_.at(which_lap);
  }
  // Returns the time interval between the latest completed lap.
  double GetLatestLap() { return GetLap(GetLapCount() - 1); };
  // Returns time elapsed since latest lap.  Retains no memory of this in Stopwatch.
  double GetElapsedOfCurrentLap() {
    Assert(is_running_,
           "Stopwatch.Elapsed() cannot be called while Stopwatch is not running.");
    time_point now = GetCurrentTime();
    double elapsed_seconds = DurationToSeconds(latest_lap_, now);
    return elapsed_seconds;
  };
  // Returns vector of all lap times.
  std::vector<double> GetLaps() { return lap_seconds_; }
  // Get number of completed laps.
  size_t GetLapCount() { return lap_seconds_.size(); }

  // Starts running stopwatch.
  void Start() {
    Assert(is_running_ == false,
           "Stopwatch.Start() cannot be called while Stopwatch is already running.");
    start_time_ = GetCurrentTime();
    latest_lap_ = start_time_;
    is_running_ = true;
  };
  // Returns time elapsed between current start/stop interval, and adds time elapsed
  // from latest lap to stop to laps.
  double Stop() {
    Assert(is_running_,
           "Stopwatch.Stop() cannot be called while Stopwatch is not running.");
    time_point now = GetCurrentTime();
    double lap_seconds = DurationToSeconds(latest_lap_, now);
    lap_seconds_.push_back(lap_seconds);
    interval_starts_.push_back(lap_seconds_.size());
    double latest_interval = GetLatestInterval();
    is_running_ = false;
    return latest_interval;
  };
  // Clears total time and all lap times from memory.
  // Option to restart watch upon clearing memory, but defaults to stopped.
  void Clear(bool restart_watch = false) {
    is_running_ = false;
    lap_seconds_.clear();
    interval_starts_.clear();
    interval_starts_.push_back(0);
    if (restart_watch) {
      Start();
    }
  };

  // Puts thread to sleep for given milliseconds.
  static void Sleep(uint64_t wait_time) {
    std::this_thread::sleep_for(std::chrono::milliseconds(wait_time));
  };

 private:
  time_point GetCurrentTime() { return clock::now(); }

  // Convert a time duration to a double in the specified time scale.
  double DurationToSeconds(time_point &start_time, time_point &end_time) {
    double seconds;
    switch (scale_) {
      case TimeScale::SecondScale:
        seconds = second_double(end_time - start_time).count();
        break;
      case TimeScale::MillisecondScale:
        seconds = millisec_double(end_time - start_time).count();
        break;
      case TimeScale::NanosecondScale:
        seconds = nanosec_double(end_time - start_time).count();
        break;
    }
    return seconds;
  };

  bool is_running_;
  TimeScale scale_;
  time_point start_time_;
  time_point latest_lap_;
  // Time intervals between each lap.
  std::vector<double> lap_seconds_;
  // Starting lap index for each time interval.
  std::vector<size_t> interval_starts_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Stopwatch") {
  Stopwatch watch(false, Stopwatch::TimeScale::MillisecondScale);
  // Functions not allowed while clock isn't running.
  CHECK_THROWS(watch.GetElapsedOfCurrentInterval());
  CHECK_THROWS(watch.GetElapsedOfCurrentLap());
  CHECK_THROWS(watch.GetLatestInterval());
  CHECK_THROWS(watch.Lap());
  CHECK_THROWS(watch.Stop());

  watch.Start();
  // Getting latest lap before first lap exists.
  CHECK_THROWS(watch.GetLatestLap());
  // Getting latest interval before first interval exists.
  CHECK_THROWS(watch.GetLatestInterval());
  Stopwatch::Sleep(3);
  watch.Stop();

  auto interval_1 = watch.GetLatestInterval();
  Stopwatch::Sleep(5);

  watch.Start();
  Stopwatch::Sleep(7);
  auto interval_1_during_next_interval = watch.GetLatestInterval();
  Stopwatch::Sleep(11);
  auto lap_1 = watch.Lap();
  Stopwatch::Sleep(13);
  auto lap_1_during_next_lap = watch.GetLatestLap();
  auto interval_2_midinterval = watch.GetElapsedOfCurrentInterval();
  Stopwatch::Sleep(17);
  watch.GetElapsedOfCurrentLap();
  Stopwatch::Sleep(19);
  auto interval_2 = watch.Stop();

  // Latest should fetch the last completed lap or interval, even if clock is running.
  CHECK_EQ(interval_1, interval_1_during_next_interval);
  CHECK_EQ(lap_1, lap_1_during_next_lap);
  // The mid-interval time should be less than the total time.
  CHECK_GT(interval_2, interval_2_midinterval);

  auto laps = watch.GetLaps();
  auto total = watch.GetTotal();
  std::vector<double> intervals = {interval_1, interval_2};
  auto sum_laps = std::accumulate(laps.begin(), laps.end(), 0.0);
  auto sum_intervals = std::accumulate(intervals.begin(), intervals.end(), 0.0);

  watch.Clear();
  CHECK_EQ(watch.GetTotal(), 0.0);
  CHECK_EQ(watch.GetLaps().size(), 0);

  watch.Start();
  // Function not allowed while clock is running.
  CHECK_THROWS(watch.Start());

  CHECK_EQ(doctest::Approx(sum_laps), total);
  CHECK_EQ(doctest::Approx(sum_intervals), total);
};

#endif  // DOCTEST_LIBRARY_INCLUDED
