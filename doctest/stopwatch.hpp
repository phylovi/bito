#pragma once

#include "../src/stopwatch.hpp"

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
