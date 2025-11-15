#include "src/infrastructure/stopWatch/stopWatch.hpp"

StopWatch::StopWatch(Parameters params) : _timeSteps(params.solver<int>("time", "timeSteps")) {
  _start = std::chrono::high_resolution_clock::now();
}

StopWatch::TimeType StopWatch::elapsed() const {
  return hhmmss(timeToDatum());
}

StopWatch::TimeType StopWatch::remaining(int currentTimeStep) const {
  auto duration = timeToDatum();
  auto averageDuration = duration / (currentTimeStep + 1);
  auto timeRemaining = (_timeSteps - (currentTimeStep + 1)) * averageDuration;
  return hhmmss(timeRemaining);
}

double StopWatch::timeToDatum() const {
  auto now = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - _start);

  // return duration as double and convert to seconds from milliseconds
  return 0.001 * duration.count();
}

StopWatch::TimeType StopWatch::hhmmss(double duration) const {
  int hh = std::floor(duration / 60.0 / 60.0);
  int mm = std::floor(duration / 60.0 - hh * 60.0);
  int ss = std::floor(duration - hh * 60.0 * 60.0 - mm * 60.0);
  return std::make_tuple(hh, mm, ss);
}