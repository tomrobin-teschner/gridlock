#pragma once

#include <chrono>
#include <cmath>
#include <tuple>

#include "src/infrastructure/parameters/parameters.hpp"

class StopWatch {
public:
  using TimeType = std::tuple<int, int, int>;
public:
  StopWatch(Parameters params);
  ~StopWatch() = default;

  TimeType elapsed() const;
  TimeType remaining(int currentTimeStep) const;

private:
  double timeToDatum() const;
  TimeType hhmmss(double ms) const;

private:
  int _timeSteps;
  std::chrono::high_resolution_clock::time_point _start;
};
