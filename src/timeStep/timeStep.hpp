#pragma once

#include <limits>
#include <memory>
#include <cmath>
#include <tuple>

#include "nlohmann/json.hpp"

#include "src/meshLooper/meshLooper.hpp"

#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArray.hpp"

class TimeStep {
public:
  using TimeStepType = std::tuple<double, double>;

public:
  TimeStep(nlohmann::json parameters, MeshLooper looper);
  ~TimeStep() = default;

public:
  TimeStepType getTimeStep(FieldArray &u, FieldArray &v, double dx, double dy, double resU, double resV,
    double resP);

private:
  MeshLooper _looper;
  double _CFL;
  double _maxCFL;
  bool _useSER;
  double _alphaSER;
};