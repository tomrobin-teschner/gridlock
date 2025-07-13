#pragma once

#include <limits>
#include <memory>
#include <cmath>
#include <tuple>

#include "nlohmann/json.hpp"

#include "src/meshLooper.hpp"

#include "src/infrastructure/utilities/data.hpp"

class TimeStep {
public:
  using TimeStepType = std::tuple<double, double>;

public:
  TimeStep(nlohmann::json parameters, MeshLooper looper);
  ~TimeStep() = default;

public:
  TimeStepType getTimeStep(std::shared_ptr<FieldType> u, std::shared_ptr<FieldType> v, std::shared_ptr<FieldType> p,
    double dx, double dy, double resU, double resV, double resP);

private:
  MeshLooper _looper;
  double _CFL;
  double _maxCFL;
  bool _useSER;
  double _alphaSER;
};