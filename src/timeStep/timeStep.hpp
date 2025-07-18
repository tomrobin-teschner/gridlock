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
  TimeStep(nlohmann::json parameters, MeshLooper looper);
  ~TimeStep() = default;

public:
  double getTimeStep(FieldArray &u, FieldArray &v, double dx, double dy);

private:
  MeshLooper _looper;
  double _CFL;
};