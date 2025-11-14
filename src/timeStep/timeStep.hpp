#pragma once

#include <limits>
#include <memory>
#include <cmath>
#include <tuple>

#include "nlohmann/json.hpp"

#include "src/mesh/mesh.hpp"

#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArrayManager.hpp"

class TimeStep {
public:
  TimeStep(nlohmann::json parameters, const Mesh& mesh, FieldArrayManager fields);
  ~TimeStep() = default;

public:
  double getTimeStep();

private:
  const Mesh& _mesh;
  double _CFL;
  FieldArrayManager _fields;
};