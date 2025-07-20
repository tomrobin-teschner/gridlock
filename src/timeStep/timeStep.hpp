#pragma once

#include <limits>
#include <memory>
#include <cmath>
#include <tuple>

#include "nlohmann/json.hpp"

#include "src/mesh/mesh.hpp"

#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArray.hpp"

class TimeStep {
public:
  TimeStep(nlohmann::json parameters, const Mesh& mesh);
  ~TimeStep() = default;

public:
  double getTimeStep(FieldArray &u, FieldArray &v);

private:
  const Mesh& _mesh;
  double _CFL;
};