#pragma once

#include <limits>
#include <memory>
#include <cmath>
#include <tuple>

#include "toml++/toml.hpp"

#include "src/mesh/mesh.hpp"

#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArrayManager.hpp"

class TimeStep {
public:
  TimeStep(toml::parse_result, const Mesh& mesh, FieldArrayManager fields);
  ~TimeStep() = default;

public:
  double getTimeStep();

private:
  const Mesh& _mesh;
  double _CFL;
  FieldArrayManager _fields;
};