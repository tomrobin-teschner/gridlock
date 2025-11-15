#pragma once

#include <limits>
#include <memory>
#include <cmath>
#include <tuple>

#include "src/mesh/mesh.hpp"

#include "src/infrastructure/utilities/data.hpp"
#include "src/infrastructure/parameters/parameters.hpp"
#include "src/fieldArray/fieldArrayManager.hpp"

class TimeInfo {
public:
  TimeInfo(Parameters params, const Mesh& mesh, FieldArrayManager fields);
  ~TimeInfo() = default;

public:
  double getTimeStep();
  double getCFL() { return _params.solver<double>("time", "CFL"); }
  int getNumTimeSteps() { return _params.solver<int>("time", "timeSteps"); }
  int getOutputFrequency() { return _params.solver<int>("output", "outputFrequency"); }

private:
  Parameters _params;
  const Mesh& _mesh;
  double _CFL;
  FieldArrayManager _fields;
};