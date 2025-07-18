#include "src/timeStep/timeStep.hpp"

TimeStep::TimeStep(nlohmann::json parameters, MeshLooper looper) : _looper(looper), _CFL(parameters["time"]["CFL"]) { }

double TimeStep::getTimeStep(FieldArray &u, FieldArray &v, double dx, double dy) {
  double dt = std::numeric_limits<double>::max();
  _looper.loopWithBoundaries([dx, dy, &dt, &u, &v, this](int i, int j) {
    auto minSpacing = std::min(dx, dy);
    auto velocityMag = std::sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j]);
    auto inviscidTimeStep = _CFL * minSpacing / velocityMag;
    if (inviscidTimeStep < dt) dt = inviscidTimeStep;
  });
  return dt;
}