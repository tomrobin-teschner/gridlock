#include "src/timeStep/timeStep.hpp"

TimeStep::TimeStep(nlohmann::json parameters, const Mesh& mesh ) : _mesh(mesh), _CFL(parameters["time"]["CFL"]) { }

double TimeStep::getTimeStep(FieldArray &u, FieldArray &v) {
  double dt = std::numeric_limits<double>::max();
  _mesh.loop().loopWithBoundaries([&dt, &u, &v, this](int i, int j) {
    auto minSpacing = std::min(_mesh.dx(), _mesh.dy());
    auto velocityMag = std::sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j]);
    auto inviscidTimeStep = _CFL * minSpacing / velocityMag;
    if (inviscidTimeStep < dt) dt = inviscidTimeStep;
  });
  return dt;
}