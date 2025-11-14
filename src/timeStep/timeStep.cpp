#include "src/timeStep/timeStep.hpp"

TimeStep::TimeStep(toml::parse_result parameters, const Mesh& mesh, FieldArrayManager fields)
  : _mesh(mesh), _CFL(parameters["time"]["CFL"].value_or(0.0)), _fields(fields) { }

double TimeStep::getTimeStep() {
  double dt = std::numeric_limits<double>::max();
  _mesh.loop().loopWithBoundaries([&dt, this](int i, int j) {
    auto minSpacing = std::min(_mesh.dx(), _mesh.dy());
    auto velocityMag = std::sqrt(std::pow(this->_fields(PV::U)[i, j], 2) + std::pow(this->_fields(PV::V)[i, j], 2));
    auto inviscidTimeStep = _CFL * minSpacing / velocityMag;
    if (inviscidTimeStep < dt) dt = inviscidTimeStep;
  });
  return dt;
}