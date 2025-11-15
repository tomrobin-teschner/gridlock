#include "src/timeInfo/timeInfo.hpp"

TimeInfo::TimeInfo(Parameters params, const Mesh& mesh, FieldArrayManager fields)
  : _params(params), _mesh(mesh), _CFL(params.solver<double>("time", "CFL")), _fields(fields) { }

double TimeInfo::getTimeStep() {
  double dt = std::numeric_limits<double>::max();
  _mesh.loop().loopWithBoundaries([&dt, this](int i, int j) {
    auto minSpacing = std::min(_mesh.dx(), _mesh.dy());
    auto velocityMag = std::sqrt(std::pow(this->_fields(PV::U)[i, j], 2) + std::pow(this->_fields(PV::V)[i, j], 2));
    auto inviscidTimeStep = _CFL * minSpacing / velocityMag;
    if (inviscidTimeStep < dt) dt = inviscidTimeStep;
  });
  return dt;
}