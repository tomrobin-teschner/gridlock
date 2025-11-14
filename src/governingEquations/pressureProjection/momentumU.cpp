#include "src/governingEquations/pressureProjection/momentumU.hpp"

MomentumU::MomentumU(FieldArrayManager fields, nlohmann::json parameters, const Mesh &mesh)
: _fields(fields), _parameters(parameters), _mesh(mesh) {
  _nu = _parameters["fluid"]["nu"];
}