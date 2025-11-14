#include "src/governingEquations/pressureProjection/momentumU.hpp"

MomentumU::MomentumU(FieldArrayManager fields, toml::parse_result parameters, const Mesh &mesh)
: _fields(fields), _parameters(parameters), _mesh(mesh) {
  _nu = _parameters["fluid"]["nu"].value_or(1.0);
}