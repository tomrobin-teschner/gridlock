#include "src/governingEquations/governingEquationsBase.hpp"

GoverningEquationsBase::GoverningEquationsBase(nlohmann::json parameters, MeshLooper &meshLooper)
  : _parameters(parameters), _meshLooper(meshLooper) {
  
}