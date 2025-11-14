#pragma once

#include <tuple>

#include "nlohmann/json.hpp"

#include "src/fieldArray/fieldArrayManager.hpp"
#include "src/mesh/mesh.hpp"

class MomentumU {
public:
  using CoefficientType = typename std::tuple<double, double, double, double, double, double>;
 
public:
  MomentumU(FieldArrayManager fields, nlohmann::json parameters, const Mesh &mesh);
  ~MomentumU() = default;

private:
  FieldArrayManager _fields;
  nlohmann::json _parameters;
  const Mesh &_mesh;
  double _nu;
};