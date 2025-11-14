#pragma once

#include <tuple>

#include "toml++/toml.hpp"

#include "src/fieldArray/fieldArrayManager.hpp"
#include "src/mesh/mesh.hpp"

class MomentumU {
public:
  using CoefficientType = typename std::tuple<double, double, double, double, double, double>;
 
public:
  MomentumU(FieldArrayManager fields, toml::parse_result parameters, const Mesh &mesh);
  ~MomentumU() = default;

private:
  FieldArrayManager _fields;
  toml::parse_result _parameters;
  const Mesh &_mesh;
  double _nu;
};