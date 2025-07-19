#pragma once

#include <Eigen/Eigen>

#include "src/governingEquations/governingEquationsBase.hpp"

class uMomentum : public GoverningEquationsBase {
public:
  uMomentum(nlohmann::json parameters, MeshLooper &meshLooper);
  virtual ~uMomentum() = default;

public:
  void solve() override;
};