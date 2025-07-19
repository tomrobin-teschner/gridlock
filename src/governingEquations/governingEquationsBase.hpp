#pragma once

#include "nlohmann/json.hpp"

#include "src/meshLooper/meshLooper.hpp"

class GoverningEquationsBase{
public:
  GoverningEquationsBase(nlohmann::json parameters, MeshLooper &meshLooper);
  virtual ~GoverningEquationsBase() = default;

public:
  virtual void solve() = 0;

protected:
  MeshLooper &_meshLooper;
  nlohmann::json _parameters;
};