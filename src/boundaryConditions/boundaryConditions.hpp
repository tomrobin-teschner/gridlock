#pragma once

#include <vector>
#include <string>
#include <map>

#include "nlohmann/json.hpp"

#include "src/meshLooper/meshLooper.hpp"
#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArray.hpp"

class BoundaryConditions {
public:
  BoundaryConditions(const FieldType &x, const FieldType &y, MeshLooper &looper, const nlohmann::json &bcParameters);

public:
  void updateGhostPoints(std::string name, FieldArray &field);

private:
  void updateEastGhostPoints(std::string name, FieldArray &field);
  void updateWestGhostPoints(std::string name, FieldArray &field);
  void updateNorthGhostPoints(std::string name, FieldArray &field);
  void updateSouthGhostPoints(std::string name, FieldArray &field);

  double dirichlet(double phiBC, double phiInternal);
  double neumann(double phiBC, double phiInternal, double spacing);

private:
  const FieldType &_x, &_y;
  MeshLooper &_looper;
  const nlohmann::json &_bcParameters;
  int _numGhostPoints;
};