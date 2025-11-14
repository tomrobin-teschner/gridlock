#pragma once

#include <vector>
#include <string>
#include <map>

#include "nlohmann/json.hpp"

#include "src/mesh/mesh.hpp"
#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArrayManager.hpp"

class BoundaryConditions {
public:
  BoundaryConditions(Mesh &mesh, FieldArrayManager fields, const nlohmann::json &bcParameters);

public:
  void updateGhostPoints(int ID);
  void updateGhostPoints(std::vector<int> IDs);

private:
  void updateEastGhostPoints(int ID);
  void updateWestGhostPoints(int ID);
  void updateNorthGhostPoints(int ID);
  void updateSouthGhostPoints(int ID);

  double dirichlet(double phiBC, double phiInternal);
  double neumann(double phiBC, double phiInternal, double spacing);

private:
  Mesh &_mesh;
  FieldArrayManager _fields;
  const nlohmann::json &_bcParameters;
};