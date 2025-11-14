#pragma once

#include <vector>
#include <string>
#include <map>

#include "toml++/toml.hpp"

#include "src/mesh/mesh.hpp"
#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArrayManager.hpp"

class BoundaryConditions {
public:
  BoundaryConditions(Mesh &mesh, FieldArrayManager fields, const toml::parse_result &bcParameters);

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
  const toml::parse_result &_bcParameters;
};