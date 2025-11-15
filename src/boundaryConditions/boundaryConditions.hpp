#pragma once

#include <vector>
#include <string>
#include <map>

#include "src/mesh/mesh.hpp"
#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArrayManager.hpp"
#include "src/infrastructure/parameters/parameters.hpp"

class BoundaryConditions {
public:
  BoundaryConditions(Parameters params, Mesh &mesh, FieldArrayManager fields);

public:
  void updateGhostPoints(int ID);
  void updateGhostPoints(std::vector<int> IDs);
  bool fullyNeumann();

private:
  void updateEastGhostPoints(int ID);
  void updateWestGhostPoints(int ID);
  void updateNorthGhostPoints(int ID);
  void updateSouthGhostPoints(int ID);

  double dirichlet(double phiBC, double phiInternal);
  double neumann(double phiBC, double phiInternal, double spacing);

private:
Parameters _params;
  Mesh &_mesh;
  FieldArrayManager _fields;
};