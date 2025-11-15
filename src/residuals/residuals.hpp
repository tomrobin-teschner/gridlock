#pragma once

#include <cmath>
#include <fstream>
#include <iomanip>
#include <tuple>
#include <iostream>
#include <string>

#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArraymanager.hpp"
#include "src/mesh/mesh.hpp"
#include "src/infrastructure/parameters/parameters.hpp"

class Residuals {
public:
  using ResidualType = std::vector<double>;
public:
  Residuals(Parameters params, FieldArrayManager pv, Mesh &mesh, std::vector<int> IDs, std::string location,
    bool writeToFile = false);
  ~Residuals();

public:
  void init();
  ResidualType getResidual(int iteration);
  double getTolerance(int ID) { return _tolerances[ID]; };

private:
  double norm(int ID);
  void write(int iteration);

private:
  FieldArrayManager _pv;
  Mesh _mesh;
  bool _writeToFile;
  std::vector<double> _start, _end, _norm, _residual;
  std::ofstream _file;
  std::vector<int> _IDs;
  std::vector<double> _tolerances;
};