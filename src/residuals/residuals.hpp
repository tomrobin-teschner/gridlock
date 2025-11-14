#pragma once

#include <cmath>
#include <fstream>
#include <iomanip>
#include <tuple>
#include <iostream>

#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArraymanager.hpp"
#include "src/mesh/mesh.hpp"

class Residuals {
public:
  using ResidualType = std::vector<double>;
public:
  Residuals(FieldArrayManager pv, Mesh &mesh, std::vector<int> IDs, bool writeToFile = false);
  ~Residuals();

public:
  void init();
  ResidualType getResidual(int iteration);

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
};