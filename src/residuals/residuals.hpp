#pragma once

#include <cmath>
#include <fstream>
#include <iomanip>
#include <tuple>
#include <iostream>

#include "src/fieldArray/fieldArray.hpp"
#include "src/mesh/mesh.hpp"

class Residuals {
public:
  using ResidualType = std::tuple<double, double, double>;
public:
  Residuals(FieldArray *u, FieldArray *v, FieldArray *p, Mesh &mesh, bool writeToFile = false);
  ~Residuals();

public:
  void init();
  ResidualType getResidual(int iteration);

private:
  double norm(FieldArray *_data);
  void write(int iteration);

private:
  FieldArray *_u, *_v, *_p;
  Mesh _mesh;
  bool _writeToFile;
  double _startU = 0.0; double _startV = 0.0; double _startP = 0.0;
  double _endU = 0.0; double _endV = 0.0; double _endP = 0.0;
  double _normU = 1.0; double _normV = 1.0; double _normP = 1.0;
  double _residualU = 0.0; double _residualV = 0.0; double _residualP = 0.0;
  std::ofstream _file;
};