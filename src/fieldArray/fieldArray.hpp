#pragma once

#include "Eigen/Eigen"

#include "src/infrastructure/utilities/data.hpp"

class FieldArray {
public:
  using FieldArrayType = typename Eigen::VectorXd;

public:
  FieldArray(int numX, int numY) : _numX(numX), _numY(numY) {};
  virtual ~FieldArray() = default;

public:
  

private:
  int _numX, _numY;
  FieldArrayType _data;
};