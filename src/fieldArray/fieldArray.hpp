#pragma once

#include "Eigen/Eigen"

#include "src/infrastructure/utilities/data.hpp"

class FieldArray {
public:
  using FieldArrayType = typename Eigen::VectorXd;

public:
  FieldArray(int numX, int numY);
  ~FieldArray() = default;

public:
  FieldArray operator=(const FieldArray &other);

  // allowing multidimensional indexing with C++23
  double& operator[](int i, int j);
  const double& operator[](int i, int j) const;

  // index the same data array with just a single index
  double& operator[](int i);
  const double& operator[](int i) const;

private:
  int map2Dto1D(int i, int j) const;

private:
  int _numX, _numY;
  FieldArrayType _data;
};