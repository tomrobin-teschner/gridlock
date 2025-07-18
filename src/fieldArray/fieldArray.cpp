#include "src/fieldArray/fieldArray.hpp"

FieldArray::FieldArray(int numX, int numY) : _numX(numX), _numY(numY) {
  _data.resize(_numX * _numY);
}

FieldArray FieldArray::operator=(const FieldArray &other) {
  _numX = other._numX;
  _numY = other._numY;
  _data = other._data;
  return *this;
}

// allowing multidimensional indexing with C++23
double& FieldArray::operator[](int i, int j) {
  return _data[map2Dto1D(i, j)];
}

const double& FieldArray::operator[](int i, int j) const {
  return _data[map2Dto1D(i, j)];
}

// index the same data array with just a single index
double& FieldArray::operator[](int i) {
  return _data[i];
}

const double& FieldArray::operator[](int i) const {
  return _data[i];
}

int FieldArray::map2Dto1D(int i, int j) const {
  return j * _numX + i;
};

