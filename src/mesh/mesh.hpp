#pragma once

#include <iostream>

#include "src/infrastructure/parameters/parameters.hpp"
#include "src/fieldArray/fieldArray.hpp"
#include "src/meshLooper/meshLooper.hpp"

class Mesh {
public:
  Mesh(Parameters params);
  ~Mesh() = default;

  void create();
  
  double x(int i, int j) const { return _x[i, j]; }
  double y(int i, int j) const { return _y[i, j]; }
  
  double Lx() const { return _lx; }
  double Ly() const { return _ly; }
  double dx() const { return _dx; }
  double dy() const { return _dy; }
  int numX() const { return _numX; }
  int numY() const { return _numY; }
  int numGhostPoints() const { return _numGhostPoints; }
  
  const MeshLooper& loop() const { return _looper; }
  
private:
  int _numGhostPoints;
  int _numX, _numY;
  double _lx, _ly;
  double _dx, _dy;
  FieldArray _x;
  FieldArray _y;
  MeshLooper _looper;
};