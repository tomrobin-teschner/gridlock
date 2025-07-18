#pragma once

class MeshLooper {
public:
  MeshLooper(int numX, int numY, int numGhostPoints) : _numX(numX), _numY(numY), _numGhostPoints(numGhostPoints) {}

public:
  template<class FunctionType>
  void loopAll(const FunctionType &f) {
    for (int i = 0; i < _numX + 2 * _numGhostPoints; ++i)
      for (int j = 0; j < _numY + 2 * _numGhostPoints; ++j)
          f(i, j);
  }

  template<class FunctionType>
  void loopWithBoundaries(const FunctionType &f) {
    for (int i = _numGhostPoints; i < _numX + _numGhostPoints; ++i)
      for (int j = _numGhostPoints; j < _numY + _numGhostPoints; ++j)
          f(i, j);
  }

  template<class FunctionType>
  void loopInterior(const FunctionType &f) {
    for (int i = _numGhostPoints + 1; i < _numX + _numGhostPoints - 1; ++i)
      for (int j = _numGhostPoints + 1; j < _numY + _numGhostPoints - 1; ++j)
          f(i, j);
  }

  template<class FunctionType>
  void eastBC(const FunctionType &f) {
    int i = _numX + _numGhostPoints - 1;
    for (int j = _numGhostPoints; j < _numY + _numGhostPoints; ++j)
      f(i, j);
  }

  template<class FunctionType>
  void westBC(const FunctionType &f) {
    int i = _numGhostPoints;
    for (int j = _numGhostPoints; j < _numY + _numGhostPoints; ++j)
      f(i, j);
  }

  template<class FunctionType>
  void northBC(const FunctionType &f) {
    int j = _numY + _numGhostPoints - 1;
    for (int i = _numGhostPoints; i < _numX + _numGhostPoints; ++i)
      f(i, j);
  }

  template<class FunctionType>
  void southBC(const FunctionType &f) {
    int j = _numGhostPoints;
    for (int i = _numGhostPoints; i < _numX + _numGhostPoints; ++i)
      f(i, j);
  }

  int getNumX() const { return _numX; }
  int getNumY() const { return _numY; }
  int getNumGhostPoints() const { return _numGhostPoints; }

private:
  int _numX, _numY, _numGhostPoints;
};
