#pragma once

class MeshLooper {
public:
  MeshLooper(int numX, int numY, int numGhostPoints) : _numX(numX), _numY(numY), _numGhostPoints(numGhostPoints) {}

public:
  template<class FunctionType>
  void loopAll(const FunctionType &f) const {
    for (int i = 0; i < _numX + 2 * _numGhostPoints; ++i)
      for (int j = 0; j < _numY + 2 * _numGhostPoints; ++j)
          f(i, j);
  }

  template<class FunctionType>
  void loopWithBoundaries(const FunctionType &f) const {
    for (int i = _numGhostPoints; i < _numX + _numGhostPoints; ++i)
      for (int j = _numGhostPoints; j < _numY + _numGhostPoints; ++j)
          f(i, j);
  }

  template<class FunctionType>
  void loopWithBoundariesReversed(const FunctionType &f) const {
    for (int j = _numGhostPoints; j < _numY + _numGhostPoints; ++j)
      for (int i = _numGhostPoints; i < _numX + _numGhostPoints; ++i)
          f(i, j);
  }

  template<class FunctionType>
  void loopInterior(const FunctionType &f) const {
    for (int i = _numGhostPoints + 1; i < _numX + _numGhostPoints - 1; ++i)
      for (int j = _numGhostPoints + 1; j < _numY + _numGhostPoints - 1; ++j)
          f(i, j);
  }

  template<class FunctionType>
  void eastBC(const FunctionType &f) const {
    int i = _numX + _numGhostPoints - 1;
    for (int j = _numGhostPoints; j < _numY + _numGhostPoints; ++j)
      f(i, j);
  }

  template<class FunctionType>
  void westBC(const FunctionType &f) const {
    int i = _numGhostPoints;
    for (int j = _numGhostPoints; j < _numY + _numGhostPoints; ++j)
      f(i, j);
  }

  template<class FunctionType>
  void northBC(const FunctionType &f) const {
    int j = _numY + _numGhostPoints - 1;
    for (int i = _numGhostPoints; i < _numX + _numGhostPoints; ++i)
      f(i, j);
  }

  template<class FunctionType>
  void southBC(const FunctionType &f) const {
    int j = _numGhostPoints;
    for (int i = _numGhostPoints; i < _numX + _numGhostPoints; ++i)
      f(i, j);
  }

  std::tuple<int, int> zeroBasedIndices(int i, int j) const {
    return {i - _numGhostPoints, j - _numGhostPoints};
  }

  std::tuple<int, int, int, int, int> getMatrixIndices(int idx, int jdx) const {
    auto ic  = map2Dto1D(idx, jdx);
    auto ip1 = map2Dto1D(idx + 1, jdx);
    auto im1 = map2Dto1D(idx - 1, jdx);
    auto jp1 = map2Dto1D(idx, jdx + 1);
    auto jm1 = map2Dto1D(idx, jdx - 1);

    return {ic, ip1, im1, jp1, jm1};
  }

  int map2Dto1D(int i, int j) const { return j * _numX + i; };

private:
  int _numX, _numY, _numGhostPoints;
};
