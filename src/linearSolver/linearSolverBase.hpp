#pragma once

#include <tuple>

template<typename MatrixType, typename VectorType, typename SolverType>
class LinearSolverBase {
public:
  using SolveType  = typename std::tuple<long long, VectorType>;

public:
  LinearSolverBase(int numX, int numY) : _numX(numX), _numY(numY) { }
  virtual ~LinearSolverBase() = default;

public:
  virtual void setZero() = 0;
  virtual void setMatrixAt(int i, int j, double value) = 0;
  virtual void setRHSAt(int i, double value) = 0;
  virtual void addRHSAt(int i, double value) = 0;
  virtual SolveType solve(int maxIterations, double tolerance) = 0;

protected:
  int _numX, _numY;
  MatrixType _A;
  VectorType _b;
  SolverType _solver;
};