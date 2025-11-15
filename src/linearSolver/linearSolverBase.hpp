#pragma once

#include <tuple>
#include <string>

#include "src/mesh/mesh.hpp"
#include "src/infrastructure/parameters/parameters.hpp"

template<typename MatrixType, typename VectorType, typename SolverType>
class LinearSolverBase {
public:
  using SolveType  = typename std::tuple<long long, VectorType>;

public:
  LinearSolverBase(Parameters params, const Mesh& mesh, std::string ID) : _numX(mesh.numX()), _numY(mesh.numY()),
    _underRelaxation(params.solver<double>("linearSolver", "underRelaxation", ID)),
    _maxIterations(params.solver<int>("linearSolver", "maxIterations", ID)),
    _tolerance(params.solver<double>("linearSolver", "tolerance", ID)) { }
  virtual ~LinearSolverBase() = default;

public:
  virtual void setZero() = 0;
  virtual void setMatrixAt(int i, int j, double value) = 0;
  virtual void setRHSAt(int i, double value) = 0;
  virtual void addRHSAt(int i, double value) = 0;
  virtual SolveType solve() = 0;
  double getUnderRelaxation() const { return _underRelaxation; }

protected:
  int _numX, _numY;
  double _underRelaxation;
  double _tolerance;
  int _maxIterations;
  MatrixType _A;
  VectorType _b;
  SolverType _solver;
};