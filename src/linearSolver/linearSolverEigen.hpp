#pragma once

#include "src/linearSolver/linearSolverBase.hpp"

template<typename MatrixType, typename VectorType, typename SolverType>
class LinearSolverEigen : public LinearSolverBase<MatrixType, VectorType, SolverType> {
public:
  using BaseClass = typename LinearSolverBase<MatrixType, VectorType, SolverType>;
public:
  LinearSolverEigen(int numX, int numY) : BaseClass(numX, numY) {
    this->_A.resize(numX * numY, numX * numY);
    this->_b.resize(numX * numY);
    setZero();
  }

public:
  virtual void setZero() override {
    this->_A.setZero();
    this->_b.setZero();
  }

  virtual void setMatrixAt(int i, int j, double value) override {
    this->_A.coeffRef(i, j) = value;
  }

  virtual void setRHSAt(int i, double value) override {
    this->_b(i) = value;
  }

  virtual void addRHSAt(int i, double value) override {
    this->_b(i) += value;
  }

  virtual typename BaseClass::SolveType solve(int maxIterations, double tolerance) override {
    this->_solver.compute(this->_A);
    this->_solver.setMaxIterations(maxIterations);
    this->_solver.setTolerance(tolerance);

    auto x = this->_solver.solve(this->_b);
    auto iterations = this->_solver.iterations();
    
    return {iterations, x};
  }
};