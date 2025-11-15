#pragma once

#include "src/linearSolver/linearSolverBase.hpp"
#include "src/mesh/mesh.hpp"

template<typename MatrixType, typename VectorType, typename SolverType>
class LinearSolverEigen : public LinearSolverBase<MatrixType, VectorType, SolverType> {
public:
  using BaseClass = typename LinearSolverBase<MatrixType, VectorType, SolverType>;
public:
  LinearSolverEigen(Parameters params, const Mesh& mesh, std::string ID) : BaseClass(params, mesh, ID) {
    this->_A.resize(this->_numX * this->_numY, this->_numX * this->_numY);
    this->_b.resize(this->_numX * this->_numY);
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

  virtual typename BaseClass::SolveType solve() override {
    this->_solver.compute(this->_A);
    this->_solver.setMaxIterations(this->_maxIterations);
    this->_solver.setTolerance(this->_tolerance);

    auto x = this->_solver.solve(this->_b);
    auto iterations = this->_solver.iterations();
    
    return {iterations, x};
  }
};