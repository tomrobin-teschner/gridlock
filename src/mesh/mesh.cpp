#include "src/mesh/mesh.hpp"

Mesh::Mesh(Parameters params) : _numGhostPoints(params.mesh<int>("mesh", "numGhostPoints")),
  _numX(params.mesh<int>("mesh", "numX")), _numY(params.mesh<int>("mesh", "numY")),
  _lx(params.mesh<double>("mesh", "Lx")), _ly(params.mesh<double>("mesh", "Ly")),
  _dx(_lx / (_numX - 1)), _dy(_ly / (_numY - 1)), _x(_numX + 2 * _numGhostPoints, _numY + 2 * _numGhostPoints),
  _y(_numX + 2 * _numGhostPoints, _numY + 2 * _numGhostPoints), _looper(_numX, _numY, _numGhostPoints) {

  // create mesh
  create();
}

void Mesh::create() {
  double alpha = 5.0;
  _looper.loopAll([this, alpha](int i, int j) {
    _x[i, j] = _dx * i - _dx * _numGhostPoints;
    _y[i, j] = _dy * j - _dy * _numGhostPoints;
    // double eta = static_cast<double>(j) / (_numY - 1);
    // _y[i, j] = 0.5 * _ly * (1 + std::tanh(alpha * (eta - 0.5)) / std::tanh(alpha / 2));
  });
}