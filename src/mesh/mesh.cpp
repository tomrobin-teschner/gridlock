#include "src/mesh/mesh.hpp"

Mesh::Mesh(toml::parse_result meshParameters, int numGhostPoints)
  : _meshParameters(meshParameters), _numGhostPoints(numGhostPoints),
  _x(_meshParameters["mesh"]["numX"].value_or(0) + 2 * _numGhostPoints,
    _meshParameters["mesh"]["numY"].value_or(0) + 2 * _numGhostPoints),
  _y(_meshParameters["mesh"]["numX"].value_or(0) + 2 * _numGhostPoints,
    _meshParameters["mesh"]["numY"].value_or(0) + 2 * _numGhostPoints),
  _looper(_meshParameters["mesh"]["numX"].value_or(0), _meshParameters["mesh"]["numY"].value_or(0), _numGhostPoints) {

    // store parameters
  _lx = _meshParameters["mesh"]["Lx"].value_or(0.0);
  _ly = _meshParameters["mesh"]["Ly"].value_or(0.0);
  _numX = _meshParameters["mesh"]["numX"].value_or(0.0);
  _numY = _meshParameters["mesh"]["numY"].value_or(0.0);
  
  // calculate dx and dy
  _dx = _lx / (_numX - 1);
  _dy = _ly / (_numY - 1);

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