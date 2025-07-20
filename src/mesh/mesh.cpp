#include "src/mesh/mesh.hpp"

Mesh::Mesh(nlohmann::json meshParameters, int numGhostPoints)
  : _meshParameters(meshParameters), _numGhostPoints(numGhostPoints),
  _x(_meshParameters["mesh"]["numX"].get<int>() + 2 * _numGhostPoints,
    _meshParameters["mesh"]["numY"].get<int>() + 2 * _numGhostPoints),
  _y(_meshParameters["mesh"]["numX"].get<int>() + 2 * _numGhostPoints,
    _meshParameters["mesh"]["numY"].get<int>() + 2 * _numGhostPoints),
  _looper(_meshParameters["mesh"]["numX"].get<int>(), _meshParameters["mesh"]["numY"].get<int>(), _numGhostPoints) {

    // store parameters
  _lx = _meshParameters["mesh"]["Lx"];
  _ly = _meshParameters["mesh"]["Ly"];
  _numX = _meshParameters["mesh"]["numX"];
  _numY = _meshParameters["mesh"]["numY"];
  
  // calculate dx and dy
  _dx = _lx / (_numX - 1);
  _dy = _ly / (_numY - 1);

  // create mesh
  create();
}

void Mesh::create() {
  _looper.loopAll([this](int i, int j) {
    _x[i, j] = _dx * i - _dx * _numGhostPoints;
    _y[i, j] = _dy * j - _dy * _numGhostPoints;
  });
}