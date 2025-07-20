#include "src/residuals/residuals.hpp"

Residuals::Residuals(FieldArray *u, FieldArray *v, FieldArray *p, Mesh &mesh, bool writeToFile)
  : _u(u), _v(v), _p(p), _mesh(mesh), _writeToFile(writeToFile) {

  if (_writeToFile) {
    _file.open("output/residual.csv");
    assert(_file.is_open() && "Could not open residual file!");
    _file << "time,u,v,p" << std::endl;
  }
}

Residuals::~Residuals() {
  if (_writeToFile)
    _file.close();
}

void Residuals::init() {
  _startU = norm(_u);
  _startV = norm(_v);
  _startP = norm(_p);
}

Residuals::ResidualType Residuals::getResidual(int iteration) {
  _endU = norm(_u);
  _endV = norm(_v);
  _endP = norm(_p);
  
  // calculate residual
  _residualU = std::fabs(_startU - _endU);
  _residualV = std::fabs(_startV - _endV);
  _residualP = std::fabs(_startP - _endP);
  
  // set normalisation factor
  if (iteration == 0) {
    if (_residualU > 0.0) _normU = _residualU;
    if (_residualV > 0.0) _normV = _residualV;
    if (_residualP > 0.0) _normP = _residualP;
  }
  
  // normalise
  _residualU /= _normU;
  _residualV /= _normV;
  _residualP /= _normP;
  
  if (_writeToFile) write(iteration);

  // normalise residual and return
  return Residuals::ResidualType{_residualU, _residualV, _residualP};
}

double Residuals::norm(FieldArray *_data) {
  double norm = 0.0;
  _mesh.loop().loopWithBoundaries([&norm, _data](int i, int j) {
    norm += std::pow((*_data)[i, j], 2);
  });
  return std::sqrt(norm);
};

void Residuals::write(int iteration) {
  _file << std::fixed << iteration + 1 << ",";
  _file << std::scientific << std::setprecision(5) << _residualU << ",";
  _file << std::scientific << std::setprecision(5) << _residualV << ",";
  _file << std::scientific << std::setprecision(5) << _residualP << std::endl;
}