#include "src/residuals/residuals.hpp"

Residuals::Residuals(Parameters params, FieldArrayManager pv, Mesh &mesh, std::vector<int> IDs, std::string location,
  bool writeToFile) : _pv(pv), _mesh(mesh), _IDs(IDs), _writeToFile(writeToFile) {

  _start.resize(IDs.size());
  _end.resize(IDs.size());
  _norm.resize(IDs.size());
  _residual.resize(IDs.size());
  
  // set up residual file writing if requested
  if (_writeToFile) {
    _file.open("output/residual.csv");
    assert(_file.is_open() && "Could not open residual file!");
    
    _file << "time";
    for (int id = 0; id < _IDs.size(); ++id)
      _file << "," << pvNames[id];
    _file << "\n";
  }

  // store convergence thresholds for residuals
  _tolerances.resize(_IDs.size());
  _tolerances[PV::U] = params.solver<double>(location, "tolerance", "u");
  _tolerances[PV::V] = params.solver<double>(location, "tolerance", "v");
  _tolerances[PV::P] = params.solver<double>(location, "tolerance", "p");
}

Residuals::~Residuals() {
  if (_writeToFile)
    _file.close();
}

void Residuals::init() {
  for (int id = 0; id < _IDs.size(); ++id)
    _start[id] = norm(_IDs[id]);
}

Residuals::ResidualType Residuals::getResidual(int iteration) {
  for (int id = 0; id < _IDs.size(); ++id)
    _end[id] = norm(_IDs[id]);
  
  for (int id = 0; id < _IDs.size(); ++id)
  _residual[id] = std::fabs(_start[id] - _end[id]);
  
  // set normalisation factor
  if (iteration == 0)
    for (int id = 0; id < _IDs.size(); ++id)
      if (_residual[id] > 0.0) _norm[id] = _residual[id];
  
  // normalise
  for (int id = 0; id < _IDs.size(); ++id)
    _residual[id] /= _norm[id];
  
  if (_writeToFile) write(iteration);

  return _residual;
}

double Residuals::norm(int ID) {
  double norm = 0.0;
  _mesh.loop().loopWithBoundaries([&norm, ID, this](int i, int j) {
    norm += std::pow(this->_pv(ID)[i, j], 2);
  });
  return std::sqrt(norm);
};

void Residuals::write(int iteration) {
  _file << std::fixed << iteration + 1;
  for (int id = 0; id < _IDs.size(); ++id)
    _file << std::scientific << std::setprecision(5) << "," <<_residual[id];
  _file << "\n";
}