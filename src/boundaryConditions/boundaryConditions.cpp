#include "src/boundaryConditions/boundaryConditions.hpp"

BoundaryConditions::BoundaryConditions(const FieldType &x, const FieldType &y, MeshLooper &looper,
  const nlohmann::json &bcParameters) : _x(x), _y(y), _looper(looper), _bcParameters(bcParameters),
  _numGhostPoints(looper.getNumGhostPoints()) { }

void BoundaryConditions::updateGhostPoints(std::string name, FieldArray &field) {
  updateEastGhostPoints(name, field);
  updateWestGhostPoints(name, field);
  updateNorthGhostPoints(name, field);
  updateSouthGhostPoints(name, field);
}

void BoundaryConditions::updateEastGhostPoints(std::string name, FieldArray &field) {
  auto eastBCs = _bcParameters["boundaries"]["east"][name];
  bool isDirichletEast = eastBCs[0] == "dirichlet" ? true : false;
  double valueEast = eastBCs[1];
  _looper.eastBC([this, isDirichletEast, valueEast, &field](int i, int j) {
    auto spacing = _x[i][j] - _x[i - 1][j];
    
    if (isDirichletEast) {
      field[i, j] = valueEast;
      for (int k = 0; k < _numGhostPoints; ++k) {
        field[i + 1 + k, j] = dirichlet(valueEast, field[i - 1 - k, j]);
      }   
    } else if (!isDirichletEast) {
      field[i, j] = field[i - 1, j];
      for (int k = 0; k < _numGhostPoints; ++k) {
        field[i + 1 + k, j] = neumann(valueEast, field[i - 1 - k, j], spacing);
      }   
    }  
  });
}

void BoundaryConditions::updateWestGhostPoints(std::string name, FieldArray &field) {
  auto westBCs = _bcParameters["boundaries"]["west"][name];
  bool isDirichletWest = westBCs[0] == "dirichlet" ? true : false;
  double valueWest = westBCs[1];  
  _looper.westBC([this, isDirichletWest, valueWest, &field](int i, int j) {
    auto spacing = _x[i + 1][j] - _x[i][j];
    
    if (isDirichletWest) {
      field[i,j] = valueWest;
      for (int k = 0; k < _numGhostPoints; ++k) {
        field[i - 1 - k, j] = dirichlet(valueWest, field[i + 1 + k, j]);
      }   
    } else if (!isDirichletWest) {
      field[i,j] = field[i + 1, j];
      for (int k = 0; k < _numGhostPoints; ++k) {
        field[i - 1 - k, j] = neumann(valueWest, field[i + 1 + k, j], spacing);
      }   
    }
  });
}

void BoundaryConditions::updateNorthGhostPoints(std::string name, FieldArray &field) {
  auto northBCs = _bcParameters["boundaries"]["north"][name];
  bool isDirichletNorth = northBCs[0] == "dirichlet" ? true : false;
  double valueNorth = northBCs[1];
  _looper.northBC([this, isDirichletNorth, valueNorth, &field](int i, int j) {
    auto spacing = _y[i][j] - _y[i][j - 1];
    
    if (isDirichletNorth) {
      field[i,j] = valueNorth;
      for (int k = 0; k < _numGhostPoints; ++k) {
        field[i, j + 1 + k] = dirichlet(valueNorth, field[i, j - 1 - k]);
      }   
    } else if (!isDirichletNorth) {
      field[i,j] = field[i, j - 1];
      for (int k = 0; k < _numGhostPoints; ++k) {
        field[i, j + 1 + k] = neumann(valueNorth, field[i, j - 1 - k], spacing);
      }   
    }    
  });
}

void BoundaryConditions::updateSouthGhostPoints(std::string name, FieldArray &field) {
  auto southBCs = _bcParameters["boundaries"]["south"][name];
  bool isDirichletSouth = southBCs[0] == "dirichlet" ? true : false;
  double valueSouth = southBCs[1];
  _looper.southBC([this, isDirichletSouth, valueSouth, &field](int i, int j) {
    auto spacing = _y[i][j + 1] - _y[i][j];
    
    if (isDirichletSouth) {
      field[i,j] = valueSouth;
      for (int k = 0; k < _numGhostPoints; ++k) {
        field[i, j - 1 - k] = dirichlet(valueSouth, field[i, j + 1 + k]);
      }   
    } else if (!isDirichletSouth) {
      field[i,j] = field[i, j + 1];
      for (int k = 0; k < _numGhostPoints; ++k) {
        field[i, j - 1 - k] = neumann(valueSouth, field[i, j + 1 + k], spacing);
      }   
    }
  });
}

double BoundaryConditions::dirichlet(double phiBC, double phiInternal) {
  return 2.0 * phiBC - phiInternal;
}

double BoundaryConditions::neumann(double phiBC, double phiInternal, double spacing) {
  return spacing * phiBC + phiInternal;
}