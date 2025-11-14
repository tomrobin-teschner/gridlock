#include "src/boundaryConditions/boundaryConditions.hpp"

BoundaryConditions::BoundaryConditions(Mesh &mesh, FieldArrayManager fields, const nlohmann::json &bcParameters)
  : _mesh(mesh), _fields(fields), _bcParameters(bcParameters) { }

void BoundaryConditions::updateGhostPoints(int ID) {
  updateEastGhostPoints(ID);
  updateWestGhostPoints(ID);
  updateNorthGhostPoints(ID);
  updateSouthGhostPoints(ID);
}

void BoundaryConditions::updateGhostPoints(std::vector<int> IDs) {
  for (const auto &ID : IDs) {
    updateEastGhostPoints(ID);
    updateWestGhostPoints(ID);
    updateNorthGhostPoints(ID);
    updateSouthGhostPoints(ID);
  }
}

void BoundaryConditions::updateEastGhostPoints(int ID) {
  auto eastBCs = _bcParameters["boundaries"]["east"][pvNames[ID]];
  bool isDirichletEast = eastBCs[0] == "dirichlet" ? true : false;
  double valueEast = eastBCs[1];
  _mesh.loop().eastBC([this, isDirichletEast, valueEast, ID](int i, int j) {
    auto spacing = _mesh.x(i, j) - _mesh.x(i - 1, j);
    
    if (isDirichletEast) {
      this->_fields(ID)[i, j] = valueEast;
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i + 1 + k, j] = dirichlet(valueEast, this->_fields(ID)[i - 1 - k, j]);
      }   
    } else if (!isDirichletEast) {
      this->_fields(ID)[i, j] = this->_fields(ID)[i - 1, j];
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i + 1 + k, j] = neumann(valueEast, this->_fields(ID)[i - 1 - k, j], spacing);
      }   
    }  
  });
}

void BoundaryConditions::updateWestGhostPoints(int ID) {
  auto westBCs = _bcParameters["boundaries"]["west"][pvNames[ID]];
  bool isDirichletWest = westBCs[0] == "dirichlet" ? true : false;
  double valueWest = westBCs[1];  
  _mesh.loop().westBC([this, isDirichletWest, valueWest, ID](int i, int j) {
    auto spacing = _mesh.x(i + 1, j) - _mesh.x(i, j);
    
    if (isDirichletWest) {
      this->_fields(ID)[i,j] = valueWest;
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i - 1 - k, j] = dirichlet(valueWest, this->_fields(ID)[i + 1 + k, j]);
      }   
    } else if (!isDirichletWest) {
      this->_fields(ID)[i,j] = this->_fields(ID)[i + 1, j];
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i - 1 - k, j] = neumann(valueWest, this->_fields(ID)[i + 1 + k, j], spacing);
      }   
    }
  });
}

void BoundaryConditions::updateNorthGhostPoints(int ID) {
  auto northBCs = _bcParameters["boundaries"]["north"][pvNames[ID]];
  bool isDirichletNorth = northBCs[0] == "dirichlet" ? true : false;
  double valueNorth = northBCs[1];
  _mesh.loop().northBC([this, isDirichletNorth, valueNorth, ID](int i, int j) {
    auto spacing = _mesh.y(i, j) - _mesh.y(i, j - 1);
    
    if (isDirichletNorth) {
      this->_fields(ID)[i,j] = valueNorth;
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i, j + 1 + k] = dirichlet(valueNorth, this->_fields(ID)[i, j - 1 - k]);
      }   
    } else if (!isDirichletNorth) {
      this->_fields(ID)[i,j] = this->_fields(ID)[i, j - 1];
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i, j + 1 + k] = neumann(valueNorth, this->_fields(ID)[i, j - 1 - k], spacing);
      }   
    }    
  });
}

void BoundaryConditions::updateSouthGhostPoints(int ID) {
  auto southBCs = _bcParameters["boundaries"]["south"][pvNames[ID]];
  bool isDirichletSouth = southBCs[0] == "dirichlet" ? true : false;
  double valueSouth = southBCs[1];
  _mesh.loop().southBC([this, isDirichletSouth, valueSouth, ID](int i, int j) {
    auto spacing = _mesh.y(i, j + 1) - _mesh.y(i, j);
    
    if (isDirichletSouth) {
      this->_fields(ID)[i,j] = valueSouth;
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i, j - 1 - k] = dirichlet(valueSouth, this->_fields(ID)[i, j + 1 + k]);
      }   
    } else if (!isDirichletSouth) {
      this->_fields(ID)[i,j] = this->_fields(ID)[i, j + 1];
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i, j - 1 - k] = neumann(valueSouth, this->_fields(ID)[i, j + 1 + k], spacing);
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