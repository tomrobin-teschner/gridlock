#include "src/boundaryConditions/boundaryConditions.hpp"

BoundaryConditions::BoundaryConditions(Mesh &mesh, FieldArrayManager fields, const toml::parse_result &bcParameters)
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
  std::string eastBCType = _bcParameters["boundaries"]["east"][pvNames[ID]][0].value_or("");
  double eastBCValue = _bcParameters["boundaries"]["east"][pvNames[ID]][1].value_or(0.0);
  bool isDirichletEast = eastBCType == "dirichlet" ? true : false;

  _mesh.loop().eastBC([this, isDirichletEast, eastBCValue, ID](int i, int j) {
    auto spacing = _mesh.x(i, j) - _mesh.x(i - 1, j);
    
    if (isDirichletEast) {
      this->_fields(ID)[i, j] = eastBCValue;
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i + 1 + k, j] = dirichlet(eastBCValue, this->_fields(ID)[i - 1 - k, j]);
      }   
    } else if (!isDirichletEast) {
      this->_fields(ID)[i, j] = this->_fields(ID)[i - 1, j];
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i + 1 + k, j] = neumann(eastBCValue, this->_fields(ID)[i - 1 - k, j], spacing);
      }   
    }  
  });
}

void BoundaryConditions::updateWestGhostPoints(int ID) {
  std::string westBCType = _bcParameters["boundaries"]["west"][pvNames[ID]][0].value_or("");
  double westBCValue = _bcParameters["boundaries"]["west"][pvNames[ID]][1].value_or(0.0);
  bool isDirichletWest = westBCType == "dirichlet" ? true : false;
 
  _mesh.loop().westBC([this, isDirichletWest, westBCValue, ID](int i, int j) {
    auto spacing = _mesh.x(i + 1, j) - _mesh.x(i, j);
    
    if (isDirichletWest) {
      this->_fields(ID)[i,j] = westBCValue;
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i - 1 - k, j] = dirichlet(westBCValue, this->_fields(ID)[i + 1 + k, j]);
      }   
    } else if (!isDirichletWest) {
      this->_fields(ID)[i,j] = this->_fields(ID)[i + 1, j];
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i - 1 - k, j] = neumann(westBCValue, this->_fields(ID)[i + 1 + k, j], spacing);
      }   
    }
  });
}

void BoundaryConditions::updateNorthGhostPoints(int ID) {
  std::string northBCType = _bcParameters["boundaries"]["north"][pvNames[ID]][0].value_or("");
  double northBCValue = _bcParameters["boundaries"]["north"][pvNames[ID]][1].value_or(0.0);
  bool isDirichletNorth = northBCType == "dirichlet" ? true : false;

  _mesh.loop().northBC([this, isDirichletNorth, northBCValue, ID](int i, int j) {
    auto spacing = _mesh.y(i, j) - _mesh.y(i, j - 1);
    
    if (isDirichletNorth) {
      this->_fields(ID)[i,j] = northBCValue;
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i, j + 1 + k] = dirichlet(northBCValue, this->_fields(ID)[i, j - 1 - k]);
      }   
    } else if (!isDirichletNorth) {
      this->_fields(ID)[i,j] = this->_fields(ID)[i, j - 1];
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i, j + 1 + k] = neumann(northBCValue, this->_fields(ID)[i, j - 1 - k], spacing);
      }   
    }    
  });
}

void BoundaryConditions::updateSouthGhostPoints(int ID) {
  std::string southBCType = _bcParameters["boundaries"]["south"][pvNames[ID]][0].value_or("");
  double southBCValue = _bcParameters["boundaries"]["south"][pvNames[ID]][1].value_or(0.0);
  bool isDirichletSouth = southBCType == "dirichlet" ? true : false;

  _mesh.loop().southBC([this, isDirichletSouth, southBCValue, ID](int i, int j) {
    auto spacing = _mesh.y(i, j + 1) - _mesh.y(i, j);
    
    if (isDirichletSouth) {
      this->_fields(ID)[i,j] = southBCValue;
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i, j - 1 - k] = dirichlet(southBCValue, this->_fields(ID)[i, j + 1 + k]);
      }   
    } else if (!isDirichletSouth) {
      this->_fields(ID)[i,j] = this->_fields(ID)[i, j + 1];
      for (int k = 0; k < _mesh.numGhostPoints(); ++k) {
        this->_fields(ID)[i, j - 1 - k] = neumann(southBCValue, this->_fields(ID)[i, j + 1 + k], spacing);
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