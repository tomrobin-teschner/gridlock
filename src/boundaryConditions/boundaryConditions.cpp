#include "src/boundaryConditions/boundaryConditions.hpp"

BoundaryConditions::BoundaryConditions(Parameters params, Mesh &mesh, FieldArrayManager fields)
  : _params(params), _mesh(mesh), _fields(fields)  { }

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

bool BoundaryConditions::fullyNeumann() {
  auto eastBCType = _params.bcs<std::string>("boundaries", "east", "p", 0);
  auto westBCType = _params.bcs<std::string>("boundaries", "west", "p", 0);
  auto northBCType = _params.bcs<std::string>("boundaries", "north", "p", 0);
  auto southBCType = _params.bcs<std::string>("boundaries", "south", "p", 0);
  return (eastBCType == "neumann" && westBCType == "neumann" && northBCType == "neumann" && southBCType == "neumann");
}

void BoundaryConditions::updateEastGhostPoints(int ID) {
  auto eastBCType = _params.bcs<std::string>("boundaries", "east", pvNames[ID], 0);
  auto eastBCValue = _params.bcs<double>("boundaries", "east", pvNames[ID], 1);
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
  auto westBCType = _params.bcs<std::string>("boundaries", "west", pvNames[ID], 0);
  auto westBCValue = _params.bcs<double>("boundaries", "west", pvNames[ID], 1);
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
  auto northBCType = _params.bcs<std::string>("boundaries", "north", pvNames[ID], 0);
  auto northBCValue = _params.bcs<double>("boundaries", "north", pvNames[ID], 1);
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
  auto southBCType = _params.bcs<std::string>("boundaries", "south", pvNames[ID], 0);
  auto southBCValue = _params.bcs<double>("boundaries", "south", pvNames[ID], 1);
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