#include "src/fieldArray/fieldArrayManager.hpp"

FieldArrayManager::FieldArrayManager(const Mesh& mesh)
: _totalSizeX(mesh.numX() + 2 * mesh.numGhostPoints()), _totalSizeY(mesh.numY() + 2 * mesh.numGhostPoints()) {
  _fields.emplace(PV::U, std::make_shared<FieldArray>(_totalSizeX, _totalSizeY));
  _fields.emplace(PV::V, std::make_shared<FieldArray>(_totalSizeX, _totalSizeY));
  _fields.emplace(PV::P, std::make_shared<FieldArray>(_totalSizeX, _totalSizeY));

  _fields.emplace(PV::U_OLD, std::make_shared<FieldArray>(_totalSizeX, _totalSizeY));
  _fields.emplace(PV::V_OLD, std::make_shared<FieldArray>(_totalSizeX, _totalSizeY));
  _fields.emplace(PV::P_OLD, std::make_shared<FieldArray>(_totalSizeX, _totalSizeY));

  _fields.emplace(PV::U_PICARD_OLD, std::make_shared<FieldArray>(_totalSizeX, _totalSizeY));
  _fields.emplace(PV::V_PICARD_OLD, std::make_shared<FieldArray>(_totalSizeX, _totalSizeY));
  _fields.emplace(PV::P_PICARD_OLD, std::make_shared<FieldArray>(_totalSizeX, _totalSizeY));
}

const FieldArray& FieldArrayManager::operator()(int ID) const {
  return *_fields.at(ID);
}

FieldArray& FieldArrayManager::operator()(int ID) {
  return *_fields.at(ID);
}

void FieldArrayManager::storeOldFields() {
  _fields.at(PV::U_OLD) = _fields.at(PV::U);
  _fields.at(PV::V_OLD) = _fields.at(PV::V);
  _fields.at(PV::P_OLD) = _fields.at(PV::P);
}

void FieldArrayManager::storePicardOldFields() {
  _fields.at(PV::U_PICARD_OLD) = _fields.at(PV::U);
  _fields.at(PV::V_PICARD_OLD) = _fields.at(PV::V);
  _fields.at(PV::P_PICARD_OLD) = _fields.at(PV::P);
}