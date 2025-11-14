#pragma once

#include <memory>
#include <unordered_map>

#include "src/fieldArray/fieldArray.hpp"
#include "src/infrastructure/utilities/data.hpp"

class FieldArrayManager {
public:
  using FieldArrayType = typename std::shared_ptr<FieldArray>;
  using FieldManagerType = typename std::unordered_map<int, FieldArrayType>;

public:
  FieldArrayManager(int totalSizeX, int totalSizeY);
  ~FieldArrayManager() = default;
  
  const FieldArray& operator()(int ID) const;
  FieldArray& operator()(int ID);

  void storeOldFields();
  void storePicardOldFields();

private:
  int _totalSizeX, _totalSizeY;
  FieldManagerType _fields;
};