#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <memory>
#include <filesystem>

#include "src/infrastructure/utilities/data.hpp"
#include "src/fieldArray/fieldArray.hpp"

class PostProcessing {
public:
  PostProcessing(std::string filename, int numX, int numY, int numGhostPoints, FieldType &x, FieldType &y);
  ~PostProcessing() = default;

public:
  void registerField(std::string name, FieldArray *field);
  void write(int iteration = -1);

private:
  std::string _filename;
  int _numX, _numY, _numZ, _numGhostPoints;
  const FieldType &_x, &_y;
  std::map<std::string, FieldArray *> _fields;
};