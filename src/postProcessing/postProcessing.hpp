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
#include "src/mesh/mesh.hpp"

class PostProcessing {
public:
  PostProcessing(std::string filename, const Mesh& mesh);
  ~PostProcessing() = default;

public:
  void registerField(std::string name, FieldArray *field);
  void write(int iteration = -1);

private:
  std::string _filename;
  const Mesh& _mesh;
  std::map<std::string, FieldArray *> _fields;
};