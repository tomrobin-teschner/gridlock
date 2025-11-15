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
#include "src/fieldArray/fieldArrayManager.hpp"
#include "src/mesh/mesh.hpp"
#include "src/infrastructure/parameters/parameters.hpp"

class PostProcessing {
public:
  PostProcessing(FieldArrayManager fields, Parameters params, const Mesh& mesh);
  ~PostProcessing() = default;

public:
  void registerFields(std::vector<int> IDs);
  void write(int iteration = -1);

private:
  FieldArrayManager _fields;
  std::string _filename;
  const Mesh& _mesh;
  std::vector<int> _IDs;
};