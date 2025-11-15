#include "src/infrastructure/parameters/parameters.hpp"

Parameters::Parameters(std::filesystem::path solverParameterFileName) {
  // read the config file first, which oints to solver specific input files to read
  _config = parseFile(solverParameterFileName);

  // now get filepaths to specific settings
  auto solver = toml::find<std::string>(_config, "inputFiles", "solverFile");
  auto mesh = toml::find<std::string>(_config, "inputFiles", "meshFile");
  auto bc = toml::find<std::string>(_config, "inputFiles", "boundaryConditionFile");
  
  // parse files
  _solverParameters = parseFile(std::filesystem::path(solver));
  _meshParameters = parseFile(std::filesystem::path(mesh));
  _boundaryConditionsParameters = parseFile(std::filesystem::path(bc));
}

Parameters::ParameterType Parameters::parseFile(std::filesystem::path filePath) {
  return toml::parse(filePath.string());
}
