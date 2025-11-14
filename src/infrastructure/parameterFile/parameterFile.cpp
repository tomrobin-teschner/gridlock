#include "src/infrastructure/parameterFile/parameterFile.hpp"

ParameterFile::ParameterFile(std::string solverParameterFileName) {
  _solverParameters = readParameterFile(solverParameterFileName);
  _meshParameters = readParameterFile(_solverParameters["config_files"]["meshFile"].value_or(""));
  _boundaryConditionsParameters = readParameterFile(_solverParameters["config_files"]
    ["boundaryConditionFile"].value_or(""));
}

ParameterFile::ParameterType ParameterFile::readParameterFile(std::string parameterFileName) {
  // reading TOML input file
  std::ifstream inputFile(parameterFileName);
  if (!inputFile.is_open())
    throw std::runtime_error("Error: Could not open input file: " + parameterFileName);

  // read content into string
  std::stringstream buffer;
  buffer << inputFile.rdbuf();
  auto content = buffer.str();

  // in case file could be read, parse it as a TOML file
  return toml::parse(content);
}

ParameterFile::ParameterType ParameterFile::getSolverParameters() {
  return _solverParameters;
};

ParameterFile::ParameterType ParameterFile::getMeshParameters() {
  return _meshParameters;
};

ParameterFile::ParameterType ParameterFile::getBoundaryConditionsParameters() {
  return _boundaryConditionsParameters;
};
