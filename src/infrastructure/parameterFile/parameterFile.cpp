#include "src/infrastructure/parameterFile/parameterFile.hpp"

ParameterFile::ParameterFile(std::string solverParameterFileName) {
  _solverParameters = readParameterFile(solverParameterFileName);
  _meshParameters = readParameterFile(_solverParameters["meshFile"]);
  _boundaryConditionsParameters = readParameterFile(_solverParameters["boundaryConditionFile"]);
}

ParameterFile::ParameterType ParameterFile::readParameterFile(std::string parameterFileName) {
  // reading JSON input file
  std::ifstream inputFile(parameterFileName);
  if (!inputFile.is_open())
    throw std::runtime_error("Error: Could not open input file: " + parameterFileName);

  // in case file could be read, parse it from JSON
  nlohmann::json tempParameterFile;
  inputFile >> tempParameterFile;
  inputFile.close();

  return tempParameterFile;
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
