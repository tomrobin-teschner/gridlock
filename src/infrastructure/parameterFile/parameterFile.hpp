#pragma once

#include <exception>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "toml++/toml.hpp"

class ParameterFile {
public:
  using ParameterType = typename toml::parse_result;

public:
  ParameterFile(std::string solverParameterFileName);
  ~ParameterFile() = default;

  ParameterType getSolverParameters();
  ParameterType getMeshParameters();
  ParameterType getBoundaryConditionsParameters();

private:
  ParameterType readParameterFile(std::string parameterFileName);

private:
  ParameterType _solverParameters;
  ParameterType _meshParameters;
  ParameterType _boundaryConditionsParameters;
};