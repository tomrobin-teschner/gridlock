#pragma once

#include <iostream>

#include <filesystem>
#include <string>

#include "toml.hpp"

class Parameters {
public:
  using ParameterType = typename toml::value;

public:
  Parameters(std::filesystem::path filePath = std::filesystem::path("input/config.toml"));
  ~Parameters() = default;

public:
  template<typename T, typename... Args>
  T solver(Args... args) { return toml::find<T>(_solverParameters, args...); }

  template<typename T, typename... Args>
  T mesh(Args... args) { return toml::find<T>(_meshParameters, args...); }

  template<typename T, typename... Args>
  T bcs(Args... args) { return toml::find<T>(_boundaryConditionsParameters, args...); }

private:
  ParameterType parseFile(std::filesystem::path filePath);

private:
  ParameterType _config;
  ParameterType _solverParameters;
  ParameterType _meshParameters;
  ParameterType _boundaryConditionsParameters;
};