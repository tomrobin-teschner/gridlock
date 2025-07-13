#pragma once

#include <iomanip>
#include <string>
#include <sstream>

class NumToStringConverter {
public:
  NumToStringConverter() = delete;
  ~NumToStringConverter() = delete;

public:
  static std::string fixedFloat(double num, int precision);
  static std::string scientificFloat(double num, int precision);
  static std::string integer(int num);
  static std::string time(int hh, int mm, int ss);
};