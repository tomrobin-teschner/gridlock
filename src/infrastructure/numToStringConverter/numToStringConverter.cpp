#include "src/infrastructure/numToStringConverter/numToStringConverter.hpp"

std::string NumToStringConverter::fixedFloat(double num, int precision) {
  std::stringstream temp;
  temp << std::fixed << std::setprecision(precision) << num;
  return temp.str();
}

std::string NumToStringConverter::scientificFloat(double num, int precision) {
  std::stringstream temp;
  temp << std::scientific << std::setprecision(precision) << num;
  return temp.str();
}

std::string NumToStringConverter::integer(int num) {
  std::stringstream temp;
  temp << num;
  return temp.str();
}

std::string NumToStringConverter::time(int hh, int mm, int ss) {
  std::stringstream temp;
  temp << std::setw(2) << std::setfill('0') << hh << ":";
  temp << std::setw(2) << std::setfill('0') << mm << ":";
  temp << std::setw(2) << std::setfill('0') << ss;
  return temp.str();
}