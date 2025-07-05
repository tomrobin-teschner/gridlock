#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <memory>

class PostProcessor {
public:
  PostProcessor(std::string filename, int numX, int numY, int numGhostPoints, FieldType &x, FieldType &y)
    : _filename(filename), _numX(numX), _numY(numY), _numGhostPoints(numGhostPoints), _x(x), _y(y) {}
  
  void registerField(std::string name, std::shared_ptr<FieldType> field) {
    _fields[name] = field;
  }

  void write(int iteration = -1) {
    
    // construct filename
    std::stringstream file;
    file << "output/";
    if (iteration != -1) {
      file << _filename << "_" << std::setw(10) << std::setfill('0') << iteration << ".dat";
    } else {
      file << _filename << ".dat";
    }

    std::ofstream solution(file.str());
    solution << "TITLE = \"Solution\"\n";
    solution << "VARIABLES = \"x\", \"y\"";
    for (const auto &fields : _fields)
      solution << ", \"" << fields.first << "\"";
    solution << "\n";

    solution << "ZONE T=\"Solution\", I=" << _numX << ", J=" << _numY << ", F=POINT\n";

    // write coordinates and data
    solution << std::scientific << std::setprecision(6);
    for (int j = _numGhostPoints; j < _numY + _numGhostPoints; ++j)
      for (int i = _numGhostPoints; i < _numX + _numGhostPoints; ++i) {
        solution << _x[i][j] << " " << _y[i][j];
        for (const auto &fields : _fields)
          solution << " " << (*fields.second)[i][j];
        solution << "\n";
      }
  }

private:
  std::string _filename;
  int _numX, _numY, _numZ, _numGhostPoints;
  const FieldType &_x, &_y;
  std::map<std::string, std::shared_ptr<FieldType>> _fields;
};